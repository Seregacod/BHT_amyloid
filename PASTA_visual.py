import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import argparse
from pathlib import Path
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages


def find_optimal_window(diff_data, mutation_pos, max_window_size=50, threshold_factor=0.3):
    """
    Находит оптимальное окно анализа на основе изменений
    """
    # Начинаем с позиции мутации
    center = mutation_pos
    n = len(diff_data)

    # Находим пиковое изменение в окрестности мутации
    search_radius = min(max_window_size // 2, center, n - center - 1)
    search_slice = slice(max(0, center - search_radius), min(n, center + search_radius + 1))

    peak_idx = np.argmax(np.abs(diff_data[search_slice]))
    peak_pos = search_slice.start + peak_idx
    peak_val = diff_data[peak_pos]

    # Определяем порог значимости
    threshold = abs(peak_val) * threshold_factor

    # Расширяем окно влево до значимых изменений
    left = peak_pos
    while left > max(0, peak_pos - max_window_size):
        if abs(diff_data[left - 1]) < threshold:
            break
        left -= 1

    # Расширяем окно вправо до значимых изменений
    right = peak_pos
    while right < min(n - 1, peak_pos + max_window_size):
        if abs(diff_data[right + 1]) < threshold:
            break
        right += 1

    # Округляем до ближайших 5 аминокислот для удобства
    left = max(0, (left // 5) * 5)
    right = min(n - 1, ((right // 5) * 5) + 5)

    # Минимальное окно, оставляем 10 аминокислот
    if right - left < 10:
        left = max(0, center - 5)
        right = min(n - 1, center + 5)
        if right - left < 10:
            right = left + 10

    return slice(left, right + 1), peak_pos, peak_val


def extract_mutation_info(filename):
    """
    Извлекает информацию о мутации из  файла
    """
    name = Path(filename).stem
    import re
    match = re.search(r'([A-Z])(\d+)([A-Z])', name)
    if match:
        wt_aa = match.group(1)
        pos = int(match.group(2))
        mut_aa = match.group(3)
        return wt_aa, pos, mut_aa, f"{wt_aa}{pos}{mut_aa}"

    return None, None, None, name


def compare_with_reference(mutant_data, ref_data, mutant_name, mutation_info, output_dir):
    """
    Сравнивает мутант с референсом и создает графики
    """
    wt_aa, mut_pos, mut_aa, mut_code = mutation_info

    # Проверяем длины
    if len(mutant_data) != len(ref_data):
        min_len = min(len(mutant_data), len(ref_data))
        mutant_data = mutant_data[:min_len]
        ref_data = ref_data[:min_len]
        print(f"  Предупреждение: разная длина, обрезано до {min_len}")

    # Вычисляем разницу
    difference = mutant_data - ref_data

    # Находим оптимальное окно анализа
    if mut_pos is not None and mut_pos < len(difference):
        analysis_slice, peak_pos, peak_val = find_optimal_window(
            difference, mut_pos - 1  # -1 т.к. позиции в файле с 1
        )
    else:
        # Если не можем определить позицию мутации, используем весь белок
        analysis_slice = slice(0, len(difference))
        peak_pos = np.argmax(np.abs(difference))
        peak_val = difference[peak_pos]

    # Создаем графики
    fig = create_comparison_plot(
        mutant_data, ref_data, difference,
        analysis_slice, mut_code, mutant_name,
        wt_aa, mut_pos, mut_aa, peak_pos, peak_val
    )

    # Сохраняем график
    plot_path = output_dir / f"comparison_{mut_code}.png"
    fig.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close(fig)

    return analysis_slice, difference[analysis_slice], peak_val


def create_comparison_plot(mutant_data, ref_data, difference,
                           analysis_slice, mut_code, mutant_name,
                           wt_aa, mut_pos, mut_aa, peak_pos, peak_val):
    """
    Создает график сравнения
    """
    positions = np.arange(analysis_slice.start, analysis_slice.stop)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'Сравнение: {mut_code} ({mutant_name})', fontsize=16, fontweight='bold')

    # 1. Основной график сравнения
    axes[0, 0].plot(positions, mutant_data[analysis_slice], 'r-', linewidth=2,
                    marker='o', markersize=5, label=f'Мутант {mut_code}')
    axes[0, 0].plot(positions, ref_data[analysis_slice], 'b-', linewidth=2,
                    marker='s', markersize=5, label='Референс', alpha=0.7)

    # Выделяем позицию мутации если известна
    if mut_pos is not None and mut_pos - 1 in positions:
        axes[0, 0].axvline(x=mut_pos, color='g', linestyle='--',
                           linewidth=1.5, alpha=0.7, label=f'Мутация {mut_code}')

    axes[0, 0].set_title(f'Профиль агрегации (окно {analysis_slice.start + 1}-{analysis_slice.stop})')
    axes[0, 0].set_xlabel('Позиция аминокислоты')
    axes[0, 0].set_ylabel('Агрегационная склонность')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)

    # 2. График разницы
    bar_colors = ['red' if d > 0 else 'blue' for d in difference[analysis_slice]]
    axes[0, 1].bar(positions, difference[analysis_slice],
                   color=bar_colors, edgecolor='black', alpha=0.7)
    axes[0, 1].axhline(y=0, color='black', linewidth=1)
    axes[0, 1].set_title('Разница (Мутант - Референс)')
    axes[0, 1].set_xlabel('Позиция аминокислоты')
    axes[0, 1].set_ylabel('ΔАгрегация')
    axes[0, 1].grid(True, alpha=0.3, axis='y')

    # Выделяем максимальную дельту
    max_diff_idx = np.argmax(np.abs(difference[analysis_slice]))
    max_diff_pos = positions[max_diff_idx]
    max_diff_val = difference[analysis_slice][max_diff_idx]

    if mut_pos is not None:
        mut_idx = mut_pos - 1 - analysis_slice.start
        if 0 <= mut_idx < len(difference[analysis_slice]):
            mut_diff = difference[analysis_slice][mut_idx]
            axes[0, 1].text(mut_pos, mut_diff,
                            f'{mut_diff:+.2e}', ha='center',
                            fontsize=8, fontweight='bold')

    # 3. Полный профиль
    axes[1, 0].plot(mutant_data, 'r-', linewidth=0.5, alpha=0.5, label='Мутант')
    axes[1, 0].plot(ref_data, 'b-', linewidth=0.5, alpha=0.5, label='Референс')
    axes[1, 0].axvspan(analysis_slice.start, analysis_slice.stop - 1,
                       alpha=0.2, color='yellow', label='Окно анализа')
    axes[1, 0].set_title('Полный профиль')
    axes[1, 0].set_xlabel('Позиция аминокислоты')
    axes[1, 0].set_ylabel('Агрегационная склонность')
    axes[1, 0].legend(loc='upper right')
    axes[1, 0].grid(True, alpha=0.3)
    plt.tight_layout()
    return fig


def main():
    parser = argparse.ArgumentParser(description='Сравнение профилей агрегации мутантов с референсом')
    parser.add_argument('input_dir', help='Директория с файлами мутантов')
    parser.add_argument('reference_file', help='Файл референса')
    parser.add_argument('--output', '-o', default='comparison_results',
                        help='Выходная директория (по умолчанию: comparison_results)')

    args = parser.parse_args()

    # Создаем выходную директорию
    output_dir = Path(args.output)
    output_dir.mkdir(exist_ok=True)

    # Загружаем референс
    print(f"Загрузка референса: {args.reference_file}")
    try:
        ref_data = np.loadtxt(args.reference_file)
        print(f"  Загружено {len(ref_data)} значений")
    except Exception as e:
        print(f"Ошибка загрузки референса: {e}")
        return

    # Находим файлы мутантов
    mutant_files = glob.glob(os.path.join(args.input_dir, "*.txt"))
    print(f"\nНайдено файлов мутантов: {len(mutant_files)}")

    if len(mutant_files) == 0:
        print("Не найдены файлы мутантов!")
        return

    # Результаты
    results = []

    # Обрабатываем каждый мутант
    for i, mutant_file in enumerate(mutant_files, 1):
        mutant_name = Path(mutant_file).stem
        print(f"\n[{i}/{len(mutant_files)}] Обработка: {mutant_name}")
        try:
            # Загружаем данные мутанта
            mutant_data = np.loadtxt(mutant_file)

            # Извлекаем информацию о мутации
            mutation_info = extract_mutation_info(mutant_file)

            # Сравниваем с референсом
            analysis_slice, diff_window, peak_diff = compare_with_reference(
                mutant_data, ref_data, mutant_name, mutation_info, output_dir
            )

            # Сохраняем результаты
            wt_aa, mut_pos, mut_aa, mut_code = mutation_info
            result = {
                'Файл': mutant_name,
                'Мутация': mut_code if mut_code else 'N/A',
                'Позиция': mut_pos if mut_pos else 'N/A',
                'Окно_начала': analysis_slice.start + 1,
                'Окно_конца': analysis_slice.stop,
                'Размер_окна': analysis_slice.stop - analysis_slice.start,
                'Макс_дельта': float(np.max(np.abs(diff_window))),
                'Позиция_макс_дельты': analysis_slice.start + np.argmax(np.abs(diff_window)) + 1,
                'Средняя_дельта': float(np.mean(diff_window)),
                'Пиковая_дельта': float(peak_diff),
                'Позиция_пика': np.argmax(np.abs(diff_window)) + analysis_slice.start + 1
            }
            results.append(result)
            print(f"  Анализ окна: {analysis_slice.start + 1}-{analysis_slice.stop}")
            print(f"  Максимальная |Δ|: {result['Макс_дельта']:.2e} (поз. {result['Позиция_макс_дельты']})")
        except Exception as e:
            print(f"  Ошибка обработки: {e}")

    # Сохраняем сводную таблицу
    if results:
        df = pd.DataFrame(results)
        summary_file = output_dir / "summary_results.csv"
        df.to_csv(summary_file, index=False, encoding='utf-8-sig')

        # Создаем сводный отчет
        print(f"\n{'=' * 60}")
        print(f"ОБРАБОТКА ЗАВЕРШЕНА!")
        print(f"{'=' * 60}")
        print(f"Всего обработано: {len(results)} мутантов")
        print(f"Результаты сохранены в: {output_dir}")
        print(f"Сводная таблица: {summary_file}")


if __name__ == "__main__":
    main()
from Bio import SeqIO
from selenium import webdriver
from selenium.webdriver.edge.service import Service
from selenium.webdriver.common.by import By
import time
import os
import urllib.request


def submit_all_fasta_at_once(fasta_file, driver_path, save_dir="PASTA_Results"):
    """
    Отправляет ВСЕ FASTA последовательности из файла в одной форме PASTA
    Возвращает один batch.tar архив со всеми результатами
    """
    os.makedirs(save_dir, exist_ok=True)
    print(f"Чтение {fasta_file}...")
    try:
        records = list(SeqIO.parse(fasta_file, "fasta"))
    except FileNotFoundError:
        print(f"Файл не найден: {fasta_file}")
        return None
    except Exception as e:
        print(f"Ошибка чтения FASTA: {e}")
        return None

    print(f"Найдено {len(records)} последовательностей")

    combined_fasta = ""
    total_length = 0

    for i, record in enumerate(records, 1):
        combined_fasta += f">{record.description}\n{record.seq}\n\n"
        total_length += len(record.seq)

        print(f"  {i}. {record.id}: {len(record.seq)} аа")

    print(f"\n📋 Итого: {len(records)} белков, {total_length} аминокислот")

    # 3. Настройка браузера
    options = webdriver.EdgeOptions()
    options.add_argument('--headless')
    options.add_argument('--start-maximized')

    driver = webdriver.Edge(
        service=Service(driver_path),
        options=options
    )

    try:
        print("\n Открываю PASTA 2.0...")
        driver.get('http://old.protein.bio.unipd.it/pasta2/')
        time.sleep(3)
        textarea = driver.find_element(By.ID, "sequence")
        print(f"Ввожу {len(records)} последовательностей...")
        textarea.clear()
        textarea.send_keys(combined_fasta)
        submit = driver.find_element(By.CSS_SELECTOR, "input[name='Submit Query']")
        submit.click()

        # 6. Умное ожидание (дольше, так как много последовательностей)

        estimated_time = max(40, 30 + len(records) * 5 + total_length * 0.2)
        estimated_time = min(estimated_time, 100)

        # Ждем с индикацией прогресса
        for i in range(int(estimated_time)):
            time.sleep(1)
            if i % 10 == 0:
                print(f"  {i} сек...")

        print("\n Ищу batch.tar архив...")
        current_url = driver.current_url
        print(f"Текущая страница: {current_url}")
        links = driver.find_elements(By.TAG_NAME, "a")
        print(f"Найдено ссылок на странице: {len(links)}")

        tar_url = None
        tar_text = ""

        for link in links:
            href = link.get_attribute("href")
            text = link.text.strip()

            if href:
                # 1. Прямая ссылка на batch.tar
                if "batch.tar" in href:
                    tar_url = href
                    tar_text = text if text else "batch.tar"
                    break

                # 2. Ссылка с текстом "batch" или "tar"
                if ("batch" in href.lower() or "tar" in href.lower()) and not href.endswith('.html'):
                    tar_url = href
                    tar_text = text if text else "архив"

        if not tar_url:
            # Если не нашли batch.tar, показываем все ссылки для отладки
            print("\n batch.tar не найден. Все ссылки на странице:")
            for i, link in enumerate(links[:20]):
                href = link.get_attribute("href")
                text = link.text[:30]
                if href:
                    print(f"  {i + 1}. '{text}' -> {href[:80]}...")

            return None

        print(f"Найден: {tar_text}")
        print(f"URL: {tar_url}")

        # Создаем имя файла на основе исходного FASTA файла
        base_name = os.path.splitext(os.path.basename(fasta_file))[0]
        filename = f"{base_name}_ALL_batch.tar"

        if '.gz' in tar_url:
            filename += '.gz'

        filepath = os.path.join(save_dir, filename)

        print(f"Скачиваю {filename}...")
        # Скачиваем файл
        urllib.request.urlretrieve(tar_url, filepath)

        # Проверяем
        if os.path.exists(filepath):
            size_kb = os.path.getsize(filepath) // 1024
            print(f"Скачан: {filename}")
            print(f"Размер: {size_kb} KB")

            if size_kb < 10:
                print("Внимание: архив очень маленький, возможно пустой")
            elif size_kb > 100:
                print(f"Архив содержит все {len(records)} результатов!")

            return filepath
        else:
            print("Ошибка: файл не скачан")
            return None

    except Exception as e:
        print(f"Ошибка: {type(e).__name__}: {e}")
        return None
    finally:
        driver.quit()


def extract_and_analyze_tar(tar_path, extract_dir=None):
    """
    Распаковывает batch.tar и показывает содержимое
    """
    import tarfile

    if extract_dir is None:
        extract_dir = os.path.splitext(tar_path)[0] + "_extracted"

    os.makedirs(extract_dir, exist_ok=True)

    try:
        print(f"\nРаспаковываю {os.path.basename(tar_path)}...")

        # Определяем формат
        if tar_path.endswith('.tar.gz') or tar_path.endswith('.tgz'):
            mode = 'r:gz'
        elif tar_path.endswith('.tar'):
            mode = 'r'
        else:
            print(f"Неизвестный формат: {tar_path}")
            return None
    except Exception as e:
        print(f"Ошибка распаковки: {e}")
        return None


if __name__ == "__main__":
    # Пути (измените под себя)
    FASTA_FILE = "RPS2_mut128.fasta"  # Ваш FASTA файл
    DRIVER_PATH = r'C:\Users\USER\msedgedriver.exe'  # Путь до exe файла поисковика
    SAVE_DIR = r'C:\Users\USER\PASTA_ALL_Results'  # Путь для сохранения

    # Запускаем
    tar_file = submit_all_fasta_at_once(
        fasta_file=FASTA_FILE,
        driver_path=DRIVER_PATH,
        save_dir=SAVE_DIR
    )

    if tar_file:
        print(f"\nУСПЕХ! Все последовательности обработаны.")
        print(f" Архив: {tar_file}")

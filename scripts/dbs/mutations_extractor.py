import pandas as pd
import numpy as np
import os
import pathlib
from argparse import ArgumentParser


def extract_mutation(file, protein):

    df = pd.read_csv(file, header=0, sep=',')
    df['aa'] = df['aa'].astype(int)
    mut_cols = df.columns[1:]
    mutations = [mut.split('_')[1] for mut in mut_cols]

    mutation_df = pd.DataFrame.from_dict({'protein': [protein]*len(mut_cols),
                   'mutation': mutations,
                   'score': [0]*len(mut_cols),
                   'average_score': [0]*len(mut_cols)}, orient='columns')
    mutation_df_new = mutation_df.copy()
    
    aa_number = int(mut_cols[1].split('_')[1][1:-1])
    start, stop = aa_number - 2, aa_number + 2
    for mut in mut_cols:
        mutation = mut.split('_')[1]
        score = df.loc[df['aa'] == aa_number, mut].values
        average_score = df.loc[start-1:stop-1, mut].mean()
        mutation_df.loc[mutation_df['mutation'] == mutation, 'score'] = score
        mutation_df.loc[mutation_df['mutation'] == mutation, 'average_score'] = average_score
        mutation_df_new = mutation_df.copy()

    return mutation_df_new

def prepare_table(directory, output):
    path = str(pathlib.Path(directory).absolute())
    final_df = pd.DataFrame()
    for file in os.listdir(path):
        protein = file.split('_')[0]

        protein_df = extract_mutation(file, protein)
        final_df = pd.concat([final_df, protein_df], axis=0)

    final_df.to_csv(output, sep='\t', header=True)

if __name__ == '__main__':

    parser = ArgumentParser(
        description='Parse mutation data for all proteins into one dataframe'
    )
    parser.add_argument('-i', '--input_dir', required=True,
                        help='Input directory with mutational data')
    parser.add_argument('-o', '--output', default='mutations_proteins.tsv',
                        help='Output TSV file')
    
    args = parser.parse_args()

    prepare_table(args.input_dir, args.output)
    print(f'Parsed all data. Saved data into {args.output} file')




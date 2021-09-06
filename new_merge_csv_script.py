#!usr/bin/env python
# -*- coding: utf-8 -*-
#  python new_merge_csv_script.py -p Data/bacteria_splitted_kmers -d Results/bacteria_splitted_kmers_cvs -sd chromosomes -e csv -f spc_names_tests.txt -k 9
import os
import sys
import time
import glob
import argparse
from collections import defaultdict
from functools import reduce
import numpy as np
import pandas as pd
from system_utils import make_me_a_folder, get_files


def get_species_name(filenames):
    species = []
    with open(filenames, 'r') as fh:
        for line in fh:
            sp = line.strip()
            species.append(sp)
    return species


def get_spc_df_merged(spc_name, filenames_dict):
    dfs = []
    for filename in filenames_dict[spc_name]:
        dfs.append(pd.read_csv(filename, dtype={'kmer': str, 'count': np.int32}))
    df_merged = pd.DataFrame(reduce(lambda left, right: pd.merge(left,
                                                                 right,
                                                                 on='kmer'),
                                    dfs).set_index('kmer').sort_index().mean(axis=1)).reset_index()
    return df_merged.rename(columns={'kmer': 'kmers', 0: 'counts'}), len(dfs)


def get_csvs_to_df_merge(spc_name, filenames_dict):
    dfs = []
    for filename in filenames_dict[spc_name]:
        dfs.append(pd.read_csv(filename, dtype={'kmer': str, 'count': np.int32}))
    dfs = [df.set_index("kmer", drop=True) for df in dfs]
    merged = pd.concat(dfs, axis=1, keys=range(len(dfs)), copy=False).sort_index().mean(axis=1).reset_index()
    return merged.rename(columns={'index': 'kmers', 0: 'counts'}), len(dfs)


def concatenation_dfs(spc_name, dict_csvs_paths):
    csvs = dict_csvs_paths[spc_name]
    concat_df = pd.DataFrame()
    for csv in csvs:
        temp_df = pd.read_csv(csv, dtype={'kmer': str,
                                          'count': np.int32}).set_index("kmer",
                                                                        drop=True)
        concat_df = pd.DataFrame(pd.concat([concat_df,
                                            temp_df],
                                           axis=1,
                                           join='outer',
                                           copy=False).sort_index().mean(axis=1))
        del temp_df
    concat_df[0] = concat_df[0].astype('float32')
    return concat_df.reset_index().rename(columns={'index': 'kmers', 0: 'counts'}), len(csvs)


def parse_arguments():
    """Parse the command line arguments to the genome_palindromes_analysis script.
    Sets up the argparse command-line parser and calls it. These options can be accessed
    using opt.option. For example opt.alp stores the alphabet provided.
    """
    parser = argparse.ArgumentParser(description='A script to analyze palindromes and kmer in genomes.')
    parser.add_argument('-p',
                        '--path',
                        metavar='path',
                        type=str,
                        required=True,
                        dest='path',
                        help='Directory of the files.')
    parser.add_argument('-sd',
                        '--sub_dir',
                        type=str,
                        dest='sub_dir',
                        help='subdirectory of the files.')
    parser.add_argument('-d',
                        '--dir_out',
                        type=str,
                        dest='dir_out',
                        help='directory name for output files.')
    parser.add_argument('-e',
                        '--ext',
                        type=str,
                        dest='extension',
                        help='Final extension for the name of the files(i.e, .csv).')
    parser.add_argument('-f',
                        '--file',
                        type=str,
                        dest='filenames',
                        help='species list  as txt.')
    parser.add_argument('-k',
                        '--kmer_len',
                        type=int,
                        action='store',
                        dest='k',
                        help='Specify the kmer/palindrome length.')
    return parser.parse_args()


def main():
    """Parses options from the command line.
    This script receive a path to the csv files, a directory name to save the
    merged csv files and a list of species names.
    The final result is a merged csv file with the mean of all kmer counts in the
    output directory.
    """
    print('Starting to process the script merge_csvs.py\n')
    cwd = os.getcwd()
    print(f'The working directory: {cwd}\n')
    start_time = time.process_time()
    opt = parse_arguments()
    dir_name = opt.path
    sub_dir = opt.sub_dir
    dir_out = opt.dir_out
    k = opt.k
    sub_sub_dir = f'kmer{k}'
    filenames = opt.filenames
    spc_names = get_species_name(filenames)
    extension = opt.extension

    cnt = 0
    all_csvs = defaultdict(list)
    for name in spc_names:
        print(f'Getting csv files from species: {name}')
        all_csvs[name] = all_csvs.get(name, [])
        full_name_dir = os.path.join(dir_name, name, sub_dir)
        for f in glob.glob(full_name_dir + f'/*.{extension}'):
            all_csvs[name].append(f)

    for name in spc_names:
        name_dir_out = os.path.join(dir_out, name, sub_dir, sub_sub_dir)
        csv_name = f'{name}_{sub_sub_dir}_merged.csv'
        full_name = os.path.join(name_dir_out, csv_name)
        dfs, len_csvs = concatenation_dfs(name, all_csvs)
        # print(f'The number of {k}-kmers in these csv file is: {dfs["kmers"].size}')
        print(f'The number of csv files for species {name} is {len_csvs}\n')
        if not os.path.exists(name_dir_out):
            os.makedirs(name_dir_out)
        dfs.to_csv(f'{full_name}', index=False)
        print(f'The final csv merged file was saved in {full_name}\n')

        cnt += 1
    end_time = time.process_time()
    total_time = end_time - start_time
    print(f'Total species: {len(spc_names)}')
    print(f'Total files: {cnt}')
    print(f'Total time for running the script: {total_time}')
    print('Script finished!')


if __name__ == "__main__":
    sys.exit(main())

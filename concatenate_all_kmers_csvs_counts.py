#!usr/bin/env python
# -*- coding: utf-8 -*-
# Usage:
# python
# concatenate_all_kmers_csvs_counts.py - d
# Data / bacteria_splitted_kmers_cvs - sd
# chromosomes - ssd
# kmer - do
# Results - sdo
# Merged - e
# csv - f
# test.txt

import os
import sys
import glob
import time
import argparse
from collections import defaultdict
import pandas as pd
import numpy as np


def generate_new_df(filename):
    return pd.read_csv(filename, dtype={'kmer': str, 'count': np.int32})


def get_df_from_csv_lists(spc_name, csv_paths):
    dfs = [generate_new_df(csv) for csv in csv_paths[spc_name]]
    return dfs


def parse_arguments():
    """Parse the command line arguments to merge csvs script.
    Sets up the argparse command-line parser and calls it. These options can be accessed
    using opt.option.
    """
    parser = argparse.ArgumentParser(description='A script to concatenate csv files.')
    parser.add_argument('-d',
                        '--dir',
                        metavar='dir',
                        type=str,
                        required=True,
                        dest='dir_in',
                        help='Directory of the files.')
    parser.add_argument('-sd',
                        '--sub_dir',
                        type=str,
                        dest='sub_dir_in',
                        help='subdirectory of the files.')
    parser.add_argument('-ssd',
                        '--ssub_dir',
                        type=str,
                        dest='ssub_dir',
                        help='subdirectory of the files.')
    parser.add_argument('-do',
                        '--dir_out',
                        type=str,
                        dest='dir_out',
                        help='directory name for output files.')
    parser.add_argument('-sdo',
                        '--subdir_out',
                        type=str,
                        dest='subdir_out',
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
    return parser.parse_args()


def main():
    """Parses options from the command line.
    This script receive a path to the csv files, a directory name to save the
    merged csv files and a list of species names.
    The final result is a concatenated csv file with the mean of all kmer counts in the
    output directory.
    """
    print('\nStarting to process the script concatenate_all_kmers_csvs_counts.py\n')
    cwd = os.getcwd()
    print(f'The working directory: {cwd}\n')
    start_time = time.process_time()
    opt = parse_arguments()
    dir_name = opt.dir_in  # 'Data/bacteria_spplited'
    sub_dir_in = opt.sub_dir_in  # 'chromosomes'
    ssub_dir = opt.ssub_dir  # 'kmer*'
    dir_out = opt.dir_out  # 'Results/test'
    subdir_out = opt.subdir_out  # Merged
    filetext = opt.filenames  # 'test_sps.txt'
    extension = opt.extension

    species = [name.strip() for name in open(filetext, 'r')]

    all_csvs = defaultdict(list)
    for name in species:
        all_csvs[name] = all_csvs.get(name, [])
        full_name_dir = os.path.join(dir_name, name, sub_dir_in, f'{ssub_dir}*')
        for f in glob.glob(full_name_dir + f'/*.{extension}'):
            all_csvs[name].append(f)

    for name, filename in all_csvs.items():
        data = get_df_from_csv_lists(name, all_csvs)
        data_frame = pd.concat(data, axis=0, ignore_index=True)
        data_frame['counts'] = data_frame['counts'].astype(int)
        fullname = os.path.join(dir_out, name, sub_dir_in, subdir_out)
        csv_name = f'{name}_allmerged.csv'
        if not os.path.exists(fullname):
            os.makedirs(fullname)
        data_frame.to_csv(f'{fullname}/{csv_name}', index=False)
        print(f'The final csv merged file was saved in {fullname}\n')

    end_time = time.process_time()
    total_time = end_time - start_time
    print(f'Total species: {len(species)}')
    print(f'Total time for running the script: {total_time}')
    print('Script finished!\n')


if __name__ == "__main__":
    sys.exit(main())

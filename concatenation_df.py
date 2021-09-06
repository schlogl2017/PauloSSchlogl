# !usr/bin/env python -*- coding: utf-8 -*- usage: concatenation_df.py [-h] -p path [-sd SUB_DIR] [-d DIR_OUT] [-ssd
# SUB_SUB_DIR] [-e EXTENSION] [-f FILENAMES] [-k NUM_CSVS] concatenation_df.py: error: the following arguments are
# required: -p/--path python concatenation_df.py -p Results/splitted_genomes_merged -sd chromosomes -d Results/test
# -ssd kmer -e csv -f spc_test.txt -k 8


import os
import sys
import time
import glob
import argparse
from collections import defaultdict
import numpy as np
import pandas as pd
from system_utils import make_me_a_folder, find_csv_filenames, get_species_name, get_csvs_to_df_concatenation


def parse_arguments():
    """Parse the command line arguments to the concatentenation_dfs script.
    Sets up the argparse command-line parser and calls it. These options can be accessed
    using opt.option. For example opt.alp stores the alphabet provided.
    """
    parser = argparse.ArgumentParser(description='A script concatenate all kmers csvs counts form genomes.')
    parser.add_argument('-p',
                        '--path',
                        metavar='path',
                        type=str,
                        required=True,
                        dest='path',
                        help='Directory of the csv files.')
    parser.add_argument('-sd',
                        '--sub_dir',
                        type=str,
                        dest='sub_dir',
                        help='subdirectory of the files. Eg. chromosomes or plasmids')
    parser.add_argument('-d',
                        '--dir_out',
                        type=str,
                        dest='dir_out',
                        help='directory name to the output file. Eg. Results')
    parser.add_argument('-ssd',
                        '--sub_sub_dir',
                        type=str,
                        dest='sub_sub_dir',
                        help='sub directory name to the output files. Eg. kmer or plasmids')
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
                        '--num_csvs',
                        type=int,
                        dest='num_csvs',
                        help='Specify the number of csv files.')
    return parser.parse_args()


def main():
    """Parses options from the command line.
    This script receive a path to the csv files, a directory name to save the
    merged csv files and a list of species names.
    The final result is a concatenated csv file with all kmers counts
    from 2-n in the output directory.
    """
    print('Starting to process the script concatenation_df.py\n')
    cwd = os.getcwd()
    print(f'The working directory: {cwd}\n')
    start_time = time.process_time()
    arg = parse_arguments()
    text_file = arg.filenames
    spc_names = get_species_name(text_file)
    suffix = arg.extension  # .csv
    dir_out = arg.dir_out  # Results/test
    dir_name = arg.path  # Results/splitted_genomes_merged
    sub_dir = arg.sub_dir  # chromosomes
    num_csvs = arg.num_csvs  # number files to concatenate
    sub_sub_dir = arg.sub_sub_dir  # kmer{i}

    num_spc = 0
    csvs = 0
    for spc in spc_names:
        filenames = find_csv_filenames(dir_name, spc, sub_dir, sub_sub_dir, num_csvs + 2, suffix)
        df, csvs_len = get_csvs_to_df_concatenation(spc, filenames)
        full_name_dir = os.path.join(dir_out, spc, sub_dir, sub_sub_dir)
        csv_name = os.path.join(full_name_dir, spc)
        if not os.path.exists(full_name_dir):
            os.makedirs(full_name_dir)
        df.to_csv(f'{csv_name}.concatenated.{suffix}', index=False)
        print(f'The final csv merged file was saved in {full_name_dir}\n')
        num_spc += 1
        csvs += csvs_len

    end_time = time.process_time()
    total_time = end_time - start_time
    print(f'Total species: {num_spc}')
    print(f'Total files: {csvs}')
    print(f'Total time for running the script: {total_time}')
    print('Script finished!')


if __name__ == "__main__":
    sys.exit(main())

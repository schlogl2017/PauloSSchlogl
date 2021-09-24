#!/usr/bin/env python
# -*- coding: utf-8 -*-
# usage: genome_length.py [-h] -di path [-do DIR_OUT] [-sd SUB_DIR]
#                         [-e EXTENSION]
# genome_length.py: error: the following arguments are required: -di/--dir_in
import os
import sys
import gzip
import glob
from time import time
from termcolor import colored
from collections import defaultdict
import argparse
import numpy as np
import pandas as pd
from system_utils import make_me_a_folder, get_all_fasta
from fasta_parser import parse_fasta


def get_genome_lenght(filenames):
    gen_len = defaultdict(dict)

    for filename in filenames:
        genus = filename.split('/')[2]
        length = 0
        for name, seq in parse_fasta(filename):
            length += len(seq)
            gen_len[genus][name] = len(seq)

    return gen_len


def get_mean_genome_lengths(data_dict):
    mean_dict = defaultdict(int, [(name, 0) for name in data_dict.keys()])
    for genus in mean_dict:
        mean_dict[genus] = int(np.array(list(data_dict[genus].values())).mean())
    return mean_dict


def parse_arguments():
    """Parse the command line arguments to the genome_length script.
    Sets up the argparse command-line parser and calls it. These args can be accessed
    using args.args.
    The resulting results are csv files from each genome contained in the genus directory
    with a list of all kmers generated from these genomes (chromosomes and plasmids).
    """
    parser = argparse.ArgumentParser(description="""A script to count all bases and GC content
    from bacterial genomes/plasmids in the fasta files.""",
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-di',
                        '--dir_in',
                        metavar='path',
                        type=str,
                        required=True,
                        dest='dir_in',
                        help='Directory root. In my case the name is conjugated with a subdir')

    parser.add_argument('-do',
                        '--dir_out',
                        type=str,
                        dest='dir_out',
                        help='directory name for output files.')

    parser.add_argument('-sd',
                        '--sub_dir',
                        type=str,
                        dest='sub_dir',
                        help='Name for a subdirectory, ex., Chromosomes.')

    parser.add_argument('-e',
                        '--extension',
                        type=str,
                        dest='extension',
                        help='Name representing the type file. Ex., gz')
    return parser.parse_args()


def main():
    # starting count the staring time of the script
    start = time()
    # checking the current directory and printing it
    cwd = os.getcwd()
    print(colored(f'\nThe working directory: {cwd}\n',
                  'green',
                  attrs=['bold']))
    # passing the arguments to the script
    args = parse_arguments()
    # name of the input directory, ex. Data/Genomes_splitted
    dir_in = args.dir_in
    # name of the sub directory to save the final result
    # Chromosomes/Plasmids
    sub_dir = args.sub_dir
    # name of the root directory to save the final result
    dir_out = args.dir_out
    # name of the root directory to save the final result
    extension = args.extension
    # get the fasta files
    filenames = glob.glob(f'{dir_in}/*/{sub_dir}/*.{extension}')
    # check if the output directory exist other wise create it
    if os.path.exists(dir_out):
        print(colored('The directory to save the files already exists!',
                      'red',
                      attrs=['bold']))
        pass
    else:
        make_me_a_folder(dir_out)
    data_len = get_genome_lenght(filenames)
    data = get_mean_genome_lengths(data_len)
    df = pd.DataFrame(data.items(), columns=['Name', 'Length'])
    file_name = 'All_plasmids_length'
    full_path = os.path.join('Results', 'Length')
    if not os.path.exists(full_path):
        os.makedirs(full_path)
    df.to_csv(f'{full_path}/{file_name}.csv', index=False)
    # the final time
    end = time()
    # print some info
    print(colored(f"Total number of genus/species analyzed: {len(data)}\n.",
                  attrs=['bold']))
    print(colored(f'Total time for the script finishes: {round(end - start, 2)}.',
                  'red',
                  attrs=['bold']))
    print(colored('Done!',
                  'green',
                  attrs=['bold']))


if __name__ == "__main__":
    sys.exit(main())

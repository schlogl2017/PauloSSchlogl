#!/usr/bin/env python
# -*- coding: utf-8 -*-
# usage: gc_stats.py [-h] -di path [-do DIR_OUT] [-sd SUB_DIR]
#                    [-ssd SUB_SUB_DIR] [-e EXTENSION]
# gc_stats.py: error: the following arguments are required: -di/--dir_in
import os
import sys
import gzip
from time import time
from termcolor import colored
from collections import Counter
import argparse
import pandas as pd
from system_utils import make_me_a_folder, get_all_fasta
from fasta_parser import parse_fasta
from basic_stats_genomes import gc_cython


def gc_count_fasta(fasta_dict, name):
    """
    Function to count all n-grams/k-mers (substrings of length n or k) in a
    big string/genome.

    Inputs:
        fasta_dict - a dictionary-like object that map a word/kmer to their value,
                    in this case a full path to the files to be analyzed.
        name - a string representing a word (key) that represent a key in a
               dictionary.

    Outputs:
        gc content - a float representing the mean of the gc content from all
                     genus/species analyzed.
    """
    # get the number of files in the names directory
    num_fastas = len(fasta_dict[name])
    # initialize the counter
    gc_tot = 0
    # iterates through the list of paths
    for filename in fasta_dict[name]:
        # reads the file and parse the content
        print(f'Reading and parsing the filename {filename}')
        for n, sequence in parse_fasta(filename):
            # add the gc content from all files
            gc_tot += gc_cython(sequence)
    # returns the mean of the gc content from all files
    return name, (gc_tot / num_fastas) * 100


def parse_arguments():
    """Parse the command line arguments to the script gc_stats.py.
    Sets up the argparse command-line parser and calls it. These args can be accessed
    using args.args.
    The resulting results are csv files from each genome contained in the genus directory
    with a list of all kmers generated from these genomes (chromosomes and plasmids).
    """
    parser = argparse.ArgumentParser(description="""A script to calculate the GC content(100%)
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

    parser.add_argument('-ssd',
                        '--sub_sub_dir',
                        type=str,
                        dest='sub_sub_dir',
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
    # extension type for fasta
    sub_sub_dir = args.sub_sub_dir
    # name of the root directory to save the final result
    extension = args.extension
    # get the list of all paths to the files in the input directory
    # ex., Data/Genomes_splitted, Chromosomes, gz
    fasta_dict = get_all_fasta(dir_in, sub_dir, extension)
    # check if the output directory exist other wise create it
    if os.path.exists(dir_out):
        print(colored('The directory to save the files already exists!',
                      'red',
                      attrs=['bold']))
        pass
    else:
        make_me_a_folder(dir_out)

    # initialize the file counter
    cnt_spc = 0
    # df list
    dfs = []
    # input the file paths and print it to show where the script is doing
    for name in fasta_dict.keys():
        print(colored(f"Start working with genus {name}\n", attrs=['bold']))
        # initialize the gc counts
        gc = gc_count_fasta(fasta_dict, name)
        # create a data frame
        df_gc = pd.DataFrame(gc).T.reset_index(drop=True)
        df_gc.columns = ['Name', 'GC_%']
        dfs.append(df_gc)
        full_path = os.path.join(dir_out, 'GC_content_all')
        file_name = 'GC_content_all'
        if not os.path.exists(full_path):
            os.makedirs(full_path)
        df_final = pd.concat(dfs, ignore_index=True)
        df_final.to_csv(f'{full_path}/{file_name}_{sub_sub_dir}.csv', index=False)
        # the number of genus/species analyzed
        cnt_spc += 1
    # the final time
    end = time()
    # print some info
    print(colored(f"Total number of genus/species analyzed: {cnt_spc}\n.",
                  attrs=['bold']))
    print(colored(f'Total time for the script finishes: {round(end - start, 2)}.',
                  'red',
                  attrs=['bold']))
    print(colored('Done!',
                  'green',
                  attrs=['bold']))


if __name__ == "__main__":
    sys.exit(main())

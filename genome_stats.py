#!/usr/bin/env python
# -*- coding: utf-8 -*-
# usage: genome_stats.py [-h] -di path [-do DIR_OUT] [-sd SUB_DIR]
#                        [-e EXTENSION]
# genome_stats.py: error: the following arguments are required: -di/--dir_in
#

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
from basic_stats_genomes import gc_cython, count_bases_cython


def gc_count_fasta(fasta_dict, name):
    """
    Function to count all n-grams/k-mers (substrings of lenght n or k) in a
    big string/genome.

    Inputs:
        fasta_dict - a dictionary-like object that map a word/kmer to their value,
                    in this case a full path to the files to be analized.
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
        for name, sequence in parse_fasta(filename):
            # add the gc content from all files
            gc_tot += gc_cython(sequence)
    # returns the mean of the gc content from all files
    return (gc_tot / num_fastas) * 100


def count_bases_fasta(fasta_dict, name):
    """
    Function to count all n-grams/k-mers (substrings of lenght n or k) in a
    big string/genome.

    Inputs:
        fasta_dict - a dictionary-like object that map a word/kmer to their value,
                    in this case a full path to the files to be analized.
        name - a string representing a word (key) that represent a key in a
               dictionary.

    Outputs:
        final_counter - a dictionary-like mapping the kmers to their calculated count
                        in the input string, from a file.
        seq_length - a integer representing the mean of the lengths from all genomes
                     files in the directory.
    """
    # get the number of files in the names directory
    num_fastas = len(fasta_dict[name])
    print(f'The number of fasta files for this genus is {num_fastas}.')
    # initialize the counter
    counter = Counter()
    # get the sequence length
    seq_len = 0
    num_files = 0
    # iterates through the list of paths
    for filename in fasta_dict[name]:
        # reads the file and parse the content
        print(f'Reading and parsing the file {filename}')
        for name, sequence in parse_fasta(filename):
            print(f'Sequence length {len(sequence)}')
            seq_len += len(sequence)
            # get the counting the kmers
            cnt = count_bases_cython(sequence)
            # add the count of the current file to the counter
            counter.update(cnt)
        num_files += 1
    # to get the mean of the kmer count for all the files
    final_counter = {k: (c // num_fastas) for k, c in counter.items()}
    return final_counter


def genome_stats_in_windows(fasta_dict, name, as_overlap=False, k=20):
    """GC Content in a DNA/RNA sub-sequence length k. In
    overlapp windows of lenght k.

    Inputs:

        sequence - a string representing a DNA sequence.
        as_overlap - boolean that represents if overlap is needed.
        k - a integer reppresenting the lengths of overlappig bases.
            Default is 20.

    Outputs:

        gc_content - an array-like object with


    """
    seq = ''
    for file in fasta_dict[name]:
        for n, seq in parse_fasta(file):
            # make sequence upper case and getting the length of it
            seq += seq.upper()
    # the array-like object to collect the data
    gc_content = []
    # non overlap sequence length
    non_overlap = range(0, len(seq) - k + 1, k)
    # overlap sequence length
    overlap = range(0, len(seq) - k + 1)
    # overlap is needed
    if as_overlap:
        # iterates to the overlap region
        for i in overlap:
            # creates the substring to count the gc_content
            subseq = seq[i:i + k]
            # count and sum up the Gs and Cs counts
            g_c = gc_cython(subseq)
            # collect the data in the array container
            gc_content.append((i, round(g_c, 4) * 100))
    # if non overlap is choosed
    else:
        # iterates to the mon overlap region
        for j in non_overlap:
            # creates the substring to count the gc_content
            subseq = seq[j:j + k]
            # count and sum up the Gs and Cs counts
            g_c = gc_cython(subseq)
            # collect the data in the array container
            gc_content.append((j, round(g_c, 4) * 100))
    return gc_content


def parse_arguments():
    """Parse the command line arguments to the the statistics.
    Sets up the argparse command-line parser and calls it. These argsions can be accessed
    using args.argsion.
    The resulting results are csv files from each genome countained in the genus directory
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
    # name of the input diretory, ex. Data/Genomes_splitted
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
    # check if the output directory existe other wise create it
    if os.path.exists(dir_out):
        print(colored('The directory to save the files already exists!',
                      'red',
                      attrs=['bold']))
        pass
    else:
        make_me_a_folder(dir_out)

    # initialize the file counter
    cnt_spc = 0
    # input the file paths and print it to show where the script is doing
    for name in fasta_dict.keys():
        print(colored(f"Start working with genus {name}\n", attrs=['bold']))
        # initialize the gc counts
        gc = gc_count_fasta(fasta_dict, name)
        # make a series
        gc_series = pd.Series(gc, index=["GC"]).reset_index()
        # create a data frame
        df_gc = pd.DataFrame(gc_series).rename(columns={'index': 'bases', 0: 'counts'})
        # count the base composition
        bases = count_bases_fasta(fasta_dict, name)
        # create a dat frame
        df = pd.DataFrame(bases.items(), columns=['bases', 'counts'])
        # concatenate the two data frames
        df_final = pd.concat([df, df_gc])
        # count the gc content in a slide window
        window = genome_stats_in_windows(fasta_dict, name, as_overlap=False, k=3000)
        # create a data frame
        dfw = pd.DataFrame(window, columns=['GC_window', 'GC_content'])
        # saving the final results as a csv file
        full_path = os.path.join(dir_out, name, 'BasicStats')
        file_name = f'{name}_basic_stats.csv'
        if not os.path.exists(full_path):
            os.makedirs(full_path)
        df_final.to_csv(f'{full_path}/{file_name}_{sub_sub_dir}', index=False)
        dfw.to_csv(f'{full_path}/{file_name}_{sub_sub_dir}_gc_window.csv', index=False)
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

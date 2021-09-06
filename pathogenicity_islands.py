#!usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
from time import time
from termcolor import colored
import argparse
import pandas as pd
from system_utils import make_me_a_folder, get_all_fasta
from fasta_parser import parse_fasta
from basic_stats_genomes import gc_cython


def gc_map(genome, k, gc_threshold):
    """
    Get the sequence and the start and end position of
    putative Pathogenic Island from prokaryotic genomes.

    """
    seq_len = len(genome) - (len(genome) % k)
    out_seq = []
    # Determine GC content of each block and change string accordingly
    for i in range(0, seq_len, k):
        gc = gc_cython(genome[i:i + k])
        if gc < gc_threshold:
            out_seq.append((genome[i:i + k], i, i + k, gc))
    return out_seq


def parse_arguments():
    """Parse the command line arguments to the GI putatie regions finder.
    Sets up the argparse command-line parser and calls it. These arguments can be accessed
    using args.option.
    The results are csv files from each genome in the genus directory
    with a list of putatie sequences, start and end positions
     generated from these genomes (chromosomes and plasmids).
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

    parser.add_argument('-pt',
                        '--pattern_file',
                        type=str,
                        dest='pattern_file',
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
    # name of the root directory to save the final result
    dir_out = args.dir_out
    # path and name of the text file with the patterns
    pattern_file = args.pattern_file
    # get the list of all paths to the files in the input directory
    # ex., Data/Genomes_splitted
    all_files = get_files(dir_in)
    # get all patterns
    all_patterns = read_patterns(pattern_file)
    # check if the output directory existe other wise create it
    if os.path.exists(dir_out):
        print(colored('The directory to save the files already exists!',
                      'red',
                      attrs=['bold']))
        pass
    else:
        make_me_a_folder(dir_out)

    # initialize the file counter
    num_files = 0
    # input the file paths and print it to show where the script is doing
    for filen in all_files:
        name = filen.split('/')[2]
        data = filen.split('/')[3]
        print(colored(f"Working with {data} from genus/species {name}", attrs=['bold']))
        # get the search done
        for n, seq in parse_fasta(filen):
            print(f'Start counting the restriction enzymes cut sites in the sequence {n}')
            cut_sites = all_re_cut_sites(seq, all_patterns)
            df = pd.DataFrame(cut_sites, columns=['site', 'positions'])
            full_path = os.path.join(dir_out, name, 'RE_cuts')
            file_name = f'{n}_{data}_re_cuts.csv'
            if not os.path.exists(full_path):
                os.makedirs(full_path)
            print(f'Saving the files in {full_path}\n')
            df.to_csv(f'{full_path}/{file_name}', index=False)
        # the number of files analyzed
        num_files += 1
    # the final time
    end = time()
    # print some info
    print(colored(f"Total number of files analyzed: {num_files}\n.",
                  attrs=['bold']))
    print(colored(f'Total time for the script finishes: {round(end - start, 2)}.',
                  'red',
                  attrs=['bold']))
    print(colored('Done!',
                  'green',
                  attrs=['bold']))


if __name__ == "__main__":
    sys.exit(main())

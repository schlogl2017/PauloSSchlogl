import os
import sys
import glob
import time
import argparse
from termcolor import colored
import numpy as np
import pandas as pd
from system_utils import make_me_a_folder
# usage: scoring_kmers.py [-h] -di path [-do DIR_OUT] [-sd SUB] [-k K]
# scoring_kmers.py: error: the following arguments are required: -di/--dir_in


def get_scores_kmer_data(filenames):
    cols = ['kmer', 'Z_score', 'Frequency', 'Rank_Z_score', 'Score']
    conc_df = pd.DataFrame()
    for filename in filenames:
        name = filename.split('/')[2]
        df_temp = pd.read_csv(filename, usecols=['kmer', 
                                                 'Observed', 
                                                 'Expected',
                                                 'Z_score'])
        df_temp.insert(0, 'Name', name)
        df_temp['Frequency'] = df_temp['Observed'] / df_temp['Observed'].sum()
        df_temp['Rank_Z_score'] = df_temp['Z_score'].rank(ascending=True)
        df_temp['Score'] = df_temp['Observed'] * np.log(df_temp['Observed']/df_temp['Expected'])
        df_temp = df_temp[['Name', 
                           'kmer', 
                           'Z_score', 
                           'Frequency', 
                           'Rank_Z_score', 
                           'Score']]
        conc_df = pd.DataFrame(pd.concat([conc_df, df_temp], copy=False))
        del df_temp
    return conc_df.reset_index(drop=True)                        


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
                        '--sub',
                        type=str,
                        dest='sub',
                        help='directory name for output files.')

    parser.add_argument('-k',
                        '--k_len',
                        type=str,
                        dest='k',
                        help='Name representing the type file. Ex., gz')
    return parser.parse_args()


def main():
    # starting count the staring time of the script
    start = time.time()
    # checking the current directory and printing it
    cwd = os.getcwd()
    print(colored(f'\nThe working directory: {cwd}\n',
                  'green',
                  attrs=['bold']))
    # passing the arguments to the script
    args = parse_arguments()
    # name of the input directory, ex. Data/Genomes_splitted
    dir_in = args.dir_in
    # name of the root directory to save the final result
    dir_out = args.dir_out
    # name of the root directory to save the final result
    sub = args.sub
    # kmer_length
    k = args.k
    # get the fasta files
    filenames = glob.glob(f'{dir_in}/*/Chromosomes/kmers_{k}/*.gz')
    # check if the output directory exist other wise create it
    if os.path.exists(dir_out):
        print(colored('The directory to save the files already exists!',
                      'red',
                      attrs=['bold']))
        pass
    else:
        make_me_a_folder(dir_out)
    print('Start calculating the scores')
    df = get_scores_kmer_data(filenames)
    file_name = f'All_kmer_{k}_{sub}_scores'
    full_path = os.path.join(f'{dir_out}', 'Scores')
    if not os.path.exists(full_path):
        os.makedirs(full_path)
    df.to_csv(f'{full_path}/{file_name}.csv', index=False)
    # the final time
    end = time.time()
    # print some info
    print(colored(f"Total number of genus/species analyzed: {len(filenames)}\n.",
                  attrs=['bold']))
    print(colored(f'Total time for the script finishes: {round(end - start, 2)}.',
                  'red',
                  attrs=['bold']))
    print(colored('Done!',
                  'green',
                  attrs=['bold']))


if __name__ == "__main__":
    sys.exit(main())

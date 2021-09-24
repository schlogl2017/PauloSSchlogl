import os
import sys
import time
import argparse
import csv
from collections import defaultdict


def parse_arguments():
    """Parse the command line arguments to the the script.
    Sets up the argparse command-line parser and calls it. These options can be accessed
    using opt.option. For example opt.alp stores the alphabet provided.
    """
    parser = argparse.ArgumentParser(description='A script to mean jellyfish count kmer in genomes.')
    parser.add_argument('-di',
                        '--dir_in',
                        metavar='path',
                        type=str,
                        required=True,
                        dest='dir_in',
                        help='Path to the files')

    parser.add_argument('-sd',
                        '--sub_dir',
                        type=str,
                        dest='sub_dir',
                        help='Subdirectory name for output files.')
    parser.add_argument('-sd2',
                        '--sub_dir2',
                        type=str,
                        dest='sub_dir2',
                        help='Subdirectory name for output files.')

    parser.add_argument('-do',
                        '--dir_out',
                        type=str,
                        dest='dir_out',
                        help='directory name for output files.')

    parser.add_argument('-k',
                        '--length',
                        type=int,
                        default=2,
                        action='store',
                        dest='k',
                        help='Specify the maximum kmer/palindrome length.')

    return parser.parse_args()


def get_number_genome_files(dir_in, sub_dir):
    """
    Return a dictionary mapping the genus/species that present a number of
    genomes greater than one.
    
    Inputs
        dir_in - string representing a root directory
        sub_dir - string representing a subdirectory
        
    Outputs
        dict - a dictionary-like object mappinf the geneus/species
               keys to their number of genomes in a dictionary.
    """
    num_gen = defaultdict(int)
    names = os.listdir(dir_in)
    for name in names:
        # get the full path to the fasta files
        # for each genus
        pwd = os.path.join(dir_in, name, sub_dir)
        # add the name and the number of fasta files for each genus
        num_gen[name] = num_gen.get(name, 0) + len(os.listdir(pwd))
    # return the dictionary filtered for genus with number genomes
    # in the director bigger tha one
    return {k: cn for (k, cn) in num_gen.items() if cn > 1}


def kmer_count_mean_genomes(dir_in, name, sub_dir, num_gen, k):
    """
    Function to get the mean of kmer counts (jellyFish) from genus/species with more
    than one genome.
    
    Inputs
        dir_in - string representing the directory of the jellyfish files (saved
                 as csv files).
        name - a string representing the genus/species.
        num_gen - a dictionary-like object mapping the genus/species to the number
                  (integer) representing the number os genomes (fasta files) fro that
                  genus.
        k -  a integer representing the length of the kmer.
    
    Outputs
        new - a dictionary-like object mapping the kmer to the calculated mean (integer)
              to be saved as a new csv file.
    """
    new = defaultdict(int)
    # get the path to the files
    path = os.path.join(dir_in, sub_dir)
    # get the full path to the csv files
    ful_path = f'{path}/{name}_k{k}.csv'
    # start reading the csv file
    # using the csv module because it reads
    # one line at time (generator)
    # could use pandas = dont have memory for big files
    with open(ful_path) as csvfile:
        data = csv.reader(csvfile)
        # iterates through the generator line by line
        for row in data:
            # read as list then get each data as needed
            kmer, cnt = row[0], int(row[1])
            # then fill the dictionary with the kmer and the calculated mean
            new[kmer] = new.get(kmer, 0) + (cnt // num_gen[name])
    return new


def main():
    """Parses options from the command line.
    Computes the k-mers count means from jellyfish software where
    bacterial genus/species with number of sequence genomes bigger than one.
    """
    # checking the wotory
    cwd = os.getcwd()
    print(f'The working directory: {cwd}')
    # counting time 
    start_time = time.process_time()
    # passing args
    arg = parse_arguments()
    dir_in = arg.dir_in
    k = arg.k
    sub_dir = arg.sub_dir
    sub_dir2 = arg.sub_dir2
    dir_out = arg.dir_out
    # checking if the output directory exist
    # if not make it
    f_pwd = os.path.join(dir_out, sub_dir2[:10])
    if os.path.exists(f_pwd):
        pass
    else:
        os.makedirs(f_pwd)
    # get the number of fasta genomic files from all genus/spcs
    num_gen = get_number_genome_files('Data/Genomes_splitted', sub_dir)
    # # get the genus names
    genus_names = num_gen.keys()
    print(f'The number of genus/species with number of sequenced genomes bigger than one is {len(genus_names)}')
    # count to check the number of files
    cnt = 0
    # starting calculating the mean
    for name in genus_names:
        km_mean_data = kmer_count_mean_genomes(dir_in, name, sub_dir2, num_gen, k)
        out_dir = f'{f_pwd}/{name}'
        csv_name = f'{name}_chr_km{k}.csv'
        if os.path.exists(out_dir):
            pass
        else:
            os.makedirs(out_dir)
        # saving the data
        with open(f'{out_dir}/{csv_name}', 'w') as fout:
            for km, mc in km_mean_data.items():
                fout.write(f'{km},{mc}\n')
        cnt += 1
    # get final time of the script
    end = time.process_time()
    total_time = end - start_time
    print(f'Where read and manipulated {cnt} files')
    print(f'The script takes {total_time} to finish!')
    # print(f'Where read and manipulated {cnt} files')
    print('Done!')


if __name__ == "__main__":
    sys.exit(main())

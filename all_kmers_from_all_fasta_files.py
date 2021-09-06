#!/usr/bin/env python
# coding: utf-8

import os
import sys
import time
import math
import argparse
from pprint import pprint
import pandas as pd
import numpy as np
from collections import defaultdict, Counter
import itertools
from alphabet import iupac_dna
from fasta_parser import fasta_item_counter, parse_fasta


def get_kmers_counts_dict_list(filenames_lst, alphabet, kmin, kmax):
    dic_list = []
    for file in filenames_lst:
        for n, seq in parse_fasta(file):
            dic_list.append(count_kmers(seq, alphabet, kmin, kmax))
    return dic_list


def merge_counts(list_dict, num_files):
    counts = Counter()
    for dic in list_dict:
        counts.update(dic)
    return {k: (cnt//num_files) for (k, cnt) in counts.items()}


def get_seq_lens_mean(filenames_lst):
    seq_len = 0
    len_files = len(filenames_lst)
    for filename in filenames_lst:
        for _, seq in parse_fasta(filename):
            seq_len += len(seq)
    total = seq_len // len_files
    return total


def kmers_frequencies(kmer_counts):
    freqs = defaultdict(float)
    total = sum(kmer_counts.values())
    for kmer, cnt in kmer_counts.items():
        freqs[kmer] = freqs.get(kmer, 0) + cnt / total
    return freqs


def get_mean_all_kmers_genomes_counts(filenames_lst, alphabet, kmin, kmax):
    """This function take as inputs the list of fasta files, an alphabet,
    the minimum and maximum kmer length. It count and update the kmer counts
    from a list of dictionaries to just one with the mean of all genomes where
    the kmers were counted.
    It returns a list of sorted tuples keys-values."""
    all_kmers = Counter()
    f_len = len(filenames_lst)
    for filename in filenames_lst:
        for name, seq in parse_fasta(filename):
            all_kmers.update(count_kmers(seq, alphabet, kmin, kmax))
    kmer_all_counts = {k: (cnt//f_len) for (k, cnt) in all_kmers.items()}
    return sorted(kmer_all_counts.items(), key=lambda k:k[0])

def count_kmers(sequence, alphabet, min_k, max_k):
    alphabet = set(alphabet)
    counts = defaultdict(int)
    for kmer in get_kmers_from_sequence(sequence, min_k, max_k):
        if set(kmer).issubset(alphabet):
            counts[kmer] = counts.get(kmer, 0) + 1
    return counts

def get_all_possible_kmers(alphabet, min_k, max_k):
    """Returns a list of all possible combinations of k-mers of
    length k from a input alphabet."""
    kmers = [''.join(letters) for n in range(min_k, max_k + 1)
             for letters in itertools.product(alphabet, repeat=n)]
    return kmers

def get_kmers_from_sequence(sequence, min_k, max_k):
    """
    Generate all DNA k-mers over the entirety of a sequence.
    Inputs:
    sequence - string where all kmers will be checked
    min_k: minimum DNA kmer length (int)
    max_k: maximum DNA kmer length (int)
    Output:
    yields all DNA kmers (str) of length min_k to max_k
    """
    limits = range(min_k, max_k + 1)
    for i in range(0, len(sequence) - max_k + 1):
        for j in limits:
            yield sequence[i:i+j]


def get_fasta_files(dir_name):
    infiles = []
    for path, subdirs, files in os.walk(dir_name):
        for name in files:
            input_files = os.path.join(path, name)
            if input_files.endswith('fa.gz') or input_files.endswith('.fa') or input_files.endswith('.fasta') or input_files.endswith('.fa.gz')                     or input_files.endswith('.fna') or input_files.endswith('.fna.gz'):
                infiles.append(input_files)
    return infiles


def get_species_name(filenames):
    species = []
    with open(filenames, 'r') as fh:
        for line in fh:
            sp = line.strip()
            species.append(sp)
    return species


def get_full_name(dir_out, sud_dir, ssub):
    return os.path.join(dir_out, sud_dir, ssub)

def get_expected_values(kmer_list, kmer_counts):
    """Calculates the expected value of all palindrome kmers in a sequence.
    Inputs:
    pal_lst = list of palindromes (str)
    counts = dictionary of kmer counts
    Output:
    expected - dictionary with expected values for all palindromes.
    The expected values are calculated as:
       E(C(W)) = C(W1-Wn-1) * C(W2-Wn) / C(W2-Wn-1)
    """
    expected = defaultdict(float)
    for kmer in kmer_list:
        suf = kmer_counts[kmer[1:]]
        pre = kmer_counts[kmer[:-1]]
        mid = kmer_counts[kmer[1:-1]]
        # to catch the divide by zero error
        if mid == 0:
            expected[kmer] = 0.0
        else:
            ex = (suf * pre) / mid
            expected[kmer] = expected.get(kmer, 0.0) + ex
    return expected


def get_z_scores(kmer_list, kmer_counts, expected_kmers, len_seq):
    """Calculates the z_score of all palindromes.
    Input:
    palindrome_lst = list of palindromes (str)
    counts = dictionary of kmer counts
    expected = dictionary of kmer expected values
    length_sequence = length of sequence (int)
    Output:
    z_score dictionary where key are palindromes and values are the calculated z_score (float)
    The z_scores are calculated as:
        Z(W) = (C(W)) - E(C(W)) / sigma(W)
    And sigma as:
        sigma(W) = sqrt(E(C(W))) * (1 - E(C(W)/N))
    """
    z_score = defaultdict(float)
    for kmer in kmer_list:
        if expected_kmers[kmer] == 0.0:
            z_score[kmer] = 0.0
        else:
            sigma = math.sqrt(expected_kmers[kmer]) * (1 - expected_kmers[kmer] / (2 * len_seq))
            z = (kmer_counts[kmer] - expected_kmers[kmer]) / sigma
            z_score[kmer] = z_score.get(kmer, 0.0) + z
    return z_score


def get_pvalues(kmer_list, z_score_kmers):
    """Calculates the p_value of all palindromes.
    Input:
    palindrome_lst - list of palindromes (str)
    z_score - a dictionary with palindromes z_scores (float)
    Output:
    a palindrome p_value dictionary.
    For probability of getting a z-value larger than t
    P(z > t) = erfc(t/sqrt(2))/2
    For probability of getting a z-value smaller than t
    P(z > t) = erfc(-t/sqrt(2))/2
    """
    p_values = defaultdict(float)
    for kmer in kmer_list:
        if z_score_kmers[kmer] < 0.0:
            under = math.erfc(-z_score_kmers[kmer] / math.sqrt(2)) / 2
            p_values[kmer] = p_values.get(kmer, 0.0) + under
        else:
            over = math.erfc(z_score_kmers[kmer] / math.sqrt(2)) / 2
            p_values[kmer] = p_values.get(kmer, 0.0) + over
    return p_values


def get_evalues(kmer_list, p_value_kmers):
    """Calculates the e_value of all palindrome kmers.
    Inputs:
    palindrome_lst - list of palindromes (str)
    p_value - a dictionary where key are the palindromes (str) and the value are p_value (float)
    Output:
    The e_value as a dictionary where the key are the palindrome (str)
    and value are e_value (float).
    """
    e_value = defaultdict(float)
    num_tests = len(kmer_list)
    for kmer in kmer_list:
        p = p_value_kmers[kmer] * num_tests
        e_value[kmer] = e_value.get(kmer, 0.0) + p
    return e_value


def get_scores(kmer_list, kmer_counts, expected_kmers):
    """Calculates de kmer escore and its expected frequency in the genome.
    Inputs:
    kmer_list - list substring of length k
    kmer_counts - a dictionary with kmer and it counts
    expected_kmers - a dictionary with kmer and it expected counts
    Output:
    scores - a dictionary with kmers and it scores, calculated as:
    S = obs - exp / obs + exp
    """
    scores = defaultdict(float)
    for kmer in kmer_list:
        if expected_kmers[kmer] == 0.0 or kmer_counts[kmer] == 0.0:
            scores[kmer] = 0.0
        else:
            scr = (kmer_counts[kmer] - expected_kmers[kmer]) / (kmer_counts[kmer] + expected_kmers[kmer])
            scores[kmer] = scores.get(kmer, 0.0) + scr
    return scores


def get_new_scores(kmer_list, kmer_counts, expected_kmers):
    """Calculates de kmer escore and its expected frequency in the genome.
    Inputs:
    kmer_list - list substring of length k
    kmer_counts - a dictionary with kmer and it counts
    expected_kmers - a dictionary with kmer and it expected counts
    Output:
    scores - a dictionary with kmers and it scores, calculated as:
    S = obs - exp / obs + exp
    """
    scores = defaultdict(float)
    for kmer in kmer_list:
        scores[kmer] = scores.get(kmer, 0.0)
        if expected_kmers[kmer] == 0.0 or kmer_counts[kmer] == 0.0:
            scores[kmer] = 0.0
        else:
            scr = kmer_counts[kmer] / (kmer_counts[kmer] + expected_kmers[kmer])
            scores[kmer] = scores.get(kmer, 0.0) + scr
    return scores


def get_odds_ratio(kmer_list, kmer_freqs):
    ors = defaultdict(float)
    for kmer in kmer_list:
        midf = kmer_freqs[kmer[1:-1]]
        pref = kmer_freqs[kmer[:-1]]
        suff = kmer_freqs[kmer[1:]]
        kmf = kmer_freqs[kmer]
        # to catch the divide by zero error
        if midf == 0.0 or kmf == 0.0 or pref == 0.0 or suff == 0.0:
            ors[kmer] = ors.get(kmer, 0.0)
        else:
            od = (kmf * midf) / (pref * suff)
            ors[kmer] = ors.get(kmer, 0.0) + od
    return ors


def get_difference(kmer_list, kmer_counts, expected_kmers):
    diff = defaultdict(float)
    for kmer in kmer_list:
        d = kmer_counts[kmer] - expected_kmers[kmer]
        diff[kmer] = diff.get(kmer, 0.0) + d
    return diff


def get_log_odds(kmer_list, kmer_counts, expected_kmers):
    log_ods = defaultdict(float)
    for kmer in kmer_list:
        log_ods[kmer] = log_ods.get(kmer, 0.0)
        if kmer_counts[kmer] == 0.0 or expected_kmers[kmer] == 0.0:
            log_ods[kmer] = 0.0
        else:
            lod = math.log(kmer_counts[kmer] / expected_kmers[kmer])
            log_ods[kmer] += lod
    return log_ods


def get_kmer_statistics(kmer_list,
                        kmer_counts,
                        kmer_expected,
                        kmer_z_scores,
                        kmer_e_values,
                        kmer_odds_ratios,
                        kmer_diffs,
                        kmer_scores,
                        kmer_new_scores,
                        kmer_log_odds):
    """"""
    stats = []
    for kmer in kmer_list:
        obs = kmer_counts[kmer]
        exp = kmer_expected[kmer]
        z_scr = kmer_z_scores[kmer]
        e_val = kmer_e_values[kmer]
        odds = kmer_odds_ratios[kmer]
        diff = kmer_diffs[kmer]
        scr = kmer_scores[kmer]
        nscr = kmer_new_scores[kmer]
        lod = kmer_log_odds[kmer]
        stats.append((kmer,
                      obs,
                      exp,
                      z_scr,
                      e_val,
                      odds,
                      diff,
                      scr,
                      nscr,
                      lod))
    return stats

def get_dataframe_from_kmer_alldata(dir_out, filename, k, kmerdata):
    df = pd.DataFrame(kmerdata, columns=["kmer",
                                         "Observed",
                                         "Expected",
                                         "Z_score",
                                         "Evalues",
                                         "Odds",
                                         "Diff",
                                         "Scores",
                                         "NScores",
                                         "Log_odds"])
    df.to_csv(f"{dir_out}/{filename}_{k}_.csv")


def parse_arguments():
    """Parse the command line arguments to the genome_kmer_analysis script.
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
    parser.add_argument('-ex',
                        '--extension',
                        type=str,
                        dest='extension',
                        help='extension of the files.')
    parser.add_argument('-do',
                        '--dir_out',
                        type=str,
                        dest='dir_out',
                        help='directory name for output files.')
    parser.add_argument('-f',
                        '--file',
                        type=str,
                        dest='filenames',
                        help='species list  as txt.')
    parser.add_argument('-e',
                        '--max_e',
                        type=float,
                        default=0.01,
                        action='store',
                        dest='max_e',
                        help='Threshold e value for report over/under represented palindromes.')
    parser.add_argument('-ma',
                        '--max_kmer_len',
                        type=int,
                        action='store',
                        dest='kmax',
                        help='Specify the maximum kmer/palindrome length.')
    parser.add_argument('-mi',
                        '--min_kmer_len',
                        type=int,
                        action='store',
                        dest='kmin',
                        help='Specify the minimum kmer/palindrome length.')
    parser.add_argument('-al',
                        '--alph',
                        type=str,
                        default=iupac_dna,
                        dest='alphabet',
                        help='The allowed characters tha compound the sequence.')
    return parser.parse_args()


def main():
    """Parses options to the command line.
    This script receive a path to the fasta files, a directory name to save the
    csv count/statistic files and a list of species names.
    The final result is a csv file with the mean of all kmer counts in the
    output directory."""
    
    print('\nStarting to process the script all_kmers_from_all_fasta_files.py\n')
    cwd = os.getcwd()
    print(f'The working directory: {cwd}\n')
    opt = parse_arguments()
    start_time = time.process_time()
    spc_file = opt.filenames
    spc_names = get_species_name(spc_file)
    alphabet = iupac_dna
    path = opt.path
    ssub = opt.sub_dir
    extension = opt.extension
    kmin = opt.kmin
    kmax = opt.kmax
    dir_out = opt.dir_out
    ssubo = f'kmer{kmax}'
    
    cnt = 0
    spc = 0
    full_dir_name = []
    for spc_name in spc_names:
        full_dir_name.append(get_full_name(path, spc_name, ssub))
        spc += 1
    
    for filename in full_dir_name:
        fastas = get_fasta_files(filename)
        num_files = len(fastas)
        km_dic_list = get_kmers_counts_dict_list(fastas, alphabet, (kmin- 2), kmax)
        kmer_counts = merge_counts(km_dic_list, num_files)
        len_seq = get_seq_lens_mean(fastas)
        kmer_freqs = kmers_frequencies(kmer_counts)
        kmer_list = get_all_possible_kmers(alphabet, kmin, kmax)
        kmer_expected = get_expected_values(kmer_list, kmer_counts)
        kmer_zscores = get_z_scores(kmer_list,
                                kmer_counts,
                                kmer_expected,
                                len_seq)
        kmer_pvalues = get_pvalues(kmer_list, kmer_zscores)
        kmer_evalues = get_evalues(kmer_list, kmer_pvalues)
        kmer_scores = get_scores(kmer_list, kmer_counts, kmer_expected)
        kmer_nscores = get_new_scores(kmer_list, kmer_counts, kmer_expected)
        kmer_odds_ratio = get_odds_ratio(kmer_list, kmer_freqs)
        kmer_diff = get_difference(kmer_list, kmer_counts, kmer_expected)
        kmer_lod = get_log_odds(kmer_list, kmer_counts, kmer_expected)
        kmer_data = get_kmer_statistics(kmer_list,
                                    kmer_counts,
                                    kmer_expected,
                                    kmer_zscores,
                                    kmer_evalues,
                                    kmer_odds_ratio,
                                    kmer_diff,
                                    kmer_scores,
                                    kmer_nscores,
                                    kmer_lod)
        df = pd.DataFrame(kmer_data, columns=["kmer", 
                                          "Observed",
                                          "Expected",
                                          "Z_score",
                                          "Evalues",
                                          "Odds",
                                          "Diff",
                                          "Scores",
                                          "NScores",
                                          "Log_odds"])
        cnt += num_files
        name_dir_out = os.path.join(dir_out, spc_name, ssub, ssubo)
        csv_name = f'{spc_name}_{ssubo}_stats.csv'
        if not os.path.exists(name_dir_out):
            os.makedirs(name_dir_out)
        df.to_csv(f"{name_dir_out}/{csv_name}")
    print(f"The final csv file was saved in {name_dir_out}/{csv_name}\n")
    end_time = time.process_time()
    total_time = end_time - start_time
    print(f'Total species: {spc}')
    print(f'Total files: {cnt}')
    print(f'Total time for running the script: {total_time}')
    print('Script finished!\n')


if __name__ == "__main__":
    sys.exit(main())


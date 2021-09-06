#!/usr/bin/env python
# coding: utf-8
# usage: python genome_kmer_analysis_final.py -p Data/bacteria_splitted -sd chromosomes \
# -do Results/test -f test_sps.txt -ma 4 -mi 4 -ex gz
import os
import sys
import glob
import argparse
import time
import math
import pandas as pd
import numpy as np
from collections import defaultdict, Counter
import itertools
from alphabet import iupac_dna
from fasta_parser import fasta_item_counter, parse_fasta


def get_species_name(filenames):
    print("Getting species names\n")
    species = []
    with open(filenames, 'r') as fh:
        for line in fh:
            sp = line.strip()
            species.append(sp)
    return species


def get_seq_lens_mean(spc_name, sequences_dict):
    seq_list = sequences_dict[spc_name]
    len_files = len(sequences_dict[spc_name])
    lens = sum([len(seq) for seq in seq_list])/len_files
    return lens

    
def get_full_name(dir_one, sub_dir, susub, extension):
    full = os.path.join(dir_one, sub_dir, susub + f"/*.{extension}")
    return full
    
    
def get_list_paths(spc_names, path, ssub, extension):
    pwd = defaultdict()
    for name in spc_names:
        pwd[name] = get_full_name(path, name, ssub, extension)
    return pwd
    
    
def get_fasta_files_paths(dict_paths):
    print('Getting the fasta files pwd\n')
    dict_fasta_paths = defaultdict(list)
    for name, path in dict_paths.items():
        dict_fasta_paths[name] = dict_fasta_paths.get(name, []) + glob.glob(path)
    return dict_fasta_paths    
    
    
def get_sequence_dict(dict_paths):
    print('Making the fasta file dictionary\n')
    all_fasta = defaultdict(list)
    for name, filenames in dict_paths.items():
        all_fasta[name] = all_fasta.get(name, [])
        for filename in dict_paths[name]:
            for n, seq in parse_fasta(filename):
                all_fasta[name].append(seq)
    return all_fasta    
    

def count_ambiguous_bases(sequence):
    sequence = sequence.upper()
    amb = ['N', 'R', 'Y', 'W', 'S', 'K', 'M']
    return sum({base: sequence.count(base) for base in amb}.values())


def list_kmer_counts(kmer_counts, kmer_list):
    kmer_count = defaultdict(int)
    for kmr in kmer_list:
        cnt = kmer_counts[kmr]
        kmer_count[kmr] = kmer_count.get(kmr, 0) + cnt
    return kmer_count


def all_kmer_counts(alphabet, min_k, max_k, *args):
    """ 
    Counts number of k-mers (kmin <= k <= kmax) in a sequence.
    Inputs:
    sequence = sequence string defined by an alphabet.
    min_k - minimum kmer length (int)
    max_k - maximum kmer length (int)
    Output:
    A dictionary with keys of k-mers and values as counts (int) and
    the length of the sequence data collected as a tuple.
    To the script work this function must have count kmers from
    mink-2 to max_k.
    """
    print('Counting all the kmers\n')
    counts = Counter()
    num_seqs = 0
    for seq in args:
        seq_len = len(seq) - count_ambiguous_bases(seq)
        num_seqs += 1
        for k in range(min_k, max_k + 1):
            for i in range(0, seq_len - k + 1):
                kmer = seq[i:i + k]
                if all(base in set(alphabet) for base in kmer):
                    counts[kmer] += 1
    return {kmr: (cnt//num_seqs) for kmr, cnt in counts.items()}    
    
    
def kmers_frequencies(kmer_counts):
    freqs = defaultdict(float)
    total = sum(kmer_counts.values())
    for kmer, cnt in kmer_counts.items():
        freqs[kmer] = freqs.get(kmer, 0) + cnt / total
    return freqs


def all_kmers_frequencies(kmer_counts):
    #dict(ChainMap(*a))
    len_kmers = set([len(keys) for keys in kmer_counts.keys()])
    freqs = []
    all_freqs = defaultdict(float)
    for lk in len_kmers:
        temp = {kmr:cnt for (kmr, cnt) in kmer_counts.items() if len(kmr) == lk}
        freqs.append(temp)
    frq = map(kmers_frequencies, freqs) # mapping method
    for d in frq:
        all_freqs.update(d)
    return all_freqs
    
    
def get_all_possible_kmers(alphabet, min_k, max_k):
    """Returns a list of all possible combinations of k-mers of
    length k from a input alphabet."""
    kmers = [''.join(letters) for n in range(min_k, max_k + 1)
             for letters in itertools.product(alphabet, repeat=n)]
    return kmers    


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
    print('Start scoring the kmers\n')
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


def get_difference(kmer_list, kmer_counts, expected_kmers):
    diff = defaultdict(float)
    for kmer in kmer_list:
        d = kmer_counts[kmer] - expected_kmers[kmer]
        diff[kmer] = diff.get(kmer, 0.0) + d
    return diff


def get_kmer_statistics(kmer_list,
                        kmer_counts,
                        kmer_expected,
                        kmer_z_scores,
                        kmer_e_values,
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
        # odds = kmer_odds_ratios[kmer]
        # diff = kmer_diffs[kmer]
        scr = kmer_scores[kmer]
        nscr = kmer_new_scores[kmer]
        lod = kmer_log_odds[kmer]
        stats.append((kmer,
                      obs,
                      exp,
                      z_scr,
                      e_val,
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
    parser.add_argument('-do',
                        '--dir_out',
                        type=str,
                        dest='dir_out',
                        help='directory name for output files.')
    parser.add_argument('-f',
                        '--file',
                        type=str,
                        dest='spc_list',
                        help='species list  as txt.')
    parser.add_argument('-ex',
                        '--extension',
                        type=str,
                        dest='extension',
                        help='extension of the fasta files.')
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
    print('\nStarting to process the script genome_kmers_analysis.py\n')
    cwd = os.getcwd()
    print(f'The working directory: {cwd}\n')
    opt = parse_arguments()
    start_time = time.process_time()
    spc_list = opt.spc_list
    spc_names = get_species_name(spc_list)
    print(f'Length of species names list: {len(spc_names)}')
    alphabet = iupac_dna
    path = opt.path
    ssub = opt.sub_dir
    kmin = opt.kmin
    kmax = opt.kmax
    dir_out = opt.dir_out
    ssubo = f'kmer{kmax}'
    extension = opt.extension

    kmer_list = get_all_possible_kmers(alphabet, kmin, kmax)

    spc = 0
    cnt = 0
    for spc_name in spc_names:
        dict_paths = get_list_paths(spc_names, path, ssub, extension)
        fasta_paths = get_fasta_files_paths(dict_paths)
        # sequences_dict = get_sequence_dict(fasta_paths)

        len_seq = get_seq_lens_mean(spc_name, sequences_dict)
        sequences = sequences_dict[spc_name]
        num_files = len(sequences)
        all_kmers_counts = all_kmer_counts(alphabet, kmin-2, kmax, *sequences)
        counts = list_kmer_counts(all_kmers_counts, kmer_list)
        # freq_total = all_kmers_frequencies(all_kmers_counts)
        expected = get_expected_values(kmer_list, all_kmers_counts)
        Z = get_z_scores(kmer_list, counts, expected, len_seq)
        P_val = get_pvalues(kmer_list, Z)
        E_val = get_evalues(kmer_list, P_val)
        SCR = get_scores(kmer_list, counts, expected)
        NSCR = get_new_scores(kmer_list, counts, expected)
        # ODR = get_odds_ratio(kmer_list, freq_total)
        DIFF = get_difference(kmer_list, counts, expected)
        LOD = get_log_odds(kmer_list, counts, expected)
        data = get_kmer_statistics(kmer_list,
                                   counts,
                                   expected,
                                   Z,
                                   E_val,
                                   ODR,
                                   DIFF,
                                   SCR,
                                   NSCR,
                                   LOD)
        name_dir_out = os.path.join(dir_out, spc_name, ssub, ssubo)
        csv_name = f'{spc_name}_{ssubo}_stats.csv'
        if not os.path.exists(name_dir_out):
            os.makedirs(name_dir_out)
        get_dataframe_from_kmer_alldata(name_dir_out, csv_name, kmax, data)
        print(f"The final csv file was saved in {name_dir_out}/{csv_name}\n")
        spc += 1
        cnt += num_files
    end_time = time.process_time()
    total_time = end_time - start_time
    print(f'Total species: {spc}')
    print(f'Total files: {cnt}')
    print(f'Total time for running the script: {total_time}')
    print('Script finished!\n')


if __name__ == "__main__":
    sys.exit(main())

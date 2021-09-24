#!usr/bin/env python

import math
from collections import defaultdict
import pandas as pd
from scipy.stats import norm



def expected_zero(kmer, len_seq, nuc_freqs):
    """
    Calculates the expected number of a substring of length k (k-mer)
    based in a zero order Markov model (zom).
    
    Inputs:
    
        kmer - a substring representing the word/k-mer (a string of length k).
        base_freqs - a dictionary-like object mapping the frequency of the
                     nucleotides/basesand their counted values normalized by
                     the length of the sequence/genome.
        len_seq - a integer representing the length of the sequence/genome where
                  the kmers were counted.
    
    Outputs:
        expected - a integer representing the kmer/substring of length k in the
                sequence/genome of length N (len_seq - len_kmer + 1).
    
    
    """
    # number positions where the kmer could be counted
    N = len_seq - len(kmer) + 1
    # get the base counts in the kmer sequence
    kbc = base_stats(kmer, 'ACGT', True, True)
    # the counter 
    expec = 1
    # iterates through the base kmer list
    # and the dictionary of nucleotides frequencies
    for bk, bs in zip(kbc, nuc_freqs):
        # multiply the nucleotides frequencies 
        # in the power of the base counts in the kmer
        expec *= nuc_freqs[bs]**kbc[bk]
    # returns the expected values fro the current kmer
    # in the genome
    return int(expec * N)


def prob_dinuc(dinuc_cnts, monuc_cnts, nuc, dinuc):
    """
    Calculates the probability of a given mononucleotide is
    if followed by mononucleotide - 1. Ex., a 'T' follow by
    an 'A'.
    
    Inputs:
    
        dinuc_cnts -  a dictionary-like object mapping the
                      dinucleotides to their counts in a
                      determineted sequence.
        monuc_cnts - a dictionary-like object mapping the
                      mononucleotides to their counts in a
                      determineted sequence.
                      
    Outputs:
        probability - a float representing the probability of
                      a given nucleotide be followed by the
                      nucleotide - 1 from a kmer of length 2.
    """
    return dinuc_cnts[dinuc] / monuc_cnts[nuc]
    

def expected_kmer_by_zom(kmer, base_freqs, len_seq):
    """
    Calculates the expected number of a substring of length k (k-mer)
    based in a zero order Markov model (zom).
    
    Inputs:
    
        kmer - a substring representing the word/k-mer (a string of length k).
        base_freqs - a dictionary-like object mapping the frequency of the
                     nucleotides/basesand their counted values normalized by
                     the length of the sequence/genome.
        len_seq - a integer representing the length of the sequence/genome where
                  the kmers were counted.
    
    Outputs:
    
        expected - a integer representing the kmer/substring of length k in the
                sequence/genome of length N (len_seq - len_kmer + 1).
    N = len(seq) - len(kmer) + 1
    E(w) = N*(nuc1*nuc2*nuc3)
    """
    # make a list of letter/kmer
    k_l = list(kmer)
    # number position in the genome where the
    # kmer was counted
    n = len_seq - len(kmer) + 1
    # counter to receive the bases values
    # to multiplicate
    cnt = 1
    # iterates in the kmer list
    # and gets the bases
    # to recover the frequencies from the dcitionary
    for base in k_l:
        # multiply all the bases values that
        # are found in the kmer
        cnt *= base_freqs[base]
    # return the expected value of the kmer
    return int(cnt * n)
        

def get_expected_higher_markov(kmer_list, kmer_counts):
    """
    Calculates the expected value for a list of kmers and their counts.
    
    Inputs:
        kmer_list - list-like object representing all kmer counted in a sequence.
                    The kmers have a length k.
        kmer_counts - dictionary-like object mapping kmer to their counts. The
                      kmer lengths must be between kmin (kmx-2) and kmax.
    
    Output:
    
        expected - dictionary-like object mapping kmer of length k to their 
                   calculated expected values.
    
    The expected values are calculated as:
    'Expected = kmer[:-1] * kmer[1:] / kmer[1:-1]'
    """
    # initialize the container
    expected = defaultdict(int)
    # iterates through the list of kmers
    for kmer in kmer_list:
        # gets the suffix, prefix and middle kmers
        suf, pref, mid = kmer_counts[kmer[1:]], kmer_counts[kmer[:-1]], kmer_counts[kmer[1:-1]]
        # deal with sero division errors
        if mid == 0:
            expected[kmer] = expected.get(kmer, 0)
        else:
            # add the kmer and it expected values 
            expected[kmer] = expected.get(kmer, 0) + int((pref * suf) / mid)
    return expected


def get_variance(kmer_list, len_seq, kmer_expected):
    """
    Calculates the variance from a list of strings of length k.
    
    Inputs:
        kmer_list - list-like object representing all kmer counted in a sequence.
                    The kmers have a length k.
        kmer_expectd - a dictionary-like object mapping kmer to their calculated
                       expected values.
        len_seq - integer representing the length of the sequence where kmers were
                  counted.
    
    Outputs:
    
        variance - a dictionary-like object mapping kmer to their calculated
                   expectd variance.
                   
    Because the model for the count is the sum of N almost independent observations, 
    each with probability P(W), it can be well modeled as a binomial distribution, 
    with varianceThe variance is calculated as:
    E(C(W)) * (1 - E(C(W))/N)    
    """
    # gets the kmer length from the list 
    # of kmers
    k = len(kmer_list[0])
    # get all possible positions to count the kmers
    N = len_seq - k + 1
    # initialize the container
    variance = defaultdict(float)
    # iterates through the kmer list
    for kmer in kmer_list:
        # gets the expected value from the dictionary
        ex_val = kmer_expected[kmer]
        # deals with zero error division
        if ex_val == 0:
            # gets the zero value if expected is zero
            variance[kmer] = variance.get(kmer, 0.0)
        else:
            # calculates the variance
            var = ex_val * (1 - ex_val / N)
            # add the kmer and the variance values in the container
            variance[kmer] = variance.get(kmer, 0.0) + var
    return variance
    
    
def get_standard_deviation(variance):
    """
    Calaculates the standard deviation from the kmers expected values.
    
    Inputs:
        variance - a dictionary-like object mapping kmer to their calculated
                   expectd variance.
    
    Outputs:
        
        std - a dictionary-like object mapping kmer to their calculated
                   expectd std.
    
    The variance is calculated as:
    sigma(W) = sqrt(Expected) * (1 - Expected/len(seq) -k + 1))
    """
    # initialize the container
    std = defaultdict(float)
    # iterates through the kmer list
    for kmer in variance:
        # deal with zero divisions errors
        if variance[kmer] == 0.0:
            std[kmer] = std.get(kmer, 0.0)
        else:
            # add the kmer and the calculated std in the
            # container
            sd = math.sqrt(variance[kmer])
            std[kmer] = std.get(kmer, 0.0) + sd
    return std    
    
    
def z_scores(kmer_exp, kmer_counts, std):
    """
    Calculates the z scores to under/over represented kmers from a sequence.
    The score is calculaated as:
    
    Z(W) = (C(W) – E(C(W))) / sigma(W), where 
    C(w) - observed values
    E(C(w)) - represents the expected value from a kmer
    sigma - represents the standard deviation
    
    Inputs:
        kmer_exp - dictionary-like object mapping kmer of length k to their 
                   calculated expected values.
        kmer_counts - dictionary-like object mapping kmer to their counts. The
                      kmer lengths must be between kmin (kmx-2) and kmax.
        std - a dictionary-like object mapping kmer to their calculated
                   expectd std.
    
    Outputs:
        z_scores - dictionary-like object mapping kmer to their z_scores.
    """
    # initialize the container
    z_scores = defaultdict(float)
    # iterates through the kmer keys
    for kmer in kmer_exp:
        # gets the kmer std value
        sd = std[kmer]
        # deals with zero error division
        if sd == 0.0:
            z_scores[kmer] = z_scores.get(kmer, 0.0)
        else:
            # calculates the z score and add 
            # the kmer and the z score values to the container
            z = (kmer_counts[kmer] - kmer_exp[kmer]) / sd
            z_scores[kmer] = z
    return z_scores   
    
    
def get_p_values(z_scores_kmers):
    """
    Calculates the p value for all kmers.
    The calculation is done as:
    over represented: P(z > t) = erfc(t/sqrt(2))/2
    under represented: P(z > t) = erfc(-t/sqrt(2))/2
    t: thresholder
    
    Inputs:
        z_scores_kmers - dictionary-like object mapping kmer to their z_scores.
        
    Outputs:
        p_vals - dictionary-like object mapping kmer to their p values.
    """
    # initialize the container
    p_vals = defaultdict(float)
    # iterates through the kmer keys
    for kmer in z_scores_kmers:
        # calculates the p values to under represented
        # kmers (negative z scores)
        # add the kmer and p values to the container
        if z_scores_kmers[kmer] < 0.0:
            under = math.erfc(-z_scores_kmers[kmer] / math.sqrt(2)) / 2
            p_vals[kmer] = p_vals.get(kmer, 0.0) + under
        else:
            # add the kmer and p values to the container to over represented
            # and all other kmers
            others = math.erfc(z_scores_kmers[kmer] / math.sqrt(2)) / 2
            p_vals[kmer] = p_vals.get(kmer, 0.0) + others
    return p_vals    


def get_scipy_p_values(z_scores_kmers):
    """
    Calculates the p value for all kmers.
    The calculation is done as:
    over represented: P(z > t) = erfc(t/sqrt(2))/2
    under represented: P(z > t) = erfc(-t/sqrt(2))/2
    t: thresholder
    
    Inputs:
        z_scores_kmers - dictionary-like object mapping kmer to their z_scores.
        
    Outputs:
        p_vals - dictionary-like object mapping kmer to their p values.
    """
    # initialize the container
    p_vals = defaultdict(float)
    # iterates through the kmer keys
    for kmer in z_scores_kmers:
        # calculates the p values to under represented
        # kmers (negative z scores)
        # add the kmer and p values to the container
        if z_scores_kmers[kmer] < 0.0:
            p_vals[kmer] = p_vals.get(kmer, 0.0) + norm.sf(abs(-z_scores_kmers[kmer]))
        else:
            # add the kmer and p values to the container to over represented
            # and all other kmers
            p_vals[kmer] = p_vals.get(kmer, 0.0) + norm.sf(abs(z_scores_kmers[kmer]))
    return p_vals
    
    
def get_e_values(kmer_list, p_vals):
    """    
    Calculates the variance from a list of strings of length k.
    
    Inputs:
        kmer_list - list-like object representing all kmer counted in a sequence.
                    The kmers have a length k.
        p_vals - a dictionary-like object mapping kmer to their calculated
                       p values.
    
    Outputs:
    
        e_values - a dictionary-like object mapping kmer to their calculated
                   e values.  
    """
    # number of tested hypoteses
    hyp_num = len(kmer_list)
    # initialize the container
    e_values = defaultdict(float)
    # iterates through the kmer list
    for kmer in kmer_list:
        # gets the p values from the input container
        p = hyp_num * p_vals[kmer]
        # calculates the e values and add the kmer 
        # and the e values to the container
        e_values[kmer] = e_values.get(kmer, 0.0) + p
    return e_values    
    
    
def gets_selected_kmers(kmer_list, 
                        kmer_counts, 
                        expected_kmers, 
                        z_scores_kmers, 
                        kmer_p_vals, 
                        kmer_e_vals, 
                        eval_cutoff=0.05):
    """
    Compile all the results from kmer analyses in a final csv for further analysis.
    
    Inputs:
        kmer_list - list-like object representing all kmer counted in a sequence.
                    The kmers have a length k.
        kmer_counts - dictionary-like object mapping kmer to their counts. The
                      kmer lengths must be between kmin (kmx-2) and kmax. 
        kmer_exp - dictionary-like object mapping kmer of length k to their 
                   calculated expected values.                   
        z_scores_kmers - dictionary-like object mapping kmer to their z_scores.
        p_vals - a dictionary-like object mapping kmer to their calculated
                       p values.
        e_values - a dictionary-like object mapping kmer to their calculated
                   e values.      
    Outputs:
    
        csv - a comma separated values with kmers and all calculated data.    
    """
    data = []
    to_check = []
    for kmer in kmer_list:
        keval = kmer_e_vals[kmer]
        if keval <= eval_cutoff:
            data.append((kmer, 
                         kmer_counts[kmer], 
                         expected_kmers[kmer], 
                         z_scores_kmers[kmer],
                         kmer_p_vals[kmer],
                         kmer_e_vals[kmer]))
        else:
            to_check.append((kmer, 
                             kmer_counts[kmer],
                             expected_kmers[kmer], 
                             z_scores_kmers[kmer],
                             kmer_p_vals[kmer],
                             kmer_e_vals[kmer]))
    df_filtered = pd.DataFrame(data, columns=['kmer', 
                                              'count',
                                              'expected',
                                              'z_score',
                                              'e_value',
                                              'p_value'
                                              ]).sort_values(by='z_score').reset_index(drop=True)
    df_tocheck = pd.DataFrame(to_check, columns=['kmer', 
                                                 'count',
                                                 'expected',
                                                 'z_score',
                                                 'e_value',
                                                 'p_value'
                                                 ]).sort_values(by='z_score').reset_index(drop=True)
    return df_filtered, df_tocheck    
    
    
def transition_matrix(base_freq, mer_freq):
    """
    Returns a transition matrix from a base_count dictionary and a
    k-mer frequency dictionary as inputs.
    
    Inputs:
        
        base_freq - a dictionary-like object mapping nucleotides to their counts.
        mer_freq - a dictionary-like object mapping kmers to their frequencies.
                   The kmers are substring with length k from a string/genome.
    
    Outputs:
        tm - a dictionary-like object mapping a substring (key) to a dictionary
             mapping th nucleotides to their probabilities. The matrix associated 
             with a change of basis for a vector space. Stochastic matrix, a 
             square matrix used to describe the transitions of a Markov chain. 
             State-transition matrix, a matrix whose product with the state vector 
             at an initial time gives at a later time .
    """
    # initialize the container
    tm = defaultdict(dict)
    # iterates to the key, value pairs
    for char1, cnt1 in mer_freq.items():
        # add the character to the container
        tm[char1] = tm.get(char1, {})
        # iterates to the key, value pairs
        for char2, cnt2 in base_freq.items():
            # add the keys and the calculated probabilities
            # to the matrix
            tm[char1][char2] = tm[char1].get(char2, cnt1 * cnt2)
    return tm    
    
    
def compute_kmer_stats(kmer_list, counts, len_genome, max_e):
    """
    This function computes the z_score to find under/over represented kmers 
    according to a cut off e-value.
    Inputs:
        kmer_list - a list of kmers
        counts - a dictionary-type with k-mers as keys and counts as values.
        len_genome - the total length of the sequence(s).
        max_e - cut off e-values to report under/over represented kmers.
    Outputs:
        results - a list of lists as [k-mer, observed count, expected count, z-score, e-value]
    """
    print(colored('Starting to compute the kmer statistics...\n',
                  'red',
                  attrs=['bold']))
    results = []
    # number of tests, used to convert p-value to e-value.
    n = len(list(kmer_list))
    for kmer in kmer_list:
        k = len(kmer)
        prefix, sufix, center = counts[kmer[:-1]], counts[kmer[1:]], counts[kmer[1:-1]]
        # avoid zero division error
        if center == 0:
            expected = 0
        else:
            expected = (prefix * sufix) // center
            observed = counts[kmer]
            sigma = math.sqrt(expected * (1 - expected / (len_genome - k + 1)))
            # avoid zero division error
            if sigma == 0.0:
                z_score = 0.0
            else:
                z_score = ((observed - expected) / sigma)
                # pvalue for all kmers/palindromes under represented
                p_value_under = (math.erfc(-z_score / math.sqrt(2)) / 2)
                # pvalue for all kmers/palindromes over represented
                p_value_over = (math.erfc(z_score / math.sqrt(2)) / 2)
                # evalue for all kmers/palindromes under represented
                e_value_under = (n * p_value_under)
                # evalue for all kmers/palindromes over represented
                e_value_over = (n * p_value_over)
                if e_value_under <= max_e:
                    results.append([kmer, observed, expected, z_score, p_value_under, e_value_under])
                elif e_value_over <= max_e:
                    results.append([kmer, observed, expected, z_score, p_value_over, e_value_over])
    return results    
    
    
    
def markov_sec_order(kmer_list, km_cnts):
    """Pode ser usada para words of 6k"""
    exp_sec_order = defaultdict(int)
    for kmer in kmer_list:
        k1, k2, k3, k4 = km_cnts[kmer[:3]], km_cnts[kmer[1:4]], km_cnts[kmer[2:5]], km_cnts[kmer[3:]]
        mid = kmer[1:-1]
        m1, m2, m3 = km_cnts[mid[:2]], km_cnts[mid[1:3]], km_cnts[mid[2:]]
        p = (k1*k2*k3*k4) / m1*m2*m3
        exp_sec_order[kmer] = exp_sec_order.get(kmer, 0) + p
    return exp_sec_order
    
    
def get_kmer_data(kmer_list,
                  kmer_counts,
                  kmer_expected,
                  kmer_z_scores,
                  kmer_e_values,
                  kmer_p_values):
    """
    Function to put together all the k-mers statistics.
    
    Inputs:
        kmer_list - a list-like objectwith all possible k-mers substrings with
                    length between kmin and kmax.                    
        kmer_counts - a dictionary-like object mapping the k-mers substrings to
                      to their counts.
        kmer_expected - a dictionary-like object mapping the k-mers substrings to
                        to their expected values according a High Order Markov Model.
        kmer_z_scores - a dictionary-like object mapping the k-mers substrings to
                        their calculated z-scores.                    
        kmer_e_values - a dictionary-like object mapping the k-mers substrings to
                        their calculated e-values.                    
        kmer_p_values - a dictionary-like object mapping the k-mers substrings to
                        their calculated p-values.
    Outputs:
        stats - a csv file that is saved in the hard disk to further analyses.
    """
    stats = []
    for kmer in kmer_list:
        obs = kmer_counts[kmer]
        exp = kmer_expected[kmer]
        z_scr = kmer_z_scores[kmer]
        e_val = kmer_e_values[kmer]
        p_val = kmer_p_values[kmer]
        stats.append((kmer,
                      obs,
                      exp,
                      z_scr,
                      e_val,
                      p_val))
    df = pd.DataFrame(stats, columns=["kmer",
                                      "Observed",
                                      "Expected",
                                      "Z_score",
                                      "E_values",
                                      "P_values"])
    
    return df    

#!usr/bin/env python


import multiprocessing as mp
from os import cpu_count
import toolz as tz
from collections import defaultdict
from scipy.stats import binom_test


def char_count(string):
    """
    Function to calculate the count of each
    character that is found in a string.
    
    Inputs:
    
        string - a string representing a word.
    
    Outputs:
    
        char_count - a dictionary-like object representing 
                     the mapping of each character in the 
                     string with their counted values
                     in the input string.
    """
    # initialize the countainer
    char_count = defaultdict(int)
    # iterates through the string
    for char in string:
        # adds the counts of each string 
        # characters to the countainer
        char_count[char] = char_count.get(char, 0) + string.count(char)
    return char_count


def get_string_frequencis(char_count):
    """
    Inputs:
    
        char_count - a dictionary-like object with the kmer (length 
                     between kmin and kmax) and their charcater
                     counts.
    """
    total = sum(char_count.values()
    freqs = defaultdict(float)
    for kmer in char_count:
        freqs[kmer] = freqs.get(kmer, 0.0) + (char_count[kmer]/total)
    return freqs


def string_probability(string, char_freqs):
    """
    Functio to calculate the Bernoulli probability of a
    random string.
    
    Inputs:
        
        string -  a string object representing a word.
        char_freqs - a dictionary-like object representing the
                     mapping of each character in the string with
                     their frequency in a text or background.
    
    Outputs:
    
        prob - a float representing the probability of the input
               string. (IID - Independent and identically distributed
               of the charcaters).
    """
    # initialyze the counts
    prob = 1
    # iterates through the string charcaters
    for i in range(len(string)):
        # multiply the initial prob with the frequency
        # of the chars in the bacground char_freqs
        prob *= char_freqs[string[i]]
    return prob


def test_hypoteshis(string_counts, prob):
    """
    Function that calculates the p-value of seeing a given string 
    in a bag of strings and checks if it is less or higher than
    the probability of a random word.
    
    Inputs:
    
        string_counts - a dictionary mapping words with 
                        their frequencies in a text or
                        background.
        prob -  a float representing the frequency of the string in a
                bacground compostion as a Bernoulli/null hypothesis.
    
    binom_test is utilized as part of the library scipy.stats.
    """
    # initialize the counter
    results = 0
    # iterates through all strings
    for kmer in string_counts:
        # add the value of the calculated p-value
        # and return it
        results += binom_test(x=string_counts[kmer],
                             n=len(string_counts),
                             p=prob,
                             alternative='two-sided')
                 
    return results


def kmer_list_probabilities(kmer_list, char_freqs):
    """
    Functio to calculate the Bernoulli probability of a
    list of strings of length k (k-mers/n-grams).
    
    Inputs:
        
        kmer_list - a list/array-like object of strings representing words
                    of length k.
        char_freqs - a dictionary-like object representing the
                     mapping of each character in the string with
                     their frequency in a text or background.
    
    Outputs:
    
        prob - a dictionary-like object representing the probability of the input
               list of string. 
    
    These probabilities are calculated by multiplying the frequency of character
    that is present in the kmer/ngram:
    'p[AATTC]' = f['A']*f['A']*f['T']*f['T']*f['C']
    """
    kmer_probs = defaultdict(float)
    for kmer in kmer_list:
        p = string_probability(kmer, char_freqs)
        kmer_probs[kmer] = kmer_probs.get(kmer, 0.0) + p
    return kmer_probs
    

def worker(arguments):
    """
    pool = mp.Pool()
    results = pool.map(worker, [('infile', n) for n in range(1, cpu_count()+1)])
    """
    n = os.cpu_count() + 1
    infile, offset = arguments
    with open(infile) as f:
        cg = 0
        totlen = 0
        count = 1
        for line in f:
            if (count % n) - offset == 0:
                if not line.startswith('>'):
                    cg += line.count('C') +
                          line.count('G')
                    totlen += len(line)
            count += 1
    return (cg, totlen)













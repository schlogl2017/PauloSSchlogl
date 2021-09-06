#!usr/bin/env python

from scipy.stats import binom_test

def binomial_kmer_probability(kmer_list, kmer_counts):
    k = len(kmer_list[0])
    prob_null = 0.25**k
    len_kmers = len(kmer_list)
    bi_probs = defaultdict(float)
    for kmer in kmer_list:
        pr = binom_test(x=kmer_counts[kmer], n=len_kmers, p=prob_null, alternative='two-sided')
        bi_probs[kmer] = bi_probs.get(kmer, 0.0) + pr
    return bi_probs

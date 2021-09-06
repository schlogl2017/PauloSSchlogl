#!usr/bin/env python
import re
from math import log
from time import process_time
from functools import reduce
from operator import mul
from pprint import pprint
from random import randint, uniform, choices
import copy
import math
from collections import Counter, defaultdict



def pattern_count(sequence, pattern):
    '''Return the count the input pattern found in to a give string.'''
    return len(re.findall(r'(?='+pattern+')', sequence))


def most_frequent_patterns(sequence, k, n=1):
    """Returns the n most frequent patterns of lenght k from a input
    sequence."""
    return Counter([sequence[i:i+k] for i in range(len(sequence) - k + 1)]).most_common(n)


def frequent_patterns(sequence, k):
    """Returns the most frequent k-mers in a string"""
    freq_pats = set()
    freq_array = compute_frequency(sequence, k)
    max_count = max(freq_array)
    for i in range(4**k -1):
        if freq_array[i] == max_count:
            pattern = number_to_pattern(i, k)
            freq_pats.add(pattern)
    return freq_pats


def reverse_complement(sequence):
    reverse = str.maketrans("ACGT","TGCA")
    return sequence.translate(reverse)[::-1]


def pattern_locations(sequence, pattern):
    """Return the pattern start positions"""
    n = len(sequence)
    m = len(pattern)
    return [i for i in range(n - m + 1) if sequence[i:i+m] == pattern]


def pattern_to_numbers(pattern):
    """Returns a integer number representing the inputed 
    pattern"""
    k = len(pattern)
    base_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    number= 0
    for i in range(k):
        number += base_index[pattern[i]]*4**(k-i-1)
    return number


def number_to_pattern(number, k):
    """Returns a pattern representing the inputed 
    number"""
    nts = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    if k == 1:
        return nts[number]
    prefix_idx = number // 4
    # remainder represents the final nucleotide of the pattern
    reminder = number % 4
    prefix_pattern = number_to_pattern(prefix_idx, k - 1)
    return prefix_pattern + nts[reminder]


def compute_frequency(sequence, k):
    """Returns the number of times that each k-mer Pattern has already 
    appeared in the sequence."""
    freq_array = [0] * (4**k)
    for i, val in enumerate(sequence[:len(sequence) - (k - 1)]):
        freq_array[pattern_to_numbers(sequence[i:i + k])] += 1
    return [f for f in freq_array]


def clump_finder(sequence, k, window_size, num_timest):
    """Returns the number of the times the clumps of patterns of k 
    length appears in a sliding window inside a inputed sequence."""
    len_seq = len(sequence)
    freq_patterns = []
    clumps = [0 for i in range(4**k)]
    window = sequence[0:window_size]
    freq = compute_frequency(window, k)
    for i in range(4 ** k):
        if freq[i] >= num_timestt:
            clumps[i] = 1
    for i in range(1, len_seq - window_size):
        first_pat = sequence[i-1:i-1+k]
        index = pattern_to_numbers(first_pat)
        freq[index] = freq[index] - 1
        last_pat = sequence[i+window_siz-k:i+window_size]
        index = pattern_to_numbers(last_pat)
        freq[index] = freq[index] + 1
        if freq[index] >= num_timest:
            clumps[index] = 1
    for i in range(4**k):
        if clumps[i] == 1:
            pat = number_to_pattern(i, k)
            freq_patterns.append(pat)
    return freq_patterns


def get_base_compostion_stats(sequence, start):
    """Returns basics statistics from a sequence in itś forward and
    revese strands. Autor: Leonard McMillan"""
    half_genome = len(sequence)//2
    ter_C = start + half_genome
    # handle genome's circular nature
    if (ter_C > len(sequence)):
        ter_C = ter_C - len(sequence) + 1
    stats = {}
    for base in "ACGT":
        total = sequence.count(base)
        # case 1: ----start========ter---->
        if (ter_C > start):                                   
            forward_count = sequence[start:ter_C].count(base)
            reverse_count = total - forward_count
        # case 2: ====ter--------start====>
        else:
            reverse_count = sequence[ter_C:start].count(base)
            forward_count = total - reverse_count
        stats[base] = (total, forward_count, reverse_count)
    return stats


def get_sequence_skew(sequence):
    """Returns the difference between the total number of 
    occurrences of G and the total number of occurrences of C in 
    the first i elements of the sequence. """
    skew = [0]
    for idx, element in enumerate(sequence):
        if sequence[idx] == 'G':
            skew.append(skew[idx] + 1)
        elif sequence[idx] == 'C':
            skew.append(skew[idx] -1)
        else:
            skew.append(skew[idx])
    return skew


def get_minimum_skew(sequence):
    """Returns a position in a sequence minimizing the skew."""
    min_skew = []
    skew = get_sequence_skew(sequence)
    m_skew = min(skew)
    return [idx for idx in range(len(sequence) + 1) if skew[idx] == m_skew]


def hamming_distance(sequence1, sequence2):
    """Return the HD form two inputed sequences"""
    assert len(sequence1 == len(sequence2)
    return len([(x,y) for x,y in zip(sequence1, sequence2) if x != y])


def hamming_pattern_count(sequence, pattern, distance):
    """Return the number of times the pattern with hamming distance d is found 
    in the inputed sequence."""
    return len([i for i in range(len(sequence) - len(pattern) + 1) if 
                hamming_distance(sequence[i:i+len(pattern)], pattern) <= distance])


def hamming_pattern_positions(sequence, pattern, distance):
    """Return list of positions of the pattern with hamming distance d from one input sequence."""
    return [i for i in range(len(sequence) - len(pattern) + 1) if 
                hamming_distance(sequence[i:i+len(pattern)], pattern) <= distance]


def get_neighbors(pattern, d):
    """Return list of all offsets of patterns with hamming distance d of `pattern"""
    if d == 0:
        return [pattern]
    if len(pattern) == 1:
        return set(['A', 'C', 'T', 'G'])
    neighborhood = set()
    neighbors = get_neighbors(pattern[1:], d)
    for kmer in neighbors:
        if hamming_distance(pattern[1:], kmer) < d:
            for char in ['A', 'C', 'T', 'G']:
                neighborhood.add(char + kmer)
        else:
            neighborhood.add(pattern[0] + kmer)
    return sorted(list(neighborhood))


def frequent_patterns_with_d_mismatches(sequence,k,distance):
   """Return list of most frequent offsets of patterns of k 
   length and with hamming distance d from one input sequence."""
    counts = {}
    for i in range(len(sequence) - k + 1):
        for kmer in get_neighbors(sequence[i:i+k], distance):
            counts[kmer] = counts.get(kmer, 0) + 1
    max_count = max(counts.values())
    return [kmer for kmer in counts if counts[kmer] == max_count]


def frequent_patterns_with_mismatches_and_rev_complements(sequence, k, distance):
    """Return list of most frequent offsets of patterns and it's reverse
    complements of k length and with hamming distance d from one input sequence."""
    counts = {}
    for i in range(len(sequence)-k+1):
        for subsequence in [sequence[i:i+k], reverse_complement(sequence[i:i+k])]:
            for kmer in get_neighbors(subsequence, distance):
                counts[kmer] = counts.get(kmer, 0) + 1
    max_count = max(counts.values())
    return [kmer for kmer in counts if counts[kmer] == max_count]


def probability(N, A, k, t, n):
    """N = lenght sequence, A = alphabet(num nucleotides), k=kmer lenght, t = times appearence and
    n = length of list dna sequences"""
    return (((N - t * (k-1))/t) / A ** (t * k)) * n


def entropy(sequence):
    p, lns = Counter(sequence), float(len(sequence))
    return -sum( count/lns * math.log(count/lns, 2) for count in p.values())


def profile_motifs(motifs):
    """Returns the frequence of the char[i] in the counts dictionary.
    It is the division of the counts of the char divide by the length
    of the motifs matrix"""
    k = len(motifs[0])
    assert all(len(x) == k for x in motifs), print('Motifs must have the same lengths')
    counts = motifs_count(motifs)
    return {char: [(count/len(motifs)) for count in counts[char]] for
           char in counts.keys()}


def motifs_count(motifs):
    """Returns the count of the nucleotides that appears in each
    columns of the motifs. Motifs are a matrix of strings(list of lists)"""
    count = defaultdict(int)
    k = len(motifs[0])
    t = len(motifs)
    for char in 'ACGT':
        count[char] = count.get(char, [0]*k)
    for i in range(t):
        for j in range(k):
            char = motifs[i][j]  # motifs[i] char[j]
            count[char][j] += 1  # ad to the count the char in j position in the cols
    return count


def motifs_scores(strings):
    score = 0
    consensus = consensus_sequence(strings)
    for i in range(len(strings)):
        for j in range(len(strings[0])):
            if strings[i][j] != consensus[j]:
                score += 1
    return score


def probability_motifs(string, profile):
    freqs = [profile[char][i] for i, char in enumerate(string)]
    return reduce(mul, freqs)


def get_most_probable_Kmer_from_profile(strings, k, profile):
    """Returns the most probable kmer of length k from a given string"""
    prob = -1  # must be -1 because you have some kmer with prob zero
    most_prob_kmer = strings[:k]
    for i in range(len(strings) - k + 1):
        seq = strings[i:i+k]
        prob_most = motif_probability(seq, profile)
        if prob_most > prob:
            prob = prob_most
            most_prob_kmer = seq
    return most_prob_kmer


def profile_motifs_pseudo_counts(motifs):
    """Returns the frequence of the char[i] in the counts dictionary.
    It is the division of the counts of the char divide by the length
    of the motifs matrix"""
    k = len(motifs[0])
    assert all(len(x) == k for x in motifs), print('Motifs must have the same lengths')
    counts = motifs_count_pseudo_counts(motifs)
    return {char: [val / (len(motifs) + 4) for val in value] for
           char, value in counts.items()}


def greedy_motif_search(string, k, t):
    best_motifs = []
    for i in range(0, t):
        best_motifs.append(string[i][0:k])
    n = len(string[0])
    for i in range(n - k + 1):
        motifs = []
        motifs.append(string[0][i:i+k])
        for j in range(1, t):
            profile = profile_motifs(motifs[0:j])
            motifs.append(get_most_probable_Kmer_from_profile(string[j], k, profile))
        if motifs_scores(motifs) < motifs_scores(best_motifs):
            best_motifs = motifs
    return (best_motifs)



def greedy_search_motif_with_pseudo_count(strings, k, t):
    best_motifs = [strings[i][0:k] for i in range(len(strings))]
    motifs = [strings[0][i:i + k] for i in range(len(strings[0]) - k + 1)]
    for motif in motifs:
        temp_motifs = [motif]
        for j in range(1, len(strings)):
            profile = profile_motifs_pseudo_counts(temp_motifs)
            temp_motifs.append(get_most_probable_Kmer_from_profile(strings[j], k, profile))
        if motifs_scores(temp_motifs) < motifs_scores(best_motifs):
            best_motifs = temp_motifs
    return best_motifs


def get_motifs(strings, k, profile):
    return [get_most_probable_Kmer_from_profile(string, k, profile) for string in strings]


def consensus_sequence_pseudo(motifs):
    """Returns the consensus sequence from a matrix of 
    motifs"""
    k = len(motifs[0])
    assert all(len(x) == k for x in motifs), print('Motifs must have the same lengths')
    counts = motifs_count_pseudo_counts(motifs)
    consensus = ''
    for i in range(k):
        count = 0
        freq_char = ''
        for char in 'ACGT':
            if counts[char][i] > count:
                count = counts[char][i]
                freq_char = char
        consensus += freq_char
    return consensus


def get_random_motifs(strings, k, t):
    random_motifs = []
    for i in range(len(strings)):
        idx = randint(0, len(strings[0]) - k)
        motif = strings[i][idx:idx+k]
        random_motifs.append(motif)
    return random_motifs


def RandomizedMotifSearch(strings, k, t):
    M = get_random_motifs(strings, k, t)
    BestMotifs = M
    score = motifs_scores(M)
    best_score = motifs_scores(BestMotifs)
    while True:
        profile = profile_motifs_pseudo_counts(M)
        M = get_motifs(strings, k, profile)
        if score < best_score:
            BestMotifs = M
        else:
            return BestMotifs


def RepeatedRandomizedMotifSearch(strings, k, t, N):
    BestScore = float('inf')
    BestMotifs = []
    for i in range(N):
        Motifs = RandomizedMotifSearch(strings, k, t)
        CurrScore = motifs_scores(Motifs)
        if CurrScore < BestScore:
            BestScore = CurrScore
            BestMotifs = Motifs
    return BestMotifs


def normalize_probabilities(probabilities):
    return {key: val / total for total in (sum(probabilities.values()),) for key, val in probabilities.items()}


def get_weighted_probabilities(probabilities):
    return choices(list(probabilities.keys()), weights=list(probabilities.values()))[0]


#random.uniform
def weighted_die(probabilities):
    s = 0.0
    for kmer, prob in probabilities.items():
        s += prob
        if uniform(0, 1) < s:
            return kmer


def profile_generated_string(strings, k, profile):
    probabilities = {}
    for i in range(len(strings) - k + 1):
        probabilities[strings[i:i+k]] = motif_probability(strings[i:i+k], profile)
    probabilities = normalize_probabilities(probabilities)
    return get_weighted_probabilities(probabilities)


# random.randint
def GibbsSampler(strings, k, t, num_searches):
    motifs = get_random_motifs(strings, k, t)
    best_motifs = motifs
    for j in range(1, num_searches):
        i = randint(0, t - 1)
        motifs.pop(i)
        profile = profile_motifs_pseudo_counts(motifs)
        motifs.insert(i, get_most_probable_Kmer_from_profile(strings[i], k, profile))
        if motifs_scores(motifs) < motifs_scores(best_motifs):
            best_motifs = motifs
    return best_motifs















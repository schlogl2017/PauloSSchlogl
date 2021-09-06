#!usr/bin/env python
# -*- coding: utf-8 -*-
import itertools
import alphabet
from sequence_utils import get_reverse_complement
from get_kmers import get_kmers_from_sequence
import fasta_parser


def iter_kmers(alphabet, k):
    """Generator function that yields every kmer (substring of length k) over an
    alphabet, which should be given as a Python set."""
    alphabets = [alphabet for i in range(k)]
    for kmer in itertools.product(*alphabets):
        yield ''.join(kmer)


def iter_palindromes(k, allowed_middle_characters, alphabet):
    """Generator function that yields every DNA reverse-complement palindrome
    of length k, including odd palindromes with center characters determined by
    allowed_middle_characters.
    allowed_middle_characters = ['G', 'T']
    """
    for kmer in iter_kmers(alphabet, k // 2):
        comp = get_reverse_complement(kmer)
        if k % 2 != 0:
            for character in allowed_middle_characters:
                yield kmer + character + comp
        else:
            yield kmer + comp


def iter_palindrome_range(k1, k2, allowed_middle_characters, alphabet):
    """Generator function that yields all DNA reverse
    complement palindromes from length k1 to k2, including
    palindromes of length k1 and k2."""
    for k in range(k1, k2 + 1):
        for palindrome in iter_palindromes(k, allowed_middle_characters, alphabet):
            yield palindrome


def gen_kmers(kmin, kmax, alphabet):
    """
    generates possible k-mers of range(kmin, kmax + 1)
    :param kmin: int, minimum kmer length
    :param kmax: int, maximum kmer length
    :param alphabet: str, accepted sequence alphabet for DNA, RNA, Amino Acids
    :return list of str, possible kmers
    """

    for n in range(kmin, kmax + 1):
        return [''.join(mer) for mer in itertools.product(alphabet, repeat=n)]


def gen_rev_palindromes(kmin, kmax, alphabet):
    """
    generate list of palindromes of length n,
    when kmin<=n<=kmax identical to their reverse complements
    :param kmin: int, min length of tested palindrome
    :param kmax:int, max length of tested palindrome
    :param bases: str, possible bases inserted in middle of odd palindrome
    :return: list, palindromes seqs identical to their reverse complements
    """
    dromes = []

    for n in range(kmin, kmax + 1):
        for left_mer in gen_kmers(n // 2, n // 2, alphabet):
            if n % 2 == 0:  # even palindrome
                dromes.append(left_mer + get_reverse_complement(left_mer))
            else:  # odd palindrome
                for midmer in alphabet:
                    dromes.append(left_mer + midmer + get_reverse_complement(left_mer))
    return dromes


def compute_stats(kmer_list, counts, N, max_e):
    """
	compute_stats computes the e-values for the supplied data.
	Pre-conditions:
	'kmer_list' - a list of kmers (for which stats will be produced)
	'counts' - any dictionary-type with k-mers as keys (min_k - 2 <= k <= max_k,
	where min_k and max_k are the bounds on the k-mer lengths in 'kmer_list')
	and counts as values.
	'N' - the total length of the sequence(s) read to produce 'counts'.
	'max_e' - the upper bound on e-values reported.
	Post-conditions:
	Retunrs a list of lists ('results') where results[i] is of the form
	[k-mer, observed count, expected count, z-score, e-value]
	"""

    # results is the list of list described in the docstring.
    results = []

    # number of tests, used to convert p-value to e-value.
    n = len(kmer_list)

    for kmer in kmer_list:

        k = len(kmer)

        observed = counts[kmer]
        expected = counts[kmer[:-1]] * counts[kmer[1:]] / counts[kmer[1:-1]]
        sigma = math.sqrt(expected * (1 - expected / (N - k + 1)))
        Z_score = (observed - expected) / sigma

        E_value_under = n * math.erfc(-Z_score / math.sqrt(2)) / 2  # E-value for under-rep
        E_value_over = n * math.erfc(Z_score / math.sqrt(2)) / 2  # E-value for over-rep

        if (E_value_under <= max_e):
            results.append([kmer, observed, expected, Z_score, E_value_under])
        elif (E_value_over <= max_e):
            results.append([kmer, observed, expected, Z_score, E_value_over])

    return results


def get_palindromes(alphabet, min_k, max_k):
    """Generates all DNA palindromes over the range from min_k to max_k.
    Inputs:
    min_k - minimum palindrome length (int)
    max_k - maximum palindrome length (int)
    Output:
    yields all possible DNA palindromes (str) of length min_k to max_k.
    Some definitions:
    A palindrome is defined as a sequence which is equal to its reverse-complement.
    Note: for odd length palindromes, the middle base does not need to be the same
    in the reverse-complement.
    Ex.: AAT is a legal palindrome even though its reverse-complement is ATT
    """
    for k in range(min_k, (max_k + 1)):
        for mer in itertools.product(alphabet, repeat=int(k / 2)):
            kmer = ''.join(mer)
            # even pal
            if k % 2 == 0:
                pal = kmer + get_reverse_complement(kmer)
                yield pal
            else:
                for base in alphabet:  # odd pal
                    pal = kmer + base + get_reverse_complement(kmer)
                    yield pal


if __name__ == '__main__':
    alphabet = alphabet.iupac_dna
    filename = "Data/test/Acidithiobacillus/chromosomes/NC_011206.1_Acidithiobacillus_ferrooxidans_ATCC_53993_complete_genome.fna.gz"
    for name, seq in fasta_parser.parse_fasta(filename):
        pal_list = list(get_palindromes(alphabet, 4, 6))
        print(len(pal_list))
        print(pal_list)

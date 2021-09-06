#!usr/bin/env python
# -*- coding: utf-8 -*-


from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict, Counter
from itertools import product
import alphabet
#from fasta_parser import fasta_item_counter, parse_fasta


def get_gc_content(sequence):
    """
    Finction to calculate the the gc content of a sequence.
    
    Inputs:
    
        sequence - a string representing a DNA sequence.
    
    Outputs:
    
        gc - a float representing the of (g + c) content of a sequence.
    
    """
    # get the sequence length and 
    # make all the sequence characters upper case
    seq_len, sequence = len(sequence), sequence.upper()
    # count all gs and cs
    c = sequence.count('C')
    g = sequence.count('G')
    # returns the gc content from a sequence
    # sum up the |Gs and Cs counts and divide 
    # by the sequence length
    return round((c + g) / seq_len, 4)


def get_at_content(gc):
    """Returns the at content of a genome.
    
    Inputs:
    
        gc - a float representing the of (g + c) content of a sequence.
        
    Outputs:
    
        at content - returns the AT content of a sequence. It is
                     defined as 1 - the GC content.
    
    """
    return 1 - gc


def count_umbiguous_bases(sequence):
    """
    Function to count all umbigous bases in a sequence.
    Ambigous bases are bases that are not in the sequence
    alphabet, ie. 'ACGT' for DNA sequences.

    Inputs:

        sequence - a string representing a DNA sequence.

    Outputs:

        integer - represent the number of ambigous bases found and
                  counted in a sequence.
    """
    sequence = sequence.upper()
    amb = ['N', 'R', 'Y', 'W', 'S', 'K', 'M']
    return sum({base: sequence.count(base) for base in amb}.values())


def get_at_gc_ratio(at, gc):
    """
    Calculates the at/gc ratio of a genome.
    
    Input:
    
       gc - a float representing the of (g + c) content of a sequence.
       at - a float representing the of (a + t) content of a sequence. 
    
    Outputs:
    
       a float representing the ratio of the gc and at content of a sequence.
    
    """
    return at / gc


def count_all_bases(sequence):
    """
    Calculates the nucleotides/base composition of a DNA sequence.
    
    Inputs:
    
        sequence - a string representing a DNA sequence.
    
    Outputs:
    
        all_bases - a dictionary-like object tha represent the count
                    of the nucleotides/bases that compound the sequence.
                    nucleotides are keys and the fruencies (as itegers) as values.
    
    """
    # create a set of bases
    bases = set(sequence)
    all_bases = defaultdict(int)
    # iterates in the base set
    for base in bases:
        # count the bases in the sequence
        all_bases[base] = sequence.count(base)
    return all_bases


def gc_content_sequence_window(sequence, as_overlap=False, k=20):
    """GC Content in a DNA/RNA sub-sequence length k. In
    overlapp windows of lenght k.
    
    Inputs:
    
        sequence - a string representing a DNA sequence.    
        as_overlap - boolean that represents if overlap is needed.
        k - a integer reppresenting the lengths of overlappig bases.
            Default is 20.
    
    Outputs:
    
        gc_content - an array-like object with 
    
    
    """
    # make sequence upper case and getting the length of it
    sequence, seq_len = sequence.upper(), len(sequence)
    # the array-like object to collect the data
    gc_content = []
    # non overlap sequence length
    non_overlap = range(0, len(sequence) - k + 1, k)
    # overlap sequence length
    overlap = range(0, seq_len - k + 1)
    # overlap is needed
    if as_overlap:
        # iterates to the overlap region
        for i in overlap:
            # creates the substring to count the gc_content
            subseq = sequence[i:i + k]
            # count and sum up the Gs and Cs counts
            g_c = subseq.count('C') + subseq.count('G')
            # collect the data in the array container
            gc_content.append(round(g_c / len(subseq), 4) * 100)
    # if non overlap is choosed
    else:
        # iterates to the mon overlap region
        for j in non_overlap:
            # creates the substring to count the gc_content
            subseq = sequence[j:j + k]
            # count and sum up the Gs and Cs counts
            g_c = subseq.count('C') + subseq.count('G')
            # collect the data in the array container
            gc_content.append(round(g_c / len(subseq), 4) * 100)
    return gc_content


def codon_frequency(sequence, codon_table):
    """
    Function to calculate the frequency of the codons in a sequence.
    
    Inputs:
        
        sequence - a string representing a DNA sequence.
        codon_table - an array-like container with all possible codon, defined
                      as trinucleotides or triplets (64 possible triplets)
    
    Outputs:
    
        counter - a dictionary-like object with codon/triplets counts from a
                      given input sequence. 
        
    """
    # initialize the counter with the list of triplets from codon_table
    counter = Counter(dict([(c, 0) for c in codon_table]))
    # create a list/array of all possible codons found in the input sequence
    triplets = [sequence.upper()[i:i + 3] for
                i in range(0, len(sequence), 3)]
    # filters the triplets list from sequences that don't have length of 3
    # nucleotides
    triplets = filter(lambda x: len(x) == 3, triplets)
    # updates counter with the triplets counts and return it
    return counter + Counter(triplets)


def gc_var(sequence, as_overlap=False, k=20):
    """
    Calculates the gc content variance in a sequence according to a 
    window of length k.
    
    Inputs:
        sequence - a string representing a DNA sequence.
        k - integer representing the length of the search window.
            default is 20.
    
    Outputs:
    
        log of the gc variantion in the sequence in a window space of
        length k.
    
    """
    # calculates the percent of gc content
    gc = get_gc_content(sequence) * 100
    # get the gc content in the window space as an array
    gc_i = np.array(gc_content_sequence_window(sequence, as_overlap, k=k))
    # get the len of the gc content in the window space
    len_gc_i = np.shape(gc_i)[0]
    # check the difference of each point 
    dif = gc_i - gc
    return np.log((1 / len_gc_i) * sum(abs(dif)))


def base_stats(sequence, alphabet, as_count=False, as_dict=False):
    """Calculates de frequency or the number of bases in a sequence.
    
    Inputs:
    
        sequence - string representing the sequence
        alphabet - a alphabet (strings characters) that compound the string sequence
        as_count - boolean set as False
        as_dict - boolean set as False
    
    Output:
    
        counts - as default returns a numpy array as frequencies (floats) or
                 as a dictionary-like object
    
    Examples:
    
    > baseFreqs(seq, 'ACGT', asCounts = False, asDict = False)
    array([0.25, 0.25, 0.25, 0.25])

    as_count - True, returns a numpy array of counts (integer)
    > baseFreqs('ACGTACGT', 'ACGT', asCounts = True, asDict = False)
    array([2, 2, 2, 2])

    as_dict - True and as_count as default (False) returns a dictionary as bases frequencies (float)
    > baseFreqs('ACGTACGT', 'ACGT', asCounts = False, asDict = True)
    {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}

    as_count True and as_dict True, returns a dictionary as base counts (integer)
    > baseFreqs('ACGTACGT', 'ACGT', asCounts = True, asDict = True)
    {'A': 2, 'C': 2, 'G': 2, 'T': 2}
    """
    # make the sequence upper case
    seq = sequence.upper()
    # count all bases in sequence and collect as an array
    counts = np.array([seq.count(i) for i in alphabet])
    # if is onle the counts
    if as_count:
        freqs = counts
    # other wise as frequencies
    else:
        freqs = counts / sum(counts * 1.0)
    # or as a dictionary like object
    if as_dict:
        return dict(zip(alphabet, freqs))
    else:
        return freqs


def get_strand_complement(sequence):
    """Returns the complement strand of the genome.
     
     Inputs:
        sequence - string representing the sequence   

    Outputs:
    
        sequence - string representing the complement of 
                   the string.    
    """
    # make the sequence upper case
    seq = sequence.upper()
    # table to change the complement characters
    change = str.maketrans('ACGT', 'TGCA')
    return seq.translate(change)


def get_reverse_complement(sequence):
    """
    Returns the reverse complement strand of the genome.

    Inputs:

        sequence - string representing the sequence.

    Outputs:

        reversed_complement_sequence - string representing the reversed
                                       sequence complement.
    """
    return get_strand_complement(sequence)[::-1]


def get_sequence_skew(sequence):
    """
    Calculates the difference between the total number of
    occurrences of G and the total number of occurrences of C in
    the first i elements of the sequence. 

    Inputs:
        sequence - string representing the sequence     
    
    Outputs:
    
        skew - an array-like object that represents the GC skew of the sequence.
    
    Ex: 
    > get_sequence_skew('ACAACGTAGCAGTAGCAGTAGT')
    [0, 0, -1, -1, -1, -2, -1, -1, -1, 0, -1, -1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 2, 2]
    
    """
    # make the sequence upper case
    sequence = sequence.upper()
    # start the array
    skew = [0]
    # iterates to the sequence elements and it indexes
    for idx, element in enumerate(sequence):
        # check if element[i] is a G
        # if so add 1
        if sequence[idx] == 'G':
            skew.append(skew[idx] + 1)
        # if the element[i] is a C
        # add to the array -1
        elif sequence[idx] == 'C':
            skew.append(skew[idx] - 1)
        else:
            # if it is not G or C add 0
            skew.append(skew[idx])
    return skew


def get_minimum_skew(sequence):
    """
    Calculates a position in a sequence minimizing the skew.
    
    Inputs:
        sequence - string representing the sequence.    
    
    Outputs:
    
        min_skew - an array-like object that represents the pistion 
                   where the GC skew is the minimized in the sequence.    
    
    Example:
    get_minimum_skew('ACAACGTAGCAGTAGCAGTAGT')
    [5]
    """
    # start the array 
    min_skew = []
    # calculates the sequence gc skew
    skew = get_sequence_skew(sequence)
    # get the minimized skew values
    m_skew = min(skew)
    # iterates to the length of the sequence
    # to get the index positions
    for idx in range(len(sequence) + 1):
        # if the position i has the same value 
        # as the minimum appende to the array
        if skew[idx] == m_skew:
            min_skew.append(idx)
    return min_skew


def plot_base_frequency_genome(x_data, y_data, x_label, y_label):
    """
    Make a plot from the base frequency distribution in a DNA sequence.
    
    Inputs:
    
        x_data - sizes of the window where the frequencies were calculated
        y_data - base_freqs obtained in the windows of the sequence.
        x_label - the labels for the x axes
        y_label - the labels for the y axes
        
    Ouputs:
    
        plot of base composition in a sequence window.
        
    """
    # color for the bases
    base_markers = {"A": "b-",
                    "C": "r-",
                    "G": "g-",
                    "T": "y-",
                    "N": "k-"}
    # drawing the plot
    fig = plt.figure(figsize=(16, 8))
    ax = fig.add_subplot(111)
    y_names = []
    for y in y_data:
        y_names.append(y)
        # adding colors to the lines representing the bases
        # the x and y data and the labels
        ax.plot(x_data, y_data[y], base_markers[y], label=y)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    ax.legend(y_names)
    plt.grid(True)


def base_content_slide_window(sequence, name, alphabet, window, step, plot=False):
    """
    Calculates the base/nucleotide frequencies in a window of size window and
    step and make a plot of the base distribution along of the sequence length.
    
     Inputs:
     
        sequence - string representing the sequence.    
        name - if plot True is a string representing the name of the plot.
        alphabet - a alphabet (strings characters) that compound the string sequence
        window - integer representing the length of the search window.
        step - integer representing the size of sequence overlap in the window.
        plot - boolean value, default is False other wise True to draw the plot.
    
    Outputs:
    
        base_freqs - a dictionary-like object mapping the bases to its frequencies.
        sizes - array-like object tha represents the size of the windows
       
    """
    # sequence as a string of upper cases characters
    # bases as a set of upper cases characters
    sequence, bases = sequence.upper(), set(alphabet.upper())
    # initialize the dictionary container and the array
    base_freqs = defaultdict(list)
    sizes = []
    # iterates to the bases and start filling the dictionary
    # with the keys and a empty array
    for base in bases:
        base_freqs[base] = base_freqs.get(base, [])
    # iterates to the sequence windows
    for i in range(0, len(sequence) - window + 1, step):
        # gets the sequence of length of the desired window
        subseq = sequence[i:i + window]
        # check if the length of the window is correct
        assert (len(subseq) == window), 'The lenght of the subsequence must have the same size of the window'
        # start calculating the frequencies
        # and feeding the containers
        for base in bases:
            freq = subseq.count(base) / len(subseq) * 100
            base_freqs[base].append(round(freq, 4))
        sizes.append((i + window))
    # if it is to plot the data
    if plot:
        plot_base_frequency_genome(sizes, base_freqs, 'Genome (kb)', 'Frequencies')
        plt.title(f"Base Distribuition in {name} genome")
        plt.savefig(f"{name}_base_freq_slidewindow_plot.pdf")
    # return the data
    return base_freqs, sizes


def strand_stats(sequence, alphabet, start):
    """
    Calculates the DNA strand base statistics over a sequence.
    
      Inputs:
     
        sequence - string representing the sequence.    
        alphabet - a alphabet (strings characters) that compound the string sequence.
        start - integer representing the initial start position where to start
                calculate the bases statistics.
    
    Outputs:
    
        strand_stat - a dictionary-like object mapping the strands statistics.         
    
    """
    # assure the characters are upper case
    alphabet = alphabet.upper()
    # assure the characters are upper case
    # get the sequence length
    seq_len, seq = len(sequence), sequence.upper()
    # get the middle position of the sequence
    half_gen = (seq_len // 2)
    # get the final position
    ter = (start + half_gen)
    # initialyze the container
    strand_stat = defaultdict(tuple)
    # for circular genomes
    if ter > seq_len:
        ter = ter - (seq_len + 1)
    # iterates through the alphabet
    # count the bases
    for base in alphabet:
        base_total = seq.count(base)
        # check the strand
        if ter > start:
            for_strand = seq[start:ter].count(base)
            rev_strand = (base_total - for_strand)
        else:
            rev_strand = seq[ter:start].count(base)
            for_strand = (base_total - rev_strand)
        # calculates the differences between strand
        dif = (for_strand - rev_strand)
        # get the data in the container
        strand_stat[base] = (base_total, for_strand, rev_strand, dif)
    return strand_stat


def print_strand_stats(strand_statistics):
    """
    Prints the strand statistics.
    
      Inputs:
     
        strand_statistics - a dictionary-like object mapping the strands statistics.    
    
    Outputs:
    
        Prints out the strands statistics.    
    
    """
    print('   Total\tFor\tRev\tDif')
    for base, count in strand_statistics.items():
        print(f'{base}: {str(count[0])}\t{str(count[1])}\t{str(count[2])}\t{str(count[3])}')


def get_counts_from_kmer_list(filenames_lst, alphabet, kmin, kmax):
    """
    Prints the strand statistics.

      Inputs:

        filenames_lst - a array-like object that represents the paths to the
                        source of the fasta files.
        alphabet - a alphabet (strings characters) that compound the string sequence.
        kmin - minimum DNA kmer length (int)
        kmax - maximum DNA kmer length (int)

    Outputs:

        dic_list - a array-like object with kmers counts or each input fasta file.

    """
    # initialize the array container
    dic_list = []
    # iterates through the file paths
    for filename in filenames_lst:
        # get the sequences and ids
        for n, seq in parse_fasta(filename):
            # append the counts to the array
            dic_list.append(count_kmers(seq, alphabet, kmin, kmax))
    return dic_list


def get_mean_all_kmers_genomes_counts(filenames_lst, alphabet, kmin, kmax):
    """
    Count all kmers from a file source.

      Inputs:

        filenames_lst - a array-like object that represents the paths to the
                        source of the fasta files.
        alphabet - a alphabet (strings characters) that compound the string sequence.
        kmin - minimum DNA kmer length (int)
        kmax - maximum DNA kmer length (int)

    Outputs:

        Returns - a list of sorted tuples keys-values.

    """
    # initialyze the counter
    all_kmers = Counter()
    # get the number of files for each genus
    f_len = len(filenames_lst)
    # get the file path
    for filename in filenames_lst:
        # get the sequences
        for name, seq in parse_fasta(filename):
            # update the counter with kmer counts from the different
            # genomes in the directory source
            all_kmers.update(count_kmers(seq, alphabet, kmin, kmax))
    # get the average count from genomes in the directory source
    # for each genus
    kmer_all_counts = {k: (cnt // f_len) for (k, cnt) in all_kmers.items()}
    return sorted(kmer_all_counts.items(), key=lambda k: k[0])


def get_all_possible_kmers(alphabet, kmin, kmax):
    """Returns a list of all possible combinations of k-mers of
    length k from a input alphabet.

    Inputs:

        alphabet - a alphabet (strings characters) that compound the string sequence
        kmin - minimum DNA kmer length (int)
        kmax - maximum DNA kmer length (int)

    Outputs:

        kmers - list of all possible combinations of k-mers of length k with length
                between kmin and kmax.

    """
    kmers = [''.join(letters) for n in range(kmin, kmax + 1)
             for letters in product(alphabet, repeat=n)]
    return kmers


def count_kmers(sequence, alphabet, kmin, kmax):
    """Generate all DNA k-mers over the entirety of a sequence, with length
    between kmin and kmax.

    Inputs:

        sequence - string representing the sequence
        alphabet - a alphabet (strings characters) that compound the string sequence
        kmin - minimum DNA kmer length (int)
        kmax - maximum DNA kmer length (int)

    Outputs:

        kmer_counts - a dictionary-like object with kmers counts from a
                      given input string with length between kmin and
                      kmax.

    """
    alphabet = set(alphabet)
    counts = defaultdict(int)
    for kmer in get_kmers_from_sequence(sequence, kmin, kmax):
        if set(kmer).issubset(alphabet):
            counts[kmer] = counts.get(kmer, 0) + 1
    return counts


def get_kmer_counts(kmer_list, kmer_counts):
    """
    Map all k-mers counts from a list of kmers of length k.

    Inputs:

        kmer_list - list or array-like object with kmers of length k.
        kmer_counts - a dictionary-like object with kmers counts from a
                      given input string with length between kmin and
                      kmax.

    Output:

        counts - a dictionary-like object with kmers counts length kma.

    counts will be utilyzed as input to calculate all kmer statistics.
    """
    counts = defaultdict(int)
    for kmer in kmer_list:
        counts[kmer] = counts.get(kmer, 0) + kmer_counts[kmer]
    return counts


def get_kmers_from_sequence(sequence, kmin, kmax):
    """
    Generate all DNA k-mers over the entirety of a sequence.

    Inputs:

        sequence - string where all kmers will be checked
        kmin: minimum DNA kmer length (int)
        kmax: maximum DNA kmer length (int)

    Output:

        yields all DNA kmers (str) of length kmin to kmax
    """
    limits = range(kmin, kmax + 1)
    seq_range = len(sequence) - kmax + 1
    for i in range(0, seq_range):
        for j in limits:
            yield sequence[i:i + j]


def kmer_positions(sequence, alphabet, k):
    """ returns the position of all k-mers in sequence as a dictionary"""
    mer_position = defaultdict(list)
    for i in range(1, len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        if all(base in set(alphabet) for base in kmer):
            mer_position[kmer] = mer_position.get(kmer, []) + [i]
    # combine kmers with their reverse complements
    pair_position = defaultdict(list)
    for kmer, pos in mer_position.items():
        krev = get_reverse_complement(kmer)
        if kmer < krev:
            pair_position[kmer] = sorted(pos + mer_position.get(krev, []))
        elif krev < kmer:
            pair_position[krev] = sorted(mer_position.get(krev, []) + pos)
        else:
            pair_position[kmer] = pos
    return pair_position


def get_kmer_count_slide_window(sequence, alphabet, window, step, kmin, kmax):
    slide_mer_count = defaultdict(Counter)
    for chunk, s, e in get_chunks(sequence, window, step):
        pos = '_'.join((str(s), str(e)))
        slide_mer_count[pos].update(count_kmers(chunk, alphabet, kmin, kmax))
    return slide_mer_count


def get_kmer_clumps(sequence, alphabet, k, window, times):
    clumps = defaultdict(list)
    kmers = kmer_positions(sequence, alphabet, k)
    for kmer, pos in kmers.items():
        clumps[kmer] = clumps.get(kmer, [])
        for i in range(len(pos) - times):
            end = i + times - 1
            while (pos[end] - pos[i]) <= window - k:
                end += 1
                if end >= len(pos):
                    break
            if end - i >= times:
                clumps[kmer].append((pos[i], end - i))
    return clumps


def sequence_cleaner(sequence, alphabet):
    """
    Clean up a sequence from not allowed characters.
    Input:
        sequence - sequence or a string
    Output:
        sequence - cleaned sequence or a string
    """
    seq = sequence.upper()
    sequence = [base for base in seq if base in alphabet]
    return ''.join(sequence)


def insert_wild_card(word, num_n=1):
    """
    Function to insert any wild card character in a word or string, generally
    the string must be a palindrome.

    Inputs:

        word - a string of characters that is best as a palindrome.
        num_n - a integer representing a the number of wild card characters
                to be inserted in the middle of the string word.

    Outputs:

        A tuple representing the word with n characters inserted in the middle and the
        inputed word without changes, if all definited condictions were achived, other
        wise return the word as it was inputed.
    """
    mid = len(word) // 2
    # to insert only one wild card character
    # with a predefinited condiction
    if num_n == 1 and is_palindrome(word) and len(word) % 2 != 0:
        return word[:mid] + 'N' + word[mid + 1:], word
    # the even words can receive two wild card chars
    # with a predefinited condiction
    elif num_n == 2 and is_palindrome(word) and len(word) % 2 == 0:
        return word[:mid - 1] + 'NN' + word[mid + 1:], word
    # only odd words can return a word with 3 chars wild cards
    # with a predefinited condiction
    elif num_n == 3 and word[:mid - 1] == get_reverse_complement(word[mid + 2:]) and len(word) % 2 != 0 and \
            len(word) >= 5:
        return word[:mid - 1] + 'NNN' + word[mid + 2:], word
    # if the condictions were not satisfied
    # return the word
    else:
        return word


def is_palindrome(string):
    """
     Function to check if a strings is palindromic or not.

     Inputs:

         string - a string of characters (a word, kmer, n-gram...)

     Outputs:

        boolean value - True if the string is palindromic other wise False.

     """
    k, mid = len(string), len(string) // 2
    # checking even palindromes
    if k % 2 == 0:
        return string[:mid] == get_reverse_complement(string[mid:])
    # checking odd palindromes
    else:
        return string[:mid] == get_reverse_complement(string[mid + 1:])


def get_counts(filename, alphabet, kmin, kmax):
    """
    Count all string countained in filename.

    Inputs:

        filename - a complete path to the file with th strings to count.
        alphabet - a alphabet (strings characters) that compound the string sequence
        min_k - minimum DNA kmer length (int)
        max_k - maximum DNA kmer length (int)

    Outputs:
        counter -  a dictionary-like object with kmers/string
                mapped to their counts

    """
    # get the list of kmers to count with length between kmin and kmax
    kmers_list = get_all_possible_kmers(alphabet, kmin, kmax)
    # initialyze the counter with all possible kmer with length
    # between kmin and kmax with zero counts
    counter = Counter(dict([(km, 0) for km in kmers_list]))
    # open and read in the kmers/string in the file
    with gzip.open(filename, 'rt') as fh:
        # iterates through the strings
        for line in fh:
            # make the adjustments int the strings
            kmer = line.replace('\n', '')
            # check if kmer/string is in the counter
            if kmer in counter:
                # if kmer is in add 1 other wise keep the zero count
                counter[kmer] += 1
    return counter


def get_counts_from_list(string_list, alphabet, kmin, kmax):
    """
    Count all string countained in filename.

    Inputs:

        filename - a complete path to the file.
        alphabet - a alphabet (strings characters) that compound the string sequence
        min_k - minimum DNA kmer length (int)
        max_k - maximum DNA kmer length (int)

    Outputs:
        counter -  a dictionary-like object with kmers/string
                mapped to their counts

    """
    # get the list of kmers to count with length between kmin and kmax
    kmers_list = get_all_possible_kmers(alphabet, kmin, kmax)
    # initialyze the counter with all possible kmer with length
    # between kmin and kmax with zero counts
    counter = Counter(dict([(km, 0) for km in kmers_list]))
    # open and read in the kmers/string in the file
    for string in string_list:
        # check if kmer/string is in the counter
        if string in counter:
            # if kmer is in add 1 other wise keep the zero count
            counter[string] += 1
    return counter


def entropy(sequence):
    """
    Calculates the entropy of a sequence.
    
    Inputs:
        sequence - a string that representes the sequence. 
    
    Outputs:
        a float number representing the sequence entropy.       
    """
    # get the bases counted and the length of the sequence
    p, lns = Counter(sequence), float(len(sequence))
    # calculates the sequence entropy
    return -sum((count/lns) * math.log(count/lns, 2) for count in p.values())

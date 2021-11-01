from collections import defaultdict

#the first 600 nucleotides from GenBank: AAHX01097212.1
rna = ("tcccccgcagcttcgggaacgtgcgggctcgggagggaggggcctggcgccgggcgcgcg"
       "cctgcgccccaccccgccccaccctggcgggtctcgcgcgcccggcccgcctcctgtcaa"
       "ccccagcgcggcggtcaggtggtccccagcccttggccccagcctccagcttcctggtcc"
       "ctcgggctctgagtcctgtctccggcagatcgcctttctgattgttctcctgcgcagctg"
       "gaggtgtatagcccctagccgagctatggtgcctcagcagatgtgaggaggtagtgggtc"
       "aggataaacccgcgcactccataataacgtgccagggctcagtgacttgggtctgcatta")

seq = rna.upper().replace('T', 'U')

#RNA codon table from http://en.wikipedia.org/wiki/Genetic_code
degenerated = (('GCU', 'GCC', 'GCA', 'GCG'),
               ('UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'),
               ('CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
               ('AAA', 'AAG'), ('AAU', 'AAC'), ('GAU', 'GAC'),
               ('UUU', 'UUC'), ('UGU', 'UGC'), ('CCU', 'CCC', 'CCA', 'CCG'),
               ('CAA', 'CAG'), ('UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'),
               ('GAA', 'GAG'), ('ACU', 'ACC', 'ACA', 'ACG'),
               ('GGU', 'GGC', 'GGA', 'GGG'), ('CAU', 'CAC'), ('UAU', 'UAC'),
               ('AUU', 'AUC', 'AUA'), ('GUU', 'GUC', 'GUA', 'GUG'),
               ('UAA', 'UGA', 'UAG'))

#prepare the dictio of degenerated codons
degen_dict = {}
for codons in degenerated:
    for codon in codons:
        degen_dict[codon] = codons

#query_codons
max_seq = len(seq)
query_codons = [seq[i:i+3] for i in range(0, max_seq, 3)]

#prepare dictio of counts:
counts = defaultdict(int)
for codon in query_codons:
    counts[codon] +=1

#actual calculation of frecuencies
data = {}
for codon in query_codons:
    if codon in  degen_dict:
        totals = sum(counts[deg] for deg in degen_dict[codon])
        frecuency = float(counts[codon]) / totals
    else:
        frecuency = 1.00

    data[codon] = frecuency

#print results
for codon, frecuency in data.iteritems():
    print "%s  -> %.2f" %(codon, frecuency)


#produces:
GUC  -> 0.57
AUA  -> 1.00
ACG  -> 0.50
AAC  -> 1.00
CCU  -> 0.25
UAU  -> 1.00
..........
GCU  -> 0.19
GAU  -> 1.00
UAG  -> 0.33
CUC  -> 0.38
UUA  -> 0.13
UGA  -> 0.33

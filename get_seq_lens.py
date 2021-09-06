#!usr/bin/env python
import os
from Bio import SeqIO


def get_sequence_lengths(fasta_filenames):
    """Returns dictionary of sequence lengths, keyed by organism.
    Biopython's SeqIO module is used to parse all sequences in the FASTA
    file corresponding to each organism, and the total base count in each
    is obtained.
    NOTE: ambiguity symbols are not discounted.
    """
    tot_lengths = {}
    for fn in fastafilenames:
        tt_lengths[os.path.splitext(os.path.split(fn)[-1])[0]] = \
            sum(o[len(s) for s in SeqIO.parse(fn, 'fasta')])
    return tot_lengths

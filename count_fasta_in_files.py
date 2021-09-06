#!usr/bin/env python
# -*- coding: utf-8 -*-
# args dir name get_fasta_files/filenames for fasta_item_counter
import sys
from fasta_parser import fasta_item_counter
from system_utils import get_fasta_files

if len(sys.argv) < 2:
    print('USAGE: count_fasta_in_files.py <dir_name>')
    sys.exit(1)

if __name__ == "__main__":
    for filename in get_fasta_files(sys.argv[1]):
        name = filename.split('/')[-1]
        num_seqs = fasta_item_counter(filename)
        print(f'The file {name} has: {num_seqs} sequencias')

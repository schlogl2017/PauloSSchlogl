#!/usr/bin/env python
import sys
import glob
import subprocess
import argparse



filenames = ['Downloads/Data/GCF_000016725.1_chr.fna']
for filename in filenames:
    name = filename.split('/')[2][:-4]
    print(filename)
    print(name)
    #subprocess.call(f'wordcount -sequence {filename} -wordsize {sys.argv[1]} -outfile {name}.wordcount.txt')
    subprocess.call(['wordcount', 
                     '-sequence', 
                     filename, 
                     '-wordsize', 
                     sys.argv[1], 
                     # minimum count is 1
                     '-mincount',
                     sys.argv[2],
                    '-outfile', 
                    name + f'.kmer_{sys.argv[1]}.tab'])



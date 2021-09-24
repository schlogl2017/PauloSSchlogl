#!/usr/bin/env bash

NAMES=$(ls Data/Test/)

for name in ${NAMES}
    do
        dsk -file Data/Test/${name}/Chromosomes/*.fna.gz -kmer-size $1 -abundance-min 0 -out ${name}_k$1 -out-dir Results/DSK_count
    done


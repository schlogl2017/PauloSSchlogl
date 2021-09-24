#!/usr/bin/env bash

NAMES=$(ls Data/Genomes_splitted/)

for name in $NAMES
    do
        find . -name ${name}_k$1.csv | xargs wc -l
    done

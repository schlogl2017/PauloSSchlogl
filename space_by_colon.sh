#!/usr/bin/env bash

NAMES=$(ls Results/kmer_count_jelly_k$1_dump/*)

for i in ${NAMES}
do
    echo "Changing space for colon in : ${i}"
    sed -e 's/\s\+/,/g' ${i} > "${i%.jf.dump}_k$1.csv"
done

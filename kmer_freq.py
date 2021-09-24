def get_kmer_frequency(csv_file):
    tot = 0
    freq = defaultdict(float)
    with open(csv_file, 'r') as fh:
        data = csv.reader(fh)
        for row in data:
            kmer, cnt = row[0], float(row[1])
            freq[kmer] = freq.get(kmer, 0.0) + cnt
            tot += cnt
    return {k: cnt/tot for k, cnt in freq.items()}

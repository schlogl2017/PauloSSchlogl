def write_fasta(seqs, fasta_file, wrap=80):
    """Write sequences to a fasta file.

    Parameters
    ----------
    seqs : dict[seq_id] -> seq
        Sequences indexed by sequence id.
    fasta_file : str
        Path to write the sequences to.
    wrap: int
        Number of AA/NT before the line is wrapped.
    """
    with open(fasta_file, 'w') as f:
        for gid, gseq in seqs.items():
            f.write('>{}\n'.format(gid))
            for i in range(0, len(gseq), wrap):
                f.write('{}\n'.format(gseq[i:i + wrap])) 

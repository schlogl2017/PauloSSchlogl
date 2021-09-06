import sys,os,subprocess
import numpy as np
from Bio import SeqIO
from pyfaidx import Fasta
import multiprocessing as mp
# rawbcf = mpileup_parallel(bam_files, ref, outpath, threads=10)

def get_chrom(filename):
    rec = list(SeqIO.parse(filename, 'fasta'))[0]
    return rec.id

def get_fasta_length(filename):
    """Get length of reference sequence"""

    refseq = Fasta(filename)
    key = list(refseq.keys())[0]
    l = len(refseq[key])
    return l



def worker(region,out,bam_files):
    """Run bcftools for a single region."""

    cmd = 'bcftools mpileup -r {reg} -O b -o {o} -f {r} {b}'.format(r=ref, reg=region, b=bam_files, o=out)
    #print (cmd)
    subprocess.check_output(cmd, shell=True)
    cmd = 'bcftools index {o}'.format(o=out)
    subprocess.check_output(cmd, shell=True)
    return

def mpileup_parallel(bam_files, ref, outpath, threads=4, callback=None):
    """Run mpileup in parallel over multiple regions, then concat vcf files.
    Assumes alignment to a bacterial reference with a single chromosome."""

    bam_files = ' '.join(bam_files)    
    rawbcf = os.path.join(outpath,'raw.bcf')
    tmpdir = 'tmp'
    chr = get_chrom(ref)  
    length = get_fasta_length(ref)

    #find regions
    bsize = int(length/(threads-1))
    x = np.linspace(1,length,threads,dtype=int)
    blocks=[]
    for i in range(len(x)):
        if i < len(x)-1:
            blocks.append((x[i],x[i+1]-1))

    pool = mp.Pool(threads)    
    outfiles = []    

    for start,end in blocks:        
        print (start, end)
        region = '{c}:{s}-{e}'.format(c=chr,s=start,e=end)
        out = '{o}/{s}.bcf'.format(o=tmpdir,s=start)
        f = pool.apply_async(worker, [region,out,bam_files])
        outfiles.append(out)

    pool.close()
    pool.join()

    #concat files
    cmd = 'bcftools concat {i} -O b -o {o}'.format(i=' '.join(outfiles),o=rawbcf)
    subprocess.check_output(cmd, shell=True)
    #remove temp files
    for f in outfiles:
        os.remove(f)
    return rawbcf



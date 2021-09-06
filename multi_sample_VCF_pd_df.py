import pandas as pd
import vcf
from gzip import open as gzopen

def vcf_to_dataframe(vcf_file):
    """
    Convert a multi sample vcf to dataframe. Records each samples FORMAT fields.
    Args:
        vcf_file: input multi sample vcf
    Returns: pandas DataFrame
    """

    ext = os.path.splitext(vcf_file)[1]
    if ext == '.gz':
        file = gzopen(vcf_file, "rt")
    else:
        file = open(vcf_file)
    vcf_reader = vcf.Reader(file,'r')
    res=[]
    cols = ['sample','REF','ALT','mut','DP','ADF','ADR','AD','chrom',
             'var_type','sub_type','start','end','QUAL']

    for rec in vcf_reader:
        x = [rec.CHROM, rec.var_type, rec.var_subtype, rec.start, rec.end, rec.QUAL]
        for sample in rec.samples:
            if sample.gt_bases == None:
                #no call
                mut=''
                row = [sample.sample, rec.REF, sample.gt_bases, mut, 0,0,0,0]
            elif rec.REF != sample.gt_bases:
                mut = str(rec.end)+rec.REF+'>'+sample.gt_bases
                cdata = sample.data
                row = [sample.sample, rec.REF, sample.gt_bases, mut, cdata[2],
                      cdata[4] ,cdata[5], cdata[6]] + x
            else:
                #call is REF
                mut = str(rec.end)+rec.REF              
                cdata = sample.data
                row = [sample.sample, rec.REF, sample.gt_bases, mut, cdata[2],
                      cdata[4] ,cdata[5], cdata[6]] + x

            res.append(row)
    res = pd.DataFrame(res,columns=cols)
    res = res[~res.start.isnull()]
    return res

df = vcf_to_dataframe('test.vcf.gz')  

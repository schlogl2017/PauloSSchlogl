import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from snpgenie import tools


#put the genome annotation into a dataframe (one row per feature)
g = tools.genbank_to_dataframe('LT708304_latest.gb')
#remove 'gene' features
g = g[g.feat_type!='gene']

#find the hypotheticals in the genome using the product field
labels=['hypothetical protein','conserved protein','conserved hypothetical',
        'unknown protein','uncharacterized','hypothetical alanine']
def find_unknown(x):
    for l in labels:
        if l in str(x).lower():           
            return True
    return False

g['hypothetical'] = g['product'].apply(find_unknown,1)
#filter
hypo = g[g.hypothetical==True]


def find_pfam_domains(g, res=None):
    """Find pfam domains"""

    new=[]
    for i,r in g.iterrows():
        if res!=None:
            if r.locus_tag in list(res.locus_tag):
                continue
        print (r.locus_tag)
        try:
            d = searchPfam(r.translation)
        except Exception as e:
            print (e)
            continue
        found = []
        for k in d:
            found = d[k]       
            found['locus_tag'] = r.locus_tag
            found['length'] = len(r.translation)
            found['product'] = r['product']
            new.append(found)    
    df=pd.DataFrame(new)
    if len(df)>0 and g is not None:
        locs = df['locations'].apply(pd.Series)
        df=pd.concat([df, locs], axis=1)
        df=df.drop(columns=['locations'])
    df = pd.concat([res,df])
    return df
    
    
res = find_pfam_domains(hypo)
res.to_csv('hypothetical_pfam.csv',index=False)

#get only the non DUF domains
known = res[~res.id.str.contains('DUF')]
#we can merge with the original dataframe if we want
cols=['locus_tag','translation','protein_id','gene']
known = known.merge(g[cols],on='locus_tag')
#get only the non DUF domains    

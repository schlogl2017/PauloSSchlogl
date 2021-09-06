import pandas as pd
from Bio import SeqIO



def features_to_dataframe(recs, cds=False):
    """Get genome records from a biopython features object into a dataframe
      returns a dataframe with a row for each cds/entry"""

    genome = recs[0]
    #preprocess features
    allfeat = []
    for (item, f) in enumerate(genome.features):
        x = f.__dict__
        q = f.qualifiers
        x.update(q)
        d = {}
        d['start'] = f.location.start
        d['end'] = f.location.end
        d['strand'] = f.location.strand
        for i in featurekeys:
            if i in x:
                if type(x[i]) is list:
                    d[i] = x[i][0]
                else:
                    d[i] = x[i]
        allfeat.append(d)
    df = pd.DataFrame(allfeat,columns=featurekeys)
    df['length'] = df.translation.astype('str').str.len()
    #print (df)
    df = check_tags(df)
    if cds == True:
        df = get_cds(df)
        df['order'] = range(1,len(df)+1)
    #print (df)
    if len(df) == 0:
        print ('ERROR: genbank file return empty data, check that the file contains protein sequences '\
               'in the translation qualifier of each protein feature.' )
    return df



def embl_to_dataframe(infile, cds=False):
    recs = list(SeqIO.parse(infile,'embl'))
    df = features_to_dataframe(recs, cds)
    return df


def index_genbank_features(gb_record, feature_type, qualifier):
    """Index features by qualifier value for easy access"""

    answer = dict()
    for (index, feature) in enumerate(gb_record.features):
        #print (index, feature)
        if feature.type==feature_type:
            if qualifier in feature.qualifiers:
                values = feature.qualifiers[qualifier]
                if not type(values) is list:
                    values = [values]
                for value in values:
                    if value in answer:
                        print ("WARNING - Duplicate key %s for %s features %i and %i" \
                           % (value, feature_type, answer[value], index))
                    else:
                        answer[value] = index
    return answer


def update_features(genome, df, field, key='locus_tag'):
    """Use a dataframe to update a genbank file with new or existing qualifier
    Returns a seqrecord object.
    """

    index = sequtils.index_genbank_features(genome, "CDS","locus_tag")    
    c=0
    for i,r in df.iterrows():        
        t=r[key]
        if t not in index:
            continue
        #print (r)
        #print (t,index[t])
        new = r[field]
        cds = genome.features[index[t]]
        if field not in cds.qualifiers:
            cds.qualifiers[field] = new
            c+=1
        else:
            curr = cds.qualifiers[field][0]
            #print (curr,new)
            if new != curr:
                cds.qualifiers[field] = new
                c+=1
    print ('updated %s features' %c)
    return genome


df = embl_to_dataframe('file.embl','embl')
#edit the dataframe in some way
feats = SeqIO.read('file.embl','embl')
new = update_features(feats, df, 'product')
SeqIO.write(new, 'new.embl', "embl")






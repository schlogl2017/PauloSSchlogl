import os, glob, time, subprocess, shutil
# get_fastani('fastani.out')

def run_test(n):
    #names are extracted from a dataframe here and used to make query
    #and reference lists
    names = list(asm.Assembly_nover)[:n]
    l=[]
    for f in glob.glob('assemblies/*.fa.gz'):
        n=os.path.basename(f).split('.')[0]
        if n in names:
            l.append(f)
    with open('query.txt', 'w') as infile:
        infile.write('\n'.join(l))
    shutil.copyfile('query.txt','reference.txt')
    #run
    cmd = 'fastANI --ql reference.txt --rl query.txt -o fastani.out -t 10'
    st=time.time()
    subprocess.check_output(cmd,shell=True)
    t = time.time()-st
    return t

times=[]
step=5
n=2
#run test with increasing increments of genomes
while n<=100:
    t=run_test(n)
    print (n,t)
    times.append((n,t))
    n+=step
    step+=2
    
    
def get_fastani(filename):
    """Get fastANI results into pairwise matrix"""
    import re
    df = pd.read_csv(filename,sep='\t',names=['query','ref','ident','x','y'])
    #these two lines are only needed to extract the label from the file names
    df['query'] = df['query'].apply(lambda x: re.split(r"[\./]+",x)[1])
    df['ref'] = df['ref'].apply(lambda x: re.split(r"[\./]+",x)[1])
    x = pd.pivot_table(df,index='query',values='ident',columns=['ref'])
    return x

    

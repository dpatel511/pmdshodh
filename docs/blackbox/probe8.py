import sys, random; sys.path.insert(0,'.')
from harness import *
random.seed(7)
def randread(i, strand):
    L=random.randint(25,90)
    ref=''.join(random.choice('ACGT') for _ in range(L))
    quals=[random.choice([2,10,20,30,37,40]) for _ in range(L)]
    # build molecule-orientation damage, then express in SAM orientation
    mref = ref if strand==0 else revcomp(ref)
    mseq=list(mref)
    for j in range(L):
        if mref[j]=='C' and random.random()<0.3: mseq[j]='T'
        if mref[j]=='G' and random.random()<0.3: mseq[j]='A'
    mseq=''.join(mseq)
    seq = mseq if strand==0 else revcomp(mseq)
    mq = quals if strand==0 else quals[::-1]
    s_ss=0; s_ds=0
    for j in range(L):
        z=j+1; y=L-j
        if mref[j]=='C' and mseq[j] in 'CT':
            s_ss+=llr_pmdtools_ss(mseq[j]=='C', z, y, mq[j])
    for j in range(L):  # ds is computed in SAM orientation (strand-symmetric)
        if ref[j]=='C' and seq[j] in 'CT': s_ds+=llr_pmdtools_ds(seq[j]=='C', j+1, quals[j])
        if ref[j]=='G' and seq[j] in 'GA': s_ds+=llr_pmdtools_ds(seq[j]=='G', L-j, quals[j])
    rec=sam_record(f'r{i}',ref,seq,flag=16*strand,qual=''.join(chr(33+q) for q in quals))
    return rec,s_ss,s_ds
data=[randread(i, i%2) for i in range(300)]
recs=[d[0] for d in data]
for mode,idx,args in (('ss',1,['--ss']),('ds',2,[])):
    rc,o,e=run(recs,['--printDS']+args)
    got=[float(l.split('\t')[-1]) for l in o.splitlines() if l.strip()]
    d=[abs(g-x[idx]) for g,x in zip(got,data)]
    print(mode,'n',len(got),'max diff',max(d), 'fwd max',max(d[0::2]),'rev max',max(d[1::2]))

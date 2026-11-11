import sys, random; sys.path.insert(0,'.')
from harness import *
# 1) multi-site additivity on random reads (forward strand, no indels)
random.seed(1)
recs=[]; preds={}
for i in range(200):
    L=random.randint(25,90)
    ref=''.join(random.choice('ACGT') for _ in range(L))
    seq=list(ref); quals=[random.choice([2,10,20,25,30,37,40]) for _ in range(L)]
    for j in range(L):
        if ref[j]=='C' and random.random()<0.25: seq[j]='T'
        if ref[j]=='G' and random.random()<0.25: seq[j]='A'
    seq=''.join(seq)
    s=0
    for j in range(L):
        if ref[j]=='C' and seq[j] in 'CT': s+=llr_pmdtools_ds(seq[j]=='C', j+1, quals[j])
        if ref[j]=='G' and seq[j] in 'GA': s+=llr_pmdtools_ds(seq[j]=='G', L-j, quals[j])
    recs.append(sam_record(f'r{i}',ref,seq,qual=''.join(chr(33+q) for q in quals))); preds[f'r{i}']=s
rc,o,e=run(recs,['--printDS'])
got=[float(l.split('\t')[-1]) for l in o.splitlines() if l.strip()]
print('records',len(recs),'lines',len(got),'rc',rc, e[:200])
diffs=[abs(g-p) for g,p in zip(got,preds.values())]
print('max |pmdtools - fitted model| over 200 random fwd reads:', max(diffs))

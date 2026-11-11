import sys; sys.path.insert(0,'.')
from harness import *
L=30
def one(pos, refb, readb, Q=40):
    ref=['A']*L; ref[pos-1]=refb; ref=''.join(ref)
    seq=list(ref); seq[pos-1]=readb; seq=''.join(seq)
    q=chr(33+Q)*L
    return sam_record(f'p{pos}{refb}{readb}Q{Q}',ref,seq,qual=q)
cases=[]
for pos in (1,2,3,5,10,15,20,25,28,29,30):
    cases+= [(pos,'C','T'),(pos,'C','C'),(pos,'G','A'),(pos,'G','G')]
out=[]
for c in cases:
    rc,o,e=run([one(*c)],['--printDS'])
    f=[x.strip() for x in o.split('\t')] if o.strip() else ['NA']
    pos,rb,qb=c
    z = pos if rb=='C' else L-pos+1
    mm = rb!=qb
    model = site_llr(not mm, p_match_ds(z,40,1), p_match_ds(z,40,0))
    print(f'pos{pos:>2} {rb}>{qb} z={z:>2}  pmdtools={f[-1]:>22}  paper={model:+.6f}  L_pmd={f[0] if len(f)>1 else ""}  L_null={f[1] if len(f)>1 else ""}')

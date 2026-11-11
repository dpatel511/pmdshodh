import sys; sys.path.insert(0,'.')
from harness import *
L=30; Q=[40]*L
def sc(ref, seq, args, qual):
    rc,o,e=run([sam_record('x',ref,seq,qual=''.join(chr(33+q) for q in qual))],['--printDS']+args); return float(o.split('\t')[-1])
# Give idx0 Q20 and idx29 Q30 so the two end sites have distinguishable terms.
qual=[20]+[40]*28+[30]
cC=lambda m: llr_pmdtools_ds(m,30,30)   # C at idx29: z from 5' = 30, Q30
gG=lambda m: llr_pmdtools_ds(m,30,20)   # G at idx0: z from 3' = 30, Q20
for mode in (['--CpG'],['--UDGhalf']):
    print('MODE',mode)
    for c in (('C','C'),('C','T'),('A','A')):
        for g in (('G','G'),('G','A'),('A','A')):
            ref=g[0]+'A'*28+c[0]; seq=g[1]+'A'*28+c[1]
            s=sc(ref,seq,mode,qual)
            opts={'none':0.0}
            if c[0]=='C': opts['C_last']=cC(c[0]==c[1])
            if g[0]=='G': opts['G_first']=gG(g[0]==g[1])
            if c[0]=='C' and g[0]=='G': opts['both']=cC(c[0]==c[1])+gG(g[0]==g[1])
            expl=[k for k,v in opts.items() if abs(v-s)<1e-9]
            print(f'  idx0 {g[0]}>{g[1]}  last {c[0]}>{c[1]}: {s:+.6f} -> {expl}')

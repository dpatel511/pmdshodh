import sys; sys.path.insert(0,'.')
from harness import *
L=30
def sc(ref, seq, args):
    rc,o,e=run([sam_record('x',ref,seq)],['--printDS']+args); return float(o.split('\t')[-1])
t=lambda m,z: llr_pmdtools_ds(m,z,40)
print('Target: C at LAST position (idx 29). Partner at idx 0 varies. Mode --CpG')
for tgt in (('C','C'),('C','T')):
    for partner in (('G','G'),('G','A'),('A','G'),('A','A')):
        ref=partner[0]+'A'*28+tgt[0]; seq=partner[1]+'A'*28+tgt[1]
        s=sc(ref,seq,['--CpG'])
        cterm=t(tgt[0]==tgt[1],30)            # C at z=30 from 5'
        gterm=t(partner[0]==partner[1],30) if partner[0]=='G' else None  # G at idx0: z=30 from 3'
        expl=[]
        for name,val in (('C_last',cterm),('G_first',gterm)):
            if val is not None and abs(s-val)<1e-9: expl.append(name)
        if gterm is not None and abs(s-(cterm+gterm))<1e-9: expl.append('both')
        if abs(s)<1e-12: expl.append('none')
        print(f' C_last {tgt[0]}>{tgt[1]}, idx0 ref {partner[0]} read {partner[1]}: {s:+.6f} -> {expl}')

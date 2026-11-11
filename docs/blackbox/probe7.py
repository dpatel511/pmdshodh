import sys; sys.path.insert(0,'.')
from harness import *
L=30; A='A'
def cols(pos,rb,qb,Q=40,Lr=L):
    ref=[A]*Lr; ref[pos-1]=rb; ref=''.join(ref); seq=list(ref); seq[pos-1]=qb
    rc,o,e=run([sam_record('x',ref,''.join(seq),qual=chr(33+Q)*Lr)],['--printDS','--ss'])
    return [float(v) for v in o.split('\t')]
pn=p_match_ds_null_pmdtools
for Q in (40,20):
  print('Q',Q)
  for pos in (1,2,3,5,10,15,16,20,29,30):
    for mm in (True,False):
        lp,ln,_,s=cols(pos,'C','T' if mm else 'C',Q)
        z=pos; y=L-pos+1
        # candidate: PMD model as S1 with both D_z and D_y; null as observed
        cand_s1=p_match_ss(z,y,Q,True)
        print(f' pos{pos:>2} {"mm" if mm else "m "} L_pmd={lp:.12f} L_null={ln:.12f}  S1 P(obs)={(1-cand_s1) if mm else cand_s1:.12f}  null_ds P(obs)={(1-pn(Q)) if mm else pn(Q):.12f}  ratio_null={ln/((1-pn(Q)) if mm else pn(Q)):.9f}')

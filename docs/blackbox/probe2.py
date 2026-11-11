import sys; sys.path.insert(0,'.')
from harness import *
L=30
def rec(Q, mm=True):
    ref='C'+'A'*(L-1); seq=('T' if mm else 'C')+'A'*(L-1)
    return sam_record('q',ref,seq,qual=chr(33+Q)+'I'*(L-1))
def get(Q,args=()):
    rc,o,e=run([rec(Q)],['--printDS']+list(args))
    f=[float(x) for x in o.split('\t')]
    return f
print('Q   L_pmd(mm,z=1)           L_null(mm)            implied e(pmd)       implied pi_null')
for Q in (0,2,5,10,13,20,30,35,40,41,42,45,50,60,70,93):
    try: lp,ln,_,s=get(Q)
    except Exception as ex: print(Q,'ERR',ex); continue
    D1=0.31; pi=0.001
    # solve e from PMD model: P(mm)=1-[(1-pi)(1-e)(1-D)+(1-pi)eD+pi e(1-D)]
    a=(1-pi)*(1-D1); b=-(1-pi)*(1-D1)+(1-pi)*D1+pi*(1-D1)
    e_p=((1-lp)-a)/b
    # null: P(mm)=pi+e-2 pi e with e=e_p
    pin=(ln-e_p)/(1-2*e_p)
    print(f'{Q:<3} {lp:.15f}  {ln:.15f}  {e_p:.6e}  {pin:.9f}   PMDS={s:.6f}')
print()
for args in (['--polymorphism_ancient','0.001'],['--polymorphism_ancient','0.005'],['--polymorphism_contamination','0.001'],['--polymorphism_contamination','0.005'],['--PMDpparam','0.5'],['--PMDconstant','0.0']):
    lp,ln,_,s=get(40,args); print(args, f'L_pmd={lp:.12f} L_null={ln:.12f} PMDS={s:.6f}')

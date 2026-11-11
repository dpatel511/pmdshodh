import sys; sys.path.insert(0,'.')
from harness import *
def sc(rec,args=()):
    rc,o,e=run([rec],['--printDS']+list(args))
    return (o.split('\t')[-1].strip()[:9] if o.strip() else 'NO-OUT') + (f' ERR:{e.strip()[-80:]}' if e.strip() else '')
def mk(ref,seq,flag=0,qual=None,mapq=37): return sam_record('x',ref,seq,flag=flag,qual=qual,mapq=mapq)
A='A'
cases={
 'C>T z1, next ref A (non-CpG)': ('CA'+A*28,'TA'+A*28),
 'C>T z1, next ref G (CpG)':     ('CG'+A*28,'TG'+A*28),
 'C>T z5 non-CpG':                (A*4+'CA'+A*24, A*4+'TA'+A*24),
 'C>T z5 CpG':                    (A*4+'CG'+A*24, A*4+'TG'+A*24),
 'C>T z5 CpG, read G>A too':      (A*4+'CG'+A*24, A*4+'TA'+A*24),
 'G>A 3p z1, prev ref A':         (A*28+'AG',A*28+'AA'),
 'G>A 3p z1, prev ref C (CpG)':   (A*28+'CG',A*28+'CA'),
 'G>A 3p z5 non-CpG':             (A*24+'AG'+A*4,A*24+'AA'+A*4),
 'G>A 3p z5 CpG':                 (A*24+'CG'+A*4,A*24+'CA'+A*4),
 'C>T at last pos (context unknown)': (A*29+'C',A*29+'T'),
 'G>A at first pos (context unknown)': ('G'+A*29,'A'+A*29),
}
modes=[[],['--CpG'],['--UDGplus'],['--UDGhalf'],['--UDGminus'],['--noCpG']]
print(f'{"case":38}'+''.join(f'{" ".join(m) or "default":>12}' for m in modes))
for k,(r,s) in cases.items():
    print(f'{k:38}'+''.join(f'{sc(mk(r,s),m):>12}' for m in modes))
print('--- CpG with ss')
for k in ('C>T z5 CpG','C>T z5 non-CpG','C>T at last pos (context unknown)'):
    r,s=cases[k]; print(f'{k:38}', sc(mk(r,s),['--ss','--CpG']), sc(mk(r,s),['--ss','--UDGhalf']))
print('--- filters')
r,s=('C'+A*29,'T'+A*29)
q10='+'+'I'*29
print('C>T z1 at Q10, -q 20:', sc(mk(r,s,qual=q10),['-q','20']), '| no -q:', sc(mk(r,s,qual=q10)))
print('mapq 37, -m 30:', sc(mk(r,s),['-m','30']), '| -m 37:', sc(mk(r,s),['-m','37']),'| -m 38:', sc(mk(r,s),['-m','38']))

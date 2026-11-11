import sys; sys.path.insert(0,'.')
from harness import *
def ds(recs,args=()):
    rc,o,e=run(recs,['--printDS']+list(args)); return rc,[l for l in o.splitlines() if l.strip()],e
def kept(recs,args):
    rc,o,e=run(recs,list(args)); return rc,[l.split('\t')[0] for l in o.splitlines() if l.strip() and not l.startswith('@')],e
L=30
base=sam_record('t1','C'+'A'*29,'T'+'A'*29)
rc,l,e=ds([base]); s=float(l[0].split('\t')[-1]); print('score',repr(s))
s_str=l[0].split('\t')[-1].strip()
print('THRESHOLD exact:', kept([base],['--threshold',s_str]))
print('THRESHOLD 3dp-rounded', kept([base],['--threshold',f'{s:.3f}']))
print('THRESHOLD slightly below', kept([base],['--threshold',str(s-1e-9)]))
print('UPPER exact:', kept([base],['--upperthreshold',s_str]))
print('THRESHOLD int 5 vs score', kept([base],['--threshold','5']))
# I=0 read: no C or G in ref
z=sam_record('noinf','A'*30,'A'*30)
print('I=0 printDS:', ds([z]))
for t in ('0','-1','-0.0001'): print(' I=0 threshold',t, kept([z],['--threshold',t]))
# other mismatch types ignored?
o1=sam_record('CtoA','C'+'A'*29,'A'*30); print('C>A at z1:',ds([o1])[1])
o2=sam_record('NinRead','C'+'A'*29,'N'+'A'*29); print('N read at C:',ds([o2])[1])
# mapq 0 and unmapped flag
print('mapq0:',ds([sam_record('mq0','C'+'A'*29,'T'+'A'*29,mapq=0)])[1])
print('unmapped flag 4:',ds([sam_record('um','C'+'A'*29,'T'+'A'*29,flag=4)])[1:])
# missing MD
r=sam_record('nomd','C'+'A'*29,'T'+'A'*29).rsplit('\t',1)[0]
print('no MD:',ds([r]))
# lowercase in MD / n in MD (eager #324)
print('MD with n:',ds([sam_record('mdn','C'+'A'*29,'T'+'A'*29,md='0n29')]))

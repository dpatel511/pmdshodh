import sys; sys.path.insert(0,'.')
from harness import *
def sc(rec,args=()):
    rc,o,e=run([rec],['--printDS']+list(args))
    return (o.split('\t')[-1].strip() if o.strip() else 'NO-OUTPUT') + (f' ERR:{e.strip()[-100:]}' if e.strip() else '')
pn=p_match_ds_null_pmdtools
def ss_llr(match,z,y,Q=40): return site_llr(match,p_match_ss(z,y,Q,True),pn(Q))
L=30; A='A'
def one(pos,rb,qb,flag=0):
    ref=[A]*L; ref[pos-1]=rb; ref=''.join(ref); seq=list(ref); seq[pos-1]=qb
    return sam_record(f'{rb}{qb}{pos}f{flag}',ref,''.join(seq),flag=flag)
print('--- ss, forward strand')
for pos,rb,qb in [(1,'C','T'),(30,'C','T'),(15,'C','T'),(2,'C','T'),(29,'C','T'),(1,'C','C'),(30,'C','C'),(30,'G','A'),(1,'G','A'),(30,'G','G')]:
    z=pos; y=L-pos+1
    print(f'{rb}>{qb} pos{pos:>2}: ss={sc(one(pos,rb,qb),["--ss"]):<28} ds={sc(one(pos,rb,qb)):<22} paperSS(z={z},y={y})={ss_llr(rb==qb,z,y):+.6f}')
print('--- ss, reverse strand (flag 16)')
for pos,rb,qb in [(1,'G','A'),(30,'G','A'),(1,'C','T'),(30,'C','T'),(15,'G','A'),(1,'G','G')]:
    print(f'{rb}>{qb} pos{pos:>2} flag16: ss={sc(one(pos,rb,qb,16),["--ss"]):<28} ds={sc(one(pos,rb,qb,16))}')

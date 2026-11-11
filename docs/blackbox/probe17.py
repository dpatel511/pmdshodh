import sys; sys.path.insert(0,'.')
from harness import *
def sc(ref, seq, args):
    rc,o,e=run([sam_record('x',ref,seq)],['--printDS']+args); return o.split('\t')[-1].strip()
for L in (30,109):
    a='A'*(L-2)
    for prev in 'AG':
        for rb,qb in (('C','C'),('C','T')):
            ref=a+prev+rb; seq=a+prev+qb
            print(f'L={L} last base {rb}>{qb} prev {prev}: default {sc(ref,seq,[]):>22} UDGhalf {sc(ref,seq,["--UDGhalf"]):>22} CpG {sc(ref,seq,["--CpG"])}')
    for nxt in 'AC':
        for rb,qb in (('G','G'),('G','A')):
            ref=rb+nxt+a; seq=qb+nxt+a
            print(f'L={L} first base {rb}>{qb} next {nxt}: default {sc(ref,seq,[]):>22} UDGhalf {sc(ref,seq,["--UDGhalf"]):>22} CpG {sc(ref,seq,["--CpG"])}')

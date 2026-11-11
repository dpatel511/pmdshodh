import sys; sys.path.insert(0,'.')
from harness import *
def sc(ref, seq, args, flag=0):
    rc,o,e=run([sam_record('x',ref,seq,flag=flag)],['--printDS']+args); return o.split('\t')[-1].strip()
a='A'*26
print('last C, first base G   :', sc('GA'+a+'AC','GA'+a+'AC',['--CpG']), '(C last kept if next wraps to first G)')
print('last C, first base A   :', sc('AA'+a+'AC','AA'+a+'AC',['--CpG']))
print('first G, last base C   :', sc('GA'+a+'AC','TA'+a+'AT',['--CpG']), '(both mismatched? no: read T at C? check)')
print('first G>A, last C>C    :', sc('GA'+a+'AC','AA'+a+'AC',['--CpG']))
print('first G>A, last A      :', sc('GA'+a+'AA','AA'+a+'AA',['--CpG']))
print('ss, last C, first G    :', sc('GA'+a+'AC','GA'+a+'AT',['--ss','--CpG']))
print('ss rev, SAM first G>A, last C:', sc('GA'+a+'AC','AA'+a+'AC',['--ss','--CpG'],flag=16))

import sys; sys.path.insert(0,'.')
from harness import *
def ok(k):
    ref='A'*(k-1)+'C'; rc,o,e=run([sam_record('x',ref,ref,qual='I'*k)],['--printDS']); return rc==0 and o.strip()!=''
lo,hi=100,1000
assert ok(lo) and not ok(hi)
while hi-lo>1:
    m=(lo+hi)//2
    if ok(m): lo=m
    else: hi=m
print('largest z for a C site (5p) that works:',lo,' first failing:',hi)
# G site from 3p: G at first base of read length k -> z=k
def okg(k):
    ref='G'+'A'*(k-1); rc,o,e=run([sam_record('x',ref,ref,qual='I'*k)],['--printDS']); return rc==0 and o.strip()!=''
print('G at z',lo,okg(lo),' z',hi,okg(hi))

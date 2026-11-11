import sys; sys.path.insert(0,'.')
from harness import *
A='A'
def sc(rec,args=()):
    rc,o,e=run([rec],['--printDS']+list(args))
    return (f"{float(o.split(chr(9))[-1]):.6f}" if o.strip() else 'NO-OUT') + (f' ERR:{e.strip()[-90:]}' if e.strip() else '')
t1=llr_pmdtools_ds(False,5,40); print('ref: C>T z5 term', round(t1,6), ' z1', round(llr_pmdtools_ds(False,1,40),6))
modes=[[],['--CpG'],['--UDGhalf']]
def row(label, rec): print(f'{label:58}', ' | '.join(f'{" ".join(m) or "default"}: {sc(rec,m)}' for m in modes))
# C at read pos5 (ref pos5) then deletion of G: ref ...C G A..., read ...T A...
ref=A*4+'C'+'G'+A*24; seq=A*4+'T'+A*24
row('C>T, next ref base G is deleted (5M1D24M)', sam_record('a',ref,seq,cigar='5M1D24M'))
ref=A*4+'C'+'A'+A*24; seq=A*4+'T'+A*24
row('C>T, next ref A deleted (5M1D24M), then A', sam_record('a2',ref,seq,cigar='5M1D24M'))
# C then insertion then ref G: ref ...C G..., read ...T [GG ins] G...
ref=A*4+'CG'+A*24; seq=A*4+'T'+'TT'+'G'+A*24
row('C>T, insertion, then ref G (5M2I25M)', sam_record('b',ref,seq,cigar='5M2I25M'))
# C before deletion of A then G: ref C A G -> read T _ G (CpG broken by deleted A)
ref=A*4+'CAG'+A*23; seq=A*4+'T'+'G'+A*23
row('C>T, deleted A, then G (5M1D24M)', sam_record('c',ref,seq,cigar='5M1D24M'))
# Soft clip: first aligned base is C>T after 3S
ref='C'+A*29; seq='GGG'+'T'+A*29
row('3S then C>T at first aligned (z=4)', sam_record('d',ref,seq,cigar='3S30M'))
ref=A*29+'G'; seq=A*29+'A'+'TTT'
row('G>A last aligned then 3S (z=4)', sam_record('e',ref,seq,cigar='30M3S'))
# CpG spanning soft clip: C last aligned, clipped next base G
ref=A*29+'C'; seq=A*29+'T'+'GGG'
row('C>T last aligned, then 3S starting G', sam_record('f',ref,seq,cigar='30M3S'))
# = and X ops
ref='C'+A*29; seq='T'+A*29
row('CIGAR 1X29=', sam_record('g',ref,seq,cigar='1X29='))
row('CIGAR 1M5N29M (ref skip)', sam_record('h',ref,seq,cigar='1M5N29M',md='0C29'))
row('lowercase read t at C', sam_record('i','C'+A*29,'t'+'a'*29,md='0C29'))
row('lowercase c match', sam_record('j','C'+A*29,'c'+'a'*29,md='30'))
print('--- ss + modes')
modes=[['--ss'],['--ss','--CpG'],['--ss','--UDGhalf']]
row('ss: C>T at first base, non-CpG', sam_record('k','C'+A*29,'T'+A*29))
row('ss: C>T at last base, non-CpG', sam_record('l',A*29+'C',A*29+'T'))
row('ss: C>T z5 non-CpG', sam_record('m',A*4+'C'+A*25,A*4+'T'+A*25))
row('ss flag16: G>A at SAM pos1 (molecule 3p C)', sam_record('n','G'+A*29,'A'+A*29,flag=16))
row('ss flag16: G>A at SAM pos30 (molecule 5p C)', sam_record('o',A*29+'G',A*29+'A',flag=16))
row('ss flag16: CpG at SAM 4-5, read CA (molecule C>T in CpG)', sam_record('p',A*3+'CG'+A*25,A*3+'CA'+A*25,flag=16))

import sys; sys.path.insert(0,'.')
from harness import *
def show(label, ref, seq, cigar):
    rec=sam_record('x',ref,seq,cigar=cigar)
    rc,o,e=run([rec],['--printalignments'])
    ls=o.splitlines()
    md=rec.split('MD:Z:')[1]
    print(f'== {label} CIGAR {cigar} MD {md}'); print('  read', ls[0] if ls else ''); print('  ref ', ls[2] if len(ls)>2 else '', e.strip()[-80:])
A='A'
# reference aligned bases (excl. clipped) ; deletion 'GG' between
ref = 'ACGTACGTAC'+'GG'+'TTCCAAGGTT'   # 10 + 2 deleted + 10
seq_al = 'ACGTACGTAC'+'TTCCAAGGTT'
show('no clip, del', ref, seq_al, '10M2D10M')
show('2S, del', ref, 'NN'+seq_al, '2S10M2D10M')
show('5S, del', ref, 'NNNNN'+seq_al, '5S10M2D10M')
show('del, 3S at end', ref, seq_al+'NNN', '10M2D10M3S')
show('2S, ins (no del)', 'ACGTACGTAC'+'TTCCAAGGTT', 'NN'+'ACGTACGTAC'+'GG'+'TTCCAAGGTT', '2S10M2I10M')
show('2S, 2 deletions', 'ACGTA'+'GG'+'CGTAC'+'G'+'TTCCAAGGTT', 'NN'+'ACGTA'+'CGTAC'+'TTCCAAGGTT', '2S5M2D5M1D10M')

import sys; sys.path.insert(0,'.')
from harness import *
def show(label, ref, seq, cigar):
    rec=sam_record('x',ref,seq,cigar=cigar); rc,o,e=run([rec],['--printalignments']); ls=o.splitlines()
    print(f'== {label}: CIGAR {cigar} MD {rec.split("MD:Z:")[1]}\n   read {ls[0] if ls else ""}\n   ref  {ls[2] if len(ls)>2 else ""} {e.strip()[-80:]}')
R1='ACGTACGTAC'; R2='TTCCAAGGTT'; R3='CAGTCAGTCA'
show('D then I', R1+'GG'+R2+R3, R1+R2+'NN'+R3, '10M2D10M2I10M')
show('I then D', R1+R2+'GG'+R3, R1+'NN'+R2+R3, '10M2I10M2D10M')
show('I only', R1+R2, R1+'NN'+R2, '10M2I10M')
show('I with =', R1+R2, R1+'NN'+R2, '10=2I10=')
show('2S with =', R1+R2, 'NN'+R1+R2, '2S20=')
show('2S with M', R1+R2, 'NN'+R1+R2, '2S20M')
show('20= then 2S', R1+R2, R1+R2+'NN', '20=2S')

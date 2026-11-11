import sys; sys.path.insert(0,'.')
from harness import *
A='A'
def show(label, recs, args=('--printDS',)):
    rc,o,e=run(recs,list(args)); print(f'{label}: rc={rc} out={o.strip()[:160]!r} err={e.strip()[-160:]!r}')
show('QUAL *', [sam_record('nq','C'+A*29,'T'+A*29,qual='*')])
show('SEQ * (secondary-style)', [sam_record('ns','C'+A*29,'*',cigar='30M',qual='*',md='0C29',flag=256)])
# heavy damage long read: 400 C>T mismatches
L=800; ref=('CA'*400); seq=('TA'*400)
show('400 C>T mismatches, 800bp', [sam_record('heavy',ref,seq,qual='I'*L)])
s=sum(llr_pmdtools_ds(False,i+1,40) for i in range(0,L,2)); print('  sum-of-logs model:', s)
# many matches: 4000 C matches
L=8000; ref=('CA'*4000)
show('4000 C matches, 8000bp', [sam_record('longm',ref,ref,qual='I'*L)])
print('  sum-of-logs model:', sum(llr_pmdtools_ds(True,i+1,40) for i in range(0,L,2)))

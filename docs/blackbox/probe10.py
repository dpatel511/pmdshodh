import sys; sys.path.insert(0,'.')
from harness import *
A='A'
def show(label, recs, args):
    rc,o,e=run(recs,args)
    print(f'== {label} {args}\n   rc={rc} out={o.strip()[:300]!r}\n   err={e.strip()[-300:]!r}')
# issue #5: MD without 0 between deletion and mismatch. ref: AAAA [CCTT deleted] TC AAAA ; read: AAAA AA AAAA? build: read mismatches T->A? Use MALT example shape.
ref='A'*4+'CCTT'+'TC'+'A'*24     # 34 ref bases
seq='A'*4+'CT'+'A'*24            # read has C at ref T (mismatch) and T at ref C (C>T!) right after deletion
good=sam_record('md_ok',ref,seq,cigar='4M4D26M')
bad=good.replace('MD:Z:'+make_md(ref,seq,'4M4D26M'),'MD:Z:'+make_md(ref,seq,'4M4D26M').replace('^CCTT0','^CCTT'))
print('good MD',good.split('MD:Z:')[1],' bad MD',bad.split('MD:Z:')[1])
show('#5 spec MD',[good],['--printDS']); show('#5 MALT-style MD',[bad],['--printDS'])
# issue #9: DS:Z: inside quality string
r=sam_record('ds_in_qual','C'+A*29,'T'+A*29,qual='DS:Z:'+'I'*25)
show('#9 DS:Z: in QUAL',[r],['--printDS']); show('#9 DS:Z: in QUAL',[r],['--threshold','0','--writesamfield'])
# existing DS tag overwrite / read
r2=sam_record('has_ds','C'+A*29,'T'+A*29,extra='DS:Z:-99')
show('existing DS:Z:-99',[r2],['--printDS']); show('existing DS:Z:-99',[r2],['--threshold','0']); show('existing DS:Z:-99',[r2],['--threshold','0','--writesamfield'])
# --first (issues #8/#14)
recs=[sam_record(f'f{i}','C'+A*29,('T' if i%3==0 else 'C')+A*29) for i in range(30)]
show('#8/#14 --first',recs,['--first']); show('--first --CpG',recs,['--first','--CpG'])
# deamination table (issue #15 sanity)
show('--deamination --range 5',recs,['--deamination','--range','5'])
# -n (issue #7) exit behaviour
show('-n 5 printDS',recs,['--printDS','-n','5'])
# eager #324: MD containing Y / n, very short reads
show('MD with Y',[sam_record('y','C'+A*29,'T'+A*29,md='0Y29')],['--printDS'])
for L in (1,2,3,5):
    show(f'short read L={L}',[sam_record(f's{L}','C'+A*(L-1),'T'+A*(L-1))],['--printDS'])
# MD inconsistent with CIGAR length
show('MD shorter than CIGAR',[sam_record('mdshort','C'+A*29,'T'+A*29,md='0C10')],['--printDS'])
show('MD longer than CIGAR',[sam_record('mdlong','C'+A*29,'T'+A*29,md='0C40')],['--printDS'])
# long read
L=5000; show('long read 5000',[sam_record('long','C'+A*(L-1),'T'+A*(L-1),qual='I'*L)],['--printDS'])
# SAM with no header line and empty input
rc,o,e=run([],['--printDS'],header=False); print('== empty input rc',rc,repr(o),repr(e[-200:]))
rc,o,e=run([],['--threshold','3','--header']); print('== header-only input rc',rc,repr(o),repr(e[-200:]))

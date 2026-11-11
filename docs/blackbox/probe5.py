import sys; sys.path.insert(0,'.')
from harness import *
def sc(rec,args=()):
    rc,o,e=run([rec],['--printDS']+list(args))
    return o.split('\t')[-1].strip() if o.strip() else f'(no output rc={rc} err={e.strip()[-120:]})'
ref_by_z={z: llr_pmdtools_ds(False,z,40) for z in range(1,40)}
def which_z(v):
    try: v=float(v)
    except: return v
    best=min(ref_by_z,key=lambda z:abs(ref_by_z[z]-v)); return f'{v:.4f} (= C>T at z={best})' if abs(ref_by_z[best]-v)<1e-9 else f'{v:.4f} (no single-z match)'
A='A'
# soft clip of 3 bases before a C>T at first aligned base
ref='C'+A*29; seq='GGG'+'T'+A*29
r=sam_record('sclip',ref,seq,cigar='3S30M'); print('3S then C>T at aligned pos1:', which_z(sc(r)), '| --noclips:', sc(r,['--noclips']))
# hard clip
r=sam_record('hclip',ref,'T'+A*29,cigar='3H30M'); print('3H then C>T:', which_z(sc(r)))
# soft clip at 3' end with G>A near 3'
ref2=A*29+'G'; seq2=A*29+'A'+'TTT'
r=sam_record('sclip3',ref2,seq2,cigar='30M3S'); print("G>A at last aligned then 3S:", which_z(sc(r)))
# insertion before C site: read = AA ins(TT) ... C at aligned ref pos 5
ref3=A*4+'C'+A*25; seq3=A*2+'TT'+A*2+'T'+A*25
r=sam_record('ins',ref3,seq3,cigar='2M2I28M'); print('2M2I then C>T at read pos7/ref pos5:', which_z(sc(r)), '| --noindels:', sc(r,['--noindels']),'| --noinsertions:',sc(r,['--noinsertions']))
# deletion before C site: ref has 2 extra bases deleted
ref4=A*2+'GG'+A*2+'C'+A*25; seq4=A*4+'T'+A*25
r=sam_record('del',ref4,seq4,cigar='2M2D28M'); print('2M2D then C>T at read pos5/ref pos7:', which_z(sc(r)), '| --nodeletions:', sc(r,['--nodeletions']))
# paired flags
for fl in (1,1+64,1+128,1+2+64,1+16+128,256,2048,512,1024):
    r=sam_record(f'f{fl}','C'+A*29,'T'+A*29,flag=fl); print('flag',fl,':',which_z(sc(r)))

import sys; sys.path.insert(0,'.')
from harness import *
A='A'
def sc(rec,args):
    rc,o,e=run([rec],['--printDS']+args); return float(o.split('\t')[-1]) if o.strip() else None
def find(v, Lseq):
    for z in range(1,Lseq+3):
        for y in range(1,Lseq+3):
            if abs(llr_pmdtools_ss(False,z,y,40)-v)<1e-9: return (z,y)
# C>T at last aligned base, then 3S: SEQ length 33
r=sam_record('a',A*29+'C',A*29+'T'+'GGG',cigar='30M3S'); v=sc(r,['--ss']); print('30M3S, C>T aligned pos30:',v,'(z,y)=',find(v,33))
r=sam_record('b','C'+A*29,'GGG'+'T'+A*29,cigar='3S30M'); v=sc(r,['--ss']); print('3S30M, C>T aligned pos1:',v,'(z,y)=',find(v,33))
ref='A'*4+'C'+'A'*25; seq='A'*2+'TT'+'A'*2+'T'+'A'*25
r=sam_record('c',ref,seq,cigar='2M2I28M'); v=sc(r,['--ss']); print('2M2I28M C>T at SEQ pos7 (len32):',v,'(z,y)=',find(v,32))
ref='A'*2+'GG'+'A'*2+'C'+'A'*25; seq='A'*4+'T'+'A'*25
r=sam_record('d',ref,seq,cigar='2M2D28M'); v=sc(r,['--ss']); print('2M2D28M C>T at SEQ pos5 (len30):',v,'(z,y)=',find(v,30))

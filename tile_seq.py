from __future__ import print_function
from string import maketrans
import sys


fasta_file = sys.argv[1] # Input fasta file
chrom = sys.argv[2] # Input interesting sequence IDs, one per line
start_pos = int(sys.argv[3])-1 # Start position (0-indexed)
stop_pos = int(sys.argv[4])-1 # Stop position (0-indexed)
window = int(sys.argv[5]) # window movement length
oligo_size = int(sys.argv[6]) # Oligo size

start_pos=start_pos-oligo_size

def read_fasta(fp, chrom):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name and name == chrom: yield (name, ''.join(seq))
            line=line.split(" ")[0]
            line=line.replace('>','')    
            name, seq = line, []
            print(name, file=sys.stderr)
        elif name == chrom:
            seq.append(line)
        else:
            pass
    if name and name == chrom: yield (name, ''.join(seq))

def revcomp(seq):
    return seq.translate(maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]

with open(fasta_file) as fp:
    for name, seq in read_fasta(fp, chrom):
        if name == chrom:
            seq_length=stop_pos-start_pos
            split_seq=list(seq)
            for i in range(0, seq_length, window):
                oligo_id=name+":"+str(start_pos+i+1)+"-"+str(start_pos+i+oligo_size)
                print(oligo_id, seq[start_pos+i:start_pos+i+oligo_size])
                seq_rc= revcomp(seq[start_pos+i:start_pos+i+oligo_size])
                print(oligo_id+"_RC", seq_rc)MLG-BH0047:oligo_design tewher$ 
MLG-BH0047:oligo_design tewher$ 
MLG-BH0047:oligo_design tewher$ cat tile_seq.py
from __future__ import print_function
from string import maketrans
import sys


fasta_file = sys.argv[1] # Input fasta file
chrom = sys.argv[2] # Input interesting sequence IDs, one per line
start_pos = int(sys.argv[3])-1 # Start position (0-indexed)
stop_pos = int(sys.argv[4])-1 # Stop position (0-indexed)
window = int(sys.argv[5]) # window movement length
oligo_size = int(sys.argv[6]) # Oligo size

start_pos=start_pos-oligo_size

def read_fasta(fp, chrom):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name and name == chrom: yield (name, ''.join(seq))
            line=line.split(" ")[0]
            line=line.replace('>','')    
            name, seq = line, []
            print(name, file=sys.stderr)
        elif name == chrom:
            seq.append(line)
        else:
            pass
    if name and name == chrom: yield (name, ''.join(seq))

def revcomp(seq):
    return seq.translate(maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]

with open(fasta_file) as fp:
    for name, seq in read_fasta(fp, chrom):
        if name == chrom:
            seq_length=stop_pos-start_pos
            split_seq=list(seq)
            for i in range(0, seq_length, window):
                oligo_id=name+":"+str(start_pos+i+1)+"-"+str(start_pos+i+oligo_size)
                print(oligo_id, seq[start_pos+i:start_pos+i+oligo_size])
                seq_rc= revcomp(seq[start_pos+i:start_pos+i+oligo_size])
                print(oligo_id+"_RC", seq_rc)

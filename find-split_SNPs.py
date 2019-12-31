from __future__ import print_function
import sys
import re

SNP_list = sys.argv[1] # Internal bed-like file with SNPs
chrom = str(sys.argv[2]) # Chromosome
start_pos = int(sys.argv[3])-1 # Start position (0-indexed)
stop_pos = int(sys.argv[4])-1 # Stop position (0-indexed)


print("Finding SNPs at "+str(chrom)+":"+str(start_pos)+"-"+str(stop_pos), file=sys.stderr)

with open(SNP_list) as fp:
    for line in fp:
        line = line.rstrip()        
        line=line.split("\t")
        
        if str(line[0]) == chrom and int(line[1]) >= start_pos and int(line[1]) <= stop_pos:
            snps=line[4].split(",")
            for i in snps:
                if(len(re.sub("[ACGTacgt]","",i)) > 0):
                    print("REMOVED: "+line[0]+"\t"+line[1]+"\t"+line[0]+":"+line[1]+":"+line[3]+":"+i+"\t"+line[3]+"\t"+i, file=sys.stderr)
                else:
                    print(line[0]+"\t"+line[1]+"\t"+line[0]+":"+line[1]+":"+line[3]+":"+i+"\t"+line[3]+"\t"+i)

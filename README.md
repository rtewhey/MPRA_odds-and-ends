# MPRA_odds-and-ends


**Create_probes.pl**
Create oligo sequences from SNP list.

INPUT SNP File (Strand = Strand for design): 
  Chr	Pos	ID	Ref	Alt Strand
  5	98073865	rs190431736	T	A	+
  5	98206082	rs145364999	T	A	+

perl Create_probes.pl -L 100 -R 100 -O 200 -P wC  Genome_Fasta_File.fa SNP_Input.txt > Output_File.txt

  -L Left sequence flank
  -R Right sequence flank
  -O Oligo Length
  -P ID postfix

**Add_adapters.pl**
Adds PCR handles to oligo sequences

perl Add_adapters.pl Create_probes_Output.txt > Oligo.out

**replace_restriciton.pl**
Optional, checks for AsiSI restriction site and alters sequence if a cut site is found.

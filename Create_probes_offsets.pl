#!/usr/bin/perl
#
#
###########
# Build Probes and introduce alleles
###########

use strict;
use warnings;
use Bio::SeqIO;
use POSIX;
use Getopt::Std;

sub reverse_complement_IUPAC;
sub combinations;
sub checkOligoOverlap;
sub checkSNPOverlap;

my %options=();
getopts('L:R:O:P:CHF', \%options);

#####
#
#-L = Length of left arm
#-R = Length of right arm
#-O = Total oligo length
#-P = Prefix to attach to ID
#-C = new seq ID used for Encode (chr:pos:ref:alt:[R|A]:window:metaInfo)
#-H = If 1 skip haplotypes (default is 0)
#-F = Variable flanking regions. Left flank length in col 7. Cannot currently work with haplotypes. Must use -H flag.
### When needed. Use $chr_rs as search for rs IDs in combination. Switch the position/allele pull for the adjacent SNP to the RS search while maintaining the internal ID for the primary variant
####
if(exists($options{L}) && exists($options{R})){die "Set length to equal left and right (-O)\n" unless($options{O}==($options{L}+$options{R}));}

my $ID_prefix = $options{P} || "";
my $c_flag = $options{C} || 0;
my $h_flag = $options{H} || 0;
my $f_flag = $options{F} || 0;

print STDERR "Using chr:pos:ref:alt:[R|A]:window:metaInfo naming\n" if($c_flag == 1);
print STDERR "Skipping haplotype building\n" if($h_flag == 1);
print STDERR "Using arm lengths from file\n" if($f_flag == 1);


my $REF = $ARGV[0];
my $ALLELES = $ARGV[1];

#######
## Global parameters
my $oligo_length = $options{O} || 200; #Length of cloned enhancer region
my $max_indel = 50; #Max length of variant 

my $left_length = $options{L} || $oligo_length/2;
my $right_length = $options{R} || $oligo_length/2;

print STDERR "\nRef:$REF\nAllele File: $ALLELES\n\nOligo Length: $oligo_length\nLeft Length: $left_length\nRight Length: $right_length\n\n";
#######


my @inline;
my $internal_ID; #SNPID:left_flank
my %allele_A; #Chr -> ID -> Base
my %allele_B;
my %SNPID;
my %pos;
my %pos_rs;
my %end_pos;
my %end_pos_rs;
my %chr;
my %chr_rs;
my %pos_strand;
my %neg_strand;
my $tmp_base;
my %OFFSET_ID;

my %var_left_length;
my %var_right_length;

#######################
###Load Allele File
#######################

open (ALLELES, $ALLELES) or die "couldn't open $ALLELES\n";
while (<ALLELES>)
{
	next if($_ =~ m/^#/ || $_ =~ m/^Chr/ || $_ =~ m/^chr/ || $_ =~ m/^CHR/);
    $_ =~ s/[\n\r]//g;
    @inline = split(/\t/);
	
	$left_length=$inline[6] if($f_flag == 1);
	$right_length=$oligo_length-$inline[6] if($f_flag == 1);
	
    $internal_ID=$inline[2].":wP".$left_length;
    $SNPID{$internal_ID}=$inline[2];
    $OFFSET_ID{$internal_ID}=":wP".$left_length;
    print STDERR "Skipping $inline[2] for being too large - ".length($inline[3])." - ".length($inline[4])."\n"  if($max_indel <= length($inline[3]) || $max_indel <= length($inline[4]));
    next if($max_indel <= length($inline[3]) || $max_indel <= length($inline[4]));
    $allele_A{$inline[0]}{$internal_ID}=$inline[3];
    $allele_B{$inline[0]}{$internal_ID}=$inline[4];
    die "ERROR: SNP IDs must be unique, duplicate found: internal_ID\n" if(exists $pos{$internal_ID});
	$pos{$internal_ID} = $inline[1];
	$pos_rs{$inline[2]} = $inline[1];

	$end_pos{$internal_ID} = $inline[1]+length($inline[3])-1;
	$end_pos_rs{$inline[2]} = $inline[1]+length($inline[3])-1;

	$chr{$internal_ID} = $inline[0];
	$chr_rs{$inline[2]} = $inline[0];

	$pos_strand{$internal_ID} = 0 unless(exists($pos_strand{$internal_ID}));
	$neg_strand{$internal_ID} = 0 unless(exists($neg_strand{$internal_ID}));
	$pos_strand{$internal_ID} = 1 if($inline[5] eq "+");
	$neg_strand{$internal_ID} = 1 if($inline[5] eq "-");


	$var_left_length{$internal_ID} = $left_length;
	$var_right_length{$internal_ID} = $right_length;

	}
close ALLELES;

#######################
###Verify reference alleles. Switch ones that don't match.
#######################
my $rs;
my $chr;
my $count;
my $tmp_value;

print STDERR "\n\n######\nFinding Reference Allele\n";
###Find Reference call
my $in  = Bio::SeqIO->new(-file => $REF , '-format' => 'Fasta');
while ( my $seq = $in->next_seq() ) 
{
    print STDERR "Chromosome: ",$seq->id," length ",$seq->length(),"\n";
    $chr = $seq->id;
    $count = keys %{$allele_A{$chr}};
    print STDERR "\t$count variants to check\n";
    
    foreach $rs (keys %{$allele_A{$chr}})
    {
    		if($allele_A{$chr}{$rs} eq $seq->subseq($pos{$rs},$end_pos{$rs}))
    			{
    			}
    		elsif($allele_A{$chr}{$rs} ne $seq->subseq($pos{$rs},$end_pos{$rs}))
    			{
   				 $tmp_value=$allele_A{$chr}{$rs};
    			 $allele_A{$chr}{$rs}=$allele_B{$chr}{$rs};
    			 $allele_B{$chr}{$rs}=$tmp_value;
				 $end_pos{$rs} = $pos{$rs}+length($allele_A{$chr}{$rs})-1;
				 $tmp_base=$seq->subseq($pos{$rs},$end_pos{$rs});
				 print STDERR "Flipped alleles for $chr,$pos{$rs},$rs: $allele_B{$chr}{$rs} -> $allele_A{$chr}{$rs} | REF $tmp_base\n";
    			}
    		$tmp_base=$seq->subseq($pos{$rs},$end_pos{$rs});
			print STDERR  "Non matching alleles\n$chr,$pos{$rs},$tmp_base\t$allele_A{$chr}{$rs}\n" if($seq->subseq($pos{$rs},$end_pos{$rs}) ne $allele_A{$chr}{$rs});
	}
}

#######################
###Find Alleles that will sit adjacent on the same oligo
#######################

my %adjacent;
my %combinations;
my %adjacent_ct;
my %adjacent_sort_key;
my $rs_2;
my $rs_3;
my $overlap_oligo;
my $overlap_snp;
my %combinations_ct;

my $i;
my $tmp1;
my $tmp2;

foreach $chr (keys %allele_A)
{
	foreach $rs (keys %{$allele_A{$chr}})
	{
		$combinations_ct{$rs}=1; ## Used later on to build unique ID names
		
		@{$adjacent{$rs}} = ();
		push(@{$adjacent{$rs}}, $rs);
		$adjacent_ct{$rs} = 1;
		
		if($h_flag != 1)  #Skip looking for pairs if flag is set to 1
		{
			foreach $rs_2 (keys %{$allele_A{$chr}})
			{
			#Find SNPs that overlap on the oligo but do not overlap with the centered SNP
			$overlap_oligo = checkOligoOverlap($rs, $rs_2);
			$overlap_snp = checkSNPOverlap($rs, $rs_2);		
			push(@{$adjacent{$rs}}, $rs_2) if($overlap_oligo == 1 && $overlap_snp == 0);
			$adjacent_ct{$rs}++ if($overlap_oligo == 1);	
			}
		}
		
		#Find all combinations of SNPs. Subroutine will produce one empty array that is signifies centered allele A.
		push(@{$combinations{$rs}},[@$_]) for combinations(@{$adjacent{$rs}});

		#Walk through all the combinations and make sure none of the adjacent SNPs overlap each other. It they do toss that combination.
		for($i=0;$i<scalar(@{$combinations{$rs}});$i++)
			{
			#print join("-",$tmp2)."\n";
			$overlap_snp=0;
			if(scalar(@{$combinations{$rs}[$i]}) > 1)
				{
				foreach $rs_2 (@{$combinations{$rs}[$i]})
					{
					foreach $rs_3 (@{$combinations{$rs}[$i]})
						{
						$overlap_snp++ if(checkSNPOverlap($rs_2, $rs_3) ==  1);
						}
					}
				}
			#Code to toss overlapping combinations followed by subtraction of the counter by one because the array is now one element shorter.
 			#print join("\t","Alt Config with overlapping SNPS (tossed): ",join(",",@{$combinations{$rs}}[$i])) if($overlap_snp > 0);		
        	splice(@{$combinations{$rs}},$i,1) if($overlap_snp > 0);
			$i-- if($overlap_snp > 0);
			}	
	}
}
print STDERR "DONE\n";

my %tmp_ct;
my %tmp_ct2;

foreach $rs (keys %adjacent_ct)
{
	$tmp_ct{$adjacent_ct{$rs}}++;
	$tmp_ct2{scalar(@{$combinations{$rs}})}++;
}

print STDERR "\n\n######\n\nMulti Allelic Probes (before pruning)\n";
my $key;
for $key ( sort {$a<=>$b} keys %tmp_ct) 
{
	print STDERR "($key)->($tmp_ct{$key})\n";
}

print STDERR "\n\n######\n\nProbe Combinations (post pruning)\n";
for $key ( sort {$a<=>$b} keys %tmp_ct2) 
{
	print STDERR "($key)->($tmp_ct2{$key})\n";
}

##################
##################
##################
#######################
###Process genome and print probes
#######################


my $oligo_left;
my $oligo_allele;
my $oligo_right;


my $l_start; 
my $l_stop;
my $r_start;
my $r_stop;
my $flank_size;
my $length;
my $spacer_R;
my $spacer_A;
my $l_seq;
my $r_seq;
my $a_ct;
my $c_ct;
my $g_ct;
my $t_ct;
my $gc;
my $tmp_seq;
my $tmp_allele;

my $tmp_rs;
my @tmp_l_flank;
my @tmp_r_flank;
my $tmp_array_pos;
my $alt_l_seq;
my $alt_r_seq;
my @alt_rs;

my $tmp_rev_5;
my $tmp_rev_3;
my $tmp_rev_allele;

my $combo_id_ct;
my $tmp_comb_id;

print STDERR "\n\n######\nBuilding Probes\n";

$in  = Bio::SeqIO->new(-file => $REF , '-format' => 'Fasta');
while ( my $seq = $in->next_seq() ) 
{
    print STDERR "Chromosome: ",$seq->id," length ",$seq->length(),"\n";
    $chr = $seq->id;
    $count = keys %{$allele_A{$chr}};
    print STDERR "\t$count Enhancer oligos to build\n";   
        
      foreach $rs (keys %{$allele_A{$chr}})
    	{
    	for($i=0;$i<scalar(@{$combinations{$rs}});$i++)
    		{
    		$spacer_R = "";
			$spacer_A = "";
			$length = ($end_pos{$rs}-$pos{$rs}+1);
			$length = $length+(length($allele_B{$chr}{$rs})-length($allele_A{$chr}{$rs})) if(length($allele_B{$chr}{$rs}) > length($allele_A{$chr}{$rs}));  ##Adjust for insertions
			#$spacer_R = '-' x (length($allele_B{$chr}{$rs})-length($allele_A{$chr}{$rs})) if(length($allele_B{$chr}{$rs}) > length($allele_A{$chr}{$rs}));
			#$spacer_A = '-' x (length($allele_A{$chr}{$rs})-length($allele_B{$chr}{$rs})) if(length($allele_A{$chr}{$rs}) > length($allele_B{$chr}{$rs}));

			print STDERR "Skipping $rs for being too large - $length\n" if($max_indel <= $length);
			next if($max_indel <= $length);
    		
			$l_start = $pos{$rs}-(floor($var_left_length{$rs}-($length/2)));
			$l_stop = $pos{$rs}-1;
			$r_stop = $end_pos{$rs}+(ceil($var_right_length{$rs}-($length/2)));
			$r_start = $end_pos{$rs}+1;
    	
			my $l_seq = $seq->subseq($l_start,$l_stop);
			my $r_seq = $seq->subseq($r_start,$r_stop);
    		my $l_seq_updated = "";
    		my $r_seq_updated = "";
    		    	
    		my $tmp_rs;
    		my @flanking_rs = ();
    		my $change_middle=0;
    		
    		foreach $tmp_rs (@{$combinations{$rs}[$i]})
    			{
    			$change_middle = 1 if($tmp_rs eq $rs);
    			push(@flanking_rs, $tmp_rs) if($tmp_rs ne $rs);
    			}
    			
    		$l_seq_updated = introduceMutations($chr,$l_start,$l_stop,$l_seq,@flanking_rs);
    		$r_seq_updated = introduceMutations($chr,$r_start,$r_stop,$r_seq,@flanking_rs); 	
    		
    		$l_seq_updated=substr($l_seq_updated, -(floor($var_left_length{$rs}-($length/2))));
    		$r_seq_updated=substr($r_seq_updated, -(ceil($var_right_length{$rs}-($length/2))));


    		my $alleleToPrint = "";
    		$alleleToPrint=$allele_A{$chr}{$rs} if($change_middle == 0);
			$alleleToPrint=$allele_B{$chr}{$rs} if($change_middle == 1);
    		
    		my $id = "";
    		if(scalar(@flanking_rs) >= 1)
    			{
    			#If there are flanking SNPs first see if that combination has been designed with alt middle allele. If true take the prior alt ID count. If it's the first time seen assign a unique id and log it.
    			$tmp_comb_id = join(",",(sort @flanking_rs));
    			if(exists($adjacent_sort_key{$rs}{$tmp_comb_id}))
    				{
    				$combo_id_ct=$adjacent_sort_key{$rs}{$tmp_comb_id};
    				}
    			else
    				{
    				$adjacent_sort_key{$rs}{$tmp_comb_id}=$combinations_ct{$rs};
    				$combo_id_ct=$adjacent_sort_key{$rs}{$tmp_comb_id};
    				$combinations_ct{$rs}++;
    				}
    			if($change_middle == 0 && $ID_prefix ne "") {$id = $rs."_A_alt-".$combo_id_ct."_".$ID_prefix;} 
    			if($change_middle == 0 && $ID_prefix eq "") {$id = $rs."_A_alt-".$combo_id_ct;}
    			if($change_middle == 1 && $ID_prefix ne "") {$id = $rs."_B_alt-".$combo_id_ct."_".$ID_prefix;}
    			if($change_middle == 1 && $ID_prefix eq "")  {$id = $rs."_B_alt-".$combo_id_ct;}
    			
    			if($change_middle == 0 && $ID_prefix ne "" && $c_flag == 1) {$id = $rs.":R:".$ID_prefix.":alt-".$combo_id_ct;} 
    			if($change_middle == 0 && $ID_prefix eq "" && $c_flag == 1) {$id = $rs.":R:NA:alt-".$combo_id_ct;}
    			if($change_middle == 1 && $ID_prefix ne "" && $c_flag == 1) {$id = $rs.":A:".$ID_prefix.":alt-".$combo_id_ct;}
    			if($change_middle == 1 && $ID_prefix eq "" && $c_flag == 1)  {$id = $rs.":A:NA:alt-".$combo_id_ct;}
    			}
    		else
    			{
    			if($change_middle == 0 && $ID_prefix ne "") {$id = $rs."_A_".$ID_prefix;} 
    			if($change_middle == 0 && $ID_prefix eq "") {$id = $rs."_A";}
    			if($change_middle == 1 && $ID_prefix ne "") {$id = $rs."_B_".$ID_prefix;}
    			if($change_middle == 1 && $ID_prefix eq "") {$id = $rs."_B";}
    			
    			if($change_middle == 0 && $ID_prefix ne "" && $c_flag == 1) {$id = $rs.":R:".$ID_prefix;} 
    			if($change_middle == 0 && $ID_prefix eq "" && $c_flag == 1) {$id = $rs.":R:";}
    			if($change_middle == 1 && $ID_prefix ne "" && $c_flag == 1) {$id = $rs.":A:".$ID_prefix;}
    			if($change_middle == 1 && $ID_prefix eq "" && $c_flag == 1) {$id = $rs.":A:";}
				if($change_middle == 0 && $ID_prefix ne "" && $c_flag == 1 && $f_flag == 1) {$id = $SNPID{$rs}.":R:".$OFFSET_ID{$rs}.":".$ID_prefix;} 
    			if($change_middle == 1 && $ID_prefix ne "" && $c_flag == 1 && $f_flag == 1) {$id = $SNPID{$rs}.":A:".$OFFSET_ID{$rs}.":".$ID_prefix;}
				if($change_middle == 0 && $ID_prefix eq "" && $c_flag == 1 && $f_flag == 1) {$id = $SNPID{$rs}.":R:".$OFFSET_ID{$rs};} 
    			if($change_middle == 1 && $ID_prefix eq "" && $c_flag == 1 && $f_flag == 1) {$id = $SNPID{$rs}.":A:".$OFFSET_ID{$rs};}

    			}
    		
    		my $list_ofAlts="-";
    		$list_ofAlts=join(",",@flanking_rs) if(scalar(@flanking_rs) >= 1);
    		
    		$tmp_allele = $allele_A{$chr}{$rs} if($change_middle == 0);
			$tmp_allele = $allele_B{$chr}{$rs} if($change_middle == 1);	
    		
    		if($pos_strand{$rs} == 1)
				{
				$tmp_seq = $l_seq_updated.$alleleToPrint.$r_seq_updated;
				


				$a_ct=0;$c_ct=0;$g_ct=0;$t_ct=0;
				$a_ct = uc($tmp_seq) =~ tr/A//;
				$c_ct = uc($tmp_seq) =~ tr/C//;
				$g_ct = uc($tmp_seq) =~ tr/G//;
				$t_ct = uc($tmp_seq) =~ tr/T//;
				$gc	= sprintf("%.1f",100*(($c_ct+$g_ct)/($a_ct+$c_ct+$g_ct+$t_ct)));	
				print  join("\t",$id,$SNPID{$rs},$chr{$rs},$pos{$rs},$gc,length($tmp_seq),$length,$allele_A{$chr}{$rs},$allele_B{$chr}{$rs},$tmp_allele); 
				print  "\t",$l_seq_updated,"\t",$alleleToPrint,"\t",$r_seq_updated;
				print  "\t",$l_seq_updated,$alleleToPrint,$r_seq_updated,"\t",$list_ofAlts,"\n";	
				}
				
			if($neg_strand{$rs} == 1)
				{
				$tmp_rev_5 = reverse_complement_IUPAC($l_seq_updated);
				$tmp_rev_allele = reverse_complement_IUPAC($alleleToPrint);
				$tmp_rev_3 = reverse_complement_IUPAC($r_seq_updated);
				$tmp_seq = $tmp_rev_5.$tmp_rev_allele.$tmp_rev_3;

				$a_ct=0;$c_ct=0;$g_ct=0;$t_ct=0;
				$a_ct = uc($tmp_seq) =~ tr/A//;
				$c_ct = uc($tmp_seq) =~ tr/C//;
				$g_ct = uc($tmp_seq) =~ tr/G//;
				$t_ct = uc($tmp_seq) =~ tr/T//;
				$gc	= sprintf("%.1f",100*(($c_ct+$g_ct)/($a_ct+$c_ct+$g_ct+$t_ct)));	
				print  join("\t",$id."_RC",$SNPID{$rs},$chr{$rs},$pos{$rs},$gc,length($tmp_seq),$length,$allele_A{$chr}{$rs},$allele_B{$chr}{$rs},$tmp_allele); 
				print  "\t",$tmp_rev_3,"\t",$tmp_rev_allele,"\t",$tmp_rev_5;
				print  "\t",$tmp_rev_3,$tmp_rev_allele,$tmp_rev_5,"\t",$list_ofAlts,"\n";
				}
    		}
		}
}

    
    
    
sub introduceMutations {
    my ($chr,$start, $stop, $seq, @snps) = @_;
    
    my $rs;
    my $i;
    my @snps_toAdd;
    my @inserts_toAdd;
    my @deletions_toAdd;
    my @MNV_toAdd;
    
    foreach $rs (@snps)
    	{
    		push(@snps_toAdd, $rs) if($pos{$rs} >= $start && $end_pos{$rs} <= $stop && length($allele_A{$chr}{$rs}) == 1 && length($allele_B{$chr}{$rs}) == 1);  
    		push(@inserts_toAdd, $rs) if($pos{$rs} >= $start && $end_pos{$rs} <= $stop && length($allele_A{$chr}{$rs}) == 1 && length($allele_B{$chr}{$rs}) > 1);    		    	
    		push(@deletions_toAdd, $rs) if($pos{$rs} >= $start && $end_pos{$rs} <= $stop && length($allele_A{$chr}{$rs}) > 1 && length($allele_B{$chr}{$rs}) == 1);    		    	
  		    print STDERR "MNV deletion - spotcheck $rs - $chr:$start-$stop\n\t".join(",",@snps)."\n" if($pos{$rs} >= $start && $end_pos{$rs} <= $stop && length($allele_A{$chr}{$rs}) > 1 && length($allele_B{$chr}{$rs}) > 1 && length($allele_A{$chr}{$rs}) > length($allele_B{$chr}{$rs}));
  		  	print STDERR "MNV insertion - spotcheck $rs - $chr:$start-$stop\n\t".join(",",@snps)."\n" if($pos{$rs} >= $start && $end_pos{$rs} <= $stop && length($allele_A{$chr}{$rs}) > 1 && length($allele_B{$chr}{$rs}) > 1 && length($allele_A{$chr}{$rs}) < length($allele_B{$chr}{$rs}));
   		    push(@deletions_toAdd, $rs) if($pos{$rs} >= $start && $end_pos{$rs} <= $stop && length($allele_A{$chr}{$rs}) > 1 && length($allele_B{$chr}{$rs}) > 1 && length($allele_A{$chr}{$rs}) > length($allele_B{$chr}{$rs}));
  		  	push(@inserts_toAdd, $rs) if($pos{$rs} >= $start && $end_pos{$rs} <= $stop && length($allele_A{$chr}{$rs}) > 1 && length($allele_B{$chr}{$rs}) > 1 && length($allele_A{$chr}{$rs}) < length($allele_B{$chr}{$rs}));
  		    print STDERR "Equal Length MNV insertion - spotcheck $rs - $chr:$start-$stop\n\t".join(",",@snps)."\n" if($pos{$rs} >= $start && $end_pos{$rs} <= $stop && length($allele_A{$chr}{$rs}) > 1 && length($allele_B{$chr}{$rs}) > 1 && length($allele_A{$chr}{$rs}) == length($allele_B{$chr}{$rs}));
            push(@MNV_toAdd, $rs) if($pos{$rs} >= $start && $end_pos{$rs} <= $stop && length($allele_A{$chr}{$rs}) > 1 && length($allele_B{$chr}{$rs}) > 1 && length($allele_A{$chr}{$rs}) == length($allele_B{$chr}{$rs}));
  		    die "One of the alleles is empty?: $rs\n" if($pos{$rs} >= $start && $end_pos{$rs} <= $stop && (length($allele_A{$chr}{$rs}) < 1 || length($allele_B{$chr}{$rs}) < 1))	

    	}
    
    my @SeqBases = split(//,$seq);
    
    foreach $rs (@snps_toAdd)
    	{
    	die "does not match: $rs $allele_A{$chr}{$rs} $SeqBases[$pos{$rs}-$start] \n" if($allele_A{$chr}{$rs} ne $SeqBases[$pos{$rs}-$start]);
    	$SeqBases[$pos{$rs}-$start] = $allele_B{$chr}{$rs};
    	}
    my $tmp_seq="";
    foreach $rs (@deletions_toAdd)
    	{
    	$tmp_seq="";
    		for($i=$pos{$rs};$i<=$end_pos{$rs};$i++)
    			{
    			$tmp_seq=$tmp_seq.$SeqBases[$i-$start];
    			$SeqBases[$i-$start]="";
    			}
    	die "does not match: $rs $allele_A{$chr}{$rs} $tmp_seq\n" if($allele_A{$chr}{$rs} ne $tmp_seq);
    	$SeqBases[$pos{$rs}-$start] = $allele_B{$chr}{$rs};
    	}
	foreach $rs (@inserts_toAdd)
		{
		$tmp_seq="";
			for($i=$pos{$rs};$i<=$end_pos{$rs};$i++)
				{
    			$tmp_seq=$tmp_seq.$SeqBases[$i-$start];
    			$SeqBases[$i-$start]="";
    			}
			die "does not match: $rs $allele_A{$chr}{$rs} $tmp_seq\n" if($allele_A{$chr}{$rs} ne $tmp_seq);
			$SeqBases[$pos{$rs}-$start] = $allele_B{$chr}{$rs};
			###Does not account for when SNP is at the -1 position. Coded so that I can easily fix this in the future.	
		}
    foreach $rs (@MNV_toAdd)
        {
        $tmp_seq="";
            for($i=$pos{$rs};$i<=$end_pos{$rs};$i++)
                {
                $tmp_seq=$tmp_seq.$SeqBases[$i-$start];
                $SeqBases[$i-$start]="";
                }
            die "does not match: $rs $allele_A{$chr}{$rs} $tmp_seq\n" if($allele_A{$chr}{$rs} ne $tmp_seq);
            $SeqBases[$pos{$rs}-$start] = $allele_B{$chr}{$rs};
            ###Does not account for when SNP is at the -1 position. Coded so that I can easily fix this in the future.  
        }

	$tmp_seq=join("",@SeqBases);
	return($tmp_seq);
}   
    
sub reverse_complement_IUPAC {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        return $revcomp;
}

sub combinations {
  return [] unless @_;
  my $first = shift;
  my @rest = combinations(@_);
  return @rest, map { [$first, @$_] } @rest;
}

sub checkOligoOverlap (@) {
	my $overlap = 0;
	my $rs_1 = $_[0];
	my $rs_2 = $_[1];
	if($pos{$rs_2} <= (floor(($pos{$rs_1}+$end_pos{$rs_1})*0.5)+($var_right_length{$rs_1})) && $pos{$rs_2} >= (floor(($pos{$rs_1}+$end_pos{$rs_1})*0.5)-($var_left_length{$rs_1})) && $end_pos{$rs_2} <= (floor(($pos{$rs_1}+$end_pos{$rs_1})*0.5)+($var_right_length{$rs_1})) && $end_pos{$rs_2} >= (floor(($pos{$rs_1}+$end_pos{$rs_1})*0.5)-($var_left_length{$rs_1}))  && $rs_1 ne $rs_2 && $chr{$rs_1} eq $chr{$rs_2})
			{
			$overlap=1
			#print STDERR "Two Oligos\t$rs_1\t$pos{$rs_1}\t$rs_2\t$pos{$rs_2}\n"
			}
	return $overlap;
}

sub checkSNPOverlap (@) {
	my $overlap = 0;	
	my $rs_1 = $_[0];
	my $rs_2 = $_[1];
				if((($pos{$rs_2} <= $pos{$rs_1} && $end_pos{$rs_2} >= $pos{$rs_1}) || ($pos{$rs_1} <= $pos{$rs_2} && $end_pos{$rs_1} >= $pos{$rs_2})) && $rs_1 ne $rs_2 && $chr{$rs_1} eq $chr{$rs_2})
				{
				$overlap=1
				#print STDERR "!!Overlapping Sites!!\t$rs_1\t$pos{$rs_1}\t$rs_2\t$pos{$rs_2}\n"
				}
		return $overlap;
}



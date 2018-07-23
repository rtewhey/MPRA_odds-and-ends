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
 
my $REF = $ARGV[0];
my $ALLELES = $ARGV[1];

#######
## Global parameters
my $oligo_length = 150; #Length of cloned enhancer region
my $max_indel = 35; #Max length of variant 

#######


my @inline;
my %allele_A; #Chr -> ID -> Base
my %allele_B;
my %pos;
my %end_pos;
my %chr;
my %pos_strand;
my %neg_strand;



#######################
###Load Allele File
#######################

open (ALLELES, $ALLELES) or die "couldn't open $ALLELES\n";
while (<ALLELES>)
{
    $_ =~ s/[\n\r]//g;
    @inline = split(/\t/);
    $allele_A{$inline[0]}{$inline[2]}=$inline[3];
    $allele_B{$inline[0]}{$inline[2]}=$inline[4];
	$pos{$inline[2]} = $inline[1];
	$end_pos{$inline[2]} = $inline[1]+length($inline[3])-1;
	$chr{$inline[2]} = $inline[0];
	$pos_strand{$inline[2]} = 0 unless(exists($pos_strand{$inline[2]}));
	$neg_strand{$inline[2]} = 0 unless(exists($neg_strand{$inline[2]}));
	$pos_strand{$inline[2]} = 1 if($inline[5] eq "+");
	$neg_strand{$inline[2]} = 1 if($inline[5] eq "-");
}
close ALLELES;



#######################
###Find Alleles that will sit adjacent on the same oligo
#######################

my %adjacent;
my %adjacent_ct;
my $chr;
my $rs;
my $rs_2;

foreach $chr (keys %allele_A)
{
	foreach $rs (keys %{$allele_A{$chr}})
	{
		$adjacent_ct{$rs} = 1;
		foreach $rs_2 (keys %{$allele_A{$chr}})
		{
			if($pos{$rs_2} <= ($pos{$rs}+($oligo_length/2)) && $pos{$rs_2} >= ($pos{$rs}-($oligo_length/2)) && $rs ne $rs_2 && $chr{$rs} eq $chr{$rs_2})
			{
				if(($pos{$rs_2} <= $pos{$rs} && $end_pos{$rs_2} >= $pos{$rs}) || ($pos{$rs} <= $pos{$rs_2} && $end_pos{$rs} >= $pos{$rs_2}))
				{
					#print STDERR "!!Overlapping Sites!!\t$rs\t$pos{$rs}\t$rs_2\t$pos{$rs_2}\n"
				}
				else
				{
					$adjacent_ct{$rs}++;
					push(@{$adjacent{$rs}},$rs_2);
						#print STDERR "\t\t$chr{$rs}\t$rs\t$pos{$rs}\t$chr{$rs_2}\t$rs_2\t$pos{$rs_2}\n"
				}
			}
		}
	}
}

my %tmp_ct;
foreach $rs (keys %adjacent_ct)
{
	$tmp_ct{$adjacent_ct{$rs}}++
}

print STDERR "\n\n######\n\nMulti allelelic probes\n";
my $key;
for $key ( sort {$a<=>$b} keys %tmp_ct) 
{
	print STDERR "($key)->($tmp_ct{$key})\n";
}

#######################
###Process genome and print probes
#######################

my $count;
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


print STDERR "\n\n######\n\nLoading genome\n";
my $in  = Bio::SeqIO->new(-file => $REF , '-format' => 'Fasta');

while ( my $seq = $in->next_seq() ) 
{
    print STDERR "Chromosome: ",$seq->id," length ",$seq->length(),"\n";
    $chr = $seq->id;
    $count = keys %{$allele_A{$chr}};
    print STDERR "\t$count Enhancer oligos to build\n";
    
    foreach $rs (keys %{$allele_A{$chr}})
    {
			$spacer_R = "";
			$spacer_A = "";
			$length = ($end_pos{$rs}-$pos{$rs}+1);
			$length = $length+(length($allele_B{$chr}{$rs})-length($allele_A{$chr}{$rs})) if(length($allele_B{$chr}{$rs}) > length($allele_A{$chr}{$rs}));  ##Adjust for insertions
			$spacer_R = '-' x (length($allele_B{$chr}{$rs})-length($allele_A{$chr}{$rs})) if(length($allele_B{$chr}{$rs}) > length($allele_A{$chr}{$rs}));
			$spacer_A = '-' x (length($allele_A{$chr}{$rs})-length($allele_B{$chr}{$rs})) if(length($allele_A{$chr}{$rs}) > length($allele_B{$chr}{$rs}));
			
			print STDERR "Skipping $rs for being too large - $length\n" if($max_indel <= $length);
			next if($max_indel <= $length);
		
			$l_start = $pos{$rs}-(floor(($oligo_length-$length)/2));
			$l_stop = $pos{$rs}-1;
			$r_stop = $end_pos{$rs}+(ceil(($oligo_length-$length)/2));
			$r_start = $end_pos{$rs}+1;

			my $l_seq = $seq->subseq($l_start,$l_stop);
			my $r_seq = $seq->subseq($r_start,$r_stop);
			
				if($pos_strand{$rs} == 1)
					{
					$tmp_seq = $l_seq.$seq->subseq($pos{$rs},$end_pos{$rs}).$r_seq;
					$a_ct=0;$c_ct=0;$g_ct=0;$t_ct=0;
					$a_ct = uc($tmp_seq) =~ tr/A//;
					$c_ct = uc($tmp_seq) =~ tr/C//;
					$g_ct = uc($tmp_seq) =~ tr/G//;
					$t_ct = uc($tmp_seq) =~ tr/T//;
					$gc	= sprintf("%.1f",100*(($c_ct+$g_ct)/($a_ct+$c_ct+$g_ct+$t_ct)));	
					print  join("\t",$rs."_A",$rs,$chr{$rs},$pos{$rs},$gc,length($tmp_seq),$length,$allele_A{$chr}{$rs},$allele_B{$chr}{$rs},$allele_A{$chr}{$rs}); 
					print  "\t",$l_seq,"\t",$seq->subseq($pos{$rs},$end_pos{$rs}),"$spacer_R\t",$r_seq;
					print  "\t",$l_seq,$seq->subseq($pos{$rs},$end_pos{$rs}),$r_seq,"\n";	
					
					
					$tmp_seq = $l_seq.$allele_B{$chr}{$rs}.$r_seq;
					$a_ct=0;$c_ct=0;$g_ct=0;$t_ct=0;
					$a_ct = uc($tmp_seq) =~ tr/A//;
					$c_ct = uc($tmp_seq) =~ tr/C//;
					$g_ct = uc($tmp_seq) =~ tr/G//;
					$t_ct = uc($tmp_seq) =~ tr/T//;
					$gc	= sprintf("%.1f",100*(($c_ct+$g_ct)/($a_ct+$c_ct+$g_ct+$t_ct)));	
					print  join("\t",$rs."_B",$rs,$chr{$rs},$pos{$rs},$gc,length($tmp_seq),$length,$allele_A{$chr}{$rs},$allele_B{$chr}{$rs},$allele_B{$chr}{$rs}); 
					print  "\t",$l_seq,"\t",$allele_B{$chr}{$rs},"$spacer_R\t",$r_seq;
					print  "\t",$l_seq,$allele_B{$chr}{$rs},$r_seq,"\n";
					}
		
				if($neg_strand{$rs} == 1)
					{
					$tmp_seq = $l_seq.$seq->subseq($pos{$rs},$end_pos{$rs}).$r_seq;
					$a_ct=0;$c_ct=0;$g_ct=0;$t_ct=0;
					$a_ct = uc($tmp_seq) =~ tr/A//;
					$c_ct = uc($tmp_seq) =~ tr/C//;
					$g_ct = uc($tmp_seq) =~ tr/G//;
					$t_ct = uc($tmp_seq) =~ tr/T//;
					
					$tmp_rev_5 = reverse_complement_IUPAC($l_seq);
					$tmp_rev_allele = reverse_complement_IUPAC($seq->subseq($pos{$rs},$end_pos{$rs}));
					$tmp_rev_3 = reverse_complement_IUPAC($r_seq);

					$gc	= sprintf("%.1f",100*(($c_ct+$g_ct)/($a_ct+$c_ct+$g_ct+$t_ct)));	
					print  join("\t",$rs."_RC_A",$rs,$chr{$rs},$pos{$rs},$gc,length($tmp_seq),$length,$allele_A{$chr}{$rs},$allele_B{$chr}{$rs},$allele_A{$chr}{$rs}); 
					print  "\t",$tmp_rev_3,"\t",$tmp_rev_allele,"$spacer_R\t",$tmp_rev_5;
					print  "\t",$tmp_rev_3,$tmp_rev_allele,$tmp_rev_5,"\n";	
					
					
					$tmp_seq = $l_seq.$allele_B{$chr}{$rs}.$r_seq;
					$a_ct=0;$c_ct=0;$g_ct=0;$t_ct=0;
					$a_ct = uc($tmp_seq) =~ tr/A//;
					$c_ct = uc($tmp_seq) =~ tr/C//;
					$g_ct = uc($tmp_seq) =~ tr/G//;
					$t_ct = uc($tmp_seq) =~ tr/T//;
					
					$tmp_rev_5 = reverse_complement_IUPAC($l_seq);
					$tmp_rev_allele = reverse_complement_IUPAC($allele_B{$chr}{$rs});
					$tmp_rev_3 = reverse_complement_IUPAC($r_seq);
					
					$gc	= sprintf("%.1f",100*(($c_ct+$g_ct)/($a_ct+$c_ct+$g_ct+$t_ct)));	
					print  join("\t",$rs."_RC_B",$rs,$chr{$rs},$pos{$rs},$gc,length($tmp_seq),$length,$allele_A{$chr}{$rs},$allele_B{$chr}{$rs},$allele_B{$chr}{$rs}); 
					print  "\t",$tmp_rev_3,"\t",$tmp_rev_allele,"$spacer_R\t",$tmp_rev_5;
					print  "\t",$tmp_rev_3,$tmp_rev_allele,$tmp_rev_5,"\n";
					
					}
			#print STDERR "\t\t",$l_seq," ",$allele_B{$chr}{$rs},"$spacer_A ",$r_seq,"\n";

			die "Non matching alleles\n$rs\t$seq->subseq($pos{$rs},$end_pos{$rs})\t$allele_A{$chr}{$rs} " if($seq->subseq($pos{$rs},$end_pos{$rs}) ne $allele_A{$chr}{$rs});

		##Build Alt haplotype
		if($adjacent_ct{$rs} > 1)
		{
			$tmp_rs = "";
			@alt_rs = ();
			@tmp_l_flank = split(//,$l_seq);
			@tmp_r_flank = split(//,$r_seq);

			foreach $tmp_rs (@{$adjacent{$rs}})
			{
				
				if(length($allele_A{$chr}{$tmp_rs})==1 && length($allele_B{$chr}{$tmp_rs})==1)
				{
					push(@alt_rs,$tmp_rs);
					if($pos{$tmp_rs} >= $l_start && $pos{$tmp_rs} <= $l_stop)
					{
						$tmp_array_pos = $pos{$tmp_rs}-$l_start;
						die "Non matching alleles during alt haplotype reconstruction\n$rs adding $tmp_rs - $tmp_l_flank[$tmp_array_pos] - $allele_A{$chr}{$tmp_rs} \n " if($tmp_l_flank[$tmp_array_pos] ne $allele_A{$chr}{$tmp_rs});
						$tmp_l_flank[$tmp_array_pos] = $allele_B{$chr}{$tmp_rs};
					}
					elsif($pos{$tmp_rs} >= $r_start && $pos{$tmp_rs} <= $r_stop)
					{
						$tmp_array_pos = $pos{$tmp_rs}-$r_start;
						die "Non matching alleles during alt haplotype reconstruction\n$rs adding $tmp_rs - $tmp_r_flank[$tmp_array_pos] - $allele_A{$chr}{$tmp_rs} \n " if($tmp_r_flank[$tmp_array_pos] ne $allele_A{$chr}{$tmp_rs});
						$tmp_r_flank[$tmp_array_pos] = $allele_B{$chr}{$tmp_rs}	;					
					}			
				}
			}
			
			##Print Alt haplotype enhancer oligos
			if(scalar(@alt_rs) >= 1)
			{
				$alt_l_seq = join('',@tmp_l_flank);
				$alt_r_seq = join('',@tmp_r_flank);
				
				
				if($pos_strand{$rs} == 1)
					{
					$tmp_seq = $alt_l_seq.$seq->subseq($pos{$rs},$end_pos{$rs}).$alt_r_seq;
					$a_ct=0;$c_ct=0;$g_ct=0;$t_ct=0;
					$a_ct = uc($tmp_seq) =~ tr/A//;
					$c_ct = uc($tmp_seq) =~ tr/C//;
					$g_ct = uc($tmp_seq) =~ tr/G//;
					$t_ct = uc($tmp_seq) =~ tr/T//;
					$gc	= sprintf("%.1f",100*(($c_ct+$g_ct)/($a_ct+$c_ct+$g_ct+$t_ct)));	
					print  join("\t",$rs."_altA",$rs,$chr{$rs},$pos{$rs},$gc,length($tmp_seq),$length,$allele_A{$chr}{$rs},$allele_B{$chr}{$rs},$allele_A{$chr}{$rs}); 
					print  "\t",$alt_l_seq,"\t",$seq->subseq($pos{$rs},$end_pos{$rs}),"$spacer_R\t",$alt_r_seq;
					print  "\t",$alt_l_seq,$seq->subseq($pos{$rs},$end_pos{$rs}),$alt_r_seq,"\t";
					print  join (",",@alt_rs)."\n";
					
					$tmp_seq = $alt_l_seq.$allele_B{$chr}{$rs}.$alt_r_seq;
					$a_ct=0;$c_ct=0;$g_ct=0;$t_ct=0;
					$a_ct = uc($tmp_seq) =~ tr/A//;
					$c_ct = uc($tmp_seq) =~ tr/C//;
					$g_ct = uc($tmp_seq) =~ tr/G//;
					$t_ct = uc($tmp_seq) =~ tr/T//;
					$gc	= sprintf("%.1f",100*(($c_ct+$g_ct)/($a_ct+$c_ct+$g_ct+$t_ct)));	
					print  join("\t",$rs."_altB",$rs,$chr{$rs},$pos{$rs},$gc,length($tmp_seq),$length,$allele_A{$chr}{$rs},$allele_B{$chr}{$rs},$allele_B{$chr}{$rs}); 
					print  "\t",$alt_l_seq,"\t",$allele_B{$chr}{$rs},"$spacer_R\t",$alt_r_seq;
					print  "\t",$alt_l_seq,$allele_B{$chr}{$rs},$alt_r_seq,"\t";
					print  join (",",@alt_rs)."\n";
					}
				if($neg_strand{$rs} == 1)
					{				
					$tmp_seq = $alt_l_seq.$seq->subseq($pos{$rs},$end_pos{$rs}).$alt_r_seq;
					$a_ct=0;$c_ct=0;$g_ct=0;$t_ct=0;
					$a_ct = uc($tmp_seq) =~ tr/A//;
					$c_ct = uc($tmp_seq) =~ tr/C//;
					$g_ct = uc($tmp_seq) =~ tr/G//;
					$t_ct = uc($tmp_seq) =~ tr/T//;
					
					$tmp_rev_5 = reverse_complement_IUPAC($alt_l_seq);
					$tmp_rev_allele = reverse_complement_IUPAC($seq->subseq($pos{$rs},$end_pos{$rs}));
					$tmp_rev_3 = reverse_complement_IUPAC($alt_r_seq);
					
					$gc	= sprintf("%.1f",100*(($c_ct+$g_ct)/($a_ct+$c_ct+$g_ct+$t_ct)));	
					print  join("\t",$rs."_RC_altA",$rs,$chr{$rs},$pos{$rs},$gc,length($tmp_seq),$length,$allele_A{$chr}{$rs},$allele_B{$chr}{$rs},$allele_A{$chr}{$rs}); 
					print  "\t",$tmp_rev_3,"\t",$tmp_rev_allele,"$spacer_R\t",$tmp_rev_5;
					print  "\t",$tmp_rev_3,$tmp_rev_allele,$tmp_rev_5,"\t";
					print  join (",",@alt_rs)."\n";
					
					$tmp_seq = $alt_l_seq.$allele_B{$chr}{$rs}.$alt_r_seq;
					$a_ct=0;$c_ct=0;$g_ct=0;$t_ct=0;
					$a_ct = uc($tmp_seq) =~ tr/A//;
					$c_ct = uc($tmp_seq) =~ tr/C//;
					$g_ct = uc($tmp_seq) =~ tr/G//;
					$t_ct = uc($tmp_seq) =~ tr/T//;
					
					$tmp_rev_5 = reverse_complement_IUPAC($alt_l_seq);
					$tmp_rev_allele = reverse_complement_IUPAC($allele_B{$chr}{$rs});
					$tmp_rev_3 = reverse_complement_IUPAC($alt_r_seq);
					
					$gc	= sprintf("%.1f",100*(($c_ct+$g_ct)/($a_ct+$c_ct+$g_ct+$t_ct)));	
					print  join("\t",$rs."_RC_altB",$rs,$chr{$rs},$pos{$rs},$gc,length($tmp_seq),$length,$allele_A{$chr}{$rs},$allele_B{$chr}{$rs},$allele_B{$chr}{$rs}); 
					print  "\t",$tmp_rev_3,"\t",$tmp_rev_allele,"$spacer_R\t",$tmp_rev_5;
					print  "\t",$tmp_rev_3,$tmp_rev_allele,$tmp_rev_5,"\t";
					print  join (",",@alt_rs)."\n";
					}
				}
					
		}    
    }
}




sub reverse_complement_IUPAC {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        return $revcomp;
}


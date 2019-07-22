#!/usr/bin/perl
#
#
###########
# Build Probes and introduce alleles
###########

use strict;
use warnings;
use POSIX;
use Getopt::Std;


my %options=();
getopts('B', \%options);

#####
#
#-B = Use flag -B if BASH format (ID Oligo)
#
#####


my $ALLELES = $ARGV[0];

#######
## Global parameters
my $oligo_length = 230; #Length of cloned enhancer region 

my $l_adapter = "ACTGGCCGCTTGACG";
my $l_pad = "GTACGGGAGGTATTGGACAGGCCGCAATAAAATATCTTTATTTTCATTACATCTGTGTGTTGGTTTTTTGTGTGAATCGATAGTACTAACATACGCTCTCCATCAAAACAAAACGAAACAAAACAAACTAGCAAAATAGGCTGTCCCCAGTGCAAGTGCAGGTGCCAGAACATTTCTCTGGCCTA";
my $r_adapter = "CACTGCGGCTCCTGC";
my $r_pad = "GATCGCGTCGACGAACCTCTAGAAAAAAAAAAAAAAAAAAAAAAGATCGGAAGAGCGTCGGCGGCCAAGCTAGTCGGGGCGGCCGGCCGCTTCGAGCAGACATGATAAGATACATTGATGAGTTTGGACAAACCACAACTAGAATGCAGTGAAAAAAATGCT";




#######


my @inline;
my $length;
my $residual;
my $left_add;
my $right_add;

my $oligo;
my $oligo_rc;
my $enhancer_rc;

my $tmp_seq;
my $a_ct;
my $c_ct;
my $g_ct;
my $t_ct;
my $gc;

#######################
###Load Allele File
#######################

if($options{B}) {
open (ALLELES, $ALLELES) or die "couldn't open $ALLELES\n";
while (<ALLELES>)
{
    $_ =~ s/[\n\r]//g;
    @inline = split(/\s+/);
	$length = length($inline[1]);
	$oligo = $l_adapter.$inline[1].$r_adapter;
	
	die "Enhancer is to large to fit an oligo of $oligo_length bp\n" if(length($oligo) > $oligo_length);
	
	if(length($oligo) < $oligo_length)
		{
		$residual = $oligo_length - length($oligo);
		
		$left_add = floor(($residual)/2);
		$right_add = ceil(($residual)/2);
		
		if($left_add < 1)
			{
			#$oligo_rc=$oligo_rc.substr($r_pad,0,$right_add) ;
			$oligo=$oligo.substr($r_pad,0,$right_add);
			}
		else
			{
			#$oligo_rc=substr($l_pad,-$left_add).$oligo_rc.substr($r_pad,0,$right_add) ;
			$oligo=substr($l_pad,-$left_add).$oligo.substr($r_pad,0,$right_add);	
			}
		}
		
	$tmp_seq = $inline[1];
	$a_ct=0;$c_ct=0;$g_ct=0;$t_ct=0;
	$a_ct = uc($tmp_seq) =~ tr/A//;
	$c_ct = uc($tmp_seq) =~ tr/C//;
	$g_ct = uc($tmp_seq) =~ tr/G//;
	$t_ct = uc($tmp_seq) =~ tr/T//;
	$gc	= sprintf("%.1f",100*(($c_ct+$g_ct)/($a_ct+$c_ct+$g_ct+$t_ct)));	
				
	print join("\t",$inline[0],"NA","NA","NA",$gc,length($inline[1]),length($oligo),uc($oligo),$inline[1])."\n";
	#print join("\t",$inline[0]."_rc",$inline[1],$inline[2],$inline[3],$inline[4],$inline[5],length($oligo_rc),$oligo_rc)."\n";	
}
close ALLELES;
}

else {
open (ALLELES, $ALLELES) or die "couldn't open $ALLELES\n";
while (<ALLELES>)
{
    $_ =~ s/[\n\r]//g;
    @inline = split(/\t/);
	$length = length($inline[13]);
	$oligo = $l_adapter.$inline[13].$r_adapter;
	
	die "Enhancer is to large to fit an oligo of $oligo_length bp\n" if(length($oligo) > $oligo_length);
	
	if(length($oligo) < $oligo_length)
		{
		$residual = $oligo_length - length($oligo);
		
		$left_add = floor(($residual)/2);
		$right_add = ceil(($residual)/2);
		
		if($left_add < 1)
			{
			#$oligo_rc=$oligo_rc.substr($r_pad,0,$right_add) ;
			$oligo=$oligo.substr($r_pad,0,$right_add);
			}
		else
			{
			#$oligo_rc=substr($l_pad,-$left_add).$oligo_rc.substr($r_pad,0,$right_add) ;
			$oligo=substr($l_pad,-$left_add).$oligo.substr($r_pad,0,$right_add);	
			}
		}
	print join("\t",$inline[0],$inline[1],$inline[2],$inline[3],$inline[4],$inline[5],length($oligo),$oligo,$inline[13])."\n";
	#print join("\t",$inline[0]."_rc",$inline[1],$inline[2],$inline[3],$inline[4],$inline[5],length($oligo_rc),$oligo_rc)."\n";	
}
close ALLELES;
}

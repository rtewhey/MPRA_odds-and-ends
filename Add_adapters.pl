#!/usr/bin/perl
#
#
###########
# Build Probes and introduce alleles
###########

use strict;
use warnings;
use POSIX;

my $ALLELES = $ARGV[0];

#######
## Global parameters
my $oligo_length = 230; #Length of cloned enhancer region 


my $l_adapter = "ACTGGCCGCTTGACG";
my $l_pad = "GCCAGAACATTTCTCTGGCCTA";
my $r_adapter = "CACTGCGGCTCCTGC";
my $r_pad = "GATCGCGTCGACGAACCTCTAGA";

#######


my @inline;
my $length;
my $residual;
my $left_add;
my $right_add;

my $oligo;
my $oligo_rc;
my $enhancer_rc;


#######################
###Load Allele File
#######################

open (ALLELES, $ALLELES) or die "couldn't open $ALLELES\n";
while (<ALLELES>)
{
    $_ =~ s/[\n\r]//g;
    @inline = split(/\t/);
	$length = length($inline[13]);
	#$enhancer_rc = reverse($inline[13]);
    #$enhancer_rc =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
	#$oligo_rc = $l_adapter.$enhancer_rc.$r_adapter;
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
	print join("\t",$inline[0],$inline[1],$inline[2],$inline[3],$inline[4],$inline[5],length($oligo),$oligo)."\n";
	#print join("\t",$inline[0]."_rc",$inline[1],$inline[2],$inline[3],$inline[4],$inline[5],length($oligo_rc),$oligo_rc)."\n";

	
	
	
}
close ALLELES;

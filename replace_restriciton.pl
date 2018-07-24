#!/usr/bin/perl
#
#
###########
# Build Probes and introduce alleles
###########

use strict;
use warnings;

my $ALLELES = $ARGV[0];


my $search = "GCGATCGC";
my $replace = "GCGATTGC";


#######################
###Load Allele File
#######################

my @inline;
my $count;
my $variant;
my $left;
my $right;


open (ALLELES, $ALLELES) or die "couldn't open $ALLELES\n";
while (<ALLELES>)
{
    $_ =~ s/[\n\r]//g;
    @inline = split(/\t/);
	if ($inline[13] =~ /$search/)
	{
		$count = 0;
		$count++ while $inline[13] =~ /$search/g;
		print STDERR "Adjusted $inline[0] - $count hits\n";
		
		$left = substr($inline[10],-(length($search)-1));
		$right = substr($inline[12],0,(length($search)-1));
		$variant = $left.$inline[11].$right;
		
		if($variant =~ /$search/)
			{
			print STDERR "WARNING: Restriction site appears to fall over variant! - $inline[0]\n";
			}
		
		$inline[10] =~ s/$search/$replace/g;
		$inline[12] =~ s/$search/$replace/g;
		$inline[13] =~ s/$search/$replace/g;
	}
	
print join("\t",@inline)."\n";			
}

close ALLELES;

#!/usr/bin/perl
#
#
###########
# Build Probes and introduce alleles
###########

use strict;
use warnings;

sub all_match_positions;

my $ALLELES = $ARGV[0];


my $search = "GCGATCGC";
my $replace = "GCGATTGC";

my $adapter_l = "CACTGCGGCTCCT";

#######################
###Load Allele File
#######################

my @inline;
my $count;
my $variant;
my $left;
my $right;
my @loc;
my $length;

my $i;

open (ALLELES, $ALLELES) or die "couldn't open $ALLELES\n";
while (<ALLELES>)
{
    $_ =~ s/[\n\r]//g;
    @inline = split(/\t/);
   	#if ($inline[13] =~ /$search/)
	
	if ($inline[7] =~ /(?<!$adapter_l)$search/)
	{
		$count = 0;
		$count++ while $inline[7] =~ /(?<!$adapter_l)$search/g;
		print STDERR "Adjusted $inline[0] - $count hits\n";
		
		$length = length($inline[7]);
		
		#$left = substr($inline[10],-(length($search)-1));
		#$right = substr($inline[12],0,(length($search)-1));
		#$variant = $left.$inline[11].$right;
		
		@loc = all_match_positions($search, $inline[7]);
		
		for($i=0;$i<scalar(@loc);$i++)
			{
			print STDERR "$loc[$i][0]\t $loc[$i][1]\n";
			if($loc[$i][0] <= ($length/2)+1 && $loc[$i][1] >= ($length/2)-1    )
				{
				print STDERR "WARNING: Restriction site appears to fall over middle of oligo! - $inline[0]\n";
				}
			}
		$inline[7] =~ s/(?<!$adapter_l)$search/$replace/g;
		$inline[8] =~ s/(?<!$adapter_l)$search/$replace/g;

	}
	
print join("\t",@inline)."\n";			
}

close ALLELES;



sub all_match_positions {
    my ($regex, $string) = @_;
    my @ret;
    while ($string =~ /($regex)/g) {
        push @ret, [pos($string)-length $1, pos($string)];
    }
    return(@ret);
    }
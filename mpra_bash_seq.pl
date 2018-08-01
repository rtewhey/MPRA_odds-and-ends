#!/usr/bin/perl

use strict;
use warnings;
use POSIX;

my $prefix = $ARGV[0];
my $size = $ARGV[1];
my $sequence = $ARGV[2]; 

my @SEQ = split(//,$sequence);

my $ct = 1;

my $length = length($sequence);
my $edge = ($length-$size)/2;
my $start = ceil($edge);
my $stop = $length-floor($edge)-1;


my $i;
my $print_seq;
my $real_pos;

=cut
for($i=1;$i<=10;$i++)
{
	my $print_seq=$sequence;
	print "$prefix$ct-$i\t$print_seq\n";
}
=cut

print "$prefix\t$sequence\n";


my $base;

for($i=$start;$i<=$stop;$i++)
	{
	$print_seq=$sequence;
	
	$base=substr($sequence,$i,1);
	$real_pos=$i+1;
	
	if(uc($base) ne "A")
		{
		$print_seq=$sequence;
		$ct++;
		substr($print_seq,$i,1,"A");
		print $prefix."_M:".$base.$real_pos."A\t".$print_seq."\n";
		}
	if(uc($base) ne "C")
		{
		$print_seq=$sequence;
		$ct++;
		substr($print_seq,$i,1,"C");
		print $prefix."_M:".$base.$real_pos."C\t".$print_seq."\n";
		}
	if(uc($base) ne "G")
		{
		$print_seq=$sequence;
		$ct++;
		substr($print_seq,$i,1,"G");
		print $prefix."_M:".$base.$real_pos."G\t".$print_seq."\n";
		}
	if(uc($base) ne "T")
		{
		$print_seq=$sequence;
		$ct++;
		substr($print_seq,$i,1,"T");
		print $prefix."_M:".$base.$real_pos."T\t".$print_seq."\n";
		}
=cut
	if(uc($base) ne "-")
		{
		$print_seq=$sequence;
		$ct++;
		substr($print_seq,$i,1,"-");
		print $prefix."_D:".$real_pos.":".$real_pos."t".$print_seq."\n";
		}
=cut
	
	}


=cut
my $del_size=5;
my $window_move=1;

for($i=$start;$i<=$stop;$i+=$window_move)
	{
		$print_seq=$sequence;
		$ct++;
		substr($print_seq,$i,$del_size,"-----");
		print "$prefix$ct\t$print_seq\n";

	}

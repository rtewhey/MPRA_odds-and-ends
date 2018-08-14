#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use Getopt::Std;


my %options=();
getopts('R:N:', \%options);

#####
#
#-R = comma separated list of elements to repeat (WT = input seq, G150C, A100G etc etc)
#-N = comma separated list of number of times to repeat element
#
#####


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


my @repeat_ID;
my @repeat_num;
my %repeat_print;

if(exists($options{R}) && exists($options{N}))
	{
	@repeat_ID = split(/,/,$options{R});
	@repeat_num = split(/,/,$options{N});
	
	die "Repeat ID and number not equal lengths:".scalar(@repeat_ID)." ".scalar(@repeat_num)."\n" if(scalar(@repeat_ID) != scalar(@repeat_num));
	
	print STDERR "Duplicating the following oligos\n";

	for($i=0;$i<scalar(@repeat_ID);$i++)
		{
		$repeat_print{$repeat_ID[$i]}=$repeat_num[$i];
		print STDERR "$repeat_ID[$i]\t$repeat_num[$i]\n";
		}
	}

print "$prefix\t$sequence\n";

print "$prefix\t$sequence\n" x ($repeat_print{WT}-1) if(exists($repeat_print{WT}));
	

my $base;
my $print_str;

for($i=$start;$i<=$stop;$i++)
	{
	$print_seq=$sequence;
	$print_str="";
	
	$base=substr($sequence,$i,1);
	$real_pos=$i+1;
	
	if(uc($base) ne "A")
		{
		$print_seq=$sequence;
		$ct++;
		substr($print_seq,$i,1,"A");
		$print_str=$prefix."_M:".$base.$real_pos."A\t".$print_seq."\n";
		print $print_str;
		print $print_str x ($repeat_print{$base.$real_pos."A"}-1) if(exists($repeat_print{$base.$real_pos."A"}));
		}
	if(uc($base) ne "C")
		{
		$print_seq=$sequence;
		$ct++;
		substr($print_seq,$i,1,"C");
		$print_str=$prefix."_M:".$base.$real_pos."C\t".$print_seq."\n";
		print $print_str;
		print $print_str x ($repeat_print{$base.$real_pos."C"}-1) if(exists($repeat_print{$base.$real_pos."C"}));
		}
	if(uc($base) ne "G")
		{
		$print_seq=$sequence;
		$ct++;
		substr($print_seq,$i,1,"G");
		$print_str=$prefix."_M:".$base.$real_pos."G\t".$print_seq."\n";
		print $print_str;
		print $print_str x ($repeat_print{$base.$real_pos."G"}-1) if(exists($repeat_print{$base.$real_pos."G"}));
		}
	if(uc($base) ne "T")
		{
		$print_seq=$sequence;
		$ct++;
		substr($print_seq,$i,1,"T");
		$print_str=$prefix."_M:".$base.$real_pos."T\t".$print_seq."\n";
		print $print_str;
		print $print_str x ($repeat_print{$base.$real_pos."T"}-1) if(exists($repeat_print{$base.$real_pos."T"}));
		}
=cut
	if(uc($base) ne "-")
		{
		$print_seq=$sequence;
		$ct++;
		substr($print_seq,$i,1,"-");
		$print_str=$prefix."_D:".$real_pos.":".$real_pos."\t".$print_seq."\n";
		print $print_str;
		print $print_str x ($repeat_print{$base.$real_pos."-"}-1) if(exists($repeat_print{$base.$real_pos."-"}));
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

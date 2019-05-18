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

my %options=();
getopts('L:R:O:P:', \%options);

#####
#
#-L = Length of left arm
#-R = Length of right arm
#-O = Total oligo length
#-P = Prefix to attach to ID
#
#####
if(exists($options{L}) && exists($options{R})){die "Set length to equal left and right (-O)\n" unless($options{O}==($options{L}+$options{R}));}

my $ID_prefix = $options{P} || "";

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
my %start; #Chr -> ID -> Base
my %stop;
my %chr;
my %strand;



#######################
###Load Allele File
#######################

open (ALLELES, $ALLELES) or die "couldn't open $ALLELES\n";
while (<ALLELES>)
{
	#next if($_ =~ m/^#/ || $_ =~ m/^Chr/ || $_ =~ m/^chr/ || $_ =~ m/^CHR/);
    $_ =~ s/[\n\r]//g;
    @inline = split(/\t/);
    $inline[0]=~s/^chr//;
    die "ERROR: SNP IDs must be unique, duplicate found: $inline[3]\n" if(exists $start{$inline[3]});
    $chr{$inline[3]}=$inline[0];
    $start{$inline[3]}=$inline[1];
    $stop{$inline[3]}=$inline[2];
    $strand{$inline[3]}=$inline[4];
}
close ALLELES;

#######################
###Verify reference alleles. Switch ones that don't match.
#######################
my $id;
my $seq_chr;
my $count;
my $tmp_value;
my $target_seq;
my $a_ct;
my $c_ct;
my $g_ct;
my $t_ct;
my $gc;
my $tmp_seq;
my $updatedID;
my $center;

print STDERR "\n\n######\nFinding Sequence\n";
###Find Reference call
my $in  = Bio::SeqIO->new(-file => $REF , '-format' => 'Fasta');
while ( my $seq = $in->next_seq() ) 
{
	$count=0;
    print STDERR "Chromosome: ",$seq->id," length ",$seq->length(),"\n";
    $seq_chr = $seq->id;
    
    my %count;
	for my $value (values %chr) {
    $count++ if($value eq $seq_chr);
	}
    print STDERR "\t$count sequences to write\n";
    
    foreach $id (keys %chr)
    {
    	if($chr{$id} eq $seq_chr) 
    	{
			$target_seq=$seq->subseq($start{$id},$stop{$id});

    		$target_seq = reverse_complement_IUPAC($target_seq) if($strand{$id} eq "-");
    		die "Strand information is not +/- :$strand{$id}:\n" unless($strand{$id} eq "-" || $strand{$id} eq "+");
    		
    		$tmp_seq=$target_seq;
    		$a_ct=0;$c_ct=0;$g_ct=0;$t_ct=0;
			$a_ct = uc($tmp_seq) =~ tr/A//;
			$c_ct = uc($tmp_seq) =~ tr/C//;
			$g_ct = uc($tmp_seq) =~ tr/G//;
			$t_ct = uc($tmp_seq) =~ tr/T//;
			$gc	= sprintf("%.1f",100*(($c_ct+$g_ct)/($a_ct+$c_ct+$g_ct+$t_ct)));	
			
			$center=$start{$id}+floor(($stop{$id}-$start{$id})/2);
			$updatedID = $chr{$id}.":".$center.":NA:NA";
			print  join("\t",$updatedID,$id,$chr{$id},$start{$id}."-".$stop{$id},$gc,length($target_seq),$strand{$id},"NA","NA","NA","NA","NA","NA",$target_seq,"-")."\n"; 
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


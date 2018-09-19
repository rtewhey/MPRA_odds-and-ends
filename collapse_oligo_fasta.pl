#!/usr/bin/perl
#
#by Ryan Tewhey, 2012
#rtewhey@broadinstitute.edu

use strict;
use warnings;
use Getopt::Std;

sub reverse_complement_IUPAC;

my %options=();
getopts('RW:', \%options);

#####
#
#-R = Collapse reverse complement sequence
#-W = Wiggle room for sequence match
#
#####

my $collapse_flag;

if(exists($options{R}))
	{
	print STDERR "Collapsing reverse complement sequences\n";
	$collapse_flag = 1;
	}
else {$collapse_flag = 0;}

my $wiggle = $options{W} || 0;

print STDERR "Allowing wiggele of $wiggle bp at either end during match\n";

my $FILE1 = $ARGV[0]; #fasta List
my $OUTFILE = $ARGV[1]; #outfile

######


my $inline;
my %by_id;
my %multi_id;
my %multi_id_rc;
my $id;
my $seq;
my $rc_seq;

open (FILE1, "$FILE1") or die("ERROR: can not open $FILE1\n");
while (<FILE1>)
	{
	$_ =~ s/[\n\r]//g;
	$_ =~ s/>//g;
	$id=$_;
 
	$inline=<FILE1>;
	$inline =~ s/[\n\r]//g;
	$seq=$inline;
	
	if(exists $by_id{$id})
		{
			die "Sequences with same id do not match!\n$id\n" unless($by_id{$id} eq $seq);
			next;
		}	
	else
		{
		$by_id{$id}=$seq;
		}
	 }


my $truncated;
my $rev_truncated;

my %match_ct;
my $search_id;
my $search_seq;
my %delete;

my @fwd_id;
my @rc_id;
my @id_list;

my $rc_txt = "RC"; #text denoting reverse complement

foreach $id (keys %by_id)
	{
	if($id =~ m/$rc_txt/)
		{
		push(@rc_id,$id);
		}
	else
		{
		push(@fwd_id,$id);
		}
	}

@id_list=(@fwd_id,@rc_id);
@id_list = sort { length $a <=> length $b } @id_list;

foreach $id (@id_list) #walk through IDs sorted for fwd strand so RCs are collapsed.
{
	$match_ct{$id}=0;
	
	if(exists($delete{$id}))
	{
	next;
	}
	
	$truncated = substr($by_id{$id}, $wiggle, -$wiggle) if($wiggle > 0);
	$truncated = $by_id{$id} if($wiggle == 0);
	die "Wiggle is an invalid value: $wiggle\n" unless($wiggle >= 0);
	
	$rev_truncated  = reverse_complement_IUPAC($truncated);
	
	$truncated=uc($truncated);
	foreach $search_id (keys %by_id)
	{
		if(uc($by_id{$search_id}) =~ m/$truncated/)
		{
		$match_ct{$id}++;
		push(@{$multi_id{$id}},$search_id) unless($search_id eq $id);
		$delete{$search_id}=1 unless($search_id eq $id);
		}
		elsif(uc($by_id{$search_id}) =~ m/$rev_truncated/ && $collapse_flag == 1)
		{
		$match_ct{$id}++;
		push(@{$multi_id_rc{$id}},$search_id) unless($search_id eq $id);
		$delete{$search_id}=1 unless($search_id eq $id);
		}
	}
	die "No matches.... it should have matched itself!\n$id\n" if ($match_ct{$id} eq 0);
}


foreach $id (keys %delete)
	{
	#print STDERR "$id\n";
	delete $by_id{$id};	
	}


my $new_ID;

open (FASTA, ">$OUTFILE.fa") or die("ERROR: can not create $OUTFILE.fa: $!\n");
open (KEYFILE, ">$OUTFILE.key") or die("ERROR: can not create $OUTFILE.key: $!\n");

@id_list = sort { $a cmp $b } (keys %by_id);

foreach $id (@id_list)
	{
	
	$new_ID = "(".$id;
	
	if(exists($multi_id{$id}))
		{
		$new_ID = $new_ID.",".join(";", @{$multi_id{$id}});
		}
	
	$new_ID = $new_ID.")";
	
	if(exists($multi_id_rc{$id}))
		{
		$new_ID = $new_ID."[".join(";", @{$multi_id_rc{$id}})."]";
		}
			
	if(exists($multi_id_rc{$id}) || exists($multi_id{$id}))
		{
		print FASTA ">".$new_ID."\n";
		print FASTA $by_id{$id}."\n";
		
		if(exists($multi_id_rc{$id}) && exists($multi_id{$id}))
			{print KEYFILE $id."\t".join(";", @{$multi_id{$id}})."\t".join(",", @{$multi_id_rc{$id}})."\n";}
		elsif(exists($multi_id_rc{$id}))
			{print KEYFILE $id."\t-\t".join(";", @{$multi_id_rc{$id}})."\n";}
		elsif(exists($multi_id{$id}))
			{print KEYFILE $id."\t".join(";", @{$multi_id{$id}})."\t-\n";}

		}
	else
		{
		print FASTA ">".$id."\n";
		print FASTA $by_id{$id}."\n";		
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

#!/usr/bin/perl
#
#by Ryan Tewhey, 2012
#rtewhey@broadinstitute.edu

use strict;
use warnings;
use Getopt::Std;

sub reverse_complement_IUPAC;

my %options=();
getopts('RFW:I:S:', \%options);

#####
#
#-R = Collapse reverse complement sequence  (Default no collapse)
#-W = Wiggle room for sequence match
#-F = full length match (will not collapse two oligos 200 and 199 even if perfect)
#-I = Column with ID (default = 1)
#-S = Column with sequence (default = 14)
#####

my $collapse_flag;
my $full_length_flag;

if(exists($options{R}))
	{
	print STDERR "Collapsing reverse complement sequences\n";
	$collapse_flag = 1;
	}
else {$collapse_flag = 0;}

my $wiggle = $options{W} || 0;

print STDERR "Allowing wiggele of $wiggle bp at either end during match\n";

if(exists($options{F}))
	{
	print STDERR "Requiring full length matches\n";
	$full_length_flag = 1;
	die "Cannot have wiggle AND Full length (-F) flag\n" if($wiggle > 0);
	}
else {$full_length_flag = 0;}

my $ID_col = $options{I} || 1;
my $SEQ_col = $options{S} || 14;

print STDERR "ID in column $ID_col\nSequencing in column $SEQ_col\n\n";



my $FILE1 = $ARGV[0]; #fasta List
my $OUTFILE = $ARGV[1]; #outfile

######


my @inline;
my %line;
my %by_id;
my %multi_id;
my %multi_id_rc;
my $id;
my $seq;
my $rc_seq;

my %seq_ct;

##Fast match first to narrow candidates for perfect matches

open (FILE1, "$FILE1") or die("ERROR: can not open $FILE1\n");
while (<FILE1>)
	{
	$_ =~ s/[\n\r]//g;
    @inline = split(/\s+/);
	$id=$inline[$ID_col-1];
	$seq=$inline[$SEQ_col-1];
	@{$line{$id}}=@inline;
	
	if(exists $by_id{$id})
		{
			die "Sequences with same id do not match!\n$id\n" unless($by_id{$id} eq $seq);
			next;
		}	
	else
		{
		$by_id{$id}=$seq;
		$seq_ct{$seq}++;
		}
	 }
	 
close FILE1;

%by_id=();
my %singles_by_id;

open (FILE1, "$FILE1") or die("ERROR: can not open $FILE1\n");
while (<FILE1>)
	{
	$_ =~ s/[\n\r]//g;
    @inline = split(/\s+/);
	$id=$inline[$ID_col-1];
	$seq=$inline[$SEQ_col-1];
	@{$line{$id}}=@inline;
	
	if($seq_ct{$seq} > 1)
		{
		$by_id{$id}=$seq;
		}
	else
		{
		$singles_by_id{$id}=$seq;
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

my $ct=1;

foreach $id (@id_list) #walk through IDs sorted for fwd strand so RCs are collapsed.
{
	$match_ct{$id}=0;
	
	$ct++;
	print STDERR "\r$ct";
	
	if(exists($delete{$id}))
	{
	next;
	}
	
	$truncated = substr($by_id{$id}, $wiggle, -$wiggle) if($wiggle > 0);
	$truncated = $by_id{$id} if($wiggle == 0);
	die "Wiggle is an invalid value: $wiggle\n" unless($wiggle >= 0);
	
	$rev_truncated  = reverse_complement_IUPAC($truncated) if($collapse_flag == 1);
	$truncated=uc($truncated);
	foreach $search_id (keys %by_id)
	{
		if($full_length_flag == 1)
			{
			if(uc($by_id{$search_id}) eq $truncated)
				{
				$match_ct{$id}++;
				push(@{$multi_id{$id}},$search_id) unless($search_id eq $id);
				$delete{$search_id}=1 unless($search_id eq $id);
				}
			elsif( ( $collapse_flag == 1) && (uc($by_id{$search_id}) eq $rev_truncated) )
				{
				$match_ct{$id}++;
				push(@{$multi_id_rc{$id}},$search_id) unless($search_id eq $id);
				$delete{$search_id}=1 unless($search_id eq $id);
				}
			
			}
		else
			{
			if(uc($by_id{$search_id}) =~ m/$truncated/)
				{
				$match_ct{$id}++;
				push(@{$multi_id{$id}},$search_id) unless($search_id eq $id);
				$delete{$search_id}=1 unless($search_id eq $id);
				}
			elsif( ( $collapse_flag == 1) && (uc($by_id{$search_id}) =~ m/$rev_truncated/) )
				{
				$match_ct{$id}++;
				push(@{$multi_id_rc{$id}},$search_id) unless($search_id eq $id);
				$delete{$search_id}=1 unless($search_id eq $id);
				}
			}
	}
#	die "No matches.... it should have matched itself!\n$id\n" if ($match_ct{$id} eq 0);
}


foreach $id (keys %delete)
	{
	#print STDERR "$id\n";
	delete $by_id{$id};	
	}

for $id (keys %singles_by_id)
	{
	if(exists $by_id{$id})
		{
		die "ID showed up in multi and single hash!\n$id\n";
		}
	else
		{
		$by_id{$id} = $singles_by_id{$id};
		}
	}

    
my $new_ID;
my $tmp_id;

open (FASTA, ">$OUTFILE.seq") or die("ERROR: can not create $OUTFILE.fa: $!\n");
open (KEYFILE, ">$OUTFILE.key") or die("ERROR: can not create $OUTFILE.key: $!\n");
open (KEYFILE_FULL, ">$OUTFILE.fullkey") or die("ERROR: can not create $OUTFILE.key: $!\n");

@id_list = sort { $a cmp $b } (keys %by_id);

foreach $id (@id_list)
	{
	
	$new_ID = "(".$id;
	
	if(exists($multi_id{$id}))
		{
		$new_ID = $new_ID.";".join(";", @{$multi_id{$id}});
		}
	
	$new_ID = $new_ID.")";
	
	if(exists($multi_id_rc{$id}))
		{
		$new_ID = $new_ID."[".join(";", @{$multi_id_rc{$id}})."]";
		}
			
	if(exists($multi_id_rc{$id}) || exists($multi_id{$id}))
		{
		if($SEQ_col==14 && $ID_col==1)
			{
			print FASTA join("\t",$new_ID,"multiple",${$line{$id}}[2],${$line{$id}}[3],${$line{$id}}[4],${$line{$id}}[5],${$line{$id}}[6],${$line{$id}}[7],${$line{$id}}[8],${$line{$id}}[9],${$line{$id}}[10],${$line{$id}}[11],${$line{$id}}[12],$by_id{$id},${$line{$id}}[14])."\n";
			}
		else
			{
			${$line{$id}}[$SEQ_col-1]=$by_id{$id};
			${$line{$id}}[$ID_col-1]=$new_ID;
			print FASTA join("\t",@{$line{$id}})."\n";
			}

		if(exists($multi_id_rc{$id}) && exists($multi_id{$id}))
			{
			print KEYFILE $id."\t".join(";", @{$multi_id{$id}})."\t".join(";", @{$multi_id_rc{$id}})."\n";
			print KEYFILE_FULL $id."\t".$new_ID."\n";
			foreach $tmp_id (@{$multi_id{$id}})
				{
				print KEYFILE_FULL $tmp_id."\t".$new_ID."\n";
				}
			}
		elsif(exists($multi_id_rc{$id}))
			{
			print KEYFILE $id."\t-\t".join(";", @{$multi_id_rc{$id}})."\n";
			print KEYFILE_FULL $id."\t".$new_ID."\n";
			foreach $tmp_id (@{$multi_id_rc{$id}})
				{
				print KEYFILE_FULL $tmp_id."\t".$new_ID."\n";
				}
			}
		elsif(exists($multi_id{$id}))
			{
			print KEYFILE $id."\t".join(";", @{$multi_id{$id}})."\t-\n";
			print KEYFILE_FULL $id."\t".$new_ID."\n";
			foreach $tmp_id (@{$multi_id{$id}})
				{
				print KEYFILE_FULL $tmp_id."\t".$new_ID."\n";
				}		
			}

		}
	else
		{
		print FASTA join("\t",@{$line{$id}})."\n";		
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

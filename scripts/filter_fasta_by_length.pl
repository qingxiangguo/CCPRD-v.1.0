#!/usr/bin/perl

use warnings;
use Bio::SeqIO;

$scriptname=$0; $scriptname =~ s/.+\///g;
if ($#ARGV != 3 || $ARGV[0] eq "-h")
        {print "\nRead in a fasta file and filter it by length\n";
        print "Usage:\t $scriptname <fasta_file> <min_length> <max_length> <output_name>\n";
        print "\n"; exit;
        }

my $file = $ARGV[0]; 
my $min = $ARGV[1];
my $max = $ARGV[2];
my $out = $ARGV[3];

open (FILE, ">>$out") or die ("Error : Cannot open file $out for writing..!\n");

my $seq_in  = Bio::SeqIO->new( -format => 'fasta',-file => $file);

while( my $seq1 = $seq_in->next_seq() ) {	
	
	my $id  = $seq1->primary_id;
	chomp $id;
	my $seq = $seq1->seq;
	chomp $seq;
	my $lseq = length($seq);
        my $desc = $seq1->desc;
	chomp $desc;
	if($lseq>=$min && $lseq <=$max){
		print FILE ">",$id,"\t", $desc, "\n",$seq,"\n";	
	}
}

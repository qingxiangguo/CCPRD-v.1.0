#!/usr/bin/perl

use warnings;
use Bio::SeqIO;


my $file = $ARGV[0];

open (FILE, ">>flated") or die ("Error : Cannot open file for writing..!\n");

my $seq_in  = Bio::SeqIO->new( -format => 'fasta',-file => $file);

while( my $seq1 = $seq_in->next_seq() ) { 

my $id  = $seq1->primary_id;
        chomp $id;
        my $seq = $seq1->seq;
        chomp $seq;
                print FILE ">",$id,"\n",$seq,"\n";
}


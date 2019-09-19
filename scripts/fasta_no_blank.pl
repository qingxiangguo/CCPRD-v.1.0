#!/usr/bin/perl
use strict;

my $usage = "USAGE:\n\tperl fasta_no_blank.pl fastaFile > noBlankFastaFile\n";
if (@ARGV == 0) { die $usage }

$_ = <>;
print;
my $total_seq_num = 1;
while (<>) {
	if (/^>/) { $total_seq_num ++; print "\n$_" }
	else      { s/\s+?$//g; print }
}
print "\n";
print STDERR "$total_seq_num Sequences in total\n";

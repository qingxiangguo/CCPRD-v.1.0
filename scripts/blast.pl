#!/usr/bin/perl
use strict;
my $usage = "Usage:\n\tperl blast.pl blastProgramType blastDB fastaFile evalue threads outPrefix outfmt\n\n\t7 parameters should be given, and the final result is outPrefix.xml or outPrefix.tab\n\nAn example:\n\tperl blast.pl blastp nr proteins.fa 1e-3 24 nr 5\n\n";
if (@ARGV != 7){die $usage}

my ($program, $db, $fasta, $evalue, $threads, $outPrefix, $outfmt) = @ARGV;
open IN, $fasta or die $!;
my ($seqID, @seqID, %seq);
while (<IN>) {
    chomp;
    if (/^>(\S+)/) { $seqID = $1; push @seqID, $seqID; }
    else           { $seq{$seqID} .= $_ }
}

mkdir "blast.tmp" unless -e "blast.tmp";
open COM, '>', "command.blast.list" or die $!;

foreach my $seq_id (@seqID) {
    $_ = $seq_id;
    s/\|/\\\|/;
    open FASTA, '>', "blast.tmp/$seq_id.fa" or die $!;
    print FASTA ">$seq_id\n$seq{$seq_id}\n";
    close FASTA;

    print COM "$program -query blast.tmp/$_.fa -db $db -evalue $evalue -num_threads 1 -outfmt $outfmt -out blast.tmp/$_.out\n";
}

system `ParaFly -c command.blast.list -CPU $threads &> command.blast.log`;

unless ( -e "FailedCommands" ) {
	open OUT, '>', "$outPrefix.tab" or die $!;
    foreach (@seqID) {
        open IN, "blast.tmp/$_.out" or die $!;
        print OUT <IN>;
    }
}

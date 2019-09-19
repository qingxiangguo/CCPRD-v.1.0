#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 genome_seq.fasta homolog_proteins.fasta homolog_gene_region.tab flanking_length CPUs

For example:
    perl $0 genome_seq.fasta homolog_proteins.fasta homolog_gene_region.tab 2000 24

    The file homolog_gene_region.tab is result of blast2geneRegion.pl;
    flanking length is better to be the average gene length.

USAGE
if (@ARGV==0){die $usage}

open IN, $ARGV[0] or die $!;
my (%genome_seq, %homolog_seq, $id);
while (<IN>) {
    chomp;
    if (/>(\S+)/) { $id = $1; }
    else { $genome_seq{$id} .= $_; }
}
close IN;

open IN, $ARGV[1] or die $!;
while (<IN>) {
    chomp;
    if (/>(\S+)/) { $id = $1; }
    else { $homolog_seq{$id} .= $_; }
}
close IN;

mkdir "geneRegion_genewise.tmp" unless -e "geneRegion_genewise.tmp";
open IN, $ARGV[2] or die $!;
open COMMAND, ">", "command.genewise.list" or die $!;
my @out;
while (<IN>) {
    @_ = split /\t/;
    my $length = $_[2] - $_[1] + 1 + $ARGV[3] + $ARGV[3];
    my $start = $_[1] - 1 - $ARGV[3];
    $start = 0 if $start < 0;
    my $seq = substr($genome_seq{$_[0]}, $start, $length);
    open OUT, ">", "geneRegion_genewise.tmp/$_[0].$_[1]\_$_[2].fasta" or die $!;
    print OUT ">$_[0].$_[1]\_$_[2]\n$seq\n";
    push @out, "$_[0]\t$_[1]\_$_[2]\t$_[3]\t$start";
    close OUT;

    my $homolog_seq = $homolog_seq{$_[3]};
    open OUT, ">", "geneRegion_genewise.tmp/$_[3].fasta" or die $!;
    print OUT ">$_[3]\n$homolog_seq\n";

    my $strand = $_[4];
    if ($strand eq "plus") {
        print COMMAND "genewise geneRegion_genewise.tmp/$_[3].fasta geneRegion_genewise.tmp/$_[0].$_[1]\_$_[2].fasta -tfor -pep -quiet -silent >> geneRegion_genewise.tmp/genome_homolog.pep\n";
    }
    else {
        print COMMAND "genewise geneRegion_genewise.tmp/$_[3].fasta geneRegion_genewise.tmp/$_[0].$_[1]\_$_[2].fasta -trev -pep -quiet -silent >> geneRegion_genewise.tmp/genome_homolog.pep\n";
    }
}
close COMMAND;

system("ParaFly -c command.genewise.list -CPU $ARGV[4] &> command.genewise.log") == 0 or die "failed to execute: $!\n";

unless (-e "FailedCommands") {
    foreach (@out) {
		@_ = split /\t/;
		my ($seq, $region, $proteinID, $pos) = ($_[0], $_[1], $_[2], $_[3]);
        open IN, "geneRegion_genewise.tmp/$seq.$region.gff" or die $!;
        while (<IN>) {
            next if m/\/\//;
            @_ = split /\t/;
            my $start = $_[3] + $pos;
            my $end = $_[4] + $pos;
			my $ID = $_[8];
			chomp($ID);
            print "$seq\t$_[1]\t$_[2]\t$start\t$end\t$_[5]\t$_[6]\t$_[7]\tID=$ID;Name=$proteinID;\n";
        }
        close IN;
        print "\n";
    }
}

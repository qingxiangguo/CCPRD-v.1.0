#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usge:
    perl $0 genewise.gff genome.fasta [filterMiddleStopCodon=yes|CompleteGeneModelsOnly=yes] > evm_protein_alignment.gff3 2> genewise_gene_models_completeness_check.txt

USAGE
if (@ARGV==0){die $usage}

my (%gene, @gene);
open IN, $ARGV[0] or die $!;
while (<IN>) {
    if (m/\tcds\t.*ID=(\S+?);/i) {
        my $gene_id = $1;
        push @gene, $gene_id unless exists $gene{$1};
        s/\tcds\t/\tprotein_match\t/;
		@_ = split /\t/;
		if ($_[6] eq "-") {
			$_ = "$_[0]\t$_[1]\t$_[2]\t$_[4]\t$_[3]\t$_[5]\t$_[6]\t$_[7]\t$_[8]";
		}
        $gene{$gene_id} .= $_;
    }
}
close IN;
#print STDERR "file $ARGV[0] was read over\n";

open IN, $ARGV[1] or die $!;
my (%seq, $seq_id);
while (<IN>) {
    chomp;
    if (/^>(\S+)/) { $seq_id = $1; }
    else { $seq{$seq_id} .= $_; }
}
close IN;
#print STDERR "file $ARGV[1] was read over\n";

my ($total_genewise_models_num, $middle_stopCodon_number, $complete_number, %middle, %complete);
foreach my $id (@gene) {
    $total_genewise_models_num ++;

    my @cds = split /\n/, $gene{$id};
    @_ = split /\t/, $cds[0];
    my ($seqID, $strand, $frame) = ($_[0], $_[6], $_[7]);

    my $cds_seq;
    my $line_number;
    if ($strand eq "+") {
        foreach (@cds) {
            $line_number ++;
            @_ = split /\t/;
            my $start = $_[3] - 1;
            my $len = $_[4] - $start;
            $len += 3 if $line_number == @cds;
            $cds_seq .= substr($seq{$seqID}, $start, $len);
        }
    }
    elsif ($strand eq "-") {
        my %sort;
        foreach (@cds) { @_ = split /\t/; $sort{$_} = $_[3]; }
        @cds = sort {$sort{$a} <=> $sort{$b}} @cds;
        foreach (@cds) {
            $line_number ++;
            @_ = split /\t/;
            my $start = $_[3] - 1;
            $start -= 3 if ($line_number == 1 && $start >= 3);
            my $len = $_[4] - $start;
            $cds_seq .= substr($seq{$seqID}, $start, $len);
        }
        $cds_seq = reverse $cds_seq;
        $cds_seq =~ tr/ATCGatcg/TAGCtagc/;
    }

    $cds_seq =~ s/^\w{$frame}//;
    #print STDERR ">$id\n";
    my $pep_seq = &cds2pep($cds_seq, $id);
    #print STDERR ">$pep_seq\n";

    if ($pep_seq =~ m/\*\w/) {
        $middle_stopCodon_number ++;
        $middle{$id} = 1;
    }
    elsif ($pep_seq =~ m/^M.*\*$/) {
        $complete_number ++;
        $complete{$id} = 1;
    }
}

print STDERR "total genewise models number: $total_genewise_models_num
gene models with middle stop codon: $middle_stopCodon_number
complete gene models number: $complete_number\n";

foreach my $gene_id (@gene) {
	my @out = split /\n/, $gene{$gene_id};
	my %sort;
	foreach (@out) { @_ = split /\t/; $sort{$_} =  $_[3]; }
	@out = sort {$sort{$a} <=> $sort{$b}} @out;
	my $out = join "\n", @out;

    if ($ARGV[2] eq "filterMiddleStopCodon=yes") {
        print "$out\n" unless exists $middle{$gene_id};
    }
    elsif ($ARGV[2] eq "CompleteGeneModelsOnly=yes") {
        print "$out\n" if exists $complete{$gene_id};
    }
    else {
        print "$out\n";
    }
}

sub cds2pep {
    my %cds2pep = (
        "TTT" => "F",
        "TTC" => "F",
        "TTA" => "L",
        "TTG" => "L",
        "TCT" => "S",
        "TCC" => "S",
        "TCA" => "S",
        "TCG" => "S",
        "TAT" => "Y",
        "TAC" => "Y",
        "TAA" => "*",
        "TAG" => "*",
        "TGT" => "C",
        "TGC" => "C",
        "TGA" => "*",
        "TGG" => "W",
        "CTT" => "L",
        "CTC" => "L",
        "CTA" => "L",
        "CTG" => "L",
        "CCT" => "P",
        "CCC" => "P",
        "CCA" => "P",
        "CCG" => "P",
        "CAT" => "H",
        "CAC" => "H",
        "CAA" => "Q",
        "CAG" => "Q",
        "CGT" => "R",
        "CGC" => "R",
        "CGA" => "R",
        "CGG" => "R",
        "ATT" => "I",
        "ATC" => "I",
        "ATA" => "I",
        "ATG" => "M",
        "ACT" => "T",
        "ACC" => "T",
        "ACA" => "T",
        "ACG" => "T",
        "AAT" => "N",
        "AAC" => "N",
        "AAA" => "K",
        "AAG" => "K",
        "AGT" => "S",
        "AGC" => "S",
        "AGA" => "R",
        "AGG" => "R",
        "GTT" => "V",
        "GTC" => "V",
        "GTA" => "V",
        "GTG" => "V",
        "GCT" => "A",
        "GCC" => "A",
        "GCA" => "A",
        "GCG" => "A",
        "GAT" => "D",
        "GAC" => "D",
        "GAA" => "E",
        "GAG" => "E",
        "GGT" => "G",
        "GGC" => "G",
        "GGA" => "G",
        "GGG" => "G",
    );
    my $seq = shift @_;
    my $gene = shift @_;
    my $pep;
    while ((length $seq) >= 3) {
        $seq =~ s/(\w{3})//;
        $pep .= $cds2pep{$1};
    }
    #print STDERR "Warning: CDS length of $gene is not multiple of 3\n" if (length $seq) > 0;
    #print STDERR "Warning: Stop Codon appear in the middle of $gene\n" if $pep =~ m/\*\w/;
    return $pep;
}

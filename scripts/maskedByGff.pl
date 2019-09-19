#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    perl $0 file.gff genome.fasta > out.masked.fasta

	--mask_type <string>    default: hardmaskN
	--mask_all

USAGE
if (@ARGV==0){die $usage}

my ($mask_type, $mask_all);
GetOptions(
	"mask_type:s" => \$mask_type,
	"mask_all" => \$mask_all,
);
$mask_type ||= "hardmaskN";

open GFF, '<', $ARGV[0];
open FASTA, '<', $ARGV[1];

my ($fasta_head, @fasta_head, %fasta);
while (<FASTA>) {
    chomp;
    if (/^>(\S+)/) {
        $fasta_head = $1;
        push @fasta_head, $fasta_head;
    }else {
        $fasta{$fasta_head} .= $_;
    }
}

while (<GFF>) {
	unless ($mask_all) {
		next if m/Low_complexity/;
		next if m/Simple_repeat/;
	}
    if (/(\S+)\t.+\t.+\t(\d+)\t(\d+)\t/) {
        my $id = $1;
        my $start = $2 - 1;
        my $length = $3 - $start;
        if ($mask_type eq "softmask") {
            substr($fasta{$id},$start,$length) =~ tr/ATCGN/atcgn/;
        }
        elsif ($mask_type eq "hardmaskX") {
            substr($fasta{$id},$start,$length) =~ tr/ATCGNatcgn/XXXXXxxxxx/;
        }
        elsif ($mask_type eq "hardmaskN") {
            substr($fasta{$id},$start,$length) =~ tr/ATCGXatcgx/NNNNNnnnnn/;
        }
        else {
            die "The right mask type should be: softmask or hardmaskX or hardmaskN!\n";
        }
    }
}

foreach (@fasta_head) {
    my $seq = $fasta{$_};
    $seq =~ s/(\w{60})/$1\n/g;
    chomp($seq);
    print ">$_\n$seq\n";
}

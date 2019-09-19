#!/usr/bin/perl

open (IN, "$ARGV[0]");
open (OUT, ">rename_all.fasta");

my $count = 1;

while (<IN>) {

	chomp;

	if ($_ =~ /^>/) {

	print OUT ">$count\n";

	$count ++ ;
}  else {

	print OUT "$_\n"

}

}

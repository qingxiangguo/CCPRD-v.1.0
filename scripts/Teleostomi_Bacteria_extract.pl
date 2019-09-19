#!/usr/bin/perl

open (IN, "$ARGV[0]");
open (OUT, ">contam_accession");

while (<IN>) {
	chomp;
	if ($_ =~ /Teleostomi|Bacteria/) {
	$_ =~ /^(\S+)(\t)/ ;
	print OUT "$1\n";

}

}


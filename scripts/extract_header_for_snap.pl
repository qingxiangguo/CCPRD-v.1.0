#!/usr/bin/perl

open (IN, "$ARGV[0]");
open (OUT, ">header");

while (<IN>) {
	chomp;
	if ($_ =/^>(\S+)/) {
	print OUT "$1\n";


}


}

#!/usr/bin/perl

my $file = $ARGV[0];

open (IN, $file);
open (OUT, ">seq_remove_by_bam0");

while (<IN>) {
	chomp;
     $_ =~ /^(\S+)(\t)(\S+)(\t)(\S+)(\t)(\S+)(\t)(\S+)(\t)(\S+)(\t)(\S+)(\t)(\S+)/; 

	if ($13 >= 200) {

	if ($_ =~ /Streptophyta|Chordata/) {

     $_ =~ /^(\S+)(\t)(\S+)(\t)(\S+)(\t)(\S+)(\t)(\S+)(\t)(\S+)(\t)(\S+)(\t)(\S+)/;

	if ($9 <=50 ) { 

	print OUT "$1\n";

}
} 

}

}

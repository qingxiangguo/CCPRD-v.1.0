#!/usr/bin/perl

my $seqfile = $ARGV[0];
open (IN, $seqfile );
open (OUT, ">duplicate_remove");

while (<IN>) {
	$hash{$_} = 1; 
	}

my @uniq = keys %hash;

foreach $value (@uniq) {
	print OUT "$value";
}
  

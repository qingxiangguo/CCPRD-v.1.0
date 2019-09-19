#!/usr/bin/perl

open (IN1, "$ARGV[0]");
open (IN2, "$ARGV[1]");
open (OUT, ">final_contam_header");

while (<IN1>) {

	chomp;
	$hash{$_} = 2;
	}



while (<IN2>) {
	chomp;
	$_ =~ /^(\S+)(\t)(\S+)/;
	if ($hash{$3}) { 
	print OUT "$1\n";

}
}

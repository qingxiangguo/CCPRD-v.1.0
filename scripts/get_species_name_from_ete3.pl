#!/usr/bin/perl

open (IN, "<result_1");
open (OUT, ">names");

while (<IN>) {
	chomp;
	our @matches = $_ =~ /u'(.+?)'/g ;

}


foreach $ele (@matches) {
	$ele =~ s/(\s)/_/gm;
	print OUT "$ele\n";

}

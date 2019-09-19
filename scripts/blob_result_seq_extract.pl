#!/usr/bin/perl

my $file = $ARGV[0];

open (IN, $file);
open (OUT, ">seq_for_blast");

while (<IN>) {
	chomp;
     $_ =~ /^(\S+)(\t)(\S+)(\t)(\S+)(\t)(\S+)(\t)(\S+)(\t)(\S+)(\t)(\S+)(\t)(\S+)/; 

	if ($13 >= 200) {

	if ($_ =~ /Firmicutes|Proteobacteria|Chordata/) {

	$_ =~ /^(\S+)/;

	print OUT "$&\n";

} 

}

}

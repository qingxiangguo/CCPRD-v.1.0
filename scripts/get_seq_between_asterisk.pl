#!/usr/bin/perl


$scriptname=$0; $scriptname =~ s/.+\///g;
 if ($#ARGV != 1 || $ARGV[0] eq "-h")
        {print "\nRead in the output of transeq and get a certain length of seq between the asterisks (stop codon)\n";
        print "Usage:\t $scriptname <fasta_file> <min_length>\n";
        print "\n";
        }

 our $min = $ARGV[1];
 our $count = 1;

 open (IN1, "$ARGV[0]" );
 open (OUT, ">between_asterisk");

while (<IN1>) {
	chomp;
	if ($_ =~ />/) {
} else {
	push (@seq_array, $_);

}
}

	foreach $content (@seq_array) {
	 @fields = split /\*/, $content;
	
	foreach $seed (@fields)  {	
	$seq_len = length ($seed);
	if ($seq_len >= $min)  { 
	if ($seed =~ /X/) {

} else {
	print OUT ">",$count, "\n",$seed,"\n";
	$count ++;
} 
}

}

}


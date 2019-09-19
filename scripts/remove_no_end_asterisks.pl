#!/usr/bin/perl

open (IN, "$ARGV[0]");

open (OUT, ">no_inner_asterisks");
	my $count =1;
	while (<IN>) {
	chomp;
     if	($_ =~ /\*/) {

	if ($_ =~ /\*$/) {

		if ($_ =~ /\*\*$/) {		
		s/\*\*/\*/  ;
		print OUT "$_\n";
		$count++;

} elsif ($_ =~ /[A-Z](\*)[A-Z]/) {
	$_ =~ s/\*/X/g;
	 print OUT "$_\n";
        $count ++;

} else {
	print OUT "$_\n";
}


} else {
	$_ =~ s/\*/X/g;
	print OUT "$_\n";
	$count ++;
}

} else {
	print OUT "$_\n";

}
}

	print "$count inner asterisks have been replaced!\n";

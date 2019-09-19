#!/usr/bin/perl

$scriptname=$0; $scriptname =~ s/.+\///g;
if ($#ARGV != 1 || $ARGV[0] eq "-h"|| $ARGV[0] eq "-help" || $ARGV[0] eq "--help")
        {print "\nChange the header of a fasta file by number, you can also designate your name.\n";
                print "Usage:\t $scriptname <fasta_file> <designated_name>\n";
              	print "example:\t $scriptname MH_MCPID.fasta WL\n";
                print "Written by Guo qingxiang, guoqing\@webmail.hzau.edu.cn. \n";
                print "Distributed without any guarantees or restrictions\n";
                        print "\n"; exit;
                                }


open (IN, "$ARGV[0]");
open (OUT, ">ordered.fa");

	my $count =1;

while (<IN>)  {
	chomp;
	if ($_ =~ /^>(.*)/) {
	s/^>(.*)/>$count/;
	print OUT "$_"."_$ARGV[1]"."\n";
	$count ++;
}else {
	print OUT "$_\n";


}

}

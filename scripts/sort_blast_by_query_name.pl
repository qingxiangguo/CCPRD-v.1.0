#!/usr/bin/perl


$scriptname=$0; $scriptname =~ s/.+\///g;
if ($ARGV[0] eq "-h")
        {print "\nSort the order of a Blast output fmt6 by query name\n";
                print "Usage:\t $scriptname <outfmt6_file\n";
                print "Written by Guo qingxiang, guoqing\@webmail.hzau.edu.cn. \n";
                print "Distributed without any guarantees or restrictions\n";
                        print "\n"; exit;
                                }


$file = $ARGV[0];

open (IN, "$file");
open (OUT, ">sorted_output");

while (<IN>) {
	push (@array, "$_");				
}

@sorted = map  { $_->[0] }                          #Sort the header of contig file using Schwartzian Transform
             sort { $a->[1] <=> $b->[1] }
             map  { [$_, $_=~/(\d+)/] }
                 @array;

foreach $key (@sorted) {                           #Output the contig result
	chomp;	
	print OUT "$key";
}



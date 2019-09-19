#!/usr/bin/perl
use List::Compare; 

$scriptname=$0; $scriptname =~ s/.+\///g;
if ($#ARGV != 1 || $ARGV[0] eq "-h")
        {print "\nRead in two file  and give you the intersection, union, unique file of them\n";
        print "Usage:\t $scriptname A B\n";
        print "Output:\n";
        print "\t A_only\t         	     file contains what only appears in A\n";
        print "\t B_only\t         	     file contains what only appears in B\n";
        print "\t inter_of_A_and_B\t         file contains the intersection of A and B\n";
        print "\t union_of_H_and_M\t         file contains the union of A and B\n";
        print "Written by Guo qingxiang, guoqing\@webmail.hzau.edu.cn. \n";
                print "Distributed without any guarantees or restrictions\n";

        print "\n"; exit;
        }


open (IN1, "$ARGV[0]");
open (IN2, "$ARGV[1]");



open (OUT1, ">inter_of_$ARGV[0]_and_$ARGV[1]");
open (OUT2, ">union_of_$ARGV[0]_and_$ARGV[1]");
open (OUT3, ">$ARGV[0]_only");
open (OUT4, ">$ARGV[1]_only");

while (<IN1>) {
        chomp;
        push (@array1, "$_");
                   }

 while (<IN2>) {
        chomp;
        push (@array2, "$_");
       }

 $lc = List::Compare->new(\@array1, \@array2); 

@array1_only = $lc->get_unique;

@array2_only = $lc->get_complement;

@inter = $lc->get_intersection;

 @union = $lc->get_union;

 foreach (@inter) {
print OUT1 "$_\n";

} 

foreach (@union) {
print OUT2 "$_\n";
}

 foreach (@array1_only) {
print OUT3 "$_\n";

}


 foreach (@array2_only) {
print OUT4 "$_\n";

}

#!/usr/bin/perl
use strict;
use warnings;

# accept SNAP format such as
# XA8	SNAP	Einit	6161	7325	-5.800	-	.	XA8r-snap.4
# XA8	SNAP	Eterm	5974	6008	5.650	-	.	XA8r-snap.4
# output standard GFF

my $current_id="";
my $real_start=0;
my $real_end=0;
my @arr=();
while(my $line=<>) {
        chomp($line);
        my ($scf,$source,$type,$start,$end,$blank1,$blank2,$blank3,$id)=split('\t',$line);
        if($type eq "Einit") {
                if($current_id ne "") {
                        print "$scf\t$source\tmRNA\t$real_start\t$real_end\t$blank1\t$blank2\t$blank3\tID=$current_id\n";
                }
                foreach (@arr) {
                        print $_."\n";
                }
                @arr=();
                $current_id=$id;
                $real_start=$start;
        }
        push(@arr, "$scf\t$source\tCDS\t$start\t$end\t$blank1\t$blank2\t$blank3\tParent=$current_id");
        $real_end=$end;
}

#!usr/bin/perl

use strict; 

my $scriptname=$0; $scriptname =~ s/.+\///g;
if ($#ARGV != 0 || $ARGV[0] eq "-h" || $ARGV[0] eq "-help")
        {
        print"##############################################################\n";
        print "#";
        print "\n#  Script name:        species_distribution.pl\n";
        print "#\n#  Version:   V 1.0\n";
        print "#\n#     Written by Guo qingxiang, guoqing\@webmail.hzau.edu.cn. \n";
        print "#\n#     Distributed without any guarantees or restrictions\n";
        print "#";
        print "\n#############################################################\n";
        print "#\n#  Function:\tReading in a list of Latin species names, \n count their frequency and then output by descending order of the frequencies\n";
        print "#\n#  Usage:\t$scriptname <species_list>\n";
        print "\n"; exit;
        }


my @inputs; 
my %in; 
my @afile= @ARGV ; 
my $filename; 
my $sucess;

open OUT, ">Species_distribution";

foreach $filename ( @afile ) { 
my $sucess = open FILE,$filename; 
if( ! $sucess ){ 
die "cannot read $filename:$!"; 
} 
while(<FILE>){ 
chomp; 
push(@inputs,$_); 
} 
} 
foreach( @inputs){ 
$in{$_}++; 
} 

sub hashValueDescendingNum {  #create a function for sort hash by value, not keys!
   $in{$b} <=> $in{$a};
}


my @k=sort hashValueDescendingNum (keys %in); 
foreach (@k){ 
	
   print OUT "$_","\t","$in{$_}\n"; 

}



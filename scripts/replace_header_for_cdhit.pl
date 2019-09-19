#!/usr/bin/perl

use Getopt::Long;

my ($cdhit_file, $help, $tag);

GetOptions ('cdhit_file=s' => \$cdhit_file,'tag=s' => \$tag, 'help!' => \$help) or die;

if ($help) {
print "\nTake the output of CDHIT and add a tag to the header for downstream hybrid\n\n";
                print "Usage:\t $scriptname -c <cd-hit_output_file> -t <The_name_tag>\n\n";
                print "Arguments:\n\n";
                print "-cdhit_file or -c: cd-hit output file.\n\n";
                print "-tag or -t: name tag, MN stand for myxo_nucl e. g. :\n\n";
                print "Written by Guo qingxiang, guoqing\@webmail.hzau.edu.cn. \n\n";
                print "Distributed without any guarantees or restrictions\n\n";
                        print "\n"; exit;


}

open (IN1, $cdhit_file);
open (OUT1, ">>${tag}_cdhit.fasta");
open (OUT2, ">countpart");

while (<IN1>) {
	chomp;
	if ($_ =~ /^>(.*)/) {
	$count +=1;
	print OUT1 ">${tag}_${count}\n";
	print OUT2 "$_ => ${tag}_${count}\n"; 

}else {
	print OUT1 "$_\n";   


}
}


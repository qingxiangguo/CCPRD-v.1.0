#!/usr/bin/perl

$scriptname=$0; $scriptname =~ s/.+\///g;
if ($#ARGV != 0 || $ARGV[0] eq "-h")
        {print "\nThe DIAMOND blastx is fast but without taxonomy information, this script can add taxon information to blastx result. \n";
                print "Dependency:\t acc2tax, see https://github.com/richardmleggett/acc2tax\n";
                print "Usage:\t $scriptname <diamond_blastx_outfmt6_result>\n";
                print "Written by Guo qingxiang, guoqing\@webmail.hzau.edu.cn. \n";
                print "Distributed without any guarantees or restrictions\n";
                        print "\n"; exit;
                                }

open IN, "$ARGV[0]";
open OUT, ">diamond_blastx_with_tax";
open TMP1, ">tmp1";
while (<IN>) {
	chomp;
	$_ =~ /gi\|(\S+?)\|(\w+?)\|(\S+?)\.(\d)\|/;
	print TMP1 "$3\n";
}
system ("acc2tax -i tmp1 -o tmp_result -d /home/train/public_database/accession2taxid/ -p");
open IN2, "tmp_result";
	while (<IN2>) {
	chomp;
	$_ =~ /(\S+)\t(.*)/;
	$hash{$1} = $2;
}

open IN, "$ARGV[0]";
while (<IN>) {
	chomp;
	$_ =~ /gi\|(\S+?)\|(\w+?)\|(\S+?)\.(\d)\|/;
	if ($hash{$3}) {
	print OUT "$_\t$hash{$3}\n";
} else {

}
}

system ("rm tmp1");
system ("rm tmp_result");

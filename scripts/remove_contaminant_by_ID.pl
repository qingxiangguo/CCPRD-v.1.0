#!/usr/bin/perl

use warnings;
use Bio::SeqIO;

$scriptname=$0; $scriptname =~ s/.+\///g;
if ($#ARGV != 1 || $ARGV[0] eq "-h")
        {print "\nFilter a fasta file according to a given list, and then producing a fasta file containing surviving sequences\n";
                print "Usage:\t $scriptname <fasta_file> <bad_list>\n";
		print "Written by Guo qingxiang, guoqing\@webmail.hzau.edu.cn. \n";
		print "Distributed without any guarantees or restrictions\n";
                        print "\n"; exit;
                                }

my $file = $ARGV[0];
my $contami = $ARGV[1];

open (IN, $contami);
open (OUT, ">survive.fasta");

while (<IN>) {
	chomp;
	$hash{$_}=1;
	}

my @uniq = keys %hash;

	my $seq_in  = Bio::SeqIO->new( -format => 'fasta',-file => $file);
	while( my $seq1 = $seq_in->next_seq() ) {
		my $id  = $seq1->primary_id;
		chomp $id;
	        my $seq = $seq1->seq;
       		 chomp $seq;
		if ($id ~~ @uniq) {
	} else {
		print OUT ">",$id,"\n",$seq,"\n";
		}
	} 
				

                                

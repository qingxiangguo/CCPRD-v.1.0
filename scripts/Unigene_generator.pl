#!/usr/bin/perl

$scriptname=$0; $scriptname =~ s/.+\///g;

use Getopt::Long;
use Bio::SeqIO;

my ($singleton, $contig, $tag);

GetOptions ('singleton=s' => \$singleton,'tag=s' => \$tag,'contig=s' => \$contig, 'help!' => \$help) or die;

open (IN1, $singleton);
open (IN2, $contig);
open (OUT1, ">>${tag}_Unigene.fasta");
open (OUT2, ">countpart");

if ($help) {
print "\nTake the output Singleton and Contig file of TGICL, and combining them into the final Unigene\n\n";
                print "Usage:\t $scriptname -s <singleton_fasta_file> -c <all_contig_fasta_file> -t <The_name_tag>\n\n";
                print "Arguments:\n\n";
                print "-singleton or -s: extracted singleton fasta file.\n\n";
                print "-contig or -c: combined contig files (If you use multi-process of TGICL)\n\n";
                print "-tag or -t: name tag, the species name e. g. :\n\n";
                print "Written by Guo qingxiang, guoqing\@webmail.hzau.edu.cn. \n\n";
                print "Distributed without any guarantees or restrictions\n\n";
                        print "\n"; exit;


}

###### Extract and sort the header of Contig file (for downstream sorting)#######

while (<IN2>) {
	if ($_ =~ />(.*)/) {
	push (@array, "$1");				
}
}

@sorted = map  { $_->[0] }                          #Sort the header of contig file using Schwartzian Transform
             sort { $a->[1] <=> $b->[1] }
             map  { [$_, $_=~/(\d+)/] }
                 @array;

my $seq_in = Bio::SeqIO->new(-format=>'fasta',-file=>$contig);
	while( my $seq1 = $seq_in->next_seq() ) {
                my $id  = $seq1->primary_id;
                my $seq = $seq1->seq;
	$data{$id} = $seq;                         #Conneect the ID and Seq by hash
	}

foreach $key (@sorted) {                           #Output the contig result
	chomp;	
	print OUT1 ">".${key},"_",$tag, "\n", $data{$key}, "\n";


}

########## Replacing the header of singleton file ##################################

while (<IN1>) {
        chomp;
        if ($_ =~ /^>(.*)/) {
        $count +=1;
        $match{$count} = $_;                      # connect the seq_id with number through hash
        print OUT1 ">Unigene${count}_$tag\n";    #change the header if header was detected
        }else{
        print OUT1 "$_\n";                       #Keep the content if no title was detected

}
} 
########## Output the countpart ###################################################

        @array_key = keys %match;
        foreach $key2 (sort {$a <=> $b} @array_key) {
        print OUT2 "$match{$key2}  => Unigene$key2\_$tag\n";
        }       
        
close (IN1);


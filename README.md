# CCPRD: A novel analytical framework for comprehensive proteomic reference database construction of non-model organisms

* Contributors
* What is CCPRD?
* Dependencies
* Pipelines
* Examples
* Citing CCPRD and software called by CCPRD
* License

# Contributors
Qingxiang (Allen) Guo  
Postdoctoral Fellow  
Northwestern University, Feinberg School of Medicine

# What is CCPRD?
Here we propose an efficient framework for constructing the comprehensive protein reference database, "customized comprehensive proteomic reference database (CCPRD)", which incorporated the draft genomes and transcriptomes. Compared with previous protocols, our method has superiorities in peptide and protein identification numbers, number of entries database size, and completeness.

# Dependencies
Perl version >= 5.10, Python version = 2.7, acc2tax version 0.4, Anvi’o version 2.0.2, Augustus version 3.2.2, BLAST version 2.2.25, Blobtools version 0.9.19, Bowtie2 version 2.2.9, CD-HIT version 4.6.6, Diamond version 0.8.31.93, EMBOSS version 6.6.0, Ete 3, EVidenceModeler version 1.1.1, GeneMarkS-T version 5.1, GeneMark-ES-ET version 2.5m, PASA version 2.2.0, RepeatMasker version 4.0.5, RepeatModeler version 1.0.4, Samtools version 1.3.1, SNAP version 2006-07-28, SOAPdenovo version 2.04, TGICL version 2.1, Tophat version 2.1.1, TransDecoder version 3.0.0, Trinity version 2.2.0, Wise version 2.4.1

# Pipelines
1.Assembly the transcriptomes and genomes

1.1 Assembly the transcriptomes from clean data

Trinity --seqType fq --left M.wulii_1_fixed.fastq --right M.wulii_2_fixed.fastq --max_memory 200G --min_contig_length 200 --min_glue 3 --group_pairs_distance 250 --path_reinforcement_distance 85 --min_kmer_cov 2 --jaccard_clip --normalize_reads --CPU 32

mv Trinity.fasta wulii.fasta

1.2 Remove redundancy in transcriptomes with CD-HIT

cd-hit-est -i wulii.fasta -o wulii.fasta_cdhit -c 0.95 -n 10 -d 0 -M 16000 -T 8 &

1.3 Further cluster the transcriptomes with TGICL

tgicl -l 40 -c 10 -v 25 -O '-repeat_stringency 0.95 -minmatch 35 -minscore 35' -F wulii.fasta_cdhit

1.4 Process the cd-hit results to get unigenes

fast_extract_seq_from_fasta.pl wulii.fasta_cdhit wulii.fasta_cdhit.singletons > extracted.fasta

cat asm_1/contigs asm_2/contigs asm_3/contigs asm_4/contigs > all.contigs

Unigene_generator.pl -s extracted.fasta  -c all.contigs  -t WL

#This will output WL_Unigene.fasta for downstream analyses.

1.5 Assembly the genomes from clean data

SOAPdenovo all -s lib.cfg -K 51 -D 1 -o WL >> soap.log

mv WL.scafSeq genome.fasta

Gapcloser -a genome.fasta -b config.txt -o gapcloser.fasta -t 32

mv gapcloser.fasta WL_genome.fasta

2.Remove potential host and bacterial contamination in transcriptome data with "conservative reciprocal best blast hit" method

2.1	Collected all available proteins and nucleotide sequences for constructing host or close-related species no-redundant database
#Host database contained (proteins or transcripts): Cyprinus carpio genome predicted proteins (GCF_000951615.1), Danio rerio genome predicted proteins (GCA_000002035.3), Carassius proteins in NCBI (4087), Carassius auratus EST (11307), Carassius auratus transcriptomes, Cyprinus carpio EST (497380); close-related no-redundant database contained (proteins or transcripts): Anemonia viridis proteins in NCBI (346)，Clytia hemisphaerica proteins in NCBI (529), Exaiptasia pallida genome predicted proteins (GCA_001417965.1), Hydra vulgaris genome predicted proteins (Hydra_RP_1.0), Hydractinia echinata proteins in NCBI (121), Metridium senile proteins in NCB (123), Nematostella vectensis genome predicted protein (GCA_000209225.1), Thelohanellus kitauei genome predicted proteins (GCA_000827895.1), Buddenbrockia plumatellae EST (765), Tetracapsuloides bryosalmonae EST (308), Polypodium hydriforme transcriptomes, Myxobolus pendula transcriptome, Myxobolus cerebralis transcriptomes, Chironex fleckeri transcriptomes.

2.2	Remove redundancy in each database

cd-hit-est -i nucl_host.fasta -o nucl_host_cdhit.fasta -c 0.9 -n 8 -T 8

cd-hit-est -i nucl_myxo.fasta -o nucl_myxo_cdhit.fasta -c 0.9 -n 8 -T 8

cd-hit -i prot_host.fasta -o prot_host_cdhit.fasta -c 0.9

d-hit -i prot_myxo.fasta -o prot_myxo_cdhit.fasta -c 0.9

2.3	Add tags in databases

replace_header_for_cdhit.pl -c nucl_host_cdhit.fasta -t HN

replace_header_for_cdhit.pl -c prot_host_cdhit.fasta -t HP

replace_header_for_cdhit.pl -c nucl_myxo_cdhit.fasta -t MN

replace_header_for_cdhit.pl -c prot_myxo_cdhit.fasta -t MP

#Sequence number: HN (45486), HP (59764), MN (103090), MP (78786).

2.4	Start hybridization

#For proteins:

tblastn -query HP_cdhit.fasta -db /home/gqx/transcriptome/assembly/M.wulii/clustering/MWdb -max_hsps 1 -out result_HP -evalue 1e-5 -outfmt 6 -num_threads 8

tblastn -query MP_cdhit.fasta -db /home/gqx/transcriptome/assembly/M.wulii/clustering/MWdb -max_hsps 1 -out result_MP -evalue 1e-5 -outfmt 6 -num_threads 8

#For nucleotides:

tblastx -query MN_cdhit.fasta -db /home/gqx/transcriptome/assembly/ M.wulii/clustering/HHdb -max_hsps 1 -out result_MN -evalue 1e-5 -outfmt 6 -num_threads 8

tblastx -query HN_cdhit.fasta -db /home/gqx/transcriptome/assembly/ M.wulii/clustering/HHdb -max_hsps 1 -out result_HN -evalue 1e-5 -outfmt 6 -num_threads 8

2.5	Process above results and remove the transcripts that only matched to host databases

#Extract non-redundant subjects for each comparison result.

cat result_MP | cut -f 2 > 1 | remove_duplicate.pl 1 | mv duplicate_remove MP

cat result_HP | cut -f 2 > 1 | remove_duplicate.pl 1 | mv duplicate_remove HP

cat result_MN | cut -f 2 > 1 | remove_duplicate.pl 1 | mv duplicate_remove MN

cat result_HN | cut -f 2 > 1 | remove_duplicate.pl 1 | mv duplicate_remove HN

#Integrating and compare non-redundant lists.

mkdir combo | cd combo/

cp ../1_host/nucl_host/hybrid_HN/HN ./

cp ../1_host/prot_host/hybrid_HP/HP ./

cp ../2_myxo/nucl_myxo/hybrid_MN/MN ./

cp ../2_myxo/prot_myxo/hybrid_MP/MP ./

cat HN HP > H_all

cat MN MP > M_all

remove_duplicate.pl H_all  | mv duplicate_remove H

remove_duplicate.pl M_all | mv duplicate_remove M

list_compare.pl M H

#Process the comparison results.

mkdir 1_only_to_host  2_only_to_myxo  3_both_match  4_neither_match  5_delete_host

cd 1_only_to_host | cp ../../3_combo/H_only ./Host_only_list

cd ../2_only_to_myxo/ | cp ../../3_combo/M_only ./Myxo_only_list

cd ../3_both_match/ | cp ../../3_combo/inter_of_M_and_H ./both_match_list

cd ../4_neither_match | cp ../../3_combo/union_of_M_and_H ./

ln -s ~/transcriptome/3_1_assembly/M.wulii/WL_Unigene.fasta ./

remove_contaminant_by_ID.pl WL_Unigene.fasta union_of_M_and_H

mv survive.fasta M.wulii_neither_match.fasta

extract_fasta_header.pl M.wulii_neither_match.fasta

mv header neither_match_list

cd 5_delete_host

remove_contaminant_by_ID.pl /home/gqx/transcriptome/assembly/WL_Unigene.fasta ../1_only_to_host/Host_only_list

mv survive.fasta WL_host_delete_Unigene.fasta

2.6	Collected proteins for constructing bacterial no-redundant database

#Bacterial no-redundant database contains: Genome predicted proteins from, Aeromonas caviae, Aeromonas hydrophila, Aeromonas sobria, Aeromonas veronii B565, Citrobacter freundii CFNIH1, Escherichia coli, Flavobacterium columnare, Klebsiella pneumoniae, Lactobacillus sakei, Pseudomonas aeruginosa, Salmonella enterica, Shewanella putrefaciens, Plesiomonas shigelloides.

2.7	Remove redundancy in bacterial database and blast the host-seq-removed transcriptomes

cd-hit -i bac_protein.fasta -o bac_protein_cdhit.fasta -c 0.95 -n 5 -T 8

#45530 proteins retained in bacterial databases after redundancy removal.

makeblastdb -in bac_protein_cdhit.fasta -out BAC -dbtype prot -parse_seqids -hash_index

blastx -query WL_host_delete_Unigene.fasta -db BAC -out result_bac_e10 -evalue 1e-10 -outfmt 6 -num_threads 32

2.8	First round removement of bacterial and confirmation

cat result_bac_e10 | cut -f1 > 1

remove_duplicate.pl 1

mv duplicate_remove bacteria_contam_list

extract_seq_from_fasta.pl WL_host_delete_Unigene.fasta bacteria_contam_list

mv extracted.fasta bacteria_contam_first.fasta

#Confirm those sequences by blasting with nr database.

blastx -query bacteria_contam_first.fasta -db nr -out bacteria_contam_first_nr_result -evalue 1e-5 -max_target_seqs 1 -num_threads 56 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle frames sscinames sskingdoms"

grep "Bacteria" bacteria_contam_first_nr_result | cut -f 1 > 2

remove_duplicate.pl 2 | rm 2

mv duplicate_remove true_bacteria_contam.list

remove_contaminant_by_ID.pl WL_host_delete_Unigene.fasta true_bacteria_contam.list

#This step will produce the final clean transcriptome with host and bacterial contamination conservatively removed (WL_all_filter_Unigene.fasta), however, if the proportion of true bacterial contamination in bacteria_contam_list is too large, a second-round decontamination is recommended.

2.9	Second round removement of bacterial contamination

#Blast transcriptomes in last step with uniref90 (with diamond) and nt database (with megablast).

diamond makedb --in uniref90.fasta --db uniref90

diamond blastx -q WL_host_delete_Unigene.fasta --sensitive -k 20 -c 1 --threads 32 --db uniref90 --out diamond_result

#Add taxonomy information.

perl -lne 'BEGIN{open UT, "<uniref90.taxlist" or die $!; while (<UT>) { $ut{$1}=$2 if /^(\S+)\t(\S+)$/ } } {print "$_\t$ut{$1}" if /^\S+\t(\S+)/ and exists $ut{$1}}' < assembly_diamond_10.out > assembly_diamond_10.out_taxid
  
blastn -task megablast -query WL_host_delete_Unigene.fasta -db nt -culling_limit 5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle staxids sscinames sskingdoms' -num_threads 48 -evalue 1e-25 -out assembly_megablast_25.out
  
#Add taxonomy information.

perl -lne 'BEGIN{open UT, "<acc2tax_nucl_all.txt" or die $!; while (<UT>) { $ut{$2}=$3 if /^(\S+)\t(\S+)\t(\S+)/ } } {print "$_\t$ut{$1}" if /^\S+\t(\S+)/ and exists $ut{$1}}' < assembly_megablast_25.out > assembly_megablast_25.out_taxid
  
#Get the species distribution in transcriptomes.

cat ../diamond/assembly_diamond_10.out_taxid ../megablast/assembly_megablast_25.out_taxid > all

export LANG=C; export LC_ALL=C; sort -k1,1 -k12,12gr -k11,11g -k3,3gr all | sort -u -k1,1 --merge > bestHits

cat bestHits | rev | cut -f 1 | rev > 2

#Transfer taxid to species name by ete3.

tax2name.py > result_1

get_species_name_from_ete3.pl

species_distribution.pl names

#Redownload host and bacterial sequences according to the species distribution results.

cat ./* > fish_protein.fasta

cd-hit -i fish_protein.fasta -o fish_protein_cdhit.fasta -c 0.95 -n 5 -T 8

makeblastdb -in fish_protein_cdhit.fasta -out FISH -dbtype prot -parse_seqids -hash_index

cat ./* > bac_protein.fasta

cd-hit -i bac_protein.fasta -o bac_protein_cdhit.fasta -c 0.95 -n 5 -T 8

makeblastdb -in bac_protein_cdhit.fasta -out BAC -dbtype prot -parse_seqids -hash_index

#Reblast with the new databases.

blastx -query WL_all_filter_Unigene.fasta -db BAC -out result_bac_1 -evalue 1 -outfmt 6 -num_threads 48

blastx -query WL_all_filter_Unigene.fasta -db FISH -out result_fish_1 -evalue 1 -outfmt 6 -num_threads 48

#Extract the potiential contaminations.

cat result_bac_1 | cut -f1 > 1

remove_duplicate.pl 1

mv duplicate_remove bacteria_contam_list

cat result_fish_1 | cut -f1 > 1

remove_duplicate.pl 1

mv duplicate_remove fish_contam_list

extract_seq_from_fasta.pl WL_all_filter_Unigene.fasta bacteria_contam_list

mv extracted.fasta bacteria_contam_first.fasta

extract_seq_from_fasta.pl WL_all_filter_Unigene.fasta fish_contam_list

mv extracted.fasta fish_contam_first.fasta

#Blast extracted potiential sequences with nr database.

diamond makedb --in nr -d nr -p 24

diamond blastx -q fish_contam_first.fasta --sensitive -k 20 -c 1 --evalue 1e-5 --threads 48 --db nr.dmnd --out fish_diamond_5.out

diamond blastx -q bacteria_contam_first.fasta --sensitive -k 20 -c 1 --evalue 1e-5 --threads 48 --db nr.dmnd --out bacteria_diamond_5.out

export LANG=C; export LC_ALL=C; sort -k1,1 -k12,12gr -k11,11g -k3,3gr bacteria_diamond_5.out | sort -u -k1,1 --merge >> diamond_bestHits

export LANG=C; export LC_ALL=C; sort -k1,1 -k12,12gr -k11,11g -k3,3gr fish_diamond_5.out | sort -u -k1,1 --merge >> diamond_bestHits

#Add taxonomy information, this step requires software acc2tax.

give_tax_2_diamond_blastx.pl diamond_bestHits

#Extract contamination keys.

grep "Bacteria" diamond_blastx_with_tax > list_2

grep " Teleostomi" diamond_blastx_with_tax > list_2

cat bacteria/list_2 fish/list_2 > bad_list

cat bad_list | cut -f1 > 1

remove_duplicate.pl 1

mv duplicate_remove true_bad_list

#Get the second round decontamination results, WL_second_all_filter_Unigene.fasta.

remove_contaminant_by_ID.pl WL_all_filter_Unigene.fasta true_bad_list

3.Remove contamination from genomes by TAGC methods

3.1	Install Blobtools, details see https://blobtools.readme.io/docs

3.2	Process the genomes and get the mapping results

filter_fasta_by_length.pl WL_genome.fasta 200 200000 WL_genome_200.fasta

bowtie2-build WL_genome_200.fasta index --threads 8

bowtie2 -p 24 -x index -1 M.wulii_1.fq -2 M.wulii_2.fq -k 1 --very-fast-local -S out.sam

samtools view -bS out.sam > out.bam

3.3	Blast against NCBI Nucleotide database using megablast and against UniRef90 using diamond BLASTX

#Megablast:

blastn -task megablast -query WL_genome_200.fasta -db nt -culling_limit 5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle staxids sscinames sskingdoms' -num_threads 48 -evalue 1e-25 -out assembly_megablast_25.out
#Add taxonomy information.

perl -lne 'BEGIN{open UT, "<acc2tax_nucl_all.txt" or die $!; while (<UT>) { $ut{$2}=$3 if /^(\S+)\t(\S+)\t(\S+)/ } } {print "$_\t$ut{$1}" if /^\S+\t(\S+)/ and exists $ut{$1}}' < assembly_megablast_25.out > assembly_megablast_25.out_taxid
  
#Change the result style for Blobtools.

awk -v OFS="\t" -F"\t" '{print $1,$17,$12}' assembly_megablast_25.out_taxid  > mega.out

#Diamond:

diamond makedb --in uniref90.aa -d uniref90

diamond blastx -q WL_genome_200.fasta --sensitive -k 20 -c 1 --evalue 1e-10 --threads 48 --db uniref90.dmnd --out assembly_diamond_10.out

#Add taxonomy information.

perl -lne 'BEGIN{open UT, "<uniref90.taxlist" or die $!; while (<UT>) { $ut{$1}=$2 if /^(\S+)\t(\S+)$/ } } {print "$_\t$ut{$1}" if /^\S+\t(\S+)/ and exists $ut{$1}}' < assembly_diamond_10.out > assembly_diamond_10.out_taxid
  
#Change the result style for Blobtools.

awk -v OFS="\t" -F"\t" '{print $1,$13,$12}' assembly_diamond_10.out_taxid  > diamond.out

3.4	Process above results by the blobtools script to annotate each scaffold

#Sort the blast results.

cat mega.out diamond.out > blast.out

sort_blast_by_query_name.pl blast.out

mv sorted_output blast.out

#Create blob database.

python2.7 blobtools create -i WL_genome_200.fasta -b out.bam -t blast.out -o M.wulii_1_blob --names names.dmp --nodes nodes.dmp

#Process the blob database.

python2.7 ../blobtools view -i M.wulii_1_blob.blobDB.json -o ./

3.5	Visualize the annotation results into Taxon-Annotated-Gc-Coverage plot (TAGC)

python2.7 blobtools blobplot -i M.wulii_1_blob.blobDB.json -o ./ --format pdf --colours colours.txt

3.6	Potential contamination is inspected manually and compared against NCBI Nucleotide database

#Reformat the result.

format.sh > result

#Inspect scaffolds with a bit-score ≥200.

blob_result_seq_extract.pl M.wulii _1_blob.blobDB.table.txt

#Extract potential contamination sequences.

extract_seq_from_fasta.pl WL_genome_200.fasta seq_for_blast

mv extracted.fasta contam_candidate.fa

#Blast nt database.

blastn -query contam_candidate.fa -db nt -evalue 1e-5 -max_target_seqs 20 -num_threads 24 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle frames sscinames sskingdoms" -out nt_result

#Add taxonomy information.

export LANG=C; export LC_ALL=C; sort -k1,1 -k12,12gr -k11,11g -k3,3gr nt_result | sort -u -k1,1 --merge > bestHits

cat bestHits | cut -f2 > acc

perl -p -i -e 's/\.(\d)//g' acc

acc2tax -i acc -o result -d accession2taxid

#Extract contamination sequences according to customized key words and get the accession number.

Teleostomi_Bacteria_extract.pl result

#Get the contamination sequence header (final_contam_header) according to the accession number.

cat bestHits | cut -f1,2 > header

perl -p -i -e 's/\.(\d)//g' header

get_contam_from_accesion.pl contam_accession header

#Also add sequences whose bam0 lower than 20 (seq_remove_by_bam0) to removal lists.

another_blob_result_seq_extract.pl M.wulii_1_blob.blobDB.table.txt

cat seq_remove_by_bam0 final_contam_header > true_bad_list

remove_duplicate.pl true_bad_list

mv duplicate_remove true_bad_list

3.7 Remove contamination in genomes

remove_contaminant_by_ID.pl WL_genome_200.fasta true_bad_list

mv survive.fasta genome.fasta

4.Genome gene prediction

4.1 GeneMark-ET provides species HMM file to Augustus for training

#Map transcriptome reads to genome.

tophat -o TophatOutput -p 4 --no-novel-juncs ~/transcriptome/myxobolus/carp_remov/index/genome ~/transcriptome/myxobolus/processing/combinedF.fastq ~/transcriptome/myxobolus/processing/combinedR.fastq

#Make hints.

bet_to_gff.pl --bed junctions.bed -gff introns.gff --label tophat2 --seq genome.fasta

gmes_petap.pl --sequence genome.fasta --ET introns.gff --et_score 10 --cores 4 --min_contig 2000

/opt/biosoft/PASApipeline-2.0.2/misc_utilities/gtf_to_gff3_format.pl genemark.gtf genome.fasta > genemark.gff3

#Filter good models from GeneMark-ET results.

filterGenemark.pl genemark.gtf introns.gff

/opt/biosoft/PASApipeline-2.0.2/misc_utilities/gtf_to_gff3_format.pl genemark.f.good.gtf genome.fasta >genemark.f.good.gff3
mv genemark.f.good.gff3 best_candidates.gff3

#Extract protein sequences.

/opt/biosoft/EVidenceModeler-1.1.1/EvmUtils/gff3_file_to_proteins.pl best_candidates.gff3 genome.fasta prot > best_candidates.fasta

remove_redundant_high_identity_genes.pl best_candidates.gff3 best_candidates.fasta 4 0.70 > best_candidates.lowIdentity.gff3 2> remove_redundant_high_identity_genes.log

#best_candidates.lowIdentity.gff3 can be used for augustus training.

4.2 Run Augustus gene prediction

#Change the style.

gff2gbSmallDNA.pl best_candidates.lowIdentity.gff3 genome.fasta 800 genes.raw.gb

#Remove the bad models in genes.raw.gb.

new_species.pl --species=for_bad_genes_removing

etraining --species=for_bad_genes_removing --stopCodonExcludedFromCDS=false genes.raw.gb 2> train.err

cat train.err | perl -pe 's/.*in sequence (\S+): .*/$1/' > badgenes.lst

filterGenes.pl badgenes.lst genes.raw.gb > genes.gb

#First training with training set.

randomSplit.pl genes.gb 100

new_species.pl --species=myxobolus_wulii

etraining --species=myxobolus_wulii genes.gb.train > train.out

#First gene prediction with test set.

augustus --species=myxobolus_wulii genes.gb.test | tee firsttest.out

#Optimize the parameters.

optimize_augustus.pl --species=myxobolus_wulii --cpus=8 genes.gb.train

#Second training with training set.

etraining --species=myxobolus_wulii genes.gb.train

#Second gene prediction with test set.

augustus --species=myxobolus_wulii genes.gb.test | tee secondtest.out

#Mask genome with RepeatMasker and RepeatModeler.

mkdir repeatMasker

cd repeatMasker

fasta_no_blank.pl genome.fasta > genome2.fasta

rm genome.fasta

mv genome2.fasta genome.fasta

RepeatMasker -pa 24 -e ncbi -species cnidaria -gff -dir repeatMasker genome.fasta

cd ..

mkdir repeatModeler

cd repeatModeler

/opt/biosoft/RepeatModeler/BuildDatabase -name wulii -engine ncbi genome.fasta

/opt/biosoft/RepeatModeler/RepeatModeler -database wulii -pa 8

/opt/biosoft/RepeatMasker/RepeatMasker -pa 4 -e ncbi -lib RM_37443.MonJan231456212017/consensi.fa.classified -dir ./ -gff genome.fasta

#Combine RepeatMasker and RepeatModeler results.

cd ..

merge_repeatMasker_out.pl repeatMasker/genome.fasta.out repeatModeler/genome.fasta.out > genome.repeat.stats

maskedByGff.pl genome.repeat.gff3 genome.fasta hardmaskN > genome.hardmaskN.fasta

#Map RNA-seq reads to masked genome.

mv genome.hardmaskN.fasta genome_db.fa

bowtie2-build genome_db.fa genome_db --threads 8

tophat2 -N 3 --read-edit-dist 3 -p 32 -i 20 -I 4000 --min-segment-intron 20 --max-segment-intron 4000 --min-coverage-intron 20 --max-coverage-intron 4000 --coverage-search --microexon-search -o result genome_db M.wulii_1_fixed.fastq M.wulii_2_fixed.fastq

#Make RNA-seq hints.

bam2hints --intronsonly --in=result/accepted_hits.bam --out=hints.gff

#Formal gene prediction.

augustus --species=myxobolus_wulii_1 --extrinsicCfgFile=extrinsic.cfg --alternatives-from-evidence=true --allow_hinted_splicesites=atac --hintsfile=hints.gff --gff3=on genome.fasta > aug.gff3

perl -p -i -e 's/\ttranscript\t/\tmRNA\t/' aug.gff3

4.3 GeneMark-ET provides species HMM file to SNAP for training 

#Use some functions in Maker.

maker2zff genemark.gff3

extract_header_for_snap.pl genome.dna

#Establish index.

fastaindex genome.fasta genome.idx

fastafetch -f genome.fasta -i genome.idx -Fq <(sort -u header) > out

mv out genome.dna

#Breakdown the genome.

fathom -categorize 1000 genome.ann genome.dna

fathom uni.ann uni.dna -export 1000 -plus

mkdir params

cd params/

forge ../export.ann ../export.dna

cd ..

hmm-assembler.pl species params/ > species.hmm

#SNAP training finish.

4.4 Run SNAP gene prediction

snap species.hmm genome.fasta -gff -quiet > snap.gff

snap2gff3.pl snap.gff > snap_ture.gff

4.5 Homology-based gene prediction by Genewise

#Prepare non-redundant close-related proteins.

rename_fasta_by_numeber.pl all_cdhit.fa

#Start the annotation.

/opt/biosoft/homolog_genewise/homolog_genewise.pl rename_all.fasta genome.hardmaskN.fasta 8 0.1 1e-9

#Filter the result.

/opt/biosoft/homolog_genewise/genewise_filter.pl genewise.gff genome.hardmaskN.fasta 15 90 1 1e-6 0.30 4 > genewise.filter.gff 2> genewise.filter.stats

#Evaluate the completeness of predicted genes and filter the sequences containing stop codons.

/opt/biosoft/homolog_genewise/genewise2EVM_input.pl genewise.filter.gff genome.hardmaskN.fasta filterMiddleStopCodon=yes > evm_protein_alignment.gff3 2> genewise_gene_models_completeness_check.txt

#Process the output style

perl -p -i -e 's/^#.*//; s/^\s*$//' genewise.gff

make_evm_recognize_gff.pl genewise.gff

/opt/biosoft/EVidenceModeler-1.1.1/EvmUtils/misc/SNAP_to_GFF3.pl evm_wise.gff3 > final_evm_wise.gff3

4.6 Gene prediction by PASA

#Prepare the transcriptome file.

perl -e 'while (<>) { print "$1\n" if />(\S+)/ }' /home/train/00.incipient_data/data_for_gene_prediction_and_RNA-seq/Trinity.fasta > tdn.accs

#End-trimming of the transcriptome.

seqclean Trinity.fasta -v /opt/biosoft/PASApipeline-2.0.2/seqclean/UniVec

#Produce config file.

cp /opt/biosoft/PASApipeline-2.0.2/pasa_conf/pasa.alignAssembly.Template.txt alignAssembly.config

DATE=`date +%Y%m%e%k%M%S | perl -pe 's/\s+//'`

echo "perl -p -i -e 's/MYSQLDB=.*/MYSQLDB=pasa_$DATE/' alignAssembly.config" | sh

#Produce mysql databases and tables.

/opt/biosoft/PASApipeline-2.0.2/scripts/create_mysql_cdnaassembly_db.dbi -r -c alignAssembly.config -S /opt/biosoft/PASApipeline-2.0.2/schema/cdna_alignment_mysqlschema

#Start mapping transcripts.

/opt/biosoft/PASApipeline-2.0.2/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -R -g genome.fasta -t Trinity.fasta.clean -T -u Trinity.fasta --ALIGNERS gmap,blat --CPU 8 --stringent_alignment_overlap 30.0 --TDN tdn.accs --MAX_INTRON_LENGTH 20000 --TRANSDECODER &> pasa.log

#This will produce pasa_*.pasa_assemblies.gff3 for downstream analyses.

4.7 Combine above gff3 file by EVM

#Augustus gff file.

/opt/biosoft/EVidenceModeler-1.1.1/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl aug.gff3 > evm_augustus.gff3

perl -p -i -e 's/^#.*//; s/^\s*$//' evm_augustus.gff3

gff3_gene_prediction_file_validator.pl evm_augustus.gff3

#Genemar-ET gff file.

ln -s ../genemark-et/genemark.gff3 evm_genemark-et.gff3

gff3_gene_prediction_file_validator.pl evm_genemark-et.gff3

#Snap gff file.

#Change the style.

/opt/biosoft/EVidenceModeler-1.1.1/EvmUtils/misc/SNAP_output_to_gff3.pl snap.zff genome.fasta > snap.gff3

gff3_gene_prediction_file_validator.pl snap.gff3

perl -p -i -e 's/^(\S+)\t(\.)\t/$1\tSNAP\t/g' snap.gff3

#PASA gff file.

cp ../new_pasa/pasa*.pasa_assemblies.gff3 ./transcript_alignments.gff3

perl -p -i -e 's/\t\S+/\tpasa_transcript_alignments/' transcript_alignments.gff3

gff3_gene_prediction_file_validator.pl transcript_alignments.gff3

#Genewise gff file.

ln -s ../genewise/wise/evm_protein_alignment.gff3 ./protein_alignments.gff3

gff3_gene_prediction_file_validator.pl protein_alignments.gff3

cat evm_augustus.gff3 evm_genemark-et.gff3 snap.gff3 | perl -pe 's/^#.*//; s/^\s*$//' > gene_predictions.gff3

#Create weights file.

echo -e "ABINITIO_PREDICTION\tAugustus\t6

ABINITIO_PREDICTION\tSNAP\t2

ABINITIO_PREDICTION\tGeneMark.hmm\t1

PROTEIN\tGeneWise\t5

TRANSCRIPT\tpasa_transcript_alignments\t10" > weights.txt

#Split the data.

partition_EVM_inputs.pl --genome genome.fasta --gene_predictions gene_predictions.gff3 --protein_alignments protein_alignments.gff3 --transcript_alignments transcript_alignments.gff3 --repeats genome.repeat.gff3 --segmentSize 500000 --overlapSize 10000 --partition_listing partitions_list.out

#Create and run commands.

write_EVM_commands.pl --genome genome.fasta --gene_predictions gene_predictions.gff3 --protein_alignments protein_alignments.gff3 --transcript_alignments transcript_alignments.gff3 --repeats genome.repeat.gff3  --weights `pwd`/weights.txt --partitions partitions_list.out --output_file_name evm.out > commands.list

ParaFly -c commands.list -CPU 4

#Combine the results and reformat to gff3.

recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out

convert_EVM_outputs_to_GFF3.pl --partitions partitions_list.out --output_file_name evm.out --genome genome.fasta

find . -regex ".*evm.out.gff3" -exec cat {} \; > EVM.all.gff3

#Extract proteins from gff3.

/opt/biosoft/EVidenceModeler-1.1.1/EvmUtils/gff3_file_to_proteins.pl EVM.all.gff3 ../evm/genome.fasta prot > genome_protein.fasta

5.Transcriptome gene prediction

5.1 De novo prediction by TransDecoder

#Predict ORF.

TransDecoder.LongOrfs -t WL_second_all_filter_Unigene.fasta -m 20

TransDecoder.Predict -t WL_second_all_filter_Unigene.fasta --cpu 12

#Get the Transdecoder result WL_second_all_filter_Unigene.fasta.transdecoder.pep.

5.2 De novo prediction by GeneMarkS-T

/opt/biosoft/GeneMarkS-T/gmst.pl --output M.wulii_gmst --fnn --faa -clean 1 WL_second_all_filter_Unigene.fasta

#Get the GeneMarkS-T result M.wulii_gmst.fasta.

5.3 Homolog-based prediction by Hercules (https://github.com/qingxiangguo/hercules-v.1.0)

#Blast transcriptome with nr, swiss, egg and kog database.

blastx -query WL_second_all_filter_Unigene.fasta -db nr -max_target_seqs 20 -out nr_result -evalue 1e-5  -num_threads 48 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle frames sscinames sskingdoms"

blastx -query WL_second_all_filter_Unigene.fasta -db kog -max_target_seqs 20 -out kog_result -evalue 1e-5  -num_threads 48 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle frames sscinames sskingdoms"

blastx -query WL_second_all_filter_Unigene.fasta -db swiss -max_target_seqs 20 -out swiss_result -evalue 1e-5  -num_threads 48 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle frames sscinames sskingdoms"

blastx -query WL_second_all_filter_Unigene.fasta -db eggnog -max_target_seqs 20 -out eggnog_result -evalue 1e-5  -num_threads 48 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle frames sscinames sskingdoms"

#Process blast results.

for next in $(cut -f1 nr_result | sort -u); do grep -w -m 20 "$next" nr_result; done > tmp

sort_blast_by_query_name.pl tmp

mv sorted_output nr | rm tmp

for next in $(cut -f1 eggnog_result | sort -u); do grep -w -m 20 "$next" eggnog_result; done > tmp

sort_blast_by_query_name.pl tmp

mv sorted_output eggnog | rm tmp

for next in $(cut -f1 swiss_result | sort -u); do grep -w -m 20 "$next" swiss_result; done > tmp

sort_blast_by_query_name.pl tmp

mv sorted_output swiss | rm tmp

for next in $(cut -f1 kog_result | sort -u); do grep -w -m 20 "$next" kog_result; done > tmp

sort_blast_by_query_name.pl tmp

mv sorted_output kog | rm tmp

hercules nr swiss egg kog

#Hercules requires the installation of anvi’o, for details, see https://github.com/merenlab/anvio/releases/v5.5.

#Hercules will combine blast results from different databases with customized priority and produce a gene-call file which can be fed into Anvi'o.

anvi-gen-contigs-database -f contigs.fa -o contigs.db --external-gene-calls gene_call

#Get the proteins.

anvi-get-aa-sequences-for-gene-calls -c contigs.db -o homolog.fa

5.4 Six-frame translation

#Those transcripts that were translated neither by de novo nor homolog-based method are translated into amino acid sequences using the Transeq script from the EMBOSS.

#Get those sequences that were translated neither by de novo nor homolog-based method.

cp ../trans_homolog/gene_call ./

grep ">" ../trans_denovo/M.wulii_gmst.faa > gmst

grep ">" ../trans_denovo/WL_second_all_filter_Unigene.fasta.transdecoder.pep > transdecoder

cat gene_call | cut -f2 > gene_call_header

perl -p -i -e "s/\s.*//" transdecoder

perl -p -i -e "s/>Gene.\d+:://" transdecoder

perl -p -i -e "s/:.*//" transdecoder

perl -p -i -e "s/\s.*//" gmst

perl -p -i -e "s/>//" gmst

cat gene_call_header gmst transdecoder > annotation_list

remove_duplicate.pl annotation_list

rm annotation_list

mv duplicate_remove annotation_list

remove_contaminant_by_ID.pl WL_second_all_filter_Unigene.fasta annotation_list

mv survive.fasta trans_left.fasta

#Six-frame translation using EMBOSS.

transeq -sequence trans_left.fasta -outseq out -frame 6

flat_the_fasta_seq.pl out

rm out

#Collect the amino acid sequence between stop codons (represented by *) and discard the sequence contained ambiguous amino acids (represented by X) or its total length is less than 30.

get_seq_between_asterisk.pl flated 30

mv between_asterisk trans_left.pep

6. Combine all the proteins predicted from genomes and transcriptomes

6.1 Collect all predicted proteins

cat genome_protein.fasta WL_second_all_filter_Unigene.fasta.transdecoder.pep M.wulii_gmst.fasta homolog.fa trans_left.pep > all.fasta

6.2 Remove redundancy and filter by length

cd-hit -i all.fasta -o all_cdhit.fasta -c 1 -T 4 -M 0

filter_fasta_by_length.pl all_cdhit.fasta 30 1000000 filtered.fasta

6.3 Give name tag to CCPRD

replace_fasta_header_by_number.pl filtered.fasta WL_MC

mv ordered.fasta CCPRD

#Get the CCPRD here.

7. Create alternative databases for comparison

7.1 Transcriptome six-frame translation

transeq -sequence WL_second_all_filter_Unigene.fasta -outseq out -frame 6 -clean

remove_no_end_asterisks.pl out

flat_the_fasta_seq.pl out

rm out

get_seq_between_asterisk.pl flated 30

cd-hit -i between_asterisk -o out -c 1 -M 160000 -T 8

replace_fasta_header_by_number.pl out WL_T6

mv ordered.fasta trans_6_frame

7.2 Genome and transcriptome six-frame translation

cat genome.fasta WL_second_all_filter_Unigene.fasta > all.fasta

transeq -sequence all.fasta -outseq out -frame 6 -clean

remove_no_end_asterisks.pl out

flat_the_fasta_seq.pl out

rm out

get_seq_between_asterisk.pl flated 30

cd-hit -i between_asterisk -o out -c 1 -M 160000 -T 8

replace_fasta_header_by_number.pl out WL_A6

7.3 CCPRD + contaminants

#Add artificial host and bacteria sequences to CCPRD, and give them the name tag WL_BAC, WL_HOST respectively.

7.4 CCPRD + sequences removed in decontamination process

#Add back sequences that were removed during the de-contamination process, and give them the name tag WL_RM.

# Examples

#The CCPRD workflow was tested in a model organism, Saccharomyces cerevisiae and provided as example file. The public data of S. cerevisiae S288C were used, which included Ensembl annotated version of S288c genome R64-1-1, GenBank GCA_000146045.2, and transcriptome shotgun assembly, GenBank GFJR01000000. Since the genome of S288c is a gold-standarded data, here we skipped the genome decontamination process by TAGC methods.

1.Download and process the transcriptomes and genomes

1.1 Download the source data

mkdir 1_source_data && cd 1_source_data

wget  https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GF/JR/GFJR01/GFJR01.1.fsa_nt.gz

gunzip GFJR01.1.fsa_nt.gz

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz

gunzip GCF_000146045.2_R64_genomic.fna.gz

1.2 Remove redundancy in transcriptomes with CD-HIT

cd ..

mkdir 2_cluster _mRNA && cd 2_cluster _mRNA

ln -s ../1_source_data/GFJR01_mRNA.fasta ./

cd-hit-est -i GFJR01_mRNA.fasta -o yeast.fasta_cdhit -c 0.95 -n 10 -d 0 -M 16000 -T 8

1.3 Further cluster the transcriptomes with TGICL

tgicl -l 40 -c 10 -v 25 -O '-repeat_stringency 0.95 -minmatch 35 -minscore 35' -F yeast.fasta_cdhit

1.4 Process the cd-hit results to get unigenes

cat asm_1/contigs asm_2/contigs asm_3/contigs asm_4/contigs asm_5/contigs asm_6/contigs asm_7/contigs asm_8/contigs asm_9/contigs asm_10/contigs > all.contigs

fast_extract_seq_from_fasta.pl yeast.fasta_cdhit yeast.fasta_cdhit.singletons > extracted.fasta

Unigene_generator.pl -s extracted.fasta -c all.contigs -t YT

#This will output YT_Unigene.fasta for downstream analyses.

2.Remove potential bacterial contamination in transcriptome data

#We use the same bacterial database used in above Pipeline section 2.6-2.9

2.1 Collected proteins for constructing bacterial no-redundant database

#Bacterial no-redundant database contains: Genome predicted proteins from, Aeromonas caviae, Aeromonas hydrophila, Aeromonas sobria, Aeromonas veronii B565, Citrobacter freundii CFNIH1, Escherichia coli, Flavobacterium columnare, Klebsiella pneumoniae, Lactobacillus sakei, Pseudomonas aeruginosa, Salmonella enterica, Shewanella putrefaciens, Plesiomonas shigelloides.

2.2 Remove redundancy in bacterial database and blast the transcriptomes

cd ..

mkdir 3_decontam_mRNA

cd 3_decontam_mRNA/

cd-hit -i bac_protein.fasta -o bac_protein_cdhit.fasta -c 0.95 -n 5 -T 8

#45530 proteins retained in bacterial databases after redundancy removal.

makeblastdb -in bac_protein_cdhit.fasta -out BAC -dbtype prot -parse_seqids -hash_index

blastx -query ../../2_cluster_mRNA/YT_Unigene.fasta -db BAC -out result_bac_e10 -evalue 1e-10 -outfmt 6 -num_threads 32

2.3 First round removement of bacterial and confirmation

cat result_bac_e10 | cut -f1 > 1

remove_duplicate.pl 1

mv duplicate_remove bacteria_contam_list

extract_seq_from_fasta.pl ../../2_cluster_mRNA/YT_Unigene.fasta bacteria_contam_list

mv extracted.fasta bacteria_contam_first.fasta

#Confirm those sequences by blasting with nr database.

blastx -query bacteria_contam_first.fasta -db nr -out bacteria_contam_first_nr_result -evalue 1e-5 -max_target_seqs 1 -num_threads 56 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle frames sscinames sskingdoms"

grep "Bacteria" bacteria_contam_first_nr_result | cut -f 1 > 2

remove_duplicate.pl 2 | rm 2

mv duplicate_remove true_bacteria_contam.list

remove_contaminant_by_ID.pl WL_host_delete_Unigene.fasta true_bacteria_contam.list

#This step will produce the final clean transcriptome with bacterial contamination conservatively removed, however, if the proportion of true bacterial contamination in bacteria_contam_list is too large, a second-round decontamination is recommended. In present case, we are ok with the yeast transcriptome, almost no bacteria contamination was detected.

3.Genome gene prediction

#Since S. cerevisiae is a model organism with good species HMM file, we used the HMM file in prediction software instead of making our own.

mkdir 4_MS_database

cd 4_MS_database

mkdir all_6_frame  CCPRD  CCPRD_contam  trans_6_frame

cd 4_MS_database/CCPRD

mkdir genome_denovo  genome_homolog  genome_pasa  trans_denovo  trans_homolog  trans_left_six_frame

3.1 GeneMark-ES prediction

cd 4_MS_database/CCPRD/ genome_denovo/GeneMark-ES

ln -s ~/examples/1_source_data/GCF_000146045.2_R64_genomic.fna ./genome.fasta

gmes_petap.pl –sequence genome.fasta –ES –fungus –cores 4

gtf_to_gff3_format.pl genemark.gtf genome.fasta > genemark.gff3

3.2 Augustus gene prediction

cd 4_MS_database/CCPRD/ genome_denovo/AUGUSTUS

ln -s ~/examples/1_source_data/GCF_000146045.2_R64_genomic.fna ./genome.fasta

augustus --species=saccharomyces_cerevisiae_S288C --gff3=on genome.fasta > aug.gff3

3.3 SNAP gene prediction

cd 4_MS_database/CCPRD/ genome_denovo/SNAP

ln -s ../GeneMark-ES/genome.fasta ./genome.fasta

ln -s ../GeneMark-ES/genemark.gff3 ./

#Use some functions in Maker.

maker2zff genemark.gff3

extract_header_for_snap.pl genome.dna

#Establish index.

fastaindex genome.fasta genome.idx

fastafetch -f genome.fasta -i genome.idx -Fq <(sort -u header) > out

mv out genome.dna

#Breakdown the genome.

fathom -categorize 1000 genome.ann genome.dna

fathom uni.ann uni.dna -export 1000 -plus

mkdir params

cd params/

forge ../export.ann ../export.dna

cd ..

hmm-assembler.pl species params/ > species.hmm

#SNAP training finish and start prediction

snap species.hmm genome.fasta -gff -quiet > snap.gff

snap2gff3.pl snap.gff > snap_true.gff

3.4 Homology-based gene prediction by Genewise

#We use the S. cerevisiae proteins in  UniProtKB/Swiss-Prot as close-related proteins.

#Start the annotation.

/opt/biosoft/homolog_genewise/homolog_genewise.pl rename_all.fasta genome.hardmaskN.fasta 8 0.1 1e-9

#Filter the result.

/opt/biosoft/homolog_genewise/genewise_filter.pl genewise.gff genome.hardmaskN.fasta 15 90 1 1e-6 0.30 4 > genewise.filter.gff 2> genewise.filter.stats

#Evaluate the completeness of predicted genes and filter the sequences containing stop codons.

/opt/biosoft/homolog_genewise/genewise2EVM_input.pl genewise.filter.gff genome.hardmaskN.fasta filterMiddleStopCodon=yes > evm_protein_alignment.gff3 2> genewise_gene_models_completeness_check.txt

#Process the output style

perl -p -i -e 's/^#.*//; s/^\s*$//' genewise.gff

make_evm_recognize_gff.pl genewise.gff

/opt/biosoft/EVidenceModeler-1.1.1/EvmUtils/misc/SNAP_to_GFF3.pl evm_wise.gff3 > final_evm_wise.gff3

3.5 Gene prediction by PASA

cd 4_MS_database/CCPRD/genome_pasa

#Prepare the transcriptome file.

perl -e 'while (<>) { print "$1\n" if />(\S+)/ }' YT_Unigene.fasta > tdn.accs

#End-trimming of the transcriptome.

seqclean YT_Unigene.fasta -v /opt/biosoft/PASApipeline-2.0.2/seqclean/UniVec

#Produce config file.

cp /opt/biosoft/PASApipeline-2.0.2/pasa_conf/pasa.alignAssembly.Template.txt alignAssembly.config

DATE=`date +%Y%m%e%k%M%S | perl -pe 's/\s+//'`

echo "perl -p -i -e 's/MYSQLDB=.*/MYSQLDB=pasa_$DATE/' alignAssembly.config" | sh

#Produce mysql databases and tables.

/opt/biosoft/PASApipeline-2.0.2/scripts/create_mysql_cdnaassembly_db.dbi -r -c alignAssembly.config -S /opt/biosoft/PASApipeline-2.0.2/schema/cdna_alignment_mysqlschema

#Start mapping transcripts.

/opt/biosoft/PASApipeline-2.0.2/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -R -g ../genome_denovo/GeneMark-ES/genome.fasta -t YT_Unigene.fasta.clean -T -u ../../../2_cluster_mRNA/YT_Unigene.fasta --ALIGNERS gmap,blat --CPU 8 --stringent_alignment_overlap 30.0 --TDN tdn.accs --MAX_INTRON_LENGTH 20000 --TRANSDECODER &> pasa.log

#This will produce pasa_*.pasa_assemblies.gff3 for downstream analyses.

3.6 Combine above gff3 file by EVM

#Augustus gff file.

/opt/biosoft/EVidenceModeler-1.1.1/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl aug.gff3 > evm_augustus.gff3

perl -p -i -e 's/^#.*//; s/^\s*$//' evm_augustus.gff3

gff3_gene_prediction_file_validator.pl evm_augustus.gff3

#Genemar-ES gff file.

ln -s 4_MS_database/CCPRD/genome_denovo/GeneMark-ES/genemark.gff3 evm_genemark-es.gff3

gff3_gene_prediction_file_validator.pl evm_genemark-es.gff3

#Snap gff file.

#Change the style.

/opt/biosoft/EVidenceModeler-1.1.1/EvmUtils/misc/SNAP_output_to_gff3.pl snap.zff genome.fasta > snap.gff3

gff3_gene_prediction_file_validator.pl snap.gff3

perl -p -i -e 's/^(\S+)\t(\.)\t/$1\tSNAP\t/g' snap.gff3

#PASA gff file.

cp ../new_pasa/pasa*.pasa_assemblies.gff3 ./transcript_alignments.gff3

perl -p -i -e 's/\t\S+/\tpasa_transcript_alignments/' transcript_alignments.gff3

gff3_gene_prediction_file_validator.pl transcript_alignments.gff3

#Genewise gff file.

ln -s ../genewise/wise/evm_protein_alignment.gff3 ./protein_alignments.gff3

gff3_gene_prediction_file_validator.pl protein_alignments.gff3

cat evm_augustus.gff3 evm_genemark-es.gff3 snap.gff3 | perl -pe 's/^#.*//; s/^\s*$//' > gene_predictions.gff3

#Create weights file.

echo -e "ABINITIO_PREDICTION\tAugustus\t6

ABINITIO_PREDICTION\tSNAP\t2

ABINITIO_PREDICTION\tGeneMark.hmm\t1

PROTEIN\tGeneWise\t5

TRANSCRIPT\tpasa_transcript_alignments\t10" > weights.txt

#Split the data.

partition_EVM_inputs.pl --genome genome.fasta --gene_predictions gene_predictions.gff3 --protein_alignments protein_alignments.gff3 --transcript_alignments transcript_alignments.gff3 --repeats genome.repeat.gff3 --segmentSize 500000 --overlapSize 10000 --partition_listing partitions_list.out

#Create and run commands.

write_EVM_commands.pl --genome genome.fasta --gene_predictions gene_predictions.gff3 --protein_alignments protein_alignments.gff3 --transcript_alignments transcript_alignments.gff3 --repeats genome.repeat.gff3  --weights `pwd`/weights.txt --partitions partitions_list.out --output_file_name evm.out > commands.list

ParaFly -c commands.list -CPU 4

#Combine the results and reformat to gff3.

recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out

convert_EVM_outputs_to_GFF3.pl --partitions partitions_list.out --output_file_name evm.out --genome genome.fasta

find . -regex ".*evm.out.gff3" -exec cat {} \; > EVM.all.gff3

#Extract proteins from gff3.

/opt/biosoft/EVidenceModeler-1.1.1/EvmUtils/gff3_file_to_proteins.pl EVM.all.gff3 ../evm/genome.fasta prot > genome_protein.fasta

4.Transcriptome gene prediction

4.1 De novo prediction by TransDecoder

#Predict ORF.

cd 4_MS_database/CCPRD/trans_denovo

TransDecoder.LongOrfs -t ../../../2_cluster_mRNA/YT_Unigene.fasta -m 20

TransDecoder.Predict -t ../../../2_cluster_mRNA/YT_Unigene.fasta --cpu 12

#Get the Transdecoder result YT_Unigene.fasta.transdecoder.pep

4.2 De novo prediction by GeneMarkS-T

gmst.pl --output YT_gmst --fnn --faa -clean 1 ../../../2_cluster_mRNA/YT_Unigene.fasta

4.3 Homolog-based prediction by Hercules (https://github.com/qingxiangguo/hercules-v.1.0)

cd 4_MS_database/CCPRD/trans_homolog

#Blast transcriptome with swiss and kog database.

blastx -query ../../../2_cluster_mRNA/YT_Unigene.fasta -db KOG -max_target_seqs 20 -out kog_result -evalue 1e-5 -num_threads 48 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle frames sscinames sskingdoms"
blastx -query ../../../2_cluster_mRNA/YT_Unigene.fasta -db swiss -max_target_seqs 20 -out swiss_result -evalue 1e-5 -num_threads 48 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle frames sscinames sskingdoms"

#Process blast results.

for next in $(cut -f1 swiss_result | sort -u); do grep -w -m 20 "$next" swiss_result; done > tmp

sort_blast_by_query_name.pl tmp

mv sorted_output swiss | rm tmp

for next in $(cut -f1 kog_result | sort -u); do grep -w -m 20 "$next" kog_result; done > tmp

sort_blast_by_query_name.pl tmp

mv sorted_output kog | rm tmp

hercules nr swiss egg kog

#Hercules requires the installation of anvi’o, for details, see https://github.com/merenlab/anvio/releases/v5.5.

#Hercules will combine blast results from different databases with customized priority and produce a gene-call file which can be fed into Anvi'o.

anvi-gen-contigs-database -f ../../../2_cluster_mRNA/YT_Unigene.fasta -o contigs.db --external-gene-calls gene_call

#Get the proteins.

anvi-get-aa-sequences-for-gene-calls -c contigs.db -o homolog.fa

4.4 Six-frame translation

#Those transcripts that were translated neither by de novo nor homolog-based method are translated into amino acid sequences using the Transeq script from the EMBOSS.

cd 4_MS_database/CCPRD/trans_left_six_frame

#Get those sequences that were translated neither by de novo nor homolog-based method.

cp ../trans_homolog/gene_call ./

grep ">" ../trans_denovo/YT_gmst.faa > gmst

grep ">" ../trans_denovo/YT_Unigene.fasta.transdecoder.pep > transdecoder

cat gene_call | cut -f2 > gene_call_header

perl -p -i -e "s/\s.*//" transdecoder

perl -p -i -e "s/>Gene.\d+:://" transdecoder

perl -p -i -e "s/:.*//" transdecoder

perl -p -i -e "s/\s.*//" gmst

perl -p -i -e "s/>//" gmst

cat gene_call_header gmst transdecoder > annotation_list

remove_duplicate.pl annotation_list

rm annotation_list

mv duplicate_remove annotation_list

remove_contaminant_by_ID.pl ../../../2_cluster_mRNA/YT_Unigene.fasta annotation_list

mv survive.fasta trans_left.fasta

#Six-frame translation using EMBOSS.

transeq -sequence trans_left.fasta -outseq out -frame 6

flat_the_fasta_seq.pl out

rm out

get_seq_between_asterisk.pl flated 30

mv between_asterisk trans_left.pep

5.Combine all the proteins predicted from genomes and transcriptomes

5.1 Collect all predicted proteins

cd 4_MS_database/CCPRD/final

cat genome_protein.fasta YT_Unigene.fasta.transdecoder.pep YT_gmst.fasta homolog.fa trans_left.pep > all.fasta

5.2 Remove redundancy and filter by length

cd-hit -i all.fasta -o all_cdhit.fasta -c 1 -T 4 -M 0

filter_fasta_by_length.pl all_cdhit.fasta 30 1000000 filtered.fasta

5.3 Give name tag to CCPRD

replace_fasta_header_by_number.pl filtered.fasta YT_CC

mv ordered.fasta CCPRD

#Get the CCPRD here.

6.Create alternative databases for comparison

6.1 Transcriptome six-frame translation

cd 4_MS_database/ trans_6_frame

transeq -sequence ../../2_cluster_mRNA/YT_Unigene.fasta -outseq out -frame 6 

remove_no_end_asterisks.pl out

flat_the_fasta_seq.pl out

rm out

get_seq_between_asterisk.pl flated 30

cd-hit -i between_asterisk -o out -c 1 -M 160000 -T 8

replace_fasta_header_by_number.pl out YT_T6

mv ordered.fa trans_6_frame

6.2 Genome and transcriptome six-frame translation

cd 4_MS_database/ all_6_frame

cat ~/examples/1_source_data/GCF_000146045.2_R64_genomic.fna ../../2_cluster_mRNA/YT_Unigene.fasta > all.fasta

transeq -sequence all.fasta -outseq out -frame 6

flat_the_fasta_seq.pl out

rm out

get_seq_between_asterisk.pl flated 30

cd-hit -i between_asterisk -o out -c 1 -M 160000 -T 8

replace_fasta_header_by_number.pl out YT_A6

mv ordered.fa all_6_frame

6.3 CCPRD + contaminants

#Add artificial bacteria sequences to CCPRD, and give them the name tag YT_BAC respectively.

cat CCPRD BAC.fa > CCPRD_contam

# Citing CCPRD and software called by CCPRD
Since CCPRD is a pipeline that depends several Bioinformatics tools, publication of results obtained by CCPRD requires that not only CCPRD is cited, but also the tools that are used by CCPRD:
Please cite:

1.	Besemer, J., Lomsadze, A., & Borodovsky, M. (2001). Genemarks: a self-training method for prediction of gene starts in microbial genomes. implications for finding sequence motifs in regulatory regions. Nucleic Acids Research, 29(12), 2607-2618.
2.	Birney, E., Clamp, M., Durbin, R. (2004). GeneWise and Genomewise. Genome Research, 14, 988-995.
3.	Eren, A. M., Özcan C. Esen, Quince, C., Vineis, J. H., Morrison, H. G., & Sogin, M. L., et al. (2015). Anvi’o: an advanced analysis and visualization platform for ‘omics data. Peerj, 3.
4.	Grabherr, M. G. , Haas, B. J. , Yassour, M. , Levin, J. Z. , Thompson, D. A. , & Amit, I. , et al. (2011). Full-length transcriptome assembly from rna-seq data without a reference genome. Nature Biotechnology, 29(7), 644-652.
5.	Haas, B. J. , Delcher, A. L. , Mount, S. M. , Wortman, J. R. , & White, O. . (2003). Improving the arabidopsis genome annotation using maximal transcript alignment assemblies. Nucleic Acids Research, 31(19), 5654-5666.
6.	Haas, B. J. , Salzberg, S. L. , Zhu, W. , & Mihaela Pertea…. (2008). Automated eukaryotic gene structure annotation using evidencemodeler and the program to assemble spliced alignments. Genome biology, 9(1).
7.	Haas, B. J., Papanicolaou, A., Yassour, M., Grabherr, M., Blood, P. D., & Bowden, J., et al. (2013). De novo transcript sequence reconstruction from rna-seq using the trinity platform for reference generation and analysis. Nature Protocols, 8(8), 1494-1512.
8.	Huerta-Cepas, J., Serra, F., & Bork, P. (2016). Ete 3: reconstruction, analysis, and visualization of phylogenomic data. Molecular Biology & Evolution, 33(6), 1635-1638. 
9.	Korf, I. (2004). Gene finding in novel genomes. Bmc Bioinformatics, 5(1), 59.
10.	Laetsch, D. R., Blaxter, M. L. (2017). Blobtools: Interrogation of genome assemblies. F1000Research, 6, 1287.
11.	Li, W., & Godzik, A. (2006). Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences. Bioinformatics, 22(13), 1658. 
12.	Lomsadze, A. , Burns, P. D. , & Borodovsky, M. . (2014). Integration of mapped rna-seq reads into automatic training of eukaryotic gene finding algorithm. Nucleic Acids Research, 42(15), e119-e119.
13.	Pertea, G., Huang, X., Liang, F., Antonescu, V., Sultana, R., & Karamycheva, S., et al. (2003). Tigr gene indices clustering tools (tgicl): a software system for fast clustering of large est datasets. Bioinformatics, 19(5), 651-652. 
14.	Stanke, M., Keller, O., Gunduz, I., Hayes, A., Waack, S., & Morgenstern, B. (2006). Augustus: ab initio prediction of alternative transcripts. Nucleic Acids Research, 34 (Web Server issue), 435-9.

# License
All source code, i.e. scripts/.pl, scripts/.sh or scripts/.py are under the MIT license.






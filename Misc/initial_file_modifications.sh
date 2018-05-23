#! /bin/sh

## Set up the features.txt fril from the StringTie assembly:

cd ~/Dropbox/RNAseq/Female_RNAseq/DnovPmRNAseq/StringTie/

echo "contig\tmin\tmax\tgene_id\ttranscript_id\tref_gene_id" > tmp.features.header.txt

awk '{ if ($3 == "transcript") print $0}' stringtie_merged.gtf | awk '{print $1"\t"$4"\t"$5"\t"$10"\t"$12"\t"$14}' | sed 's/"//g' | sed 's/;//g' | gsed '/MSTRG.*MSTRG/ s/$/NA/g' | cat tmp.features.header.txt - > features.txt

rm tmp.features.header.txt

###-------------------------------------------------------------###

## extract gene lengths from stringtie gtf file:
cd ~/Dropbox/RNAseq/Female_RNAseq/DnovPmRNAseq/StringTie/

gtf2gff3.py stringtie_merged.gtf > stringtie_merged.gff3

python3 ~/ybin/transcript_length_from_gff.py stringtie_merged.gff3 > tmp.transcript_lengths.txt
echo "transcript_id length exon_number" > tmp.tLength.header.txt
cat tmp.tLength.header.txt tmp.transcript_lengths.txt > transcript_lengths.txt
rm tmp.tLength.header.txt tmp.transcript_lengths.txt

###-------------------------------------------------------------###

## append a summary description column from the swissprot blastX hit to the Trinotate databases

cd ~/Dropbox/RNAseq/Female_RNAseq/DnovPmRNAseq/Trinotate/Genome

## genome
cat Trinotate_report.xls | awk 'BEGIN{FS=OFS="\t"}{print $3}' | sed 's/.*Full=//g' | sed 's/;.*//g' | sed 's/sprot_Top_BLASTX_hit/sprot_Top_BLASTX_hit_description/g' | paste -d"\t" Trinotate_report.xls - > Trinotate_report_with_description_column.xls

sed -i.bak 's/\#gene_id/gene_id/g' Trinotate_report_with_description_column.xls
rm Trinotate_report_with_description_column.xls.bak


## pasa
cd ~/Dropbox/RNAseq/Female_RNAseq/DnovPmRNAseq/Trinotate/PASA

cat Trinotate_report.xls | awk 'BEGIN{FS=OFS="\t"}{print $3}' | sed 's/.*Full=//g' | sed 's/;.*//g' | sed 's/sprot_Top_BLASTX_hit/sprot_Top_BLASTX_hit_description/g' | paste -d"\t" Trinotate_report.xls - > tmp_Trinotate_report_with_description_column.xls

awk 'BEGIN{FS=OFS="\t"}{print $8}' tmp_Trinotate_report_with_description_column.xls | sed 's/\^.*//g' | sed 's/dvir.proteins.fasta_BLASTX/dvir1.06.BLASTX/g' | paste -d"\t" tmp_Trinotate_report_with_description_column.xls - > Trinotate_report_with_description_column.xls 

rm tmp_Trinotate_report_with_description_column.xls

sed -i.bak 's/\#gene_id/gene_id/g' Trinotate_report_with_description_column.xls
rm Trinotate_report_with_description_column.xls.bak

###-------------------------------------------------------------###

## extract transcript lengths from the pasa transcriptome

cd ~/Dropbox/RNAseq/Female_RNAseq/DnovPmRNAseq/Trinity/pasa/bt2_reference/

echo "transcript_id\tlength" > header

get_transcript_lengths.pl compreh_init_build.fasta | cat header - > pasa_transcript_lengths.txt


###-------------------------------------------------------------###

## Generate the samples.txt file

cd ~/Dropbox/RNAseq/Female_RNAseq/DnovPmRNAseq/Genome/eXpress

head -1 genome.gene.counts.matrix | tr '\t' '\n' | sed '/^$/d' > reps

head -1 genome.gene.counts.matrix | tr '\t' '\n' | sed '/^$/d' | sed 's/_.$//g' | paste -d"\t" - reps  > ../../Misc/samples.txt

rm reps

###-------------------------------------------------------------###

## Add a header column to the pasa vs. genome blast result

cd ~/Dropbox/RNAseq/Female_RNAseq/DnovPmRNAseq/Misc

echo "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" | tr ' ' '\t' | cat - BLASTn.pasa_vs_genome.outfmt6 > BLASTn.pasa_vs_genome.outfmt6.txt


###-------------------------------------------------------------###

## create counts by minimum TPM tables

cd ~/Dropbox/RNAseq/Female_RNAseq/DnovPmRNAseq/Genome/eXpress
~/Programs/trinityrnaseq-2.1.1/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl genome.gene.TPM.not_cross_norm > genome.gene.TPM.not_cross_norm.counts_by_min_TPM

cd ~/Dropbox/RNAseq/Female_RNAseq/DnovPmRNAseq/Genome/annotated_eXpress
~/Programs/trinityrnaseq-2.1.1/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl annotated.gene.TPM.not_cross_norm > annotated.gene.TPM.not_cross_norm.counts_by_min_TPM

cd ~/Dropbox/RNAseq/Female_RNAseq/DnovPmRNAseq/Trinity/pasa/eXpress
~/Programs/trinityrnaseq-2.1.1/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl pasa.gene.TPM.not_cross_norm > pasa.gene.TPM.not_cross_norm.counts_by_min_TPM

###-------------------------------------------------------------###

## extract one-to-one blast hits for pasa vs. stringtie transcripts table

cd ~/Dropbox/RNAseq/Female_RNAseq/DnovPmRNAseq/Misc

awk '{print $1"\t"$2}' BLASTn.pasa_vs_genome.outfmt6 | sort -u > BLASTn.pasa_vs_genome.outfmt6_one_to_one.txt

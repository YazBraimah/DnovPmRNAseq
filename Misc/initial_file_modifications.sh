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



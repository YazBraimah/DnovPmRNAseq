#!/bin/bash
#$ -S /bin/bash
#$ -q regular.q
#$ -j y
#$ -o mRNA_transfer_hetSnip.pasa.mrg.$JOB_ID
#$ -N mRNA_transfer_hetSnip.pasa.mrg
#$ -M ya76@cornell.edu
#$ -m be
#$ -pe bscb 1
#$ -binding linear:1
#$ -l h_vmem=16G
#$ -l h_rt=06:00:00
#$ -cwd

date
d1=$(date +%s)

mkdir -p /workdir/$USER/$JOB_ID
cd /workdir/$USER/$JOB_ID

## copy bam files
cp /fs/cbsufsrv5/data1/ya76/Snakemake.Output.Folders/D.novamexicana_Tuxedo2_ggTrinity_RNAseq/Trinity/pasa/eXpress/Female_H3_RT_1/bowtie2.bam /workdir/$USER/$JOB_ID/H3_RT_1.bowtie2.bam
cp /fs/cbsufsrv5/data1/ya76/Snakemake.Output.Folders/D.novamexicana_Tuxedo2_ggTrinity_RNAseq/Trinity/pasa/eXpress/Female_H3_RT_2/bowtie2.bam /workdir/$USER/$JOB_ID/H3_RT_2.bowtie2.bam
cp /fs/cbsufsrv5/data1/ya76/Snakemake.Output.Folders/D.novamexicana_Tuxedo2_ggTrinity_RNAseq/Trinity/pasa/eXpress/Female_H3_RT_3/bowtie2.bam /workdir/$USER/$JOB_ID/H3_RT_3.bowtie2.bam
cp /fs/cbsufsrv5/data1/ya76/Snakemake.Output.Folders/D.novamexicana_Tuxedo2_ggTrinity_RNAseq/Trinity/pasa/eXpress/Female_C3_RT_1/bowtie2.bam /workdir/$USER/$JOB_ID/C3_RT_1.bowtie2.bam
cp /fs/cbsufsrv5/data1/ya76/Snakemake.Output.Folders/D.novamexicana_Tuxedo2_ggTrinity_RNAseq/Trinity/pasa/eXpress/Female_C3_RT_2/bowtie2.bam /workdir/$USER/$JOB_ID/C3_RT_2.bowtie2.bam
cp /fs/cbsufsrv5/data1/ya76/Snakemake.Output.Folders/D.novamexicana_Tuxedo2_ggTrinity_RNAseq/Trinity/pasa/eXpress/Female_C3_RT_3/bowtie2.bam /workdir/$USER/$JOB_ID/C3_RT_3.bowtie2.bam

samtools sort -o H3_RT_1.bowtie2.sorted.bam H3_RT_1.bowtie2.bam
samtools sort -o H3_RT_2.bowtie2.sorted.bam H3_RT_2.bowtie2.bam
samtools sort -o H3_RT_3.bowtie2.sorted.bam H3_RT_3.bowtie2.bam
samtools sort -o C3_RT_1.bowtie2.sorted.bam C3_RT_1.bowtie2.bam
samtools sort -o C3_RT_2.bowtie2.sorted.bam C3_RT_2.bowtie2.bam
samtools sort -o C3_RT_3.bowtie2.sorted.bam C3_RT_3.bowtie2.bam

ls -1 C3*bowtie2.sorted.bam > c_mrg.list
ls -1 H3*bowtie2.sorted.bam > h_mrg.list

samtools merge -r -b c_mrg.list C3_RT.bam
samtools merge -r -b h_mrg.list H3_RT.bam

samtools index C3_RT.bam
samtools index H3_RT.bam

## copy reference
cp /fs/cbsufsrv5/data1/ya76/Snakemake.Output.Folders/D.novamexicana_Tuxedo2_ggTrinity_RNAseq/Trinity/pasa/bt2_reference/compreh_init_build.fasta .
samtools faidx compreh_init_build.fasta

## copy gene list:
cp /home/ya76/GitHub_Repositories/D.novamexicana_Tuxedo2_ggTrinity_RNAseq/transfered.mRNAs/transferred.mRNA.candidates.IDs.pasa.list  .

## make bam list
ls -1 *RT.bam > bam.list

hetSnipper.pl transferred.mRNA.candidates.IDs.pasa.list bam.list compreh_init_build.fasta candidate.male-transferred.mRNAs.pasa.mrg

### run vcftools stuff:
cd candidate.male-transferred.mRNAs.pasa.mrg

## make the bzipped/indexed VCFs
for file in *.output/*bam.vcf; do bgzip -c $file > $file.gz; tabix -p vcf $file.gz; done

## merge files for each gene
for folder in *output; do cd $folder; ls -1 *bam.vcf.gz > file.list; bcftools merge -l file.list -o $folder.merged.vcf; cd ..; done
for folder in *output; do cd $folder; ls -1 *C3*bam.vcf.gz > con_file.list; bcftools merge -l con_file.list -o $folder.con_merged.vcf; cd ..; done
for folder in *output; do cd $folder; ls -1 *H3*bam.vcf.gz > het_file.list; bcftools merge -l het_file.list -o $folder.het_merged.vcf; cd ..; done

## Make list of all merged files:
ls -1 *output/*output.merged.vcf > merged_files.list
ls -1 *output/*output.con_merged.vcf > con_merged_files.list
ls -1 *output/*output.het_merged.vcf > het_merged_files.list

## concatenate all merged files:
vcf-concat -f merged_files.list > all_cat_merged.vcf
vcf-concat -f con_merged_files.list > con_cat_merged.vcf
vcf-concat -f het_merged_files.list > het_cat_merged.vcf

## filter by a minimum read-depth of 20
# /programs/vcflib/bin/vcffilter -f "DP > 20" all_cat_merged.vcf > all_cat_merged.filtered.vcf
# /programs/vcflib/bin/vcffilter -f "DP > 20" con_cat_merged.vcf > con_cat_merged.filtered.vcf
# /programs/vcflib/bin/vcffilter -f "DP > 20" het_cat_merged.vcf > het_cat_merged.filtered.vcf

vcf-simplify.py SimplifyVCF -toType table -inVCF all_cat_merged.vcf -out all_cat_merged_simplified.vcf -unphased yes
# vcf-simplify.py SimplifyVCF -toType table -inVCF all_cat_merged.filtered.vcf -out all_cat_merged_simplified.filtered.vcf -unphased yes

cp -r /workdir/$USER/$JOB_ID/candidate.male-transferred.mRNAs.pasa.mrg /fs/cbsufsrv5/data1/ya76/Snakemake.Output.Folders/D.novamexicana_Tuxedo2_ggTrinity_RNAseq/Trinity/pasa/transferred_transcripts


rm -r /workdir/$USER/$JOB_ID

date
d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)


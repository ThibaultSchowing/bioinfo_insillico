#!/bin/bash

#SBATCH -o yeast-output.txt
#SBATCH -e yeast-error.txt

#SBATCH --job-name="SNP and vcf annotation LRAT - T5-8"
#SBATCH --time=24:00:00
#SBATCH --mem=25G
#SBATCH --partition=pcourse80
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1


#load modules
source /data/users/lfalquet/BC7107_19/scripts/module.sh
#create and go to the TP directory

#mkdir /data/users/student11/snp
cd /data/users/student11/snp
#link the reference genome and the reads locally


#
#ln -s /data/users/lfalquet/BC7107_17/reference/Saccharomyces_cerevisiae.R64-1-1.92.fa R64-1-1.92.fa
#ln -s /data/users/lfalquet/BC7107_17/reference/Saccharomyces_cerevisiae.R64-1-1.92.gtf R64-1-1.92.gtf
#ln -s /data/users/student11/bwa/*.bam .
#


# call SNPs for each strain (each bam file) and index the vcf.gz
#
#samtools mpileup -m 10 -ugf R64-1-1.92.fa T5.bam | bcftools call -vcO z --threads 4 -o T5.vcf.gz
#samtools mpileup -m 10 -ugf R64-1-1.92.fa T6.bam | bcftools call -vcO z --threads 4 -o T6.vcf.gz
#samtools mpileup -m 10 -ugf R64-1-1.92.fa T7.bam | bcftools call -vcO z --threads 4 -o T7.vcf.gz
#samtools mpileup -m 10 -ugf R64-1-1.92.fa T8.bam | bcftools call -vcO z --threads 4 -o T8.vcf.gz

#
#tabix T5.vcf.gz
#tabix T6.vcf.gz
#tabix T7.vcf.gz
#tabix T8.vcf.gz
#


#time to have a look at the individual vcf file (zmore)
#what do you see? many variants!
#ideally should be done for all mutants together…
#samtools mpileup -ugf R64-1-1.92.fa ${LIST_OF_BAM_FILES} | bcftools call -vcO v --threads 4 -o {XX}.vcf
#or much faster merged with bcftools once you have all vcf files and their indexes


#
#bcftools merge -O v `ls *.vcf.gz` > T5678_all.vcf
#


#annotate vcf file using snpEff
#first download and install snpEff locally
#wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
#unzip snpEff_latest_core.zip
# get the database for the yeast genome (warning in snpEff 86=90).


#java -Xmx4g -jar ./snpEff/snpEff.jar download R64-1-1.86


#annotate the VCF file
java -Xmx4g -jar ./snpEff/snpEff.jar -no-upstream -no-downstream R64-1-1.86 T5678_all.vcf > T5678_all.annot.vcf

#remove synonymous and intergenic variants
cat T5678_all.annot.vcf | java -jar ./snpEff/SnpSift.jar filter "(( ANN[*].EFFECT != 'synonymous_variant') & ( ANN[*].EFFECT != 'intergenic_region'))" > T5678_all_coding.vcf
# keep only variants that are found in less than 4 strains
cat T5678_all_coding.vcf | java -jar ./snpEff/SnpSift.jar filter "((countVariant() < 4))" > T5678_all_coding_max3.vcf

# Tip: modify before importing into excel
perl -pe 's/;ANN=/;\tANN=/g' T5678_all_coding_max3.vcf > T5678_all_coding_max3_forexcel.vcf
perl -i -pe 's/FORMAT\t/FORMAT\t\t/' T5678_all_coding_max3_forexcel.vcf
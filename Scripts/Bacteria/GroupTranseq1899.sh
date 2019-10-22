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

cd pan_genome_sequences
transeq -trim Y group_1899.fa.aln > group_1899.prot.fa
mafft --auto --clustalout group_1899.prot.fa > group_1899.prot.aln
DE NOVO ASSEMBLY ANNOTATION


1. FASTQ, TRIM, SOAPDENOVO, SPADES
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
USER=student11
NAMED=FAM1213FAM8102
NAME1=FAM1213
NAME2=FAM8102

#create and go to the TP directory
cd /data/users/${USER}
# go to the TP directory
mkdir asm_${NAMED}
cd asm_${NAMED}
#depending on your group create a link to the raw files assigned to you
# replace with your correct file, here I create an alias for your strain
ln -s /data/users/lfalquet/BC7107_19/reads/asm/${NAME1}_R1.fastq.gz ${NAME1}_R1.fastq.gz
ln -s /data/users/lfalquet/BC7107_19/reads/asm/${NAME1}_R2.fastq.gz ${NAME1}_R2.fastq.gz

ln -s /data/users/lfalquet/BC7107_19/reads/asm/${NAME2}_R1.fastq.gz ${NAME2}_R1.fastq.gz
ln -s /data/users/lfalquet/BC7107_19/reads/asm/${NAME2}_R2.fastq.gz ${NAME2}_R2.fastq.gz
#verify that your files are in your directory (ls -l)
#check quality of your data with fastqc (Using sbatch)
fastqc -t 2 ${NAME1}_R*.fastq.gz
fastqc -t 2 ${NAME2}_R*.fastq.gz



#################################### FROM HERE just read the content of script "asmdoitall.sh" until END########################################

#if needed apply sickle to filter out bad quality reads
#sickle pe -f ${NAME}_R1.fastq.gz -r ${NAME}_R2.fastq.gz -t sanger -l 127 -n -o ${NAME}_R1.fq -p ${NAME}_R2.fq -s ${NAME}_S.fq
#with trimmotatic (using sbatch)
trimmomatic PE -threads 8 -phred33 ${NAME}_R1.fastq.gz ${NAME}_R2.fastq.gz ${NAME}_1trim.fastq.gz ${NAME}_1unpaired.fastq.gz ${NAME}_2trim.fastq.gz ${NAME}_2unpaired.fastq.gz
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:8 MINLEN:100
gunzip -c ${NAME}_1trim.fastq.gz > ${NAME}_R1.fq
gunzip -c ${NAME}_2trim.fastq.gz > ${NAME}_R2.fq
#SOAPdenovo << Short-read Assembly Method (Made especially for Illumina GA short reads
create directory:
mkdir asm${NAME}_SOAP
cd asm${NAME}_SOAP
CONFIG="max_rd_len=101\n[LIB\navg_ins=300\nreverse_seq=0\nasm_flags=3\nrank=1\nq1=/data/users/${USER}/asm_${NAME}/${NAME}_R1.fq\nq2=/data/users/${USER}/asm_${NAME}/${NAME}_R2.fq\n"

echo -e $CONFIG >> ${NAME}.config
for k in 95 85 75 65; do
SOAPdenovo-127mer all -s ${NAME}.config -K $k -d 2 -o ${NAME}_$k -p 8;
/data/users/lfalquet/BC7107_19/scripts/GapCloser -b ${NAME}.config -a ${NAME}_$k.scafSeq -p 31 -t 8 -o ${NAME}_$k.closed;
done
#The contigs output file is called *.contigs
#The scaffolds output file is called *.scafSeq, *.closed is the name chosen for the gap closed scaffolds
# back to upper directory
cd ..
spades.py --careful -m 50 -t 8 -k 21,33,55,77,99 -1 ${NAME}_R1.fq -2 ${NAME}_R2.fq -o asm${NAME}_spades

########################################################### END OF script ####################################################################

#launch the assembly with
sbatch /data/users/lfalquet/BC7107_19/scripts/asmdoitall.slurm FAM${Our_Files}

#this will take at least 1hour or more.
#Compare the different kmer assembly of SOAPdenovo with the SPAdes assembly using
/data/users/lfalquet/BC7107_19/scripts/abyss-fac-all -t 1000 *.fasta
or
/data/users/lfalquet/BC7107_19/scripts/abyss-fac-all -t 1000 *.closed
Which one is the best according to you?
Try to Blastn a large contig vs nr to identify the closest homologous genome.
http://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome

#compare to reference genome with MAUVE










2. PROKKA (ANNOTATION AND GENERATION OF GFFS FOR ROARY)
#load modules
source /data/users/lfalquet/BC7107_19/scripts/module.sh

#Annotate the best assemblies with PROKKA (for each strain):
cd /data/users/student11/final_asm/
#if not already done, clean the fasta from small contigs
#sort_contigs.pl -b -m 1000 -p -z FAM${1}.fasta FAM${1}_1K.fasta

#rename the contig names (modifies SPAdes assemblies only)
perl -i -pe 's/(>NODE_\d+)_(.*)/$1 $2/' FAM${1}_1K.fasta
#example of PROKKA command used in sbatch script below
#prokka --force --outdir out_${NAME} --prefix ${NAME} --compliant --proteins /data/users/lfalquet/BC7107_19/reference/PGHs_9prot.fa --evalue 1e-30 --species "L.helveticus" --strain "${NAME}" --locustag ${NAME} --cpus 8 ${NAME}.fasta

#for each strain do (Take all the 6 best asm's from Falquets folder:
sbatch /data/users/lfalquet/BC7107_19/scripts/prokka.slurm FAM${1}_1K

# split sbatch because we need to wait for the output of the first sbatch to continue...

# these steps were lunched in the terminal directly (no script)

# Use prokka for the reference genome (this needs to be done to have a look at the reference in MAUVE):
cp /data/users/lfalquet/BC7107_19/reference/NC_010080.fasta .

sbatch /data/users/lfalquet/BC7107_19/scripts/prokka.slurm NC_010080

# The output is in the directory called out_FAMnnnn_1K or out_NC_010080


#Search for the genes of interest in *.gff or *.gbk or *.tbl
#Visualize with MAUVE

#create a directory for the GFFs and copy all of them there
mkdir GFFs
cp FAM*1K/*.gff GFFs/



3. ROARY
#Roary: Generating the pan-genome of several annotated strains
# put GFFs of the genomes to be analyzed into an array variable
GFFs=$(ls GFFs/*gff)
echo ${GFFs[@]}
# example of roary command used in the sbatch script below:
#roary -p 8 -r -f $PWD/ROARY -e -n -v -z ${GFFs[@]} 
    # -p (number of threads) / -r create r plots // -f output directory 
    # -e create multiFASTA aligment of core genes using PRANK
    # -n fast core aligment with mafft, use with -e
    # -v verbose output to SDTOUT
    # -z dont delete intermediate files
#Try to understand the parameters used by looking at the software webpage https://sanger-pathogens.github.io/Roary/
#run roary. It will run in less than 1h.
sbatch /data/users/lfalquet/BC7107_19/scripts/roary.slurm

            # See roary.slurm!

# copy the results to your local computer with scp or winscp
# Go through the files and try the understand what is the information encoded.
# How many core genes and accessory genes are detected? Hint: Look at the file summary_statistics.txt
# Visualize roary results in phandango. http://jameshadfield.github.io/phandango/

# Drag and drop the following files:
# - gene_presence_absence.csv
# - accessory_binary_genes.fa.newick
# - Zoom using ctrl+mouse scroll
#try to extract the genes that are specific for each phenotypes:
#from within ROARY directory:
query_pan_genome -a difference --input_set_one ../GFFs/FAM19191_1K.gff,../GFFs/FAM23285_1K.gff,../GFFs/FAM8102_1K.gff --input_set_two ../GFFs/FAM1213_1K.gff,../GFFs/FAM1450_1K.gff,../GFFs/FAM22076_1K.gff
#this outputs several files that can be queried for the "Lhv_" annotation:
grep Lhv_ set_difference_unique_set_one_statistics.csv
grep Lhv_ set_difference_unique_set_two_statistics.csv
grep Lhv_ set_difference_common_set_statistics.csv

#Compare the results within Excel and try to understand which PGH is involved in the activity seen on the gel and which PGH family contains pseudogene in one or more strain
#for each group of interest do:

cd pan_genome_sequences
transeq -trim Y group_NNNN.fa.aln > group_NNNN.prot.
mafft --auto --clustalout group_NNNN.prot.fa > group_NNNN.prot.aln 



ROARY.SLURM
#!/bin/bash

#SBATCH --mail-type=end,fail
#SBATCH --job-name="roary"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=25G
#SBATCH --partition=pcourse80

source /data/users/lfalquet/BC7107_19/scripts/module.sh
export PATH=/data/users/lfalquet/BC7107_19/scripts/:$PATH

GFFs=$(ls GFFs/*gff)

roary -p 8 -r -f ${PWD}/ROARY -e -n -v -z ${GFFs[@]}

# example of roary command used in the sbatch script below:
#roary -p 8 -r -f $PWD/ROARY -e -n -v -z ${GFFs[@]} 
    # -p (number of threads) / -r create r plots // -f output directory 
    # -e create multiFASTA aligment of core genes using PRANK
    # -n fast core aligment with mafft, use with -e
    # -v verbose output to SDTOUT
    # -z dont delete intermediate files

DE-NOVO ASSEMBLY ANNOTATION
#IMPORTANT: all your best assemblies should be located in "/data/users/$USER/final_asm/" directory and have the same name structure FAM${1}.fasta


#Start Script 2
#load modules
source /data/users/lfalquet/BC7107_19/scripts/module.sh

# Annotate the best assemblies with PROKKA (for each strain):
cd /data/users/$USER/final_asm/

#if not already done, clean the fasta from small contigs << in case we use our own files that have not been cleaned
#sort_contigs.pl -b -m 1000 -p -z FAM${1}.fasta FAM${1}_1K.fasta

#rename the contig names (modifies SPAdes assemblies only)
perl -i -pe 's/(>NODE_\d+)_(.*)/$1 $2/' FAM${1}_1K.fasta

#example of PROKKA command used in sbatch script below

#prokka --force --outdir out_${NAME} --prefix ${NAME} --compliant --proteins /data/users/lfalquet/BC7107_19/
reference/PGHs_9prot.fa --evalue 1e-30 --species "L.helveticus" --strain "${NAME}" --locustag ${NAME} --cpus 8
${NAME}.fasta

#for each strain do:
sbatch /data/users/lfalquet/BC7107_19/scripts/prokka.slurm FAM${1}_1K.fasta
# End Script 2


#and for reference strain do (done in terminal):
cp /data/users/lfalquet/BC7107_19/reference/NC_010080.fasta .
sbatch /data/users/lfalquet/BC7107_19/scripts/prokka.slurm NC_010080


#The output is in the directory called out_FAMnnnn_1K or out_NC_010080
#Search for the genes of interest in *.gff or *.gbk or *.tbl
#Visualize with MAUVE

#create a directory for the GFFs and copy all of them there
mkdir GFFs
cp FAM*1K/*.gff GFFs/

#Roary: Generating the pan-genome of several annotated strains
# put GFFs of the genomes to be analyzed into an array variable
GFFs=$(ls GFFs/*gff)
echo ${GFFs[@]}

# example of roary command used in the sbatch script below:
#roary -p 8 -r -f $PWD/ROARY -e -n -v -z ${GFFs[@]} 
    # -p (number of threads) / -r create r plots // -f output directory 
    # -e create multiFASTA aligment of core genes using PRANK
    # -n fast core aligment with mafft, use with -e
    # -v verbose output to SDTOUT
    # -z dont delete intermediate files
#Try to understand the parameters used by looking at the software webpage >> https://sanger-pathogens.github.io/Roary/
#run roary. It will run in less than 1h.
sbatch /data/users/lfalquet/BC7107_19/scripts/roary.slurm

# copy the results to your local computer with scp or winscp
# Go through the files and try the understand what is the information encoded.
# How many core genes and accessory genes are detected?
Hint: Look at the file summary_statistics.txt

# Visualize roary results in phandango. >> http://jameshadfield.github.io/phandango/
# Drag and drop the following files:
# - gene_presence_absence.csv
# - accessory_binary_genes.fa.newick
# - Zoom using ctrl+mouse scroll
#try to extract the genes that are specific for each phenotypes:
#from within ROARY directory:
query_pan_genome -a difference --input_set_one ../GFFs/FAM19191_1K.gff,../GFFs/FAM23285_1K.gff,../GFFs/FAM8102_1K.gff --input_set_two ../GFFs/FAM1213_1K.gff,../GFFs/FAM1450_1K.gff,../GFFs/FAM22076_1K.gff
#this outputs several files that can be queried for the "Lhv_" annotation:
grep Lhv_ set_difference_unique_set_one_statistics.csv
grep Lhv_ set_difference_unique_set_two_statistics.csv
grep Lhv_ set_difference_common_set_statistics.csv
#Compare the results within Excel and try to understand which PGH is involved in the activity seen on the gel and which PGH family contains
pseudogene in one or more strain
#for each group of interest do:
cd pan_genome_sequences
transeq -trim Y group_{GROUP NAME}.fa.aln > group_{GROUP NAME}.prot.fa
mafft --auto --clustalout group_NNNN.prot.fa > group_{GROUP NAME}.prot.aln

# MAFFT is a multiple sequence alignment program for unix-like operating systems. It offers a range of multiple alignment methods. 
#what are the common PGHs?
#Do you find some pseudo genes?



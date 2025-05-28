# Copying the read files from Dr. Lukens' directory to my directory as several files were gzipped and could not be unzipped in the residing directory:

cp /scratch/lukens/Assignment_2_Seqs/*.fastq . /home/vbularca/scratch/Assignment_2/Sequences
cp /scratch/lukens/Assignment_2_Seqs/*.fastq.gzip /home/vbularca/scratch/Assignment_2/Sequences

# Loading necessary modules:

module load StdEnv/2020 module load star/2.7.9a
module load subread/2.0.3
module load samtools

# Creating a star_job.sh file to allocate resources and recursively align reads from each file to the reference genome:

#!/bin/bash
#SBATCH --job-name=STAR_Alignment #SBATCH --output=STAR_%j.out #SBATCH --error=STAR_%j.err #SBATCH --cpus-per-task=8
#SBATCH --mem=32000
#SBATCH --time=18:00:00

for i in ./Sequences/*.fastq; do j=$(basename $i)
STAR --runMode alignReads \
--runThreadN 8 \
--genomeDir /scratch/lukens/Assignment_2_Genome \ --readFilesIn $i \
--outFileNamePrefix $j
Done

# Running the job using sbatch: 

sbatch star_job.sh

# Converting sam files to bam files:

for sam_file in *.sam; do 
bam_file="${sam_file%.sam}.bam"
samtools view -Sb "$sam_file" > "$bam_file"
Done

# Running featureCounts to get gene counts for each sample and saving it to one text file for further processing:

featureCounts -a /scratch/lukens/Assignment_2_Genome/genomic.gtf -o gene_counts.txt -T 8 -t exon -g gene_id -s 2 SRR10551657_1.fastqAligned.out.bam SRR10551658_1.fastqAligned.out.bam SRR10551659_1.fastqAligned.out.bam SRR10551660.fastqAligned.out.bam SRR10551661.fastqAligned.out.bam SRR10551662.fastqAligned.out.bam SRR10551663_1.fastqAligned.out.bam SRR10551664_1.fastqAligned.out.bam SRR10551665_1.fastqAligned.out.bam

# gene_id to ensure counts are assigned to gene IDs from the .gtf file
# -s 2 to define strandedness of RNA-seq data (NEBNext mRNA Library Prep Kit for Illumina results in reverse strand) 

head gene_counts.txt.summary

# The summary shows that there were no unmapped reads. However, there are a lot of unassigned multimapping reads.

# Dowloading gene_counts.txt to my local drive:

scp vbularca@graham.computecanada.ca:/home/vbularca/scratch/Assignment_2/gene_counts.txt /Users/vladbularca/Desktop
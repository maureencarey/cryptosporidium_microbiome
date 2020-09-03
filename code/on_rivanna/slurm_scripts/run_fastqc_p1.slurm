#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --account=tumi
#SBATCH --output=fastqc_p1.out

module purge
module load fastqc/0.11.5

sra_genome_path="/project/tumi/medlock/crypto_microbiome/data/sequence_backup/DataPool1/fastqs"
fastqc_output_path1="/home/mac9jc/crypto_micro/fastqc_1/run1"
cd $sra_genome_path

FILES1=$sra_genome_path/*.fastq.gz
for f in $FILES1; do
   fastqc $f
done

mv ./*_fastqc.* $fastqc_output_path1


sra_genome_path="/project/tumi/medlock/crypto_microbiome/data/sequence_backup/DataPool2/fastqs"
fastqc_output_path1="/home/mac9jc/crypto_micro/fastqc_1/run2"
cd $sra_genome_path

FILES1=$sra_genome_path/*.fastq.gz
for f in $FILES1; do
   fastqc $f
done

mv ./*_fastqc.* $fastqc_output_path1

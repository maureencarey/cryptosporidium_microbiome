#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --account=tumi
#SBATCH --output=fastqc_p2.out

module purge
module load fastqc/0.11.5

analyze_genome_path="/scratch/mac9jc/crypto_micro/trimmed_seqs/run1"
fastqc_output_path2="/home/mac9jc/crypto_micro/fastqc_2/run1"
cd $analyze_genome_path

FILES1=$analyze_genome_path/*.fq
for f in $FILES1; do
   fastqc $f
done

mv ./*_fastqc.* $fastqc_output_path2

analyze_genome_path="/scratch/mac9jc/crypto_micro/trimmed_seqs/run2"
fastqc_output_path2="/home/mac9jc/crypto_micro/fastqc_2/run2"
cd $analyze_genome_path

FILES1=$analyze_genome_path/*.fq
for f in $FILES1; do
   fastqc $f
done

mv ./*_fastqc.* $fastqc_output_path2

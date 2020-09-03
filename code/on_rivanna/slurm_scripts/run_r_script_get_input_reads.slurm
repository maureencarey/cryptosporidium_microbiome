#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --account=tumi
#SBATCH --mem-per-cpu=6000
#SBATCH --output=/project/tumi/carey/crypto_micro/get_reads.out

module purge
module load gcc openmpi R/4.0.0
Rscript /project/tumi/carey/crypto_micro/R_scripts/get_input_reads.R 

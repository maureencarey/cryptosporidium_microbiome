#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=standard
#SBATCH --cpus-per-task=20
#SBATCH --account=tumi
#SBATCH --mem-per-cpu=6000
#SBATCH --time=48:00:00
#SBATCH --output=/project/tumi/carey/crypto_micro/dada2_run2.out

module purge
module load gcc openmpi R/4.0.0
Rscript /project/tumi/carey/crypto_micro/R_scripts/dada2_run2.R

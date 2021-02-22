#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --account=tumi
#SBATCH --mem-per-cpu=6000
#SBATCH --time=48:00:00
#SBATCH --output=/project/tumi/carey/crypto_micro/unifrac.out

module purge
#module load gcc/7.1.0  openmpi/3.1.4 intel/18.0  intelmpi/18.0 
#module load R/4.0.0
module load intel/18.0 intelmpi/18.0
module load goolf/7.1.0_3.1.4 
module load R/4.0.0
Rscript /project/tumi/carey/crypto_micro/R_scripts/unifrac_phylo.R

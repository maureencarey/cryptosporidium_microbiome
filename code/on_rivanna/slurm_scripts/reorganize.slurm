#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=dev
#SBATCH --cpus-per-task=1
#SBATCH --account=tumi
#SBATCH --output=reorganize.out

module purge

og_path_run1="/project/tumi/medlock/crypto_microbiome/data/sequence_backup/DataPool1/fastqs"
og_path_run2="/project/tumi/medlock/crypto_microbiome/data/sequence_backup/DataPool2/fastqs"
og_trimmed_path1="/scratch/mac9jc/crypto_micro/trimmed_seqs/run1"
og_trimmed_path2="/scratch/mac9jc/crypto_micro/trimmed_seqs/run2"
path_run1="/project/tumi/carey/crypto_microbiome/raw_seqs/run1"
path_run2="/project/tumi/carey/crypto_microbiome/raw_seqs/run2"
trimmed_path_run1="/project/tumi/carey/crypto_microbiome/trimmed_seqs/run1"
trimmed_path_run2="/project/tumi/carey/crypto_microbiome/trimmed_seqs/run2"

# copy original files
FILES=$og_path_run1/*
for f in $FILES; do
    cp $f $path_run1
done
FILES=$og_path_run2/*
for f in $FILES; do
    cp $f $path_run2
done

# move trimmed files
FILES=$og_trimmed_path1/*
for f in $FILES; do
    cp $f $trimmed_path_run1 
done
FILES=$og_trimmed_path2/*
for f in $FILES; do
    cp $f $trimmed_path_run2
done



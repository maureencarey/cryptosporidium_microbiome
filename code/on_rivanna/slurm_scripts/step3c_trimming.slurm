#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --account=tumi
#SBATCH --output=trimming3c_132.out

module purge
module load gcc bbmap
module load cutadapt/2.5

sra_genome_path="/project/tumi/medlock/crypto_microbiome/data/sequence_backup"
analyze_genome_path="/scratch/mac9jc/crypto_micro/trimmed_seqs/run2"
cd $sra_genome_path

FILES1=$sra_genome_path/DataPool2/fastqs/132*.fastq.gz
for f in $FILES1; do

   if [[ $f == *_R1_* ]]
   then
      pre_string=$sra_genome_path/DataPool2/fastqs/

      echo '_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _'
      file_without_ext="${f%.*}"
      file_without_ext="${file_without_ext%.*}"
      echo $f
      f_without_path=$(basename $f)
      f_string=$(echo $(basename $file_without_ext)| cut -d'R' -f 1)
      f_without_num=${f_string%"R1_001"}
      f_rev=$f_without_num"R2_001.fastq.gz"
      f_rev=$pre_string$f_rev

      ## Adapter trimming # always recommended prior to merging
      output_file1=$analyze_genome_path/$f_without_num"_1_Adapter_trimmed.fq"
      output_file2=$analyze_genome_path/$f_without_num"_2_Adapter_trimmed.fq"
#      bbduk.sh in1=$f in2=$f_rev out1=$output_file1 out2=$output_file2 ref=~/cparvum_genomes/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo

      ## primer trimming
      output_file3=$analyze_genome_path/$f_without_num"_1_PT.fq.gz"
      output_file4=$analyze_genome_path/$f_without_num"_2_PT.fq.gz"
      output_file5=$analyze_genome_path/$f_without_num"_1_PT2.fq.gz"
      output_file6=$analyze_genome_path/$f_without_num"_2_PT2.fq.gz"
      cutadapt -a GTGCCAGCAGCCGCGGTAA...ATTAGATACCCTGGTAGTCC -A GGACTACCAGGGTATCTAAT...TTACCGCGGCTGCTGGCAC -o $output_file3 -p $output_file4 $output_file1 $output_file2 --cores=4 > $analyze_genome_path/"cutadapt_log_"$f_without_num".txt"
      cutadapt -a GGACTACCAGGGTATCTAAT...TTACCGCGGCTGCTGGCAC -A GTGCCAGCAGCCGCGGTAA...ATTAGATACCCTGGTAGTCC --minimum-length 50 -o $output_file5 -p $output_file6 $output_file3 $output_file4 --cores=4 > $analyze_genome_path/"cutadapt_log2_"$f_without_num".txt"

      # evaluate quality of reads   
   fi
done


#!/bin/bash
# This script iterates through a given directory and runs the program kallisto on each pair of paired fastq files, writing the output 
# to a specified output directory. The path to the kallisto executable file and a necessary index file are hardcoded into this script.
# Kallisto produces 3 files for each subject(pair of paired files), which are written to a directory named with the subject info. 
# When the script is completed, there will be an output directory (specified by second argument) that contains many subdirectories each named
# for the subject that kallisto ran on.
#
# Example Usage
# run_kallisto_on_trimmed_raw_files.sh /path/to/fastq_files/ kallisto_output_dir_name

#Input:
#First argument is directory where paired fastq files are stored, with file names in the following format:
#paired example files:
#Mar2015_L1_TCGACCA_L001_R1_001.trimmed.fastq.gz
#Mar2015_L1_TCGACCA_L001_R2_001.trimmed.fastq.gz
#the second argument is the name of a directory where the output will be written. The directory name given will be
#created within the input directory.

KALLISTO_PATH="/Users/burtonigenomics/Downloads/kallisto/build/src/kallisto"
OUT_DIR=$2
KALLISTO_INDEX="/Volumes/Seagate Expansion Hard Drive/_RNAseq/kallisto_index/H_burtoni_rna_092915.fa.kallisto.idx"
#KALLISTO_INDEX="/Volumes/Seagate Expansion Hard Drive/_RNAseq/kallisto_index/H_burtoni_rna.fa.kallisto.idx"
#KALLISTO_INDEX="/Volumes/Seagate\ Expansion\ Hard\ Drive/_RNAseq/kallisto_index/H_burtoni_rna.fa.kallisto.idx"
#old location on fish studies server:
#KALLISTO_INDEX="/Volumes/fishstudies/_Burtoni_genome_files/H_burtoni_rna.fa.kallisto.idx"
cd "$1"
mkdir "$OUT_DIR"
#get unique subject info for each pair, ex: Mar2015_L1_TCGACCA_L001
SUBJECTS=`ls | grep fastq | awk '{split($1,a,".");print a[1]}' |  awk '{split($1,a,"_R");print a[1]}' | uniq`
#iterate through all pairs, calling kallisto:
for pair_info in $SUBJECTS
do
  READS1=`ls | grep $pair_info | grep R1 | grep fastq`
  READS2=`ls | grep $pair_info | grep R2 | grep fastq`
    echo "----------------------------------------------------------------------------"
	echo "-------------- Working on "$pair_info" ------------------------------------------"
	echo "----------------------------------------------------------------------------"
	echo "---------------------- "$READS1
	echo "---------------------- "$READS2
	echo "---------------------- writing to "$OUT_DIR"/"$pair_info
  #call Kallisto:
  $KALLISTO_PATH quant \
  -i "$KALLISTO_INDEX" \
  -o "${OUT_DIR}/"$pair_info \
  -b 0 \
  -t 4 \
  $READS1 $READS2
done



  

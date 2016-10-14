#!/bin/bash
#Input is directory where paired fastq files are stored, with file names in the following format:
#paired example files:
#Mar2015_L1_TCGACCA_L001_R1_001.trimmed.fastq.gz
#Mar2015_L1_TCGACCA_L001_R2_001.trimmed.fastq.gz

KALLISTO_PATH="/Users/burtonigenomics/Downloads/kallisto/build/src/kallisto"
OUT_DIR=$2
KALLISTO_INDEX="/Volumes/fishstudies/_Burtoni_genome_files/H_burtoni_rna.fa.kallisto.idx"
cd $1
mkdir $OUT_DIR
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
  -i $KALLISTO_INDEX \
  -o $OUT_DIR"/"$pair_info \
  -b 0 \
  -t 4 \
  $READS1 $READS2
done



  

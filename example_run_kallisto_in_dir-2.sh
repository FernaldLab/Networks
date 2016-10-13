#!/bin/bash

### Command line arg $1 should be name of directory containing fastq files named with format:
###  L[0-9]_AGTCAA_*_R1_*.fastq
###	 L[0-9]_AGTCAA_*_R2_*.fastq
###  where "AGTCAA" is index for a subject 
### Fastq files should be only files in directory
###
### Output written to directory OUT_DIR 

OUT_DIR="/Volumes/fishstudies-1/_Elim_RNAseq/4_expression/NCBI/kallisto/Mar2015"		# should make this command line arg $2
KALLISTO_INDEX="/Volumes/fishstudies/_Burtoni_genome_files/H_burtoni_rna.fa.kallisto.idx"
echo "Analyzing files in:"
echo -e "  "$1"\n"
echo "Using kallisto index:"
echo -e "  "$KALLISTO_INDEX"\n"
echo "Writing output to:"
echo -e "  "$OUT_DIR"\n"

cd $1
SUBJECTS=`ls | awk '{split($1,a,"_");print a[2]}' | uniq`
for s in $SUBJECTS
do
	READS1=`ls *$s*R1*fastq`
	READS2=`ls *$s*R2*fastq`
	echo "----------------------------------------------------------------------------"
 	echo "-------------- Working on "$s" ------------------------------------------"
 	echo "----------------------------------------------------------------------------"
 	echo "---------------------- "$READS1
 	echo "---------------------- "$READS2
 	echo "---------------------- writing to "$OUT_DIR"/"$s
#  	kallisto quant \
#  	-i $KALLISTO_INDEX \
#  	-o $OUT_DIR"/"$s \
#  	-b 0 \
#  	-t 4 \
#  	$READS1 $READS2
done
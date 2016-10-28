#!/bin/bash
KALLISTO_PATH="/Users/burtonigenomics/Downloads/kallisto/build/src/kallisto"
KALLISTO_INDEX_1="/Volumes/Seagate Expansion Hard Drive/_RNAseq/kallisto_index/H_burtoni_rna_092915.fa.kallisto.idx"
KALLISTO_INDEX_2="/Volumes/Seagate Expansion Hard Drive/_RNAseq/kallisto_index/H_burtoni_rna_and_pseudo_092915.fa_kallisto.idx"
#KALLISTO_INDEX="/Volumes/Seagate\ Expansion\ Hard\ Drive/_RNAseq/kallisto_index/H_burtoni_rna.fa.kallisto.idx"
#old location on fish studies server:
#KALLISTO_INDEX="/Volumes/fishstudies/_Burtoni_genome_files/H_burtoni_rna.fa.kallisto.idx"
cd "$1"
mkdir "${KALLISTO_INDEX_1}_dir"
mkdir "${KALLISTO_INDEX_2}_dir"

#get unique subject info for each pair, ex: Mar2015_L1_TCGACCA_L001
SUBJECTS=`ls | grep fastq | awk '{split($1,a,".");print a[1]}' |  awk '{split($1,a,"_R");print a[1]}' | uniq`
#iterate through all pairs, calling kallisto:
pair_info=`echo $SUBJECTS | cut -d" " -f2`
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
		-i "$KALLISTO_INDEX_1" \
		-o "${KALLISTO_INDEX_1}_dir/"$pair_info \
		-b 0 \
		-t 4 \
		$READS1 $READS2
$KALLISTO_PATH quant \
		-i "$KALLISTO_INDEX_2" \
		-o "${KALLISTO_INDEX_2}_dir/"$pair_info \
		-b 0 \
		-t 4 \
		$READS1 $READS2






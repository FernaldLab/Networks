rm(list=ls());
options(stringsAsFactors=F);
library(JASPAR2014);
library(TFBSTools);
library(Rsamtools);


# load burtoni annotation info
#  load workspace from 'annotations_prep_for_dmrs_wgcna.R'
#   should include:
#    handannos, sl0, slGenesStats, zcne, abHsMap, an, an2, annoCombo, sl
#     zcne is a data frame with the highly conserved non-coding regions
load('/Volumes/fishstudies/_Burtoni_annotations/WORKSPACE_annotations_prep_for_dmrs_wgcna.RData');

# create a reference to genome fasta file
genomeFasta = FaFile('H_burtoni_v1.assembly.fa');

# make a GRanges object from zcne, strand info is irrelevant, 
#  fourth column is an id number that needs to be stored in mcols(zGR)
zGR = GRanges(seqnames=,
              ranges=IRanges(),
              strand='*',
              mcols= );

# make DNAStringSet object for zcnes
zGRseqs = getSeq(genomeFasta, zGR);
names(zGRseqs) = as.character(mcols(zGR)[,1]);

# search for one transcription factor at a time to benchmark 
#  use androgen receptor ids for testing
ARids = c('MA0007.1','MA0007.2');
min.score= '90%'

system.time( { 
  ARhits1 = searchSeq(ARids[1], zGRseqs, min.score=min.score)
} )

system.time( { 
  ARhits2 = searchSeq(ARids[2], zGRseqs, min.score=min.score)
} )

# 
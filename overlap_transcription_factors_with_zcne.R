#Overlapping transcription factors with zcne locations;
load('/Volumes/fishstudies/_Burtoni_annotations/WORKSPACE_annotations_prep_for_dmrs_wgcna.RData');
names(zcne) = c("seqnames", "start", "stop", "ID_num", "score", "strand");
zcne_GRanges = GRanges(seqnames= zcne$seqnames, ranges=IRanges(start=zcne$start, end=zcne$stop), strand='*', mcols = zcne[, 4:6] );
load('/Volumes//Seagate Expansion Hard Drive/_RNAseq/_R_data/JASPAR_hormone_receptor_IDs_wholeGenomeHits_dataFrame.RData');
receptor_ID_hits_GRanges = GRanges(seqnames= tfhitsMultiHormoneByScaffoldDf$seqnames, ranges = IRanges(start = tfhitsMultiHormoneByScaffoldDf$start,
                                                                                                       end = tfhitsMultiHormoneByScaffoldDf$end), strand = '*', mcols = tfhitsMultiHormoneByScaffoldDf[, c(-1, -4, -5)]);
receptor_ID_zcne_overlap_hits = findOverlaps(query=zcne_GRanges, subject=receptor_ID_hits_GRanges, ignore.strand=TRUE);


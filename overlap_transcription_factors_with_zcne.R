#Overlapping transcription factors with zcne locations;
rm(list=ls());
options(stringsAsFactors=F);
library(GenomicRanges)
load('/Volumes/fishstudies/_Burtoni_annotations/WORKSPACE_annotations_prep_for_dmrs_wgcna.RData');
names(zcne) = c("seqnames", "start", "stop", "ID_num", "score", "strand");
#Load zcne data into GRange Object
zcne_GRanges = GRanges(seqnames= zcne$seqnames, ranges=IRanges(start=zcne$start, end=zcne$stop), strand='*', mcols = zcne[, 4:6] );

#Load transcription_factor location data into GRange Object
load('/Volumes//Seagate Expansion Hard Drive/_RNAseq/_R_data/JASPAR_hormone_receptor_IDs_wholeGenomeHits_dataFrame.RData');
receptor_ID_hits_GRanges = GRanges(seqnames= tfhitsMultiHormoneByScaffoldDf$seqnames, ranges = IRanges(start = tfhitsMultiHormoneByScaffoldDf$start,
                                                                                                       end = tfhitsMultiHormoneByScaffoldDf$end), strand = '*', mcols = tfhitsMultiHormoneByScaffoldDf[, c(-1, -4, -5)]);

#Overlap the two GRange Objects, producing a GRange object representing transcription factor locations that were heavily conserved
#produces a hits object which is just dataframe of where overlaps occured
receptor_ID_zcne_overlap_hits = findOverlaps(query=zcne_GRanges, subject=receptor_ID_hits_GRanges, ignore.strand=TRUE);
#subset transcription factor locations with subjectHits column of hits object that was just returned by findOverlaps:
tf_and_zcne_df = tfhitsMultiHormoneByScaffoldDf[(subjectHits(receptor_ID_zcne_overlap_hits)),];

tf_and_zcne_df$zcneRowHit = queryHits(receptor_ID_zcne_overlap_hits);

# put zcne info into data frame that now has a row for each zcne-transcript_factor match
tf_and_zcne_df$zcne_start = zcne[tf_and_zcne_df$zcneRowHit,]$start;
tf_and_zcne_df$zcne_stop = zcne[tf_and_zcne_df$zcneRowHit,]$stop;
tf_and_zcne_df$zcne_ID_num = zcne[tf_and_zcne_df$zcneRowHit,]$ID_num;
tf_and_zcne_df$zcne_score = zcne[tf_and_zcne_df$zcneRowHit,]$score;
tf_and_zcne_df$zcne_strand = zcne[tf_and_zcne_df$zcneRowHit,]$strand;
save(file = "/Volumes/Seagate Expansion Hard Drive/_RNAseq/Networks/tf_zcne_df.RData", tf_and_zcne_df)
#make GRanges object out of tf_and_zcne_df
tf_and_zcne_GRanges = GRanges(
  seqnames = tf_and_zcne_df$seqnames,
  ranges = IRanges(
    start = tf_and_zcne_df$start,
    end = tf_and_zcne_df$end
  ),
  strand = tf_and_zcne_df$strand,
  mcols = tf_and_zcne_df
);


#find closest gene on either side of transcription factor: 
nearest_gene_index_vec = nearest(x=tf_and_zcne_GRanges, subject=an$gffGenesGR, ignore.strand = FALSE); 
nearest_gene_cols_to_add = as.data.frame(an$gffGenesGR)[nearest_gene_index_vec,];
names(nearest_gene_cols_to_add) = paste(names(nearest_gene_cols_to_add), "nearest", sep="_");
tf_and_zcne_df= cbind(tf_and_zcne_df, nearest_gene_cols_to_add);

#preceding:
preceding_gene_index_vec = precede(x=tf_and_zcne_GRanges, subject=an$gffGenesGR, ignore.strand = FALSE); 
preceding_gene_cols_to_add = as.data.frame(an$gffGenesGR)[preceding_gene_index_vec,];
names(preceding_gene_cols_to_add) = paste(names(preceding_gene_cols_to_add), "preceding", sep="_");
tf_and_zcne_df= cbind(tf_and_zcne_df, preceding_gene_cols_to_add);

#following:
following_gene_index_vec = follow(x=tf_and_zcne_GRanges, subject=an$gffGenesGR, ignore.strand = FALSE); 
following_gene_cols_to_add = as.data.frame(an$gffGenesGR)[following_gene_index_vec,];
names(following_gene_cols_to_add) = paste(names(following_gene_cols_to_add), "following", sep="_");
tf_and_zcne_df= cbind(tf_and_zcne_df, following_gene_cols_to_add);

#want distance to preceding and distance to following 
#get data frame with rows that HAVE a preceding gene:







#distance to preceding:
tf_and_zcne_df_no_na = tf_and_zcne_df[!is.na(tf_and_zcne_df$mcols.geneSym_preceding),];
na_tf_and_zcne_df = tf_and_zcne_df[is.na(tf_and_zcne_df$mcols.geneSym_preceding),];
preceding_GRanges = GRanges(seqnames= tf_and_zcne_df_no_na$seqnames_preceding, ranges= IRanges(start= tf_and_zcne_df_no_na$start_preceding,
                                                                                         end = tf_and_zcne_df_no_na$end_preceding), strand = tf_and_zcne_df_no_na$strand_preceding);

tf_and_zcne_no_na_GRanges = GRanges(seqnames=tf_and_zcne_df_no_na$seqnames, ranges = IRanges(start=tf_and_zcne_df_no_na$start,
                                                                                             end = tf_and_zcne_df_no_na$end), strand = tf_and_zcne_df_no_na$strand);
preceding_distance_vec = distance(x=tf_and_zcne_no_na_GRanges, y= preceding_GRanges);
tf_and_zcne_df_no_na$distance_to_preceding = preceding_distance_vec;
na_tf_and_zcne_df$distance_to_preceding = rep(NA, length(na_tf_and_zcne_df$seqnames));
tf_and_zcne_df_dist_precede_added = rbind(tf_and_zcne_df_no_na, na_tf_and_zcne_df);



#distance to following:
tf_and_zcne_df_no_na_following = tf_and_zcne_df_dist_precede_added[!is.na(tf_and_zcne_df_dist_precede_added$mcols.geneSym_following),];
na_tf_and_zcne_df_following = tf_and_zcne_df_dist_precede_added[is.na(tf_and_zcne_df_dist_precede_added$mcols.geneSym_following),];
following_GRanges = GRanges(seqnames= tf_and_zcne_df_no_na_following$seqnames_following, ranges= IRanges(start= tf_and_zcne_df_no_na_following$start_following,
                                                                                               end = tf_and_zcne_df_no_na_following$end_following), strand = tf_and_zcne_df_no_na_following$strand_following);

tf_and_zcne_no_na_GRanges_following = GRanges(seqnames=tf_and_zcne_df_no_na_following$seqnames, ranges = IRanges(start=tf_and_zcne_df_no_na_following$start,
                                                                                             end = tf_and_zcne_df_no_na_following$end), strand = tf_and_zcne_df_no_na_following$strand);
following_distance_vec = distance(x=tf_and_zcne_no_na_GRanges_following, y= following_GRanges);
tf_and_zcne_df_no_na_following$distance_to_following = following_distance_vec;
na_tf_and_zcne_df_following$distance_to_following = rep(NA, length(na_tf_and_zcne_df_following$seqnames));
tf_and_zcne_df_dist_follow_precede_added = rbind(tf_and_zcne_df_no_na_following, na_tf_and_zcne_df_following);

#join them back together:

#as.list on GRANge

#get Number of overlapping genes:
#make GRanges Object for 

tf_and_zcne_df_dist_follow_precede_added_GR = GRanges(seqnames= tf_and_zcne_df_dist_follow_precede_added$seqnames, ranges= IRanges(start= tf_and_zcne_df_dist_follow_precede_added$start,
                                                                                                                                   end = tf_and_zcne_df_dist_follow_precede_added$end), strand = tf_and_zcne_df_dist_follow_precede_added$strand);

num_genes_overlapped_with_tf_vec = countOverlaps(query=tf_and_zcne_df_dist_follow_precede_added_GR, subject= an$gffGenesGR, ignore.strand = FALSE);
tf_and_zcne_df_dist_follow_precede_added$num_genes_overlapped = num_genes_overlapped_with_tf_vec;




#========================================================================================================
#distance to nearest:
tf_and_zcne_df_no_na_nearest = tf_and_zcne_df_dist_follow_precede_added[!is.na(tf_and_zcne_df_dist_follow_precede_added$mcols.geneSym_nearest),];
na_tf_and_zcne_df_nearest = tf_and_zcne_df_dist_follow_precede_added[is.na(tf_and_zcne_df_dist_follow_precede_added$mcols.geneSym_nearest),];
nearest_GRanges = GRanges(seqnames= tf_and_zcne_df_no_na_nearest$seqnames_nearest, ranges= IRanges(start= tf_and_zcne_df_no_na_nearest$start_nearest,
                                                                                                         end = tf_and_zcne_df_no_na_nearest$end_nearest), strand = tf_and_zcne_df_no_na_nearest$strand_nearest);

tf_and_zcne_no_na_GRanges_nearest = GRanges(seqnames=tf_and_zcne_df_no_na_nearest$seqnames, ranges = IRanges(start=tf_and_zcne_df_no_na_nearest$start,
                                                                                                                 end = tf_and_zcne_df_no_na_nearest$end), strand = tf_and_zcne_df_no_na_nearest$strand);
nearest_distance_vec = distance(x=tf_and_zcne_no_na_GRanges_nearest, y= nearest_GRanges);
tf_and_zcne_df_no_na_nearest$distance_to_nearest = nearest_distance_vec;
na_tf_and_zcne_df_nearest$distance_to_nearest = rep(NA, length(na_tf_and_zcne_df_nearest$seqnames));
tf_and_zcne_df_dist_follow_precede_nearest_added = rbind(tf_and_zcne_df_no_na_nearest, na_tf_and_zcne_df_nearest);

save(tf_and_zcne_df_dist_follow_precede_nearest_added, file="/Volumes/Seagate Expansion Hard Drive//_RNAseq/Networks/transcription_factor_zcne_regions_gene_df.RData")


#Figure out upstream and downstream genes:

# isUpstream function determines whether a transcription factor is upstream or downstream from 
# a gene. It takes as argument, 
# strand :(represented as - or +) which is the strand that both 
# the gene and the transcription factor on.
#
#tf_start: start index of transcription factor
#tf_end: end index of transcription factor
#gene_start: start index of gene
#gene_end end index of gene
#
#Returns: a string, either "upstream" or "downstream", or NA if one or more input values are NA
  
isUpstream = function(strand, tf_start, tf_end, gene_start, gene_end) {
  isNA_vec = c(is.na(strand), is.na(tf_start), is.na(tf_end), is.na(gene_start), is.na(gene_end))
  if(!all(!isNA_vec)) {
    return(NA)
  }
  #assuming no overlaps, this is just to get whether it should be upstream or downstream
  #returns True for upstream, false for downstream
  if(strand == "+") {
    dist_between = gene_start - tf_start
  } else if(strand == "-") {
    dist_between = tf_end - gene_end
  } else {
    print(paste("Strand not plus or minus, is: ", strand, sep=""))
  }

  if(dist_between >0) {
    return("upstream")
  }
  
  return("downstream")
  
  
}  


#

#for each transcription factor: 
#determine whether the transcription factor is upstream or downstream from the nearest gene
#use apply, call isUpstream function on each row

#make column representing if tf is upstream or downstream from nearest gene:
tf_and_zcne_df_dist_follow_precede_nearest_added$nearest_up_or_down = 
  apply(X=tf_and_zcne_df_dist_follow_precede_nearest_added, 
        MARGIN=1, 
        FUN=function(row) {
          isUpstream(row["strand"],
                     as.numeric(row["start"]),
                     as.numeric(row["end"]),
                     as.numeric(row["start_nearest"]),
                     as.numeric(row["end_nearest"]))
        }
  )

#make column representing if tf is upstream or downstream from preceding gene:
tf_and_zcne_df_dist_follow_precede_nearest_added$preceding_up_or_down = 
  apply(X=tf_and_zcne_df_dist_follow_precede_nearest_added, 
        MARGIN=1, 
        FUN=function(row) {
          isUpstream(row["strand"],
                     as.numeric(row["start"]),
                     as.numeric(row["end"]),
                     as.numeric(row["start_preceding"]),
                     as.numeric(row["end_preceding"]))
        }
  );

#make column representing if tf is upstream or downstream from following gene:
tf_and_zcne_df_dist_follow_precede_nearest_added$following_up_or_down = 
  apply(X=tf_and_zcne_df_dist_follow_precede_nearest_added, 
        MARGIN=1, 
        FUN=function(row) {
          isUpstream(row["strand"],
                     as.numeric(row["start"]),
                     as.numeric(row["end"]),
                     as.numeric(row["start_following"]),
                     as.numeric(row["end_following"]))
        }
  )


#distance to preceding < 5000  or overlap
#make 2 smaller versions of dataframe, one with preceding hit within 5000 bp
preceding_gene_within_5000_tf_df = 
  tf_and_zcne_df_dist_follow_precede_nearest_added[
  tf_and_zcne_df_dist_follow_precede_nearest_added$distance_to_preceding < 5000,]

#other with overlaps

#check to see how many tfs have a preceding gene within 5000 nucleotides AND have an overlapping gene (overlapping gene is never preceding), its 731
has_overlap = 
  tf_and_zcne_df_dist_follow_precede_nearest_added[
  tf_and_zcne_df_dist_follow_precede_nearest_added$num_genes_overlapped > 0,]


#match precede df with recipBL gene row and cbind:
recipBL_matched_to_precede_gene = an$recipBL[match(preceding_gene_within_5000_tf_df$mcols.geneSym_preceding, an$recipBL$gene),]
names(recipBL_matched_to_precede_gene) = paste(names(recipBL_matched_to_precede_gene), "preceding_gene", sep = "_")
preceding_gene_within_5000_tf_df = as.data.frame(cbind(preceding_gene_within_5000_tf_df, recipBL_matched_to_precede_gene))
preceding_gene_within_5000_tf_df_no_na = preceding_gene_within_5000_tf_df[!is.na(preceding_gene_within_5000_tf_df$mcols.geneSym_preceding),]

#match overlap df with recipBL gene row and cbind:
recipBL_matched_to_overlap_gene = an$recipBL[match(has_overlap$mcols.geneSym_nearest, an$recipBL$gene),]
names(recipBL_matched_to_overlap_gene) = paste(names(recipBL_matched_to_overlap_gene), "overlap_gene", sep = "_")
has_overlap_tf_df = as.data.frame(cbind(has_overlap, recipBL_matched_to_overlap_gene))

#remove NA rows:
has_overlap_tf_df_no_na = has_overlap_tf_df[!is.na(has_overlap_tf_df$num_genes_overlapped),]


#ADD column in has_overlap_tf_df_no_na of distance of transcription factor start to start of Gene
#function used to find distance of overlapping tf and gene. This returns difference of start of
# gene and tf(or end if negative strand), 
get_dist_within_overlap = function(gene_start, gene_end, tf_start, tf_end, strand) {
  dist = 0
  if(strand == "+") {
    dist = as.numeric(tf_start) - as.numeric(gene_start)
  } else if(strand == "-") {
    #print(paste(gene_end, tf_end, sep = " __ "))
    dist = as.numeric(gene_end) - as.numeric(tf_end)
  }
  
  if(dist < 0) {
    dist = 0;
  }
  return(dist);
}

tf_nearest_dist_within_overlap_vec = apply(X = has_overlap_tf_df_no_na, MARGIN = 1, 
                                           FUN = function(row) {
                                             return(get_dist_within_overlap(row["start_nearest"], row["end_nearest"], 
                                                                            row["start"], row["end"], row["strand"]))
                                           })

has_overlap_tf_df_no_na$dist_within_overlap = tf_nearest_dist_within_overlap_vec
has_overlap_tf_df_no_na$dist_within_overlap_divided_by_gene_width = has_overlap_tf_df_no_na$dist_within_overlap / has_overlap_tf_df_no_na$width_nearest

save(file = "/Volumes/Seagate Expansion Hard Drive/_RNAseq/Networks/tf_zcne_info_overlapping_and_preceding_hits.RDATA", ... = preceding_gene_within_5000_tf_df_no_na, has_overlap_tf_df_no_na)
#subset for androgen receptor transcription factors:

AR_subset_has_overlap_tf_df_no_na = subset(has_overlap_tf_df_no_na, TF =="AR" | TF == "Ar" )
AR_subset_preceding_gene_within_5000_tf_df_no_na  = subset(preceding_gene_within_5000_tf_df_no_na, TF =="AR" | TF == "Ar" )


save(file = "/Volumes/Seagate Expansion Hard Drive/_RNAseq/Networks/AR_subset_tf_zcne_info.RDATA", ... = AR_subset_preceding_gene_within_5000_tf_df_no_na, AR_subset_has_overlap_tf_df_no_na)



#NOTE: AT this point we have 2 data frames, each subsets of the original 26,000 row dataframe containing info about tf, zcne, surrounding genes, overlaps, etc.
#preceding_gene_within_5000_tf_df is a subset of the original, with preceding genes within 5000 bp of the tf. Matching info about the preceding gene has been added in from an$recipbl
#has_overlap_tf_df is a subset of the original, with tf overlapping a gene. matching info about the overlapped gene has also been added in.

# The 2 dataframes will have duplicate rows, because theyre are some tf factors with a prev gene < 5000 AND an overlapping gene:
#these cases should result in an identical row:

#to combine, neglect rows in prev that have overlaps, this is done to avoid having duplicate rows:

#overlapped_gene_OR_preceding_gene_within_5000_tf_df = as.data.frame(rbind(has_overlap_tf_df, 
 #                                                                   preceding_gene_within_5000_tf_df[preceding_gene_within_5000_tf_df$num_genes_overlapped == 0,]))


#to-do: GET RID OF NAs, we dont need them:
#add recipBL_ to names from recipBL to avoid confusion

#then subset for ar receptors and email to Austin


#add columns with info from recipBL dataframe:
#recipBL_info_to_add_to_preceding_gene_within_5000_or_overlap_tf_list =
#apply(X=preceding_gene_within_5000_or_overlap_tf_df, 
#      MARGIN=1,
#      FUN=get_recipBL_row)
#recipBL_info_to_add_to_preceding_gene_within_5000_or_overlap_tf_list = lapply(X = recipBL_info_to_add_to_preceding_gene_within_5000_or_overlap_tf_list, as.data.frame)
#make above data into dataframe, it is currently list
#recipBL_info_to_add_to_preceding_gene_within_5000_or_overlap_tf_df = do.call(args=recipBL_info_to_add_to_preceding_gene_within_5000_or_overlap_tf_list, what="rbind")
#install.packages("RDocumentation")
#library(data.table)
#recipBL_do_call_info_to_add_to_preceding_gene_within_5000_or_overlap_tf_df = do.call("rbind", recipBL_info_to_add_to_preceding_gene_within_5000_or_overlap_tf_list)
#recipBL_info_to_add_to_preceding_gene_within_5000_or_overlap_tf_df = rbindlist(l = recipBL_info_to_add_to_preceding_gene_within_5000_or_overlap_tf_list, fill = TRUE)


#make subset with ar hits, sort it by the score of hit
#preceding_gene_within_5000_or_overlap_tf_df_recipBL_info_added_all_data = cbind(preceding_gene_within_5000_or_overlap_tf_df, preceding_gene_within_5000_or_overlap_tf_df_recipBL_info_added)



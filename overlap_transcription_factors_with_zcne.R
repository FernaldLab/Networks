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

#make GRanges object out of tf_and_zcne_df
tf_and_zcne_GRanges = GRanges(seqnames = tf_and_zcne_df$seqnames, ranges = IRanges(start = tf_and_zcne_df$start,
                                                                                   end = tf_and_zcne_df$end), strand = tf_and_zcne_df$strand, mcols = tf_and_zcne_df);


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


#has NAs where no other genes exist on same scaffold of_ tf_and_zcne_GRanges range


#sapply(X=nearest_gene_index_vec, FUN= helper)






#nearest_gene_df= an$gffGenesDF[nearest_gene_index_vec,]
#names(nearest_gene_df) = paste("nearest_gene", names(an$gffGenesDF), sep="_")
#tf_and_zcne_df$nearest_gene_index = nearest_gene_index_vec;

#tf_and_znce_merged_with_nearest_gene_df = merge(tf_and_zcne_df, nearest_gene_df)
#names(an$gffGenesDF) = paste("gene", names(an$gffGenesDF), sep="_");

following_gene_index_vec = follow(x=tf_and_zcne_GRanges, subject=an$gffGenesGR, ignore.strand = FALSE);

tf_and_zcne_df$preceding_gene_sym = an$gffGenesDF[preceding_gene_index_vec,]$geneSym;
tf_and_zcne_df$following_gene_sym = an$gffGenesDF[following_gene_index_vec,]$geneSym;

#column for how many it overlaps, distance to preceding, distance to following, names of preceding, following, and nearest
#NEXT STEP:
#use match() function to get more data about each gene in corresponding column.














#tf_and_zcne_GRanges



                   
                        
                        
#

#table(tfhitsMultiHormoneByScaffoldDf$TF)
#receptor_ID_to_zcne_hit_list <- vector(mode="list", length=length((subjectHits(receptor_ID_zcne_overlap_hits))))
#names(receptor_ID_to_zcne_hit_list) = receptor_ID_hits_GRanges[unique(subjectHits(receptor_ID_zcne_overlap_hits))]$mcols.TF;
#apply(subjectHits(receptor_ID_zcne_overlap_hits), function x

#      table(receptor_ID_hits_GRanges[subjectHits(receptor_ID_zcne_overlap_hits)]$mcols.TF)
      
      
      
      
#add_to_list = function(row_of_receptor_ID_hit) {
#   transcript_factor_name = receptor_ID_hits_GRanges[row_of_receptor_ID_hit]$mcols.TF;
#   #add entry to list:
#   receptor_ID_to_zcne_hit_list[transcript_factor_name] = 
#} 
      
    #  lapply()


#minus strand - want to see
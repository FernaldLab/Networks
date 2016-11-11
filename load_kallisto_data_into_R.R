#this code loads data returned by kallisto into R
rm(list=ls());
options(stringsAsFactors=F);

# ---------------------------------------------------------------------------------------------------
# load tpm (transcripts per million) data from kallisto for every subject 
# ---------------------------------------------------------------------------------------------------

#constants used for creating dataframe
NUM_OF_SUBJECTS = 41;
NUM_OF_TRANSCRIPT_IDS = 47807;
# get name of directory that contains subdirectories for each subject
outer_dir = "/Volumes/Seagate Expansion Hard Drive/_RNAseq/kallisto_ouput_new_index_file/";

  
# get names of subject directories, full.names is true so list.dirs returns whole path to each directory
subject_dirs = list.dirs(path=outer_dir, recursive= FALSE, full.names= TRUE);
subject_names = list.dirs(path=outer_dir, recursive=FALSE, full.names= FALSE);
  
# get filenames to load by appending file name to path
files = paste0(subject_dirs, '/abundance.tsv' );

# initialize data frame where rows are transcripts and columns are subjects
#  column names should be subject_dirs 
#  row names should be transcript ids
#   need to get from the row names in the abundance.tsv files 
#    use strsplit() to parse out just the part that starts XM_ or XR_
data_raw = as.data.frame(matrix(ncol= NUM_OF_SUBJECTS, nrow=NUM_OF_TRANSCRIPT_IDS));
colnames(data_raw) = subject_names;

#get target id names from 1st file
target_ids = read.delim(file=files[1], sep="\t")$target_id;
#parse out transcript_id, make it onto vector:
transcript_id_vec = unlist(lapply(strsplit(target_ids, split="[|]"), function(x) x[4]));
#set rownames equal to transcript ids:
rownames(data_raw) = transcript_id_vec;

# loop through files
#  for each file, load the data in the "tpm" column into the appropriate column of data_raw
counter=1;
for (f in files) {
  #print(paste(f, subject_names[counter], sep="__"));
  data_raw[, subject_names[counter]] = read.delim(file=f, sep="\t")$tpm
  counter = counter + 1
}

# load burtoni annotation info
#  load workspace from 'annotations_prep_for_dmrs_wgcna.R'
#   should include:
#    handannos, sl0, slGenesStats, zcne, abHsMap, an, an2, annoCombo, sl
load('/Volumes/fishstudies/_Burtoni_annotations/WORKSPACE_annotations_prep_for_dmrs_wgcna.RData');
rm(sl, sl0, handannos); gc();

# ---------------------------------------------------------------------------------------------------
# convert transcript level tpm values to gene level
# ---------------------------------------------------------------------------------------------------

# use match() to make sure data_raw rows match the rows in an$lookup

#get subset that matches for testing purposes
#data_raw_subset = data_raw[rownames(data_raw) %in% an$lookup$transcript_id,];
#an_lookup_subset = an$lookup[an$lookup$transcript_id %in% rownames(data_raw),];
data_raw_ordered = data_raw[order(rownames(data_raw)),];
an_lookup_ordered = an$lookup[order(an$lookup$transcript_id),]
#get gene
#data_raw_subset$gene = NULL
#apply(data_raw, function(x) an$lookup$transcript_id == x)
#data_raw = data_raw[match( ), ];

# use split() to make a list where each element is a data frame holding the transcripts for a given gene
data_raw_list = split(data_raw_ordered, an_lookup_ordered$gene);

# collapse each data frame in data_raw_list into one summed row 
#  create a new data frame where rows are genes instead of transcripts
#  could write nested for loops but much faster to use a combination of apply() and sapply()
tpm_transcript_df = data_raw;
gene_transcript_tpm_list = data_raw_list;
gene_tpm_df = as.data.frame(t(sapply(data_raw_list, function(x) colSums(as.data.frame(x)))));
save(tpm_transcript_df, gene_transcript_tpm_list, gene_tpm_df, file="/Volumes/Seagate Expansion Hard Drive/_RNAseq/Networks/kallisto_output_after_cutAdapt.RData");

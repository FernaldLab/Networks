rm(list=ls());
options(stringsAsFactors=F);

# ---------------------------------------------------------------------------------------------------
# load tpm (transcripts per million) data from kallisto for every subject 
# ---------------------------------------------------------------------------------------------------

# get name of directory that contains subdirectories for each subject
outer_dir = 

# get names of subject directories
subject_dirs = 
  
# get filenames to load
files = paste0(outer_dir, subject_dirs, '/abundance.tsv' );

# initialize data frame where rows are transcripts and columns are subjects
#  column names should be subject_dirs
#  row names should be transcript ids
#   need to get from the row names in the abundance.tsv files 
#    use strsplit() to parse out just the part that starts XM_ or XR_
data_raw = as.data.frame(matrix(ncol=, nrow=));

# loop through files
#  for each file, load the data in the "tpm" column into the appropriate column of data_raw
for (f in ) {
  data_raw[, ] = ;
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
data_raw = data_raw[match( ), ];

# use split() to make a list where each element is a data frame holding the transcripts for a given gene
data_raw_list = split(data_raw, );

# collapse each data frame in data_raw_list into one summed row 
#  create a new data frame where rows are genes instead of transcripts
#  could write nested for loops but much faster to use a combination of apply() and sapply()
data = as.data.frame(sapply(data_raw_list, 
                            function(x) apply(x, , )));




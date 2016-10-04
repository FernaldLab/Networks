#################################################################################
### load raw TPM data
#################################################################################
# rm(list=ls());
# source('/Volumes/fishstudies-1/_code/preProcATH-for_web_noVSN.R');

# expressionDIR = '~/Documents/kallisto_expr_cutadaptJan27/'  
# setwd(expressionDIR);

# tpm0 = read.table(paste0(expressionDIR, 'TPM_raw.tsv'), header=T, row.names=1);
# rg0 = read.table('~/Documents/rg_withFirst30Beh.tsv', header=T, row.names=1, sep='\t');

# # make sure the samples match up before doing this
# names(tpm0) = rownames(rg0);

##########################
# set up options
init = function() {
	options(stringsAsFactors=F);
}
init();


################################################################################
# setUpForDevelopment
#
# setwd('~/Desktop/Kai Fall 2016/Networks');
# setUpForDevelopment();
#
# Description:
#	Loads the data in test_data and sets other global variables as needed to
#	create a nice, tidy development environment.
#	Before calling this function, set the working directory to the equivalent
#	of the GitHub repo home folder.
#
# Arguments:
#	None.
#
# Returns:
#	None.
################################################################################
setUpForDevelopment = function() {
	source('utils.R');
	transcription_data <<- read.table('test_data/transcription_data.tsv', header = T, row.names = 1);
	readgroup_data <<- read.table('test_data/readgroupsdata.tsv', header = T, row.names = 1, sep = '\t');
}


################################################################################
# .matchSubjectNames
#
# Description:
#	Makes transcription data and read group data have matching subject names
#
# Arguments:
#	Required:
#	transcription_data - A data frame where each row is a gene and each column
#	                     is a subject, giving the transcripts per million (TPM)
#	                     expression level of each gene
#	readgroup_data     - A data frame with rownames equal to the column
#	                     names of transcription_data, containing columns with
#	                     read group data and (optionally) behavior data
#
#	Optional:
#	None.
#
# Returns:
#	A list containing:
#	transcription_data - A data frame where each row is a gene and each column
#	                     is a subject, giving the transcripts per million (TPM)
#	                     expression level of each gene
#	readgroup_data     - A data frame with rownames equal to the column
#	                     names of transcription_data, containing columns with
#	                     read group data and (optionally) behavior data
################################################################################
.matchSubjectNames = function(transcription_data, readgroup_data) {
	# we should really assert that all the names are really the same, but for now,
	# we'll content ourselves with a hacky version that just checks the barcodes
	# are the same.
	transcription_barcodes = regmatches(names(transcription_data), regexpr('[ACTG]*$', names(transcription_data)));
	readgroup_barcodes = regmatches(rownames(readgroup_data), regexpr('[ACTG]*$', rownames(readgroup_data)));
	stopifnot(sum(transcription_barcodes != readgroup_barcodes) == 0);

	names(transcription_data) = rownames(readgroup_data);
	return(list(transcription_data=transcription_data, readgroup_data=readgroup_data));
}

################################################################################
# getAndValidateData
#
# dataList <- getAndValidateData();
# transcription_data = dataList$transcription_data;
# readgroup_data = dataList$readgroup_data;
#
# Description:
#	Gets transcription data from parameter transcription_tsv or
#	a user-selected file, gets read group data from parameter
#	readgroup_tsv or a user-selected file, then validates both data frames
#	and matches up their names.
#
# Arguments:
#	Required:
#	None.
#
#	Optional:
#	transcription_tsv - the path of a .tsv file containing the transcription
#	                    data (in transcripts per million)
#	readgroup_tsv     - the path of a .tsv file containing the read group data
#
# Returns:
#	A list containing:
#	transcription_data - A data frame where each row is a gene and each column
#	                     is a subject, giving the transcripts per million (TPM)
#	                     expression level of each gene
#	readgroup_data     - A data frame with rownames equal to the column
#	                     names of transcription_data, containing columns with
#	                     read group data and (optionally) behavior data
################################################################################
getAndValidateData = function(transcription_tsv = NA, readgroup_tsv = NA) {
	if (is.na(transcription_tsv)) {
		transcription_tsv = .getFile('Please select the file with the TPM data (as a .tsv).');
	}
	cat('Loading transcription data from file \'', transcription_tsv, '\'\n');
	transcription_data = read.table(transcription_tsv, header = T, row.names = 1);
	
	if (is.na(transcription_tsv)) {
		readgroup_tsv = .getFile('Please select the file with read group and behavior data.');
	}
	cat('Loading read group data from file \'', readgroup_tsv, '\'\n');
	readgroup_data = read.table(readgroup_tsv, header = T, row.names = 1, sep = '\t');

	# this would be a great place to assert that transcription_data and readgroup_data are valid.

	dataList = .matchSubjectNames(transcription_data=transcription_data, readgroup_data=readgroup_data);

	return(dataList);
}

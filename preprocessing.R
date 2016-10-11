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
# source('preprocessing.R');
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
	trans_and_rg_data <- getAndValidateData(transcription_tsv = 'test_data/transcription_data.tsv',
	                               readgroup_tsv = 'test_data/readgroupsdata.tsv');
	transcription_data <<- trans_and_rg_data$transcription_data;
	readgroup_data <<- trans_and_rg_data$readgroup_data;
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
	# are the same. BACKLOG
	transcription_barcodes = regmatches(names(transcription_data), regexpr('[ACTG]*$', names(transcription_data)));
	readgroup_barcodes = regmatches(rownames(readgroup_data), regexpr('[ACTG]*$', rownames(readgroup_data)));
	stopifnot(sum(transcription_barcodes != readgroup_barcodes) == 0);

	names(transcription_data) = rownames(readgroup_data);
	return(list(transcription_data=transcription_data, readgroup_data=readgroup_data));
}

################################################################################
# getAndValidateData
#
# trans_and_rg_data <- getAndValidateData();
# transcription_data = trans_and_rg_data$transcription_data;
# readgroup_data = trans_and_rg_data$readgroup_data;
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

	# this would be a great place to assert that transcription_data and readgroup_data are valid. BACKLOG

	trans_and_rg_data = .matchSubjectNames(transcription_data=transcription_data, readgroup_data=readgroup_data);

	return(trans_and_rg_data);
}


################################################################################
# getSubset
#
# Description:
#	Should either return a List of subsets (all, all but ND, females, etc.) and/or
#	fetch the data for a specified subset (from transcription_data and readgroup_data)
#	TODO implement
#	see ref lines 18-43
################################################################################
getSubset = function() {
	stop('getSubset is unimplemented.')
}


################################################################################
# removeExcludedGenes
# 
#
#
# Description:
#
# Arguments:
#
# Returns:
################################################################################
removeExcludedGenes = function(trans_and_rg_data, geneNamesToExclude = c()) {
	transcription_data = trans_and_rg_data$transcription_data;

	# remove given genes
	rowsToRemove = match(geneNamesToExclude, rownames(transcription_data));
	transcription_data = transcription_data[-rowsToRemove,];

	# remove genes with zero variance
	transcription_data = transcription_data[apply(transcription_data, 1, var) > 0,]


	return(trans_and_rg_data);
}




################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# HERE IS THE SCRIPT PART
# trans_and_rg_data = getAndValidateData()
# removeExcludedGenes() TODO implement
# removeExcludedSamples() TODO implement

















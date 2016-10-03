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
options(stringsAsFactors=F);


setUpForDevelopment = function() {
	setwd('~/Desktop/Kai Fall 2016/Networks/');
	transcription_data <<- read.table('test_data/transcription_data.tsv', header = T, row.names = 1);
	readgroup_data <<- read.table('test_data/readgroupsdata.tsv', header = T, row.names = 1, sep = '\t');
}


################################################################################
# getAndValidateData
#
# list[transcription_data, readgroup_data] <- getAndValidateData();
#
# Description:
#	Gets transcription data from environment variable transcription_data or
#	a user-selected file, gets read group data from environment variable
#	readgroup_data or a user-selected file, then validates both data frames
#	and matches up their names.
#
# Arguments:
#	None
#
# Returns:
#	transcription_data - A data frame containing (something?) TODO
#	readgroup_data     - A data frame with rownames equal to the column
#	                     names of transcription_data, containing (something?) TODO
################################################################################
getAndValidateData = function() {
	if (exists('transcription_data')) {
		cat('Transcription data found in variable transcription_data\n');
	} else {
		transcripts_filepath = .getFile('Please select the file with the TPM data (as a .tsv).');
		transcription_data = read.table(transcripts_filepath, header = T, row.names = 1);
	}
	
	if (exists('readgroup_data')) {
		cat('Read group data found in variable readgroup_data\n');
	} else {
		readgroup_filepath = .getFile('Please select the file with read group and behavior data.');
		readgroup_data = read.table(readgroup_filepath, header = T, row.names = 1, sep = '\t');
	}

	# TODO implement these!!!
	#assertValidTranscriptionData(transcription_data);
	#assertValidReadGroupData(readgroup_data);

	#matchRowNames(transcription_data, readgroup_data);

	return(list(transcription_data, readgroup_data));
}

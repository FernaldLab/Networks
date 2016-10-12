testPreprocessing = function() {
	source('preprocessing.R');
	setUpForDevelopment();
	trans_and_rg_data <<- removeExcludedGenes(trans_and_rg_data=trans_and_rg_data);
}

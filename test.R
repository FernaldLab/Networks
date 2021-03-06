.runtest = function(condition, testname) {
    if (!condition) {
        stop(paste('ASSERT FAILED:', testname));
    } else {
		cat('Test passed:', testname, '\n');
	}
}



testPreprocessing = function(verbosity = 5) {
	source('preprocessing.R');
	setUpForDevelopment();
	setVerbosity(verbosity);
	trans_and_rg_data <- removeExcludedGenesAndNormalize(trans_and_rg_data=trans_and_rg_data);
	trans_and_rg_data <- removeRarelyExpressedGenes(
		trans_and_rg_data = trans_and_rg_data,
		max_fraction_zeroes = 0,
		only_one_group_can_have_zeroes = TRUE
	);
	trans_and_rg_data <<- trans_and_rg_data;
	setVerbosity(); # unsilence printing
}

test_removeRarelyExpressedGenes = function() {
	source('preprocessing.R');
	setUpForDevelopment();
	setVerbosity(-1); # silence printing
	zero_data = list(transcription_data = read.table('test_data/zero_count_test.tsv', header = T, row.names = 1, sep = '\t'),
	                 readgroup_data = trans_and_rg_data$readgroup_data);

	result1 = removeRarelyExpressedGenes(trans_and_rg_data = zero_data, max_fraction_zeroes = 1);
	.runtest(nrow(result1$transcription_data) == nrow(zero_data$transcription_data), 'No rows removed if max is 1');

	result2 = removeRarelyExpressedGenes(trans_and_rg_data = zero_data, max_fraction_zeroes = 0);
	.runtest(nrow(result2$transcription_data) == 57, 'All rows with >= 1 zero removed');

	result3 = removeRarelyExpressedGenes(trans_and_rg_data = zero_data, max_fraction_zeroes = 0.25);
	.runtest(nrow(result3$transcription_data) == 63, 'All rows with >= 25% zeroes removed');

	result4 = removeRarelyExpressedGenes(
		trans_and_rg_data = zero_data,
		max_fraction_zeroes = 0,
		only_one_group_can_have_zeroes = TRUE
	);
	.runtest(nrow(result4$transcription_data) == 57, 'All rows with >= 1 zero removed (by group)');

	result5 = removeRarelyExpressedGenes(
		trans_and_rg_data = zero_data,
		max_fraction_zeroes = 1,
		only_one_group_can_have_zeroes = TRUE
	);
	.runtest(nrow(result5$transcription_data) == 60, 'Three rows where only one condition had a zero were kept');

	setVerbosity(); # unsilence printing
}


test_runComBat = function() {
	source('preprocessing.R');
	source('test_data/combat_test.R');
	.runComBat(
		transcription_data = combat_test$dat,
		readgroup_data = combat_test$factors,
		batch_factor = 'batch',
		condition_matrix = combat_test$factors[,2],
		impute = FALSE
	);
}

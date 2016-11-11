##########################
# set up options
init = function() {
	options(stringsAsFactors=F);

	# try to source libraries.
	tryCatch({
		library(impute);
	},
	error = function(e) {
		warning(paste('Some dependencies of preprocessing.R could not be located. Please install',
		              'the required libraries.'));
		stop(e);
	});

	# try to install ComBat if necessary
	#tryCatch({
	#	stopifnot(is.function(ComBat));
	#},
	#error = function(e) {
	#	# try to install sva
	#	cat('Installing "sva" package from BioConductor');
	#	tryCatch({
	#		biocLite("sva");
	#	},
	#	error = function(e) {
	#		source("https://bioconductor.org/biocLite.R");
	#		biocLite("sva");
	#	});
	#});

	# fail gracefully if the user failed to setwd()
	tryCatch({
		source('utils.R');
		source('preProcess_interactive.R');
		source('ComBat.R');
	},
	warning = function(w) {
		stop(paste('Some dependencies of preprocessing.R could not be located. Please set the',
		           'working directory to be the directory containing the code files (using setwd()).'));
	});
}
init();

# BACKLOG: write a function that's the equivalent of --help


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
	init();
	trans_and_rg_data <<- getAndValidateData(transcription_tsv = 'test_data/transcription_data.tsv',
	                                         readgroup_tsv = 'test_data/readgroupsdata.tsv');
	#transcription_data <<- trans_and_rg_data$transcription_data;
	#readgroup_data <<- trans_and_rg_data$readgroup_data;
}



################################################################################
# Data Accessor Methods
#
# These are intended to abstract away the specific column names in the readgroup
# data table, so they can be changed without a major code refactor.
################################################################################
.getCondition = function(readgroup_data) {
	# get experimental group (e.g. dominant, female)
	return(readgroup_data$Condition)
}

.getReadgroupColumn = function(readgroup_data, colname) {
	if (colname %in% colnames(readgroup_data)) {
		return(readgroup_data[,colname])
	}
	else {
		stop(paste('No column', colname, 'in read group matrix.'));
	}
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
	transcription_barcodes = regmatches(names(transcription_data), regexpr('[ACTG]{6}', names(transcription_data)));
	readgroup_barcodes = regmatches(rownames(readgroup_data), regexpr('[ACTG]{6}', rownames(readgroup_data)));
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
		transcription_tsv = .getFile('Please select the file with the TPM data (as a .tsv).\n');
	}
	.catlog('Loading transcription data from file \'', transcription_tsv, '\'\n', sep = '', importance = 1);
	transcription_data = read.table(transcription_tsv, header = T, row.names = 1);
	
	if (is.na(readgroup_tsv)) {
		readgroup_tsv = .getFile('Please select the file with read group and behavior data.\n');
	}
	.catlog('Loading read group data from file \'', readgroup_tsv, '\'\n', sep = '', importance = 1);
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
#	BACKLOG implement
#	see ref lines 18-43
################################################################################
getSubset = function() {
	stop('getSubset is unimplemented.')
}


################################################################################
# removeExcludedGenesAndNormalize
# 
# trans_and_rg_data <- removeExcludedGenesAndNormalize(trans_and_rg_data=trans_and_rg_data)
#
# Description:
#   Removes the genes named in genes_to_exclude, eliminates genes with zero
#   variance, and normalizes the other genes to be log(expression + 1).
#
# Arguments:
#	Required:
#	trans_and_rg_data - A list containing:
#		transcription_data - A data frame where each row is a gene and each column
#		                     is a subject, giving the transcripts per million (TPM)
#		                     expression level of each gene
#		readgroup_data     - A data frame with rownames equal to the column
#		                     names of transcription_data, containing columns with
#		                     read group data and (optionally) behavior data
#
#	Optional:
#	genes_to_exclude  - A character vector containing names of genes to exclude
#	log_addition      - The constant to add to each expression level before taking
#                       the log (necessary because log(0) = -infinity)
#                       Note that later, zeroes are removed, but not negative numbers,
#                       so a value of 1 guarantees that all values are positive.
#
# Returns:
#	A list containing:
#	transcription_data - A data frame where each row is a gene and each column
#	                     is a subject, giving the log transcripts per million (TPM)
#	                     expression level of each gene
#	readgroup_data     - A data frame with rownames equal to the column
#	                     names of transcription_data, containing columns with
#	                     read group data and (optionally) behavior data
################################################################################
removeExcludedGenesAndNormalize = function(
	trans_and_rg_data,
	genes_to_exclude = c(),
	log_addition = 1
) {
	transcription_data = trans_and_rg_data$transcription_data;

	# remove given genes
	if (length(genes_to_exclude)) {
		rowsToRemove = match(genes_to_exclude, rownames(transcription_data));
		transcription_data = transcription_data[-rowsToRemove,];
		.catlog('Removed ', length(genes_to_exclude), ' user-provided genes.\n' , sep = '', importance = 1)
	}

	# remove genes with zero variance
	.catlog('Removed ', sum(!(apply(transcription_data, 1, var) > 0)), ' genes with zero variance.\n' , sep = '', importance = 1)
	transcription_data = transcription_data[apply(transcription_data, 1, var) > 0,];

	# take log of gene expression
	transcription_data = log2(transcription_data + log_addition);

	trans_and_rg_data$transcription_data = transcription_data;
	return(trans_and_rg_data);
}


################################################################################
# .sliceByCondition
# 
# experimental_conditions = .getConditions(trans_and_rg_data$readgroup_data);
# transcription_data_by_condition = sliceByCondition(
# 	transcription_data = trans_and_rg_data$transcription_data,
# 	condition_vector = experimental_conditions
# )
#
# Description:
# 	Separates the columns of a transcription data frame into separate data frames
# 	for each condition.
#
# Arguments:
#	Required:
#	transcription_data - A data frame where each row is a gene and each column
#	                     is a subject, giving the expression level of each gene
#	condition_vector   - A vector of length ncol(transcription_data), giving the
#	                     experimental condition for each group
#
#	Optional:
#	None.
#
# Returns:
#	A list where each entry is a data frame like transcription_data but where all
#	columns are from the same experimental condition (the name of the entry)
################################################################################
.sliceByCondition = function(transcription_data, condition_vector) {
	condition_names = names(table(condition_vector));
	sliced_data = list();
	for (condition in condition_names) {
		indices = which(condition_vector == condition);
		sliced_data[[condition]] = transcription_data[,indices]
	}
	return(sliced_data);
}

		


################################################################################
# removeRarelyExpressedGenes
# 
# trans_and_rg_data <- removeRarelyExpressedGenes(trans_and_rg_data=trans_and_rg_data)
#
# Description:
#	Removes genes with too many zero-expression values in the dataset, as determined
#	by the parameters.
#
# Arguments:
#	Required:
#	trans_and_rg_data - A list containing:
#		transcription_data - A data frame where each row is a gene and each column
#		                     is a subject, giving the log transcripts per million (TPM)
#		                     expression level of each gene
#		readgroup_data     - A data frame with rownames equal to the column
#		                     names of transcription_data, containing columns with
#		                     read group data and (optionally) behavior data
#
#	Optional:
#	max_fraction_zeroes            - Numeric; the maximum fraction of subjects
#	                                 that can have a zero value (between 0 and 1)
#	                                 This is a TEST_PARAMETER
#	only_one_group_can_have_zeroes - Logical; should genes be removed if they have
#	                                 *any* zero-expression values in more than one group?
#	                                 This is a TEST_PARAMETER
#
# Returns:
#	A list containing:
#	transcription_data - A data frame where each row is a gene and each column
#	                     is a subject, giving the log transcripts per million (TPM)
#	                     expression level of each gene
#	readgroup_data     - A data frame with rownames equal to the column
#	                     names of transcription_data, containing columns with
#	                     read group data and (optionally) behavior data
################################################################################
removeRarelyExpressedGenes = function(
	trans_and_rg_data,
	max_fraction_zeroes = NULL,
	only_one_group_can_have_zeroes = NULL
) {
	# default parameters
	if (is.null(max_fraction_zeroes)) max_fraction_zeroes = 0.3;
	if (is.null(only_one_group_can_have_zeroes)) only_one_group_can_have_zeroes = FALSE;

	transcription_data = trans_and_rg_data$transcription_data
	zero_counts = apply(transcription_data == 0, 1, sum);
	
	if (only_one_group_can_have_zeroes) {
		min_num_subjects_per_condition = floor(min(table(.getCondition(trans_and_rg_data$readgroup_data))));
		zlim = min_num_subjects_per_condition * max_fraction_zeroes;
		# any gene with more than zlim zeroes will be eliminated

		# generally, we're going to try to get rid of all genes that have any zeroes, unless
		# ALL the zeroes are in the same condition (interesting!)
		# this little block is to figure out for each gene whether there is more than one
		# condition with at least one zero-expression value.
		transcription_data_by_condition = .sliceByCondition(transcription_data, .getCondition(trans_and_rg_data$readgroup_data));
		zeros_by_condition = list();
		for (i in 1:length(transcription_data_by_condition)) {
			has_at_least_one_zero = apply(transcription_data_by_condition[[i]] == 0, 1, sum) > 0;
			zeros_by_condition[[names(transcription_data_by_condition)[i]]] = has_at_least_one_zero;
		}
		zeros_by_condition = data.frame(zeros_by_condition);
		num_conditions_with_at_least_one_zero = apply(zeros_by_condition, 1, sum);
		# this is now a vector giving for each gene, the number of conditions where the expression
		# level for at least one subject in that condition was 0.

		rows_to_remove = which(zero_counts > zlim | num_conditions_with_at_least_one_zero > 1)
	} else {
		num_subjects = ncol(trans_and_rg_data$transcription_data);
		zlim = num_subjects * max_fraction_zeroes;

		rows_to_remove = which(zero_counts > zlim);
	}

	if (length(rows_to_remove)) {
		.catlog('Removed', length(rows_to_remove), 'genes with too many zero-expression values.\n', importance = 1)
		trans_and_rg_data$transcription_data = transcription_data[-rows_to_remove, ];
	} else {
		.catlog('No genes have too many zero-expression values.\n', importance = 1)
	}
	return(trans_and_rg_data);
}



################################################################################
# removeLowVarianceGenes
# 
# trans_and_rg_data <- removeLowVarianceGenes(trans_and_rg_data = trans_and_rg_data)
#
# Description:
# 	Removes the genes with coefficient of variance less than the threshold.
#
# Arguments:
#	Required:
#	trans_and_rg_data - A list containing:
#		transcription_data - A data frame where each row is a gene and each column
#		                     is a subject, giving the log transcripts per million (TPM)
#		                     expression level of each gene
#		readgroup_data     - A data frame with rownames equal to the column
#		                     names of transcription_data, containing columns with
#		                     read group data and (optionally) behavior data
#
#	Optional:
#		make_plot          - Logical; should the distribution of coefficients of
#		                     variance be plotted, and the summaries output to the
#		                     console?
#		cv_threshold       - Numeric; the minimum coefficient of variance required
#		                     to remain in the dataset. Genes below this are excluded.
#
# Returns:
#	A list containing:
#	transcription_data - A data frame where each row is a gene and each column
#	                     is a subject, giving the log transcripts per million (TPM)
#	                     expression level of each gene
#	readgroup_data     - A data frame with rownames equal to the column
#	                     names of transcription_data, containing columns with
#	                     read group data and (optionally) behavior data
################################################################################
removeLowVarianceGenes = function(
	trans_and_rg_data,
	make_plot = FALSE,
	cv_threshold = NULL
) {
	# default parameters
	if (is.null(cv_threshold)) cv_threshold = 0.01;

	transcription_data = trans_and_rg_data$transcription_data;
	
	coefficients_of_variance = apply(transcription_data, 1, .coefficientOfVariance);
	if (make_plot) {
		#dev.off(); 
		par(mfrow=c(1,2));
		hist(coefficients_of_variance, xlab = 'Coefficient of variance', main = 'Before cut');
		summary(coefficients_of_variance);
	}

	cvcut = quantile(coefficients_of_variance, cv_threshold);
	.catlog('Removed ', sum(!(coefficients_of_variance > cvcut)), ' genes with coefficient of variance',
	        ' below threshold ', cv_threshold, '.\n' , sep = '', importance = 1)
	transcription_data = transcription_data[coefficients_of_variance > cvcut,];
	if (make_plot) {
		new_cv = apply(transcription_data, 1, .coefficientOfVariance);
		hist(new_cv, xlab = 'Coefficient of variance', main = 'After cut');
		summary(new_cv);
	}
	trans_and_rg_data$transcription_data <- transcription_data;
	return(trans_and_rg_data);
}


################################################################################
# runPreprocessInteractive
# 
# preprocess_out = runPreprocessInteractive(trans_and_rg_data=trans_and_rg_data)
#
# Description:
# 	Runs the interactive preprocessing script found in preProcess_interactive.R
#
# Arguments:
#	Required:
#	trans_and_rg_data - A list containing:
#		transcription_data - A data frame where each row is a gene and each column
#		                     is a subject, giving the log transcripts per million (TPM)
#		                     expression level of each gene
#		readgroup_data     - A data frame with rownames equal to the column
#		                     names of transcription_data, containing columns with
#		                     read group data and (optionally) behavior data
#
#	Optional:
#		deviate               - Numeric; when removing outlier "probes"/tpm values, what is the maximum number
#		                        of standard deviations from the mean to allow?
#		                        TEST_PARAMETER: try 2, 2.5, 3
#		probe_thresh          - Numeric; how many probes/tpm values can be NA for a given gene? NULL uses the
#		                        default value.
#		                        TEST_PARAMETER: try 0, 1/3.
#		sample_thresh         - Numeric; how many probes/tpm values can be NA for a given subject? NULL uses the
#		                        default value.
#		                        TEST_PARAMETER: try 0, 1/3.
#		                        BACKLOG try allowing NAs iff they are all from the same condition, similar to
#		                        the function for removing zeroes.
#		IACthresh             - Numeric; when removing outlier samples, what is the maximum number of standard
#		                        deviations away from the mean that you are willing to tolerate? (IAC means
#		                        inter-subject correlation)
#								TEST_PARAMETER: this is the interactive part of the function; not necessarily a hard threshold.
#								Should be no greater than 3; after 3, don't even ask a human. Test 2 and 2.5 as well.
#
# Returns:
# 	A list including:
# 		data_removedOutlierSamples  - transcription_data but with outlier samples removed
# 		data_Qnorm                  - transctiption_data Qnormalized (end of this step's pipeline)
# 		and other entries that are irrelevant.
################################################################################
runPreprocessInteractive = function(trans_and_rg_data,
						            deviate=NULL,
						            probe_thresh=NULL,
						            sample_thresh=NULL,
						            IACthresh=NULL,
									interactive = TRUE
) {
	# default parameters
	if (is.null(deviate)) deviate = 2.5;
	if (is.null(IACthresh)) IACthresh = 2;

	preprocess_out = preProcess(
		datIN = trans_and_rg_data$transcription_data,
		removeOutlierProbes=TRUE,
		deviate=deviate,
		removeTooManyNAs=TRUE,
		probe_thresh=probe_thresh,
		sample_thresh=sample_thresh,
		removeOutlierSamples=TRUE,
		IACthresh=IACthresh,
		Qnorm=TRUE,
		showplots = .verbosity >= 2,
		interactive = interactive
	);
	
	return(preprocess_out);
}




################################################################################
# .computeAndPlotFactorEffectsOnMeanExpr
#
# Description:
# 	Computes the ANOVA significance of each factor in factor_indices on the 
# 	transcription values, then makes a plot of these.
#
# Arguments:
# 	Required:
# 	transcription_data  - the transcription_data
# 	readgroup_data      - a readgroup_data dataframe corresponding to transcription_data
# 	subject_indices     - numeric; which subjects should be considered? (vector of indices)
# 	factor_indices      - numeric; column indices of INFO that give the factors whose effects
# 	                      on mean expression should be considered.
# 	
# 	Optional:
# 	make_plot - boolean; should a plot be produced? Defaults to TRUE.
# 	main      - the title of the plot
# 	...       - any graphical parameters, to be passed on to the barplot function
#
# Returns:
# 	A numeric vector of p values, where the name of each entry is the name of the
# 	factor column it corresponds to.
################################################################################
.computeAndPlotFactorEffectsOnMeanExpr = function (
	transcription_data,
	readgroup_data,
	subject_indices,
	factor_indices,
	make_plot = TRUE,
	main = '',
	...
) {
	transcription_data = transcription_data[, subject_indices];
	readgroup_data = readgroup_data[subject_indices, ];
	transcription_data = transcription_data[apply(transcription_data,1,var) > 0, ];


	mean_expression = apply(transcription_data, 2, mean, na.rm=T);
	factors = paste0(paste0('factor(',paste0('readgroup_data[,',factor_indices,']' ),')' ), collapse='+');
	lmformula = paste0('mean_expression~', factors);

	if (.verbosity >= 2) {
		print(lmformula);
	}
	p_values = as.vector(na.omit(anova(lm(as.formula(lmformula)))$"Pr(>F)"));
	if (.verbosity >= 2) {
		print(anova(lm(as.formula(lmformula))));
	}

	names(p_values) = names(readgroup_data)[factor_indices];
	if (make_plot) {
		barplot(
			-log10(p_values),
			ylab='-log10(pval)',
			main=paste(main,' (n=',ncol(transcription_data),')',sep=''),
			...
		);
		abline(h=c(-log10(.05), -log10(.01)), col='red');
	}

	return(p_values);
}



################################################################################
# .runComBat
# 
# Description:
# 	Runs the ComBat program from ComBat.R
#
# Arguments:
# 	Required:
#	transcription_data - A data frame where each row is a gene and each column
#	                     is a subject, giving the log transcripts per million (TPM)
#	                     expression level of each gene
#	readgroup_data     - A data frame with rownames equal to the column
#	                     names of transcription_data, containing columns with
#	                     read group data and (optionally) behavior data
#	batch_factor       - Character; name of the column in readgroup_data whose
#	                     effect should be filtered out.
#	condition_matrix   - A subset of the columns in readgroup_data (possibly just one
#	                     column) representing the covariates like condition whose effect
#	                     should be preserved
#	...                - Additional parameters to pass to ComBat()
#
# Returns:
# 	A matrix with the same dimensions as transcription_data, where the effects
# 	of batch_factor have been factored out.
################################################################################
.runComBat = function(
	transcription_data,
	readgroup_data,
	batch_factor, # name of a column in readgroup_data
	condition_matrix, # nsubjects x ncovariates (usually 1)
	...
) {
	combat_data_filename = paste('.transcriptionDataForComBat', '_', Sys.time(), '.tsv', sep = '')
	combat_info_filename = paste('.infoForComBat', '_', Sys.time(), '.tsv', sep = '')


	.catlog('Running ComBat for batch', batch_factor, '\n', importance = 1.5)
	# format expression data and write to file
	transcription_data = cbind(gene = rownames(transcription_data), transcription_data);
	write.table(transcription_data, file = combat_data_filename, quote = FALSE, sep = '\t', row.names = FALSE);


	# format group data and write to a file
	batch_vector = .getReadgroupColumn(readgroup_data, batch_factor)
	info = as.data.frame(cbind(rownames(readgroup_data), rownames(readgroup_data), batch_vector, condition_matrix));
	names(info)[1:3] = c('Array name', 'Sample name', 'Batch');
	if (length(names(info)) > 3) {
		names(info)[4:ncol(info)] = paste('Covariate', 4:ncol(info) - 3);	
	}
	write.table(info, file=combat_info_filename, quote=FALSE, sep='\t', row.names=FALSE);

	# run ComBat
	combatout = ComBat(
		expression_xls=combat_data_filename, 
	    sample_info_file=combat_info_filename,
		filter=FALSE,
		write=FALSE,
		skip=1,
		prior.plots = FALSE,
		...
	);
	file.remove(combat_data_filename);
	file.remove(combat_info_filename);
	
	# ComBat returns a character error message if it failed. Throw an actual error instead.
	if (is.character(combatout)) {
		stop(combatout);
	}
	return(combatout[,-1]); #-1 is to get rid of gene name column.
}



################################################################################
# .plotFactorEffects
# 
# Description:
# 	Makes length(factors_to_plot) plots, where the [[i]]th plot shows the relative
# 	contributions of the factors named in factors_to_plot[[i]].
#
# Arguments:
#	Required:
#	transcription_data - A data frame where each row is a gene and each column
#	                     is a subject, giving the log transcripts per million (TPM)
#	                     expression level of each gene
#	readgroup_data     - A data frame with rownames equal to the column
#	                     names of transcription_data, containing columns with
#	                     read group data and (optionally) behavior data
#	factors_to_plot    - A list where each entry is a character vector of factors
#	                     (column names in readgroup_data) giving the x-values for
#	                     that plot.
#	...                - Any additional parameters to pass to
#	                     .computeAndPlotFactorEffectsOnMeanExpr()
#
#	Optional:
#	None.
#
# Returns:
# 	The p_values vector returned by the first call to .computeAndPlotFactorEffectsOnMeanExpr()
# 	(first element of factors_to_plot)
################################################################################
.plotFactorEffects = function(
	transcription_data,
	readgroup_data,
	factors_to_plot,
	num_blank_plots = 0,
	...
) {
	p_values = NA
	no_p_values = TRUE
	num_subjects = ncol(transcription_data);
	for (factors in factors_to_plot) {
		readgroup_columns = match(factors, names(readgroup_data));
		pvs = .computeAndPlotFactorEffectsOnMeanExpr(
			transcription_data = transcription_data,
			readgroup_data = readgroup_data,
			subject_indices = 1:num_subjects,
			factor_indices = readgroup_columns,
			...
		);
		if (no_p_values) {
			p_values = pvs;
			no_p_values = FALSE;
		}
		else {
			for (factor in names(pvs)) {
				if (!(factor %in% p_values) && factor != 'LibSeqTank') {
					p_values[factor] <- pvs[factor]
				}
			}
		}
	}
	while (num_blank_plots > 0) {
		plot.new();
		num_blank_plots = num_blank_plots - 1;
	}
	return(p_values)
}


################################################################################
# .runCombatForFactors
# 
# Description:
# 	Iteratively run ComBat on the factors in combat_factors_sequence, make a plot
# 	for each factor, and output the results in a list where each entry is the
# 	TPM matrix after one ComBat run.
#
# Arguments:
#	Required:
#	transcription_data       - A data frame where each row is a gene and each column
#	                           is a subject, giving the log transcripts per million (TPM)
#	                           expression level of each gene
#	readgroup_data           - A data frame with rownames equal to the column
#	                           names of transcription_data, containing columns with
#	                           read group data and (optionally) behavior data
#	combat_factors_sequence  - the character vector giving the order in which factors should
#	                           be factored out by ComBat. 
#	factors_to_plot          - A character vector of factors (column names in readgroup_data)
#	                           giving the x-values for the plots.
#
#	Optional:
#	None.
#
# Returns:
# 	A list containing the transcription_data matrix [[1]] after imputation but before
# 	ComBat runs, [[2]] after the first ComBat run, [[3]] after the third, etc.
################################################################################
.runCombatForFactors = function(
	transcription_data,
	readgroup_data,
	combat_factors_sequence,
	factors_to_plot
) {
	if (.verbosity <= 2) {
		sink('.routput')
	}
	imputed_data = as.data.frame(impute.knn(as.matrix(transcription_data))$data);
	if (.verbosity <= 2) {
		sink()
		file.remove('.routput')
	}

	combat_outputs = list(preprocessed=imputed_data);
	condition_matrix = as.matrix(.getCondition(readgroup_data));
	for (factor in combat_factors_sequence) {
		new_transcription_data = combat_outputs[[length(combat_outputs)]];
		
		out = .runComBat(
			transcription_data = new_transcription_data,
			readgroup_data = readgroup_data,
			batch_factor = factor,
			condition_matrix = condition_matrix
		);
		pvals = .plotFactorEffects(
			transcription_data = out,
			readgroup_data = readgroup_data,
			factors_to_plot = factors_to_plot,
			main = paste('After ComBat on', factor)
		);
		if (pvals['Condition'] < 0.05 && sum(pvals[names(pvals) != 'Condition'] < 0.05) == 0) {
			break;
		}
		combat_outputs[[factor]] = out;
	}
	
	return(combat_outputs);
}



################################################################################
# plotFactorsAndRunComBat
# 
# plotFactorsAndRunComBat(trans_and_rg_data=trans_and_rg_data, preprocess_out=preprocess_out)
#
# Description:
# 	Makes a matrix of plots giving the significance of various factors to the
# 	gene expression values (factors_to_plot); then iteratively runs ComBat in the
# 	order of combat_factors_sequence and plots the results.
#
# Arguments:
#	Required:
#	trans_and_rg_data - A list containing:
#		transcription_data   - A data frame where each row is a gene and each column
#		                       is a subject, giving the log transcripts per million (TPM)
#		                       expression level of each gene
#		readgroup_data       - A data frame with rownames equal to the column
#		                       names of transcription_data, containing columns with
#		                       read group data and (optionally) behavior data
#	preprocess_out           - The value returned by runPreprocessInteractive
#
#	Optional:
#	combat_factors_sequence  - the character vector giving the order in which factors should
#	                           be factored out by ComBat. 
#
# Returns:
# 	A list containing the transcription_data matrix [[1]] after imputation but before
# 	ComBat runs, [[2]] after the first ComBat run, [[3]] after the third, etc.
################################################################################
plotFactorsAndRunComBat = function(
	trans_and_rg_data,
	preprocess_out,
	factors_to_plot = NULL,
	combat_factors_sequence = NULL,
	outfile = NULL
) {
	# default parameters
	if (is.null(combat_factors_sequence)) {
		combat_factors_sequence = c('Lib.constr.date', 'Tank', 'RNAseq.date');
	}
	if (is.null(factors_to_plot)) {
		factors_to_plot = list(
			c('Condition','Tank','Lib.constr.date','RNAseq.date'),
			c('Condition','Tank','LibSeq'),
			c('Condition','LibSeqTank')
		);
	}

	transcription_data = trans_and_rg_data$transcription_data;
	readgroup_data = trans_and_rg_data$readgroup_data;

	ncol = max(length(factors_to_plot), length(combat_factors_sequence))
	if (!is.null(outfile)) {
		png(filename = outfile, width = 400 * ncol, height = 300 * 3);
	}
	par(mfrow=c(3,ncol));
	num_blank_plots = ncol - length(factors_to_plot)


	# 1. plot for unfilitered data
	.plotFactorEffects(
		transcription_data = transcription_data,
		readgroup_data = readgroup_data,
		factors_to_plot = factors_to_plot,
		num_blank_plots = num_blank_plots,
		main = 'All subjects'
	);

	# 2. plot for preprocess()ed data
	new_readgroup_data = readgroup_data[match(names(preprocess_out$data_removedOutlierSamples), rownames(readgroup_data)), ];
	factor_effects = .plotFactorEffects(
		transcription_data = preprocess_out$data_removedOutlierSamples,
		readgroup_data = new_readgroup_data,
		factors_to_plot = factors_to_plot,
		num_blank_plots = num_blank_plots,
		main = 'Outlier subjects removed'
	);

	if (.verbosity >= 1) print(factor_effects);
	factor_effects = factor_effects[order(factor_effects)]
	factor_effects = factor_effects[names(factor_effects) != 'Condition']
	separate = factor_effects[names(factor_effects) != 'LibSeq'];
	together = factor_effects[!(names(factor_effects) %in% c('Lib.constr.date', 'RNAseq.date'))]
	if (separate[1] < together[1]) {
		factor_effects = separate;
	} else if (separate[1] > together[1]) {
		factor_effects = together;
	} else if (separate[2] <= together[2]) {
		factor_effects = separate;
	} else {
		factor_effects = together;
	}
	factor_effects = names(factor_effects)[1:min(3, length(factor_effects))]
	
	# 3. print some things
	factors_to_print_about = c('Lib.constr.date', 'RNAseq.date', 'Tank', 'LibSeq', 'LibSeqTank');
	for (factor in factors_to_print_about) {
		sort(table(paste(.getCondition(new_readgroup_data), .getReadgroupColumn(new_readgroup_data, factor))));
	}

	# 4. impute missing values and iteratively run ComBat
	combat_outputs = .runCombatForFactors(
		transcription_data = preprocess_out$data_Qnorm,
		readgroup_data = new_readgroup_data,
		combat_factors_sequence = factor_effects,
		factors_to_plot = factors_to_plot[1]
	);

	if (!is.null(outfile)) {
		dev.off();
	}
	
	return(list(combat_outputs = combat_outputs, readgroup_data = new_readgroup_data));
}



################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


################################################################################
# runPipeline
# 
# runPipeline(
# 	trans_and_rg_data = trans_and_rg_data,
# 	parameters = list(
# 		removeLowVarianceGenes.cv_threshold = 0.02,
#		plotFactorsAndRunComBat.factors_to_plot = list(
#			c('Condition','Tank','Lib.constr.date','RNAseq.date')
#		)
#	)
# );
#
# Description:
# 	Runs through the entire preprocessing pipeline, using the parameters in parameters.
#
# Arguments:
#	Required:
#	trans_and_rg_data   - A list containing:
#		transcription_data   - A data frame where each row is a gene and each column
#		                       is a subject, giving the raw transcripts per million (TPM)
#		                       expression level of each gene
#		readgroup_data       - A data frame with rownames equal to the column
#		                       names of transcription_data, containing columns with
#		                       read group data and (optionally) behavior data
#
#	Optional:
#	plot_outfile        - character; the file name to which the factor effects plots should
#	                      be saved.
#	parameters          - a list() of parameters for the substeps of the preprocessing pipeline.
#						  Parameters take the form functionName.parameter_name. The default
#						  value is used for any parameters that are not provided.
#		List of parameters you can provide:
#		removeRarelyExpressedGenes.max_fraction_zeroes
#		removeRarelyExpressedGenes.only_one_group_can_have_zeroes
#		removeLowVarianceGenes.cv_threshold
#		runPreprocessInteractive.deviate
#		runPreprocessInteractive.probe_thresh
#		runPreprocessInteractive.sample_thresh
#		runPreprocessInteractive.IACthresh
#		plotFactorsAndRunComBat.factors_to_plot
#		plotFactorsAndRunComBat.combat_factors_sequence
#	interactivePreproc  - Boolean; should the interactive preprocessing script be run interactively?
#
# Returns:
# 	Nothing.
################################################################################
runPipeline = function(
	trans_and_rg_data,
	plot_outfile = NULL,
	parameters = list(
	),
	interactivePreproc = FALSE
) {
	run_features = list();

	# BACKLOG getSubsets() - only female, male, etc.

	original_num_subjects = ncol(trans_and_rg_data)
	trans_and_rg_data <- removeExcludedGenesAndNormalize(trans_and_rg_data = trans_and_rg_data);
	run_features$n_genes_after_removeExcludedGenesAndNormalized = nrow(trans_and_rg_data[[1]])
	
	# BACKLOG remove excluded samples, or high TPM genes
	
	trans_and_rg_data <- removeRarelyExpressedGenes(
		trans_and_rg_data = trans_and_rg_data,
		max_fraction_zeroes = parameters$removeRarelyExpressedGenes.max_fraction_zeroes,
		only_one_group_can_have_zeroes = parameters$removeRarelyExpressedGenes.only_one_group_can_have_zeroes
	);
	run_features$n_genes_after_removeRarelyExpressedGenes = nrow(trans_and_rg_data[[1]])
	
	trans_and_rg_data <- removeLowVarianceGenes(
		trans_and_rg_data=trans_and_rg_data,
		cv_threshold = parameters$removeLowVarianceGenes.cv_threshold
	);
	run_features$n_genes_after_removeLowVarianceGenes = nrow(trans_and_rg_data[[1]])
	
	preprocess_out <- runPreprocessInteractive(
		trans_and_rg_data=trans_and_rg_data,
		deviate = parameters$runPreprocessInteractive.deviate,
		probe_thresh = parameters$runPreprocessInteractive.probe_thresh,
		sample_thresh = parameters$runPreprocessInteractive.sample_thresh,
		IACthresh = parameters$runPreprocessInteractive.IAC_thresh,
		interactive = interactivePreproc
	);
	# TEST - num subjects
	run_features$n_subjects = ncol(preprocess_out$data_Qnorm);
	stopifnot(run_features$n_subjects >= 0.66 * original_num_subjects)
	run_features$n_genes = nrow(preprocess_out$data_Qnorm);
	# TEST - meanIAC
	IACs = cor(preprocess_out$data_Qnorm, method="p", use="complete.obs");
	run_features$meanIAC = mean(IACs[upper.tri(IACs)]);
	stopifnot(run_features$meanIAC > 0.9);

	filtered_output <- plotFactorsAndRunComBat(
		trans_and_rg_data=trans_and_rg_data,
		preprocess_out=preprocess_out,
		factors_to_plot = parameters$plotFactorsAndRunComBat.factors_to_plot,
		combat_factors_sequence = parameters$plotFactorsAndRunComBat.combat_factors_sequence,
		outfile = plot_outfile
	);
	combat_out <- filtered_output$combat_output;
	filtered_readgroup_data <- filtered_output$readgroup_data;


	transcription_data = combat_out[[length(combat_out)]];
	readgroup_data = filtered_readgroup_data;
	p_values = .plotFactorEffects(
		transcription_data = transcription_data,
		readgroup_data = readgroup_data,
		factors_to_plot = c('Condition', 'Lib.constr.date', 'Tank', 'RNAseq.date'),
		num_blank_plots = 0,
		make_plot = FALSE
	);
	run_features$p_condition = p_values['Condition'];
	stopifnot(run_features$p_condition <= 0.1)
	p_others = p_values[-which(names(p_values) == 'Condition')];
	names(p_others) <- paste('p', names(p_others), sep = '_');
	if (sum(p_others[names(p_others) != 'p_Tank'] < 0.05) > 0) {
		first_sig_thing = which(p_others < 0.05)[1];
		stop(paste(names(p_others)[first_sig_thing], 'is < 0.05'));
	}
	run_features = c(run_features, p_others)
	return(list(run_features = run_features, transcription_data = transcription_data, readgroup_data = readgroup_data));
}



.writeParametersToLog = function(filename, params) {
	write('\n\n\n======================', filename, append = TRUE)
	for (i in 1:length(params)) {
		write(paste0(names(params)[i], ': ', paste(params[[i]], collapse = ',')), filename, append = TRUE)
	}
}


findIdealParameters = function(trans_and_rg_data, parameter_sets) {
	parametersAsString = function(params, long = FALSE) {
		string = paste(if (long) names(params) else '', as.character(params), sep = '_', collapse = '__');
		# get rid of parens
		string = gsub('[()]', '-', string);
		# get rid of commas
		string = gsub(',', '.', string);
		# get rid of anything weird
		string = gsub('[^-a-zA-Z_0-9.]', '', string);
		return(paste0('params__', string));
	}
	setVerbosity(1.6);

	log_file_name = paste0('find_ideal_parameters_log_', Sys.time());
	dir_name = paste0('parameter_results', Sys.time(), '/');
	dir.create(dir_name);
	write('Log\n', log_file_name);

	metrics = NA;
	for (params in parameter_sets) {
		outfile_pref = paste(dir_name, parametersAsString(params), sep = '');
		.writeParametersToLog(log_file_name, c(params, file_prefix = outfile_pref));

		cat('\n\n\n--------------\n', parametersAsString(params, long = TRUE), '\n--------------\n', sep = '');

		tryCatch({
			output = runPipeline(
				trans_and_rg_data,
				plot_outfile = paste(outfile_pref, 'png', sep = '.'),
				parameters = params,	
				interactivePreproc = FALSE
			);
			save(output, file = paste(outfile_pref, 'RData', sep = '.'));
			.writeParametersToLog(log_file_name, output$run_features);
			print(output$run_features)
			if (is.na(is.data.frame(metrics))) {
				metrics = data.frame(output$run_features, row.names = parametersAsString(params)) 
			} else {
				metrics = data.frame(rbind(metrics, output$run_features))
				rownames(metrics)[nrow(metrics)] = parametersAsString(params)
			}
		}, error = function(e) {
			print(e);
			write(as.character(e), file = log_file_name, append = TRUE)
			if (file.exists(paste(outfile_pref, 'png', sep = '.'))) {
				file.remove(paste(outfile_pref, 'png', sep = '.'));
			}
		});
	}
	print(metrics);
	return(metrics);
}


################################################################################
# .generateAllParameters
# 
# .generateAllParameters(
# 	params_left = list(
# 		removeLowVarianceGenes.cv_threshold = c(0.05, 0.1, 0.02),
#		plotFactorsAndRunComBat.factors_to_plot = list(
#			c('Condition','Tank','Lib.constr.date','RNAseq.date')
#		),
#		plotFactorsAndRunComBat.combat_factors_sequence = list(
#			c('Tank', 'Lib.constr.date', 'RNAseq.date'),
#			c('Tank', 'LibSeq'),
#			c('Lib.constr.date', 'RNAseq.date')
#		)
#	)
# );
#
# Description:
# 	Generates a list of parameter sets to iterate over
#
# Arguments:
#	Required:
#	params_left  - A list giving the parameters that should be added to current_set. The
#	               name of each entry is the parameter name, and the value is a list or vector of
#	               possible parameter values.
#
#	Optional:
#	current_set  - The current list of parameters (starts as list())
#
# Returns:
# 	A list() of list()s, where each list is a complete set of parameters, and all
# 	combinations of parameters are represented.
################################################################################
.generateAllParameters = function(params_left, current_set = list()) {
	if (length(params_left) == 0) {
		return(list(current_set));
	}
	nextParameterName = names(params_left)[1];
	nextParameterValues = params_left[[1]];
	params_left = params_left[-1];

	set_of_all_parameters = list();
	for (val in nextParameterValues) {
		current_set[[nextParameterName]] = val;
		set_of_all_parameters = c(set_of_all_parameters, .generateAllParameters(params_left, current_set));
	}
	return(set_of_all_parameters);
}










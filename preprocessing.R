##########################
# set up options
init = function() {
	options(stringsAsFactors=F);
	library(impute); # TODO ask AUSTIN - what happens if this isn't installed?

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
#	TODO implement
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
#                       This is a TEST_PARAMETER
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
removeExcludedGenesAndNormalize = function(trans_and_rg_data,
                                           genes_to_exclude = c(),
							               log_addition = 1) {
	transcription_data = trans_and_rg_data$transcription_data;

	# remove given genes
	if (length(genes_to_exclude)) {
		rowsToRemove = match(genes_to_exclude, rownames(transcription_data));
		transcription_data = transcription_data[-rowsToRemove,];
		.catlog('Removed ', length(genes_to_exclude), ' user-provided genes.\n' , sep = '', importance = 3)
	}

	# remove genes with zero variance
	.catlog('Removed ', sum(!(apply(transcription_data, 1, var) > 0)), ' genes with zero variance.\n' , sep = '', importance = 3)
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
removeRarelyExpressedGenes = function(trans_and_rg_data,
									  max_fraction_zeroes = 0.3,
									  only_one_group_can_have_zeroes = FALSE) {
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
		.catlog('Removing', length(rows_to_remove), 'genes with too many zero-expression values.\n', importance = 3)
		trans_and_rg_data$transcription_data = transcription_data[-rows_to_remove, ];
	} else {
		.catlog('No genes have too many zero-expression values.\n', importance = 3)
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
removeLowVarianceGenes = function(trans_and_rg_data,
								  make_plot = TRUE,
								  cv_threshold = 0.01
) {
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
	        ' below threshold ', cv_threshold, '.\n' , sep = '', importance = 3)
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
#		removeOutlierProbes   - TODO ask Austin
#		deviate               - TODO ask Austin
#		removeTooManyNAs      - TODO ask Austin
#		probe_thresh          - TODO ask Austin
#		sample_thresh         - TODO ask Austin
#		removeOutlierSamples  - TODO ask Austin
#		IACthresh             - TODO ask Austin
#		Qnorm                 - TODO ask Austin
#		All of these are TEST_PARAMETERS
#
# Returns: TODO better description - ask AUSTIN.
# 	A list including:
# 		outlierProbesOUTPUT, checkMissingDataOUTPUT, outlierSamplesOUTPUT; the outputs
# 			of several intermediate steps
# 		data_removedOutlierSample 
# 		data_Qnorm
# 		and others that appear not to be used.
################################################################################
runPreprocessInteractive = function(trans_and_rg_data,
	                                removeOutlierProbes=T,
						            deviate=2.5,
						            removeTooManyNAs=T,
						            probe_thresh=NULL,
						            sample_thresh=NULL,
						            removeOutlierSamples=T,
						            IACthresh=2,
						            Qnorm=T
) {
	preprocess_out = preProcess(datIN = trans_and_rg_data$transcription_data,
	                    removeOutlierProbes=removeOutlierProbes,
						deviate=deviate,
						removeTooManyNAs=removeTooManyNAs,
						probe_thresh=probe_thresh,
						sample_thresh=sample_thresh,
						removeOutlierSamples=removeOutlierSamples,
						IACthresh=IACthresh,
						Qnorm=Qnorm
	);
	# TODO this is painfully slow.
	return(preprocess_out);
}




################################################################################
# .computeAndPlotFactorEffectsOnMeanExpr
#
# Description:
# 	helper function for makeFactorEffectsPlots. Directly copy/pasted from
# 	Austin's old code file; if it works, do we really care what it does?
# 	TODO write description of plots produced
#
# Arguments:
# 	Required:
# 	TPM       - the transcription_data
# 	INFO      - a readgroup_info dataframe corresponding to TPM
# 	SUBSET    - numeric; which subjects should be considered? (vector of indices)
# 	INFOcols  - numeric; column indices of INFO that give the factors whose effects
# 	on mean expression should be considered.
# 	
# 	Optional:
# 	MAIN      - the title of the plot
# 	...       - any graphical parameters, to be passed on to the barplot function
#
# Returns:
# 	None.
################################################################################
.computeAndPlotFactorEffectsOnMeanExpr = function (TPM, INFO, SUBSET, INFOcols, MAIN='', ...) {
	DAT = TPM[, SUBSET];
	INFO = INFO[SUBSET, ];
	DAT = DAT[apply(DAT,1,var) > 0, ];
	meanexpr = apply(DAT, 2, mean, na.rm=T);
	factors = paste0(paste0('factor(',paste0('INFO[,',INFOcols,']' ),')' ), collapse='+');
	lmformula = paste0('meanexpr~', factors);print(lmformula)
	pvals = as.vector(na.omit(anova(lm(as.formula(lmformula)))$"Pr(>F)"));
	print(anova(lm(as.formula(lmformula))));
	names(pvals) = names(INFO)[INFOcols];
	barplot(-log10(pvals), ylab='-log10(pval)',main=paste(MAIN,' (n=',ncol(DAT),')',sep=''), ...);
	abline(h=c(-log10(.05), -log10(.01)), col='red');
}



################################################################################
# .runCombat
# 
# transcription_data = .runComBat(TODO parameters something something)
#
# Description:
# 	Runs the ComBat program from ComBat.R
#
# Arguments:
# 	TODO figure this out
#
# Returns: TODO better description - ask AUSTIN.
# 	I have literally no idea TODO
################################################################################
.runComBat = function (DAT, INFO, BATCHcol, COVcol=NULL, impute=F, ...) {
	if (impute) { DAT = as.data.frame(impute.knn(as.matrix(DAT))$data) }
	# format expression data and write to file
	datToWrite = as.data.frame(DAT);
	datToWrite = as.data.frame(cbind(gene=rownames(datToWrite), datToWrite));
	write.table(datToWrite, file='datForComBat.txt', quote=F, sep='\t', row.names=F);
	
	if (is.null(COVcol)) {
		inds = c(BATCHcol);
	} else if (is.numeric(COVcol)) {
		inds = c(BATCHcol, COVcol);
	} else {
		stop('arg COVcol must be NULL or numeric');
	}
	print(inds);
	infoToWrite = cbind(gsub('-', '.', rownames(INFO)), INFO[, inds]);
	infoToWrite = as.data.frame(cbind(names(DAT), infoToWrite));
	names(infoToWrite)[1:3] = c('Array name', 'Sample name', 'Batch');
	if (!is.null(COVcol) & is.numeric(COVcol)) {
		names(infoToWrite)[4:ncol(infoToWrite)] = paste('Covariate ', 1:length(COVcol), sep='');	
	}
	write.table(infoToWrite, file='infoForComBat.txt', quote=F, sep='\t', row.names=F);
	
	combatout = ComBat(expression_xls='datForComBat.txt', 
	                   sample_info_file='infoForComBat.txt',
					   filter=F,
					   write=F,
					   skip=1,
					   ...
	);
	test_combatout <<- combatout;
	if (is.character(combatout)) {
		stop(combatout);
	}
	return(combatout[,-1]);
}


################################################################################
# makeFactorEffectsPlots
# 
# makeFactorEffectsPlots(trans_and_rg_data=trans_and_rg_data, preprocess_out=preprocess_out)
#
# Description:
# 	No side effects; returns nothing.
# 	TODO write description of plots produced
# 	TODO refactor readgroup_data$* to use new getters
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
#	preprocess_out         - The value returned by runPreprocessInteractive
#
#	Optional:
#		TODO unknown.
#
# Returns:
# 	None.
################################################################################
makeFactorEffectsPlots = function(trans_and_rg_data, preprocess_out) {
	transcription_data = trans_and_rg_data$transcription_data;
	readgroup_data = trans_and_rg_data$readgroup_data;
	

	#dev.off(); TODO resolve this - dev.off() should belong with whoever dev.on()ed
	par(mfrow=c(3,3));

	# TODO make this a parameter
	factname_groups = list(c('Condition','Tank','Lib.constr.date','RNAseq.date'),
	                       c('Condition','Tank','LibSeq'),
						   c('Condition','LibSeqTank')
	);

	# 1. plot for unfilitered data
	num_subjects = ncol(transcription_data);
	
	for (facts_to_plot in factname_groups) {
		readgroup_columns = match(facts_to_plot, names(readgroup_data));
		.computeAndPlotFactorEffectsOnMeanExpr(TPM = transcription_data,
											   INFO = readgroup_data,
											   SUBSET = 1:num_subjects,
											   INFOcols = readgroup_columns,
											   MAIN = 'All subjects'
		);
	}

	# 2. plot for preprocess()ed data
	preNormPostOutliersDAT = preprocess_out$data_removedOutlierSamples;
	preNormPostOutliersINFO = readgroup_data[match(names(preprocess_out$data_removedOutlierSamples), rownames(readgroup_data)), ];
	num_subjects = ncol(preNormPostOutliersDAT);
	for (facts_to_plot in factname_groups) {
		readgroup_columns = match(facts_to_plot, names(preNormPostOutliersINFO));
		.computeAndPlotFactorEffectsOnMeanExpr(TPM = preNormPostOutliersDAT,
											   INFO = preNormPostOutliersINFO,
											   SUBSET = 1:num_subjects,
											   INFOcols = readgroup_columns,
											   MAIN = 'After preProcess() interactive script'
		);
	}

	# 3. print some things
	sort(table(paste(preNormPostOutliersINFO$Condition, preNormPostOutliersINFO$Lib.constr.date)));
	sort(table(paste(preNormPostOutliersINFO$Condition, preNormPostOutliersINFO$RNAseq.date)));
	sort(table(paste(preNormPostOutliersINFO$Condition, preNormPostOutliersINFO$Tank)));
	sort(table(paste(preNormPostOutliersINFO$Condition, preNormPostOutliersINFO$LibSeq)));
	sort(table(paste(preNormPostOutliersINFO$Condition, preNormPostOutliersINFO$LibSeqTank)));

	# 4. TODO there are still 3 empty plots in the plotting window!!!
	warning('There are three empty spots still in this plotting window. You\'re using an incomplete function!!');


	combat_outputs = list(preprocessed=preprocess_out$data_Qnorm);
	new_readgroup_data = preNormPostOutliersINFO; # TODO refactor
	column_indices_for_plot = match(c('Condition','Tank','Lib.constr.date','RNAseq.date'), names(new_readgroup_data));
	combat_factors_sequence = c('Lib.constr.date', 'Tank', 'RNAseq.date');
	covar_i = match('Condition', names(new_readgroup_data));
	for (factor in combat_factors_sequence) {
		new_transcription_data = combat_outputs[length(combat_outputs)];
		factor_i = match(factor, names(new_readgroup_data));
		out = .runComBat(DAT = new_transcription_data,
		                 INFO = new_readgroup_data,
						 BATCHcol = factor_i,
						 COVcol = covar_i,
						 impute = F, # TODO in original script this is T on first iteration and F after - why?
						 prior.plots = F
		);
		.computeAndPlotFactorEffectsOnMeanExpr(TPM = out,
											   INFO = new_readgroup_data,
											   SUBSET = 1:ncol(out),
											   INFOcols = column_indices_for_plot,
											   MAIN = paste('After ComBat on', factor)
		);
		combat_outputs[factor] = out;
	}



}



################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
runPipeline = function() {
	trans_and_rg_data = getAndValidateData();
	# BACKLOG getSubsets() - only female, male, etc.
	trans_and_rg_data <- removeExcludedGenesAndNormalize(trans_and_rg_data = trans_and_rg_data);
	# BACKLOG remove excluded samples, or high TPM genes
    trans_and_rg_data <- removeRarelyExpressedGenes(
        trans_and_rg_data = trans_and_rg_data,
        max_fraction_zeroes = 0,
        only_one_group_can_have_zeroes = TRUE
    );
	trans_and_rg_data <- removeLowVarianceGenes(trans_and_rg_data=trans_and_rg_data);
	# BACKLOG sanity check with prostaglandin
	preprocess_out <- runPreprocessInteractive(trans_and_rg_data=trans_and_rg_data);
	makeFactorEffectsPlots(trans_and_rg_data=trans_and_rg_data, preprocess_out=preprocess_out);
}
















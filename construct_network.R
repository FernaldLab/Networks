##########################
# set up options
init = function() {
	options(stringsAsFactors=F);

	# try to source libraries.
	tryCatch({
		library(WGCNA);
	},
	error = function(e) {
		warning(paste('Some dependencies of preprocessing.R could not be located. Please install',
		              'the required libraries.'));
		stop(e);
	});

	# fail gracefully if the user failed to setwd()
	tryCatch({
		source('utils.R');
	},
	warning = function(w) {
		stop(paste('Some dependencies of preprocessing.R could not be located. Please set the',
		           'working directory to be the directory containing the code files (using setwd()).'));
	});

	allowWGCNAThreads(nThreads = 4);
}
init();

# BACKLOG: write a function that's the equivalent of --help


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

setUpForDevelopment = function() {
	init();
	load('test_data/post_preproc_data.RData');
	pp_trans_and_rg_data <<- pp_trans_and_rg_data;
}

################################################################################
# runPickSoftThreshold
#
# Description:
#   Runs pickSoftThreshold from the WGCNA library
#
# Arguments:
#	Required:
#
#	Optional:
#	network_type        - what type of network should be constructed? One of
#	                      'unsigned', 'signed', and 'signed hybrid'.
#
# Returns:
################################################################################
runPickSoftThreshold = function(
	trans_and_rg_data,
	network_type = 'unsigned',
	verbose = 5
) {
	transcription_data = trans_and_rg_data$transcription_data;
	readgroup_data = trans_and_rg_data$readgroup_data;

	output = WGCNA::pickSoftThreshold(
		data = t(transcription_data),
		dataIsExpr = TRUE,
		networkType = network_type,
		moreNetworkConcepts = TRUE,
		verbose = verbose,
		RsquaredCut = 0.65
		# blockSize = NULL - NO, try smaller for faster computation. BUT, harder 2 deal with output.
	);

	.catlog('Power estimate is ', output$powerEstimate, '(NA means no power found)\n', sep = '');
	print(output$fitIndices);
	return(output$fitIndices);
}

# This function copied from https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-auto.pdf
plotThresholdResults = function(
	rpst_output # returned by runPickSoftThreshold
) {
	par(mfrow = c(1,2));
	powers = rpst_output$Power;
	cex1 = 0.9;
	sft = list(fitIndices = rpst_output);

	plot(
		sft$fitIndices[,1],
		-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
		xlab="Soft Threshold (power)",
		ylab="Scale Free Topology Model Fit,signed R^2",
		type="n",
		main = paste("Scale independence")
	);
	text(
		sft$fitIndices[,1],
		-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
		labels=powers,
		cex=cex1,
		col="red"
	);

	plot(
		sft$fitIndices[,1],
		sft$fitIndices[,5],
		xlab="Soft Threshold (power)",
		ylab="Mean Connectivity",
		type="n",
		main = paste("Mean connectivity")
	);
	text(
		sft$fitIndices[,1],
		sft$fitIndices[,5],
		labels=powers,
		cex=cex1,
		col="red"
	);
}


constructNetworkAndModules = function(
	trans_and_rg_data,
	beta
) {
	net = blockwiseModules(
		t(trans_and_rg_data$transcription_data),
		power = beta,
		TOMType = "unsigned", # TOM = Topological Overlap Matrix
		minModuleSize = 30,
		reassignThreshold = 0, # reassign genes that are closer to another cluster
		mergeCutHeight = 0.25,
		numericLabels = TRUE,
		pamRespectsDendro = FALSE,
		saveTOMs = FALSE,
		verbose = 3
	);
	return(net);
}


plotBasicDendrogram = function(
	constructNetworkAndModules_output
) {
	net = constructNetworkAndModules_output;
	#sizeGrWindow(12, 9);
	mergedColors = labels2colors(net$colors);
	plotDendroAndColors(
		net$dendrograms[[1]],
		mergedColors[net$blockGenes[[1]]],
		"Module colors",
		dendroLabels = FALSE,
		hang = 0.03,
		addGuide = TRUE,
		guideHang = 0.05
	);
}


save_module_data = function (
	constructNetworkAndModules_output,
	outfile = 'network_modules'
) {
	net = constructNetworkAndModules_output;
	moduleLabels = net$colors;
	moduleColors = labels2colors(net$colors);
	MEs = net$MEs;
	geneTree = net$dendrograms[[1]];
	save(MEs, moduleLabels, moduleColors, geneTree,
		file = paste0(outfile, ".RData")
	)
}


# Target:
# r^2 > 0.8 (at least 0.7)
# slope between -1 and -2
# mean connectivity at least 15 to 20
# lower is better




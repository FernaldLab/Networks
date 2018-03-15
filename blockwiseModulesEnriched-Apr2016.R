# functions for iteratively removing background genes from network using WGNCA blockwiseModules() function
# Austin Hilliard, White Lab, UCLA, September 2011
#  revised in Fernald lab, Stanford, February 2013
#  blockwiseModulesEnrichedIterate added in Fernald lab, September 2014
#  extensive comments added in Fernald lab, April 2016

# =======================================================================================================================================================
# This function is a wrapper for blockwiseModulesEnriched, which is itself a wrapper for WGCNA::blockwiseModules
# -------------------------------------------------------------------------------------------------------------------------------------------------------
# Dependencies are the same as blockwiseModulesEnriched
# Nothing is returned within R. Network, data, dendrograms, density testing results are saved out
# =======================================================================================================================================================

blockwiseModulesEnrichedIterate = function( DATA, 
                                            networkType='signed', 
                                            power=14, 
                                            maxBlockSize=ncol(DATA)+1, 
					    deepSplitVec=c(2,4), 
					    mergeCutHeight=0.1, 
					    minModuleSizeVec=c(10,20,40,80,100), 
				            minKMEtoStayVec=c(0.3, 0.5),
		                            minCoreKMEVec=c(0.5, 0.7),
					    densityPermTest=F, 
					    skipThresh=300, 
					    ...) {
	
	for (DS in deepSplitVec) {
		for (MM in minModuleSizeVec) {
			for (MkME in minKMEtoStayVec) {
				for (MCkME in minCoreKMEVec) {
					saveFileBase = paste0(deparse(substitute(DATA)), '_',networkType,
					                      '_p',power,'_ds',DS,'_mm',MM,'_mch',mergeCutHeight,'_mKME',MkME,'_mCoreKME',MCkME);
					cat('Building network: deepSplit=', DS,
					    ', mergeCutHeight=', mergeCutHeight,
					    ', minModuleSize=', MM, 
					    ', minKMEtoStay=', MkME,
					    ', minCoreKME=', MCkME, 
					    '\n  filebase: ', saveFileBase,
					    '\n',
					    sep=''
					    );
					net = blockwiseModulesEnriched(DATA=DATA, maxBlockSize=maxBlockSize, power=power, networkType=networkType, minModuleSize=MM,
											   	   deepSplit=DS, minKMEtoStay=MkME, minCoreKME=MCkME, mergeCutHeight=mergeCutHeight,
											   	   densityPermTest=densityPermTest, skipThresh=skipThresh, verbose=0,
											  	   saveFileBase=saveFileBase, ...
											  	   );
				}
			}
		}
	}	
}

# =======================================================================================================================================================
# This function is a wrapper for WGCNA::blockwiseModules with optional module density permutation testing
# -------------------------------------------------------------------------------------------------------------------------------------------------------
# It will build an initial network using blockwiseModules and identify any background genes
#  (grey genes and genes from modules that failed the optional density testing)
# Background genes are removed and another network is built, and this repeats until no background genes remain
# If requested, a module density permutation test will be performed for each network iteration
#  the test essentially resamples the TOM to build a null distribution of density for modules of a certain size
#
# The density test is kludgy and could be improved, both in method and runtime
#  Possible strategies:
#   - resample the expression data and build smaller TOMs rather than dealing with one giant TOM or trying to keep track of blocks
#      would lose the actual TOM structure from the entire network since TO's not just a function of pairwise relationships like adjacency
#   - use the adjacency matrix instead of the TOM
#      would be much faster and results would probably be pretty much the same
#   - rather than running nPerm permutations for each module, do nPerm permuations overall but resample for every module during each run
#      could greatly increase speed but resampling would be less randomized
# -------------------------------------------------------------------------------------------------------------------------------------------------------
# Depends: From WGCNA library- blockwiseModules, plotDendroAndColors, TOMsimilarityFromExpr, collectGarbage
#          From this document- getModDensitiesFromTOM, modDensityPerm
# 
# Returns: Final list of network objects output by WGCNA::blockwiseModules
#  Also optionally saves .RData network and data objects and jpgs of dendrograms
# =======================================================================================================================================================

blockwiseModulesEnriched = function( DATA,                    # data frame of expression values, columns should be genes        
								                     maxBlockSize=6000,       # maxBlockSize value for blockwiseModules
								                     fixedBlockSize=TRUE,     # ignore for now (keep set to TRUE), experimental for adjusting maxBlockSize automatically
								                     networkType='signed',    # string defining signed/unsigned network type for blockwiseModules
								                     power=14,                # integer power value for blockwiseModules
								                     deepSplit=2,             # integer deepSplit value for blockwiseModules
								                     minModuleSize=10,        # integer minModuleSize value for blockwiseModules
								                     smartMinModThresh=2,     # ignore
								                     smartMinModMultiplier=2, # ignore
								                     densityPermTest=TRUE,    # boolean indicating whether to perform module density permutation testing
								                     plotModDensities=TRUE,   # if densityPermTest = T, plot module densities as a function of module size after each run
								                     nPerm=1000,              # if densityPermTest = T, integer to set number of permutations 
								                     permTestPvalThresh=.01,  # if densityPermTest = T, numeric p-value threshold for passing the density test
								                     skipThresh=500,          # if densityPermTest = T, integer defining module size required to skip density test
								                     skipGrey=TRUE,           # if densityPermTest = T, boolean defining whether to test density of grey genes
								                     onlyGreyThresh=3,        # if densityPermTest = T, integer defining number of runs with only grey background genes before density test is canceled 
								                     preTOM=NULL,             # if densityPermTest = T, user can provide pre-computed TOM matrix 
								                     verbose=1,               # integer verbose value for blockwiseModules
								                     saveNets=TRUE,           # boolean indicating whether to save network object from each run into a .RData file
								                     saveDendros=TRUE,        # boolean indicating whether to save jpgs of dendrograms from each run, one per block
								                     saveDATA=TRUE,           # boolean indicating whether to save input data from each run into a .RData file
								                     savePermOut=TRUE,        # if densityPermTest = T, boolean indicating whether to save density test results from each run into a .RData file
								                     saveFileBase='',         # string to serve as prefix for any output file names
								                     setRun=NULL,             # integer to set initial run number (mostly used if function had to be quit and restarted)
								                     graphicsType='quartz',   # if saveDendros = T, string defining graphics device 
								                     res=150,                 # if saveDendros = T, integer defining jpg resolution
								                     ...                      # additional arguments to blockwiseModules
								                    ) {
  
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Prepare for running blockwiseModules in a while-loop 
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # make sure columns are genes and rows are samples
	if (ncol(DATA) < nrow(DATA)) {
		DATA = as.data.frame( t(DATA) );
		cat('\n...transposing input data so genes are in columns\n');
	}
	
  # initialize values that control the while-loop 
  #  bgGenes is number of background genes to be removed
  #   always includes grey genes, and if densityPermTest = T, also includes genes from modules that fail density test
  #  onlyGrey is number of previous consecutive runs where no modules failed the density testing and only grey genes were removed
	bgGenes = 1;	
	onlyGrey = 0;
	
	# if a run number is provided via setRun, adjust value of run variable for console messages and output filenames 
	if (is.numeric(setRun)) { 
		run = setRun - 1;
	} else {
		run = 0;
	}
	
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Iteratively build/rebuild network until no background genes remain 
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	while (sum(bgGenes) > 0) {
		
	  # get number of genes in the data for this run, only relevant for experimental maxBlockSize adjustments
		numProbes = ncol(DATA);
		
		# update run number and print to console
		run = run + 1;
		cat('\nRUN ', run, '...\n', sep='');
		
		# if enough previous runs had all modules pass the density test then cancel density testing, provides massive speed increase
		if (onlyGrey >= onlyGreyThresh) {
			cat('   ....skipping perm tests in remaining runs\n');
			densityPermTest = FALSE;
		}
		
	# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
	# Experimental maxBlockSize adjustments, ignore for now
	# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
		if (numProbes > maxBlockSize & !fixedBlockSize) {
			adjustedBlockSize = (ceiling((numProbes/1000)) * 1000) / 2;
			cat('   ...maxBlockSize adjusted to ', adjustedBlockSize, '\n', sep='');
		} else {
			adjustedBlockSize = maxBlockSize;
			cat('   ...using maxBlockSize = ', adjustedBlockSize, '\n', sep='');
		}
		
	# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
	# Run blockwiseModules after determining how to handle TOM, potentially save data and network out to .RData files
	# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
		# if it's possible to build the network in one block, density test is on, and no TOM was provided,
		#  then have blockwiseModules save the TOM temporarily to use in the density test
		# otherwise don't save the TOM
		#  if the density test is on but no TOM is provided then it will be recomputed in one block
		#  I had issues with reloading multiple saved TOMs from different blocks so this is the not-so-great workaround
		#  I recommend building the network in a single block if possible, or if mulitple blocks are required not running the density test
		if ( (numProbes < adjustedBlockSize) & densityPermTest & is.null(preTOM) ) {
			saveTOMs = TRUE;
			saveTOMFileBase = 'tmpTOM';
		} else {
			saveTOMs = FALSE;
			saveTOMFileBase = NULL;
		}
		
		cat('...constructing network\n');
		
		# save the expression data from this run if requested
		if (saveDATA) { save(DATA, file=paste0(saveFileBase, 'run', run, 'DATA.RData')) }
		
		# run blockwiseModules 
		net = blockwiseModules(DATA,
							             maxBlockSize=adjustedBlockSize,
							             networkType=networkType,
							             power=power,
							             deepSplit=deepSplit,
							             minModuleSize=minModuleSize,
							             saveTOMs=saveTOMs,
							             saveTOMFileBase=saveTOMFileBase,
							             verbose=verbose, 
							             ...);
		
		# save the network object from this run if requested
		if (saveNets) { save(net, file=paste0(saveFileBase, 'run', run, 'NET.RData')) }
		
		# get module sizes for printing to console
		modSizes = sort(table(net$colors), decreasing=T);
		cat('\nmodules:\n', sep='');
		print(modSizes);
		
		# if ( is.numeric(smartMinModThresh) ) {
			
		# }
		
		collectGarbage();
		
	# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
	# Save jpgs of dendrogram(s) if requested, will be one file per block
	# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
		if (saveDendros) {
			cat('   ...saving dendrogram(s) as .jpgs\n');
			for (block in 1:length(net$blockGenes)) {
				jpeg(file=paste0(saveFileBase, 'run', run, 'dendro-block', block, '.jpg'), 
					   width=15, height=6, units='in', quality=100, 
					   type=graphicsType, 
					   res=res);
				plotDendroAndColors(net$dendrograms[[block]],
					                  net$colors[ net$blockGenes[[block]] ],
					                  groupLabels='module',
					                  rowText=net$colors[ net$blockGenes[[block]] ],
					                  main=paste0('block ', block, ': ', length(net$blockGenes[[block]]), ' genes'),
					                  dendroLabels=FALSE, hang=.05, addGuide=TRUE, guideHang=.05);
				dev.off();
			}
		}
		
  # ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Run module density testing if requested
  # ----------------------------------------------------------------------------------------------------------------------------------------------------------------------  
		if (densityPermTest) {
		  
		# --------------------------------------------------------------------------------------------------------------------------------------------------------------------
	  # Get the TOM
	  # --------------------------------------------------------------------------------------------------------------------------------------------------------------------  
		  # if a single-block TOM from blockwiseModules was saved, reload it
			if (saveTOMs) {
				cat('   ...loading TOM\n');
				load('tmpTOM-block.1.RData');
				TOM = as.matrix(TOM);
				
			# otherwise if a TOM was provided, use it
			#  there should be some checking here to ensure the TOM is valid 
			} else if ( !(is.null(preTOM)) ) {
				cat('   ...using preTOM\n');
				TOM = preTOM;	
				
			# otherwise recompute the TOM 
			#  this could be massive if the network had to be built in multiple blocks, try to avoid
			} else {
				cat('   ...re-computing TOM\n');
				TOM = TOMsimilarityFromExpr(DATA, networkType=networkType, power=power);
			}
			collectGarbage();
			
		# --------------------------------------------------------------------------------------------------------------------------------------------------------------------
		# Run the test and remove any background genes
		# --------------------------------------------------------------------------------------------------------------------------------------------------------------------  
			# get actual density of each module
			modDensities = getModDensitiesFromTOM(TOM, net$colors, plot=plotModDensities, skipGrey=skipGrey);
			
			# run permutation tests
			cat('\n   ...working on density perm test\n');
			permTest = modDensityPerm(TOM, net$colors, modDensities, nPerm=nPerm, skipThresh=skipThresh, skipGrey=skipGrey);
			
			# get names of any modules that failed the density test
			weakTOM = names( table(net$colors)[permTest$pvals > permTestPvalThresh] );
			cat('      ...weakTOM: ', weakTOM, '\n', sep='');                             # need to fix printing, no spaces between module names, use paste0 with collapse=' '
			
			# remove background genes from data for next iteration, includes genes from failed modules and any grey genes
			bgGenes = net$colors %in% weakTOM;
			cat('   ......removing ', sum(bgGenes), ' background genes\n', sep='');
			DATA = DATA[, !bgGenes];
			
			# save results of permutation test from this run if requested
			if (savePermOut) {
				permOut = list(modDensities=modDensities, permTest=permTest, weakTOM=weakTOM);
				save(permOut, file = paste0(saveFileBase, 'run', run, 'PermOut.RData'));
			}
			
			# check if background genes were all grey
			#  if yes increase value of onlyGrey by 1
			#  otherwise reset onlyGrey to 0
			if ( length(weakTOM)==1 ) { 
				if (weakTOM=='grey') { onlyGrey = onlyGrey + 1 }
			} else {
				cat('      ....re-setting onlyGrey...\n');
				onlyGrey = 0;
			}
			cat('      onlyGrey = ', onlyGrey, '\n', sep='');
			
	# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
	# If module density testing is not requested then just remove grey genes
	# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------  
		} else {
			bgGenes = net$colors == 'grey';
			cat('\n   skipping permutation test...\n');
			cat('    ...removing ', sum(bgGenes), ' grey module genes\n', sep='');
			DATA = DATA[, !bgGenes];
		}
		
		collectGarbage();
		
	# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
	# If there are no background genes print message to console, while-loop is broken
	# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------  
		if (sum(bgGenes) == 0) { cat('\nno more background genes, all done\n') }
	}
	
	# return blockwiseModules output
	return(net);				
	
}

# =======================================================================================================================================================
# This function uses a TOM and a vector of module assignments to compute average TO for each module
# -------------------------------------------------------------------------------------------------------------------------------------------------------
# There's no reason the TOM couldn't be an adjacency matrix, or any other kind of symmetric similarity matrix, instead
# -------------------------------------------------------------------------------------------------------------------------------------------------------
# Depends: From WGCNA library- vectorizeMatrix, collectGarbage
# Returns: Named numeric vector of module densities
# =======================================================================================================================================================

getModDensitiesFromTOM = function( TOM,          # numeric topological overlap matrix
                                   colors,       # character vector containing module assignments for the genes in TOM
                                   plot=TRUE,    # boolean indicating whether to plot module densities as a function of module size
                                   diag=FALSE,   # boolean indicating whether to include the values in the diagonal of the TOM
                                   skipGrey=TRUE # boolean indicating whether to ignore grey genes and just assign them density=0
                                  ) {
	
  # get the numbers of modules and genes in each module
	modTable = table(colors);
	nMods = length(modTable);
	
	# loop through modules to build vector of module densities
	modDensities = c();
	for (m in 1:nMods) {
	  
	  # if the current module is grey and skipGrey=T, assign density of 0 and move to next module
	  checkGrey = names(modTable)[m] == 'grey';
		if (skipGrey & checkGrey) {
			modDensities = c(modDensities, 0);
			next;
		}
	  
	  # otherwise get the indices for genes in the current module
		modGenes = colors == names(modTable)[m];
		
		# subset the TOM for genes in the current module and compute the density as mean TO
		modTOM = TOM[modGenes, modGenes];
		modDensities = c(modDensities, mean(vectorizeMatrix(modTOM, diag=diag)));
	}
	names(modDensities) = names(modTable);
	
	# plot module densities as a function of module size using default graphics device
	#  e.g. in RStudio it will show up in the plot window
	if (plot) {
		plot(as.vector(modTable), modDensities, type='n', xlab='module size', ylab='module density');
		text(as.vector(modTable), modDensities, labels=names(modDensities));
	}
	
	collectGarbage();
	return(modDensities);
	
} 
	
# =======================================================================================================================================================
# This function resamples a TOM to build randomized distributions of average density for every module based on their sizes 
# -------------------------------------------------------------------------------------------------------------------------------------------------------
# See the comments of blockwiseModulesEnriched above for possible ways to improve this procedure
#
# Towards the end of the function is a block of code to deal with NA pvals
#  I don't remember why I put this in and therefore find it disconcerting
#  Resampling is done without replacement so it can't be a too-small TOM issue
#  No values in the TOM should ever be 0, so it can't be a divide-by-0 issue
#  Probably is something trivial, from some very specific siutation, or I'm not reading my own code correctly, 
#   just need to test with some print statements inserted to figure out when it happens
#
# -------------------------------------------------------------------------------------------------------------------------------------------------------
# Depends: From WGCNA library- collectGarbage (this could just be replaced with gc())
# Returns: 2-element list containing a numeric vector of pvalues and a list of pseudo-density distributions for each module
# =======================================================================================================================================================

modDensityPerm = function( TOM,             # numeric topological overlap matrix
                           colors,          # character vector containing module assignments for the genes in TOM
                           modDensities,    # named numeric vector of module densities, usually output from getModDensitiesFromTOM
                           nPerm=1000,      # number of permutations to run for each module
                           diag=FALSE,      # boolean indicating whether to include the values in the diagonal of the TOM
                           skipThresh=5000, # integer defining module size required to skip density test
                           skipGrey=TRUE    # boolean indicating whether to ignore grey genes and just assign them density=0
                          ) {
	
  # get the numbers of modules and genes in each module
	modTable = table(colors);
	nMods = length(modTable);
	
	# check that the colors and modDensities arguments contain compatible information
	test = sum( names(modTable)==names(modDensities) ) == nMods;
	if(!test) { error('module names don\'t match colors') }
	
	# loop through modules to perform permutation test
	#  null distributions will be stored in the list modPseudoDensities   
	#  p-values will be stored in the vector pvals
	modPseudoDensities = list();
	pvals = c();
	for (m in 1:nMods) {
	  
	  # get current module
		mod = names(modTable)[m];
		
		# if the current module is grey and skipGrey=T, set modPseudoDensities$grey = NA, set pval = 1, then move on
		if (skipGrey & (mod == 'grey')) {
		  cat('skipping ', mod, ': ', modTable[m], ' probes\n', sep = '');
		  modPseudoDensities[[m]] = rep(NA, nPerm);
		  pvals = c(pvals, 1);
		  cat(' pval = ', 1, '\n', sep = '');
		  next;
		}
		
		# if current module is larger than skipThresh, set modPseudoDensities[[m]] = NA, set pval = 0, then move on
		if (modTable[m] > skipThresh) {
		  cat('skipping ', mod, ': ', modTable[m], ' probes\n', sep = '');
		  modPseudoDensities[[m]] = rep(NA, nPerm);
		  pvals = c(pvals, 0);
		  cat(' pval = ', 0, '\n', sep = '');
		  next;
		}
		
		# run permutation test loop to build a vector of randomized average module densities
		cat('working on ', mod, ': ', modTable[m], ' probes\n', sep = '');
		pseudoDensities = c();
		for (p in 1:nPerm) {
		  
		  # print update to console, should change so 1000 isn't hard-coded here
			if (p %% 1000 == 0) { cat('  permutation ', p, '\n', sep = '') }   
		  
		  # randomly select n genes from TOM, where n is size of the current module, then get mean TO
			pseudoGenes = sample(1:ncol(TOM), modTable[m]);
			pseudoTOM = TOM[pseudoGenes, pseudoGenes];
			pseudoDensities = c(pseudoDensities, mean(vectorizeMatrix(pseudoTOM, diag=diag)));
		}
		
		# add pseudo-densities for current module to the modPseudoDensities list
		modPseudoDensities[[m]] = pseudoDensities;
		
		# count number of pseudo-densities greater than the actual density of current module to get pvalue
		p = (sum(pseudoDensities > modDensities[m])) / nPerm;
		pvals = c(pvals, p);
		cat(' pval = ', p, '\n', sep = '');
		collectGarbage();
	}
	names(modPseudoDensities) = names(modTable);
	names(pvals) = names(modTable);
	
	# check for any NAs in the pvalue vector and hunt them down in the modPseudoDensities list
	#  Could just add na.rm=T to the sum(pseudoDensities > modDensities[m]) calculation but...
	#  I honestly don't remember why this would happen, i.e. where the NAs are coming from 
	#   could be from the skipping of large/grey modules, but then it wouldn't make sense to replace NAs with the median 
	#   so it seems like there are individual permutations that return NA as a pseudo-density, 
	#   maybe there are NAs in the TOM for some reason?
	permNA = which(is.na(pvals));
	if (length(permNA) > 0) {
		for (p in 1:length(permNA)) {
			permNAind = which( is.na( modPseudoDensities[[ permNA[p] ]] ) );
			permMed = median( modPseudoDensities[[ permNA[p] ]], na.rm=T );
			modPseudoDensities[[ permNA[p] ]][permNAind] = permMed;
			greater = sum( modPseudoDensities[[permNA[p]]] > modDensities[permNA[p]] );
			pvals[permNA[p]] = greater / nPerm;
		}
	}
	
	# return list containing a numeric vector of pvalues and a list of pseudo-density distributions for each module
	out = list(pvals=pvals, modPseudoDensities=modPseudoDensities);
	return(out);
	
}
	

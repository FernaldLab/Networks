# Networks

## Parameter Search Directions
- Load preprocessing code into R
- Necessary files: preprocessing.R, utils.R, preProcess_interactive.R, ComBat.R (all in same directory)
~~~
source('preprocessing.R')
~~~

- Load transcription data and readgroup data (see `getAndValidateData()`)
  - `getAndValidateData()` does this for .tsv files with the exact same format as 'test_data/transcription_data.tsv' and 'test_data/readgroupsdata.tsv'
  - Or you can just make your own `list` where
    - `list$transcription_data` is a data frame where each row is a gene and each column is a subject, giving the transcripts per million (TPM) expression level of each gene/transcript
    - `list$readgroup_data` is a data frame where each row is a subject and each column is a piece of readgroup data or behavior data.
    - `rownames(list$readgroup_data)` should be identical to `colnames(list$transcription_data)`; see `.matchSubjectNames()` for a function that does this
    - `list$readgroup_data` should have a column `Condition` giving social status of each fish

- Load different parameter options
  - To use the same ones as before: `source('parameter_sets.R')` (they end up in variable `parameter_sets_fixed`)
  - To generate your own set:
~~~
my_parameter_options <- .generateAllParameters(
  params_left = list(
 		removeLowVarianceGenes.cv_threshold = c(0.05, 0.1, 0.2),
		plotFactorsAndRunComBat.factors_to_plot = list(
			c('Condition','Tank','Lib.constr.date','RNAseq.date')
		),
		plotFactorsAndRunComBat.combat_factors_sequence = list(
			c('Tank', 'Lib.constr.date', 'RNAseq.date'),
			c('Tank', 'LibSeq'),
			c('Lib.constr.date', 'RNAseq.date')
		)
	)
);
~~~
  - Of note here:
    - Every entry of the list `params_left` represents all of the values of a parameter that should be tested. The parameter `removeLowVarianceGenes.cv_threshold` can be either 0.05, 0.1, or 0.2; `plotFactorsAndRunComBat.combat_factors_sequence` can be either `c('Tank', 'Lib.constr.date', 'RNAseq.date')`,  `c('Tank', 'LibSeq')`, or `c('Lib.constr.date', 'RNAseq.date')`.
    - If you don't provide options for a parameter, the default value for that parameter will be used.
    - Parameters that can be tested are marked in preprocessing.R with a comment that says `TEST_PARAMETER` (edit-find for these). They all have descriptions explaining exactly what they control.

- Run code. Make sure you're in a directory where you feel comfortable writing a folder full of data and a log file.
~~~
results_matrix = findIdealParameters(trans_and_rg_data=trans_and_rg_data, parameter_sets=parameter_sets_fixed)
dump('results_matrix', file='backup_file.R')
~~~
  - `results_matrix` is a data frame giving, for each successful run, the number of genes left after each preprocessing step, the number of subjects left, the mean IAC, and the p-values of each of the ComBat factors and of Condition.
  - For each successful run, a plot is written to the output folder, and an .RData file. To load the results of a run, `load()` the file, and the results will be in variable `output`. `output` contains `run_features` (the data frame row for this run), and `transcription_data` and `readgroup_data` (the preprocessed data).
  - "Success" of a run is determined in `runPipeline()` - see all of the `stopifnot()`s.


- The next coding tasks are marked with `# BACKLOG`

- A start on network construction code is in construct_network.R



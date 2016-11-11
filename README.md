# Networks
0.
~~~
source('preprocessing.R')
~~~

1. Load transcription data and readgroup data (see `getAndValidateData()`)
  - subject names should match exactly (`.matchSubjectNames()`)
  - return `list()` called `trans_and_rg_data`:
    \#   A list containing:
    \#   transcription_data - A data frame where each row is a gene and each column
    \#                        is a subject, giving the transcripts per million (TPM)
    \#                        expression level of each gene
    \#   readgroup_data     - A data frame with rownames equal to the column
    \#                        names of transcription_data, containing columns with
    \#                        read group data and (optionally) behavior data


2. Load parameters
  `source('parameter_sets.R')`
  (they end up in variable `parameter_sets_fixed`)


3. Make sure you're in a directory where you feel comfortable writing a folder full of data and a log file

4.
~~~
results_matrix = findIdealParameters(trans_and_rg_data, parameter_sets_fixed)
dump('results_matrix', file='backup_file.R')
~~~



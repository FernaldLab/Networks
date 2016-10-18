# Returns a file name selected interactively.
.getFile = function(prompt) {
	cat(prompt);
	return(file.choose());
}



# Asserts that exactly one of the arguments is provided, and the rest are NA.
# Throws a warning if any of the arguments is a function; just be aware!!
.assertOne = function(...,
                      error_message_none = 'No arguments provided',
					  error_message_two = 'Too many arguments provided') {
	argument_list <- list(...);
	num_given = sum(lapply(argument_list, !is.na));
	if (num_given == 0) {
		stop(error_message_none);
	} else if (num_given > 1) {
		stop(error_message_two);
	}
}


.coefficientOfVariance = function(x) {
	stdev = sd(as.numeric(x), na.rm=TRUE);
	avg = mean(as.numeric(x), na.rm=TRUE);
	return(stdev / avg);
}


# catlog cats a progress message
# defined here so that ALL messages can be redirected to an error/log file, or be
# silenced, or prioritized
# importance - 1 will always get printed, 3 is for only super verbose
# 0 is so important that it ALWAYS is copied to stdout.
setVerbosity = function(threshold = 5) {
	.catlog <<- function(..., importance = 1) {
		if (importance < threshold) {
			cat(...)
		}
	}
}
setVerbosity();


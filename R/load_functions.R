#' Compile sequence lengths from a FASTA format file
#' 
#' @param file The name of the file containing the information to import
#' 
#' @export
#' 
#' @return A database containing the lengths of each sequence.
#' 
load_lengths <- function(file) {
	# check line ending type
	first_line <- readLines(file, n = 1)
	if (length(first_line) == 0)	
		stop("The database file appears to be empty. Are you sure this is the correct file?", call. = FALSE)

	len_first_line <- nchar(first_line)
	if (grepl("\r$", readChar(file, len_first_line + 1)))
		stop("The line endings of your database file are incompatible with supeRbaits. You can convert your database using convert_line_endings()\n", call. = FALSE)
	# --
	# extract sequence lengths
	getlengths.time <- system.time({
		the.lengths <- callr::r(function(getChromLengths, path) 
			{
				getChromLengths(path = path)
			}, 
			args = list(getChromLengths = getChromLengths,
									path = file),
			spinner = TRUE,
			show = TRUE)
	})
	if (getOption("supeRbaits_show_times", default = FALSE)) {
		attributes(the.lengths)$time_elapsed <- getlengths.time["elapsed"]
		print(getlengths.time)
	}

	# failsafe for r-oldrel
	to.convert <- which(unlist(lapply(the.lengths, class)) == "factor")

	if (length(to.convert) > 0) {
		for (i in to.convert) {
			the.lengths[, i] <- as.character(the.lengths[, i])
		}
	}
	# ---

	return(the.lengths)
}

#' Load file with areas of the chromosome to avoid
#' 
#' @param file The name of the file containing the information to upload
#' 
#' @keywords internal
#' 
#' @return The exclusions dataframe
#' 
load_exclusions <- function(file){
	output <- as.data.frame(data.table::fread(file, showProgress = FALSE), stringsAsFactors = FALSE)
	if (ncol(output) != 3)
		stop("The exclusions file does not appear to contain exactly three columns\n", call. = FALSE)
	colnames(output) <- c("chr", "start", "stop")
	if (!is.integer(output$start))
		stop("The second column in the exclusions file does not appear to contain only numeric data.\n", call. = FALSE)
	if (!is.integer(output$stop))
		stop("The third column in the exclusions file does not appear to contain only numeric data.\n", call. = FALSE)
	if (any(is.na(output)))
		stop("There appear to be NA's in the exclusions file.\n", call. = FALSE)
	if (any(output[, 2:3] < 0))
		stop("Some data in the exclusions file appears to be negative.\n", call. = FALSE)
	if(any(output$start > output$stop))
		stop("Not all starting points are greater than the ending points in the exclusions file.\n", call. = FALSE)
	return(output)
}

#' Load file with areas of the chromosome to target
#' 
#' @inheritParams load_exclusions
#' 
#' @keywords internal
#' 
#' @return The regions dataframe
#' 
load_regions <- function(file){
	output <- as.data.frame(data.table::fread(file, showProgress = FALSE), stringsAsFactors = FALSE)
	if (ncol(output) != 3)
		stop("The regions file does not appear to contain exactly three columns\n", call. = FALSE)
	colnames(output) <- c("chr", "start", "stop")
	if (!is.integer(output$start))
		stop("The second column in the regions of interest file does not appear to contain only numeric data.\n", call. = FALSE)
	if (!is.integer(output$stop))
		stop("The third column in the regions of interest file does not appear to contain only numeric data.\n", call. = FALSE)
	if (any(is.na(output)))
		stop("There appear to be NA's in the regions of interest file.\n", call. = FALSE)
	if (any(output[, 2:3] < 0))
		stop("Some data in the regions of interest file appears to be negative.\n", call. = FALSE)
	if(any(output$start > output$stop))
		stop("Not all starting points are greater than the ending points in the regions of interest file.\n", call. = FALSE)
	return(output)
}

#' Load file with points in the chromosome to target
#' 
#' @inheritParams load_exclusions
#' 
#' @keywords internal
#' 
#' @return The targets dataframe
#' 
load_targets <- function(file){
	output <- as.data.frame(data.table::fread(file, showProgress = FALSE), stringsAsFactors = FALSE)
	if (ncol(output) != 2)
		stop("The exclusions file does not appear to contain exactly two columns\n", call. = FALSE)
	colnames(output) <- c("chr", "target")
	if (!is.integer(output$target))
		stop("The second column in the targets file does not appear to contain only numeric data.\n", call. = FALSE)
	if (any(is.na(output)))
		stop("There appear to be NA's in the targets file.\n", call. = FALSE)
	if (any(output$targets < 0))
		stop("Some data in the targets file appears to be negative.\n", call. = FALSE)
	return(output)
}

#' Load file with areas of the chromosome to avoid
#' 
#' @param file The name of the file containing the information to upload
#' 
#' @keywords internal
#' 
#' @return The exclusions dataframe
#' 
load_exclusions <- function(file){
	output <- as.data.frame(data.table::fread(file, showProgress = FALSE))
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
	output <- as.data.frame(data.table::fread(file, showProgress = FALSE))
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
	output <- as.data.frame(data.table::fread(file, showProgress = FALSE))
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

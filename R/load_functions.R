#' Load file with areas of the chromosome to avoid
#' 
#' @param file The name of the file containing the information to upload
#' 
#' @keywords internal
#' 
#' @return The exclusions dataframe
#' 
load_exclusions <- function(file){
 #	cat("debug: load_exclusions\n"); flush.console()
	output <- data.table::fread(file, sep = "\t")
	if (!is.integer(output[, 2]))
		stop("The second column in the exclusions file does not appear to contain only numeric data.\n")
	if (!is.integer(output[, 3]))
		stop("The third column in the exclusions file does not appear to contain only numeric data.\n")
	if (any(is.na(output)))
		stop("There appear to be NA's in the exclusions file.\n")
	if (any(output[, 2:3] < 0))
		stop("Some data in the exclusions file appears to be negative.\n")
	for (i in 1:nrow(output)) {
		if(output[i, 2] > output[i, 3])
			stop("The starting poing in line ", i, " is greater than the ending point in the exclusions file.\n")
	}
	colnames(output) <- c("chr", "start", "stop")
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
 #	cat("debug: load_regions\n"); flush.console()
	output <- data.table::fread(file, sep = "\t")
	if (!is.integer(output[, 2]))
		stop("The second column in the regions of interest file does not appear to contain only numeric data.\n")
	if (!is.integer(output[, 3]))
		stop("The third column in the regions of interest file does not appear to contain only numeric data.\n")
	if (any(is.na(output)))
		stop("There appear to be NA's in the regions of interest file.\n")
	if (any(output[, 2:3] < 0))
		stop("Some data in the regions of interest file appears to be negative.\n")
	for (i in 1:nrow(output)) {
		if(output[i, 2] > output[i, 3])
			stop("The starting poing in line ", i, " is greater than the ending point in the regions of interest file.\n")
	}	
	colnames(output) <- c("chr", "start", "stop")
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
 #	cat("debug: load_targets\n"); flush.console()
	output <- data.table::fread(file, sep = "\t")
	if (!is.integer(output[, 2]))
		stop("The second column in the targets file does not appear to contain only numeric data.\n")
	if (any(is.na(output)))
		stop("There appear to be NA's in the targets file.\n")
	if (any(output[, 2] < 0))
		stop("Some data in the targets file appears to be negative.\n")
	colnames(output) <- c("chr", "target")
	return(output)
}

#' Load file with the chromosomes' length
#' 
#' @inheritParams load_exclusions
#' 
#' @keywords internal
#' 
#' @return The lengths dataframe
#' 
load_lengths <- function(file){
	cat("debug: load_lengths\n"); flush.console()
	output <- read.table(file, sep = "\t")
	if (!is.integer(output[, 2]))
		stop("The second column in the chromosomes' length file does not appear to contain only numeric data.\n")
	if (any(is.na(output)))
		stop("There appear to be NA's in the chromosomes' length file.\n")
	if (any(output[, 2] < 0))
		stop("Some data in the chromosomes' length file appears to be negative.\n")
	if (any(output[, 2] == 0))
		stop("Some chromosomes' length is 0.\n")
	if (any(table(output[, 1]) > 1))
		stop("There are repeated chromosome names in the chromosomes' length file.\n")
	colnames(output) <- c("chr", "length")
	return(output)
}

#' Check chromosome name consistency
#' 
#' @param exclusions The exclusions dataframe
#' @param regions The regions dataframe
#' @param targets The targets dataframe
#' @param the.lengths The lengths dataframe
#' 
#' @keywords internal
#' 
#' @return The lengths dataframe
#' 
check_chr_names <- function(exclusions = NULL, regions = NULL, targets = NULL, the.lengths) {
 #	cat("debug: check_chr_names\n"); flush.console()
	if (!is.null(exclusions)) {
		if (any(is.na(match(exclusions[, 1], the.lengths[, 1]))))
			stop("Not all chromosomes' names in the exclusions match the names of the listed chromosomes.\n")
	}
	if (!is.null(regions)) {
		if (any(is.na(match(regions[, 1], the.lengths[, 1]))))
			stop("Not all chromosomes' names in the regions of interest match the names of the listed chromosomes.\n")
	}
	if (!is.null(targets)) {
		if (any(is.na(match(targets[, 1], the.lengths[, 1]))))
			stop("Not all chromosomes' names in the targets match the names of the listed chromosomes.\n")
	}
}

#' Check that start and end points fall whithin chromosome boundaries
#' 
#' @inheritParams check_chr_names
#' 
#' @keywords internal
#' 
#' @return The lengths dataframe
#' 
check_chr_boundaries <- function(exclusions = NULL, regions = NULL, targets = NULL, the.lengths) {
 #	cat("debug: check_chr_boundaries\n"); flush.console()
	if (!is.null(exclusions)) {
		for(i in unique(exclusions[, 1])) {
			link <- match(i, the.lengths[, 1])
			if (any(exclusions[exclusions[ ,1] == i, 2:3] > the.lengths[link, 2]))
				stop("Some data in the exclusions is off-boundaries.\n")
		}
	}
	if (!is.null(regions)) {
		for(i in unique(regions[, 1])) {
			link <- match(i, the.lengths[, 1])
			if (any(regions[regions[ ,1] == i, 2:3] > the.lengths[link, 2]))
				stop("Some data in the regions if interest is off-boundaries.\n")
		}
	}
	if (!is.null(targets)) {
		for(i in unique(targets[, 1])) {
			link <- match(i, the.lengths[, 1])
			if (any(targets[targets[ ,1] == i, 2] > the.lengths[link, 2]))
				stop("Some data in the targets is off-boundaries.\n")
		}
	}
}
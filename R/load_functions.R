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
	if (!is.null(exclusions)) {
		if (any(link <- is.na(match(exclusions$chr, the.lengths$name)))) {
			warning("Not all of the sequences' names in the exclusions match the names listed in the database. Removing orphan exclusions.", immediate. = TRUE, call. = FALSE)
			exclusions <- exclusions[!link, ]
		}
	}
	if (!is.null(regions)) {
		if (any(link <- is.na(match(regions$chr, the.lengths$name)))) {
			warning("Not all of the sequences' names in the regions of interest match the names listed in the database. Removing orphan regions.", immediate. = TRUE, call. = FALSE)
			regions <- regions[!link, ]
		}
	}
	if (!is.null(targets)) {
		if (any(link <- is.na(match(targets$chr, the.lengths$name)))) {
			warning("Not all of the sequences' names in the targets match the names listed in the database. Removing orphan targets.", immediate. = TRUE, call. = FALSE)
			targets <- targets[!link, ]
		}
	}
	return(list(exclusions = exclusions, regions = regions, targets = targets))
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
	if (!is.null(exclusions)) {
		capture <- lapply(unique(exclusions$chr), function(i) {
                        # cat(i)
			link <- match(i, the.lengths$name)
			if (any(exclusions[exclusions$chr == i, 2:3] > the.lengths$size[link]))
				stop("Exclusion data for sequence ", the.lengths$name[link], " is off-boundaries.\n", call = FALSE)
		})
	}
	if (!is.null(regions)) {
		capture <- lapply(unique(regions$chr), function(i) {
			link <- match(i, the.lengths$name)
			if (any(regions[regions[ ,1] == i, 2:3] > the.lengths$size[link]))
				stop("Region data for sequence ", the.lengths$name[link], " is off-boundaries.\n", call = FALSE)
		})
	}
	if (!is.null(targets)) {
		capture <- lapply(unique(targets$chr), function(i) {
			link <- match(i, the.lengths$name)
			if (any(targets[targets$chr == i, 2] > the.lengths$size[link]))
				stop("Target data for sequence ", the.lengths$name[link], " is off-boundaries.\n", call = FALSE)
		})
	}
}

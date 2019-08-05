#' Find valid ranges within the regions
#' 
#' @param length The length of the chromosome being analysed
#' @param n number of baits to generate
#' @param size The size of each bait
#' @param seed A number to fix the randomization process, for reproducibility
#' @param regions A subset of the regions dataframe for the target chromosome
#' @param exclusions A subset of the exclusions dataframe for the target chromosome.
#' @param chr The chromosome name.
#' 
#' @return a dataframe with starting and ending positions
#' 
#' @keywords internal
#' 
region_baits <- function(length, n, size, seed = NULL, regions, exclusions = NULL, chr) {
	temp_ranges <- regions
	if(temp_ranges$stop[nrow(temp_ranges)] > length)
		temp_ranges$stop[nrow(temp_ranges)] <- length
	if(!is.null(exclusions))
		temp_ranges <- trim_ranges(ranges = ranges, exclusions = exclusions)
	temp_ranges <- check_ranges(ranges = temp_ranges, size = size, chr = chr)
	valid_ranges <- expand_ranges(ranges = temp_ranges, size = size)
	if (length(valid_ranges) > n) {
		cat("Warning: The maximum possible number of individual targetted baits for chromosome", chr, "is lower than the desired n.\n")
		n <- length(valid.ranges)
	}
	return(get_bait_positions(valid_ranges = valid_ranges, seed = seed, size = size, n = n))
}

#' Find valid ranges within the targets
#' 
#' @inheritParams region_baits
#' @param targets A subset of the targets dataframe for the target chromosome
#' 
#' @return a dataframe with starting and ending positions
#' 
#' @keywords internal
#' 
target_baits <- function(length, n, size, seed = NULL, targets, exclusions = NULL, chr) {
	temp_ranges <- find_target_ranges(targets = targets, size = size, length = length)
	if(!is.null(exclusions))
		temp_ranges <- trim_ranges(ranges = ranges, exclusions = exclusions)
	temp_ranges <- check_ranges(ranges = temp_ranges, size = size, chr = chr)
	valid_ranges <- expand_ranges(ranges = temp_ranges, size = size)
	if (length(valid_ranges) > n) {
		cat("Warning: The maximum possible number of individual targetted baits for chromosome", chr, "is lower than the desired n.\n")
		n <- length(valid.ranges)
	}
	recipient <- get_bait_positions(valid_ranges = valid_ranges, seed = seed, size = size, n = n)
	return(list(table = as.data.frame(recipient), n = n))
}

#' Extract n number of baits randomly from parts of the chromosome length
#' 
#' @inheritParams region_baits
#' 
#' @return a dataframe with starting and ending positions
#' 
#' @keywords internal
#' 
random_baits <- function(length, n, size, seed = NULL, exclusions = NULL, chr) {
	if(!is.null(exclusions)) {
		# Check if the start is excluded
		if (exclusions[1, 2] == 1) {
			starting.point <- exclusions[1, 3] + 1
			exclusions <- exclusions[-1, ]
		} else {
			starting.point <- 1
		}
		# If there are more exclusions, check if length is excluded
		if (nrow(exclusions) > 0 && exclusions[nrow(exclusions), 3] == length) {
			length <- exclusions[nrow(exclusions), 2] - 1
			exclusions <- exclusions[-nrow(exclusions), ] 
		}
		# Find valid ranges
		if (nrow(exclusions) > 1) {
			temp_ranges <- find_random_ranges(starting.point = starting.point, length = length, exclusions = exclusions)
			temp_ranges <- check_ranges(ranges = temp_ranges, size = size, chr = chr)
			valid_ranges <- expand_ranges(ranges = temp_ranges, size = size)
		} else {
			valid_ranges <- seq(from = starting.point, to = length, by = 1)
		}
	} else {
		valid_ranges <- seq_len(length - size)
	}
	if (length(valid_ranges) > n) {
		cat("Warning: The maximum possible number of random baits for chromosome", chr, "is lower than the desired n.\n")
		n <- length(valid_ranges)
	}
	return(getBaitPositions(valid_ranges = valid_ranges, seed = seed, size = size, n = n))
}

#' Make a table of valid ranges within the chromosome
#' 
#' @param starting.point the first valid point in the chromosome
#' @param length the last valid point in the chromosome
#' @inheritParams region_baits
#' 
#' @return A table with valid ranges
#' 
#' @keywords internal
#' 
find_random_ranges <- function (starting.point, length, exclusions) {
	temp_ranges <- data.frame(
		start = c(starting.point, rep(NA, nrow(exclusions))),
		stop = c(rep(NA, nrow(exclusions)), length)
		)
	for(i in 1:nrow(exclusions)) {
		temp_ranges[i, 2] <- exclusions[i, 2] - 1
		temp_ranges[i + 1, 1] <- exclusions[i, 3] + 1
	}
	return(temp_ranges)
}

#' Make a table of valid ranges within the chromosome
#' 
#' @inheritParams region_baits
#' @inheritParams target_baits
#' 
#' @return A table with valid ranges
#' 
#' @keywords internal
#' 
find_target_ranges <- function(targets, size, length) {
	x <- data.frame(
		start = c(targets[, 2] - size),
		stop = c(targets[, 2])
		)
	if(x$stop[nrow(x)] > length)
		x$stop[nrow(x)] <- length
	return(x)
}

#' Extract the bait start and ending points from the valid ranges
#' 
#' @inheritParams region_baits
#' @param valid_ranges a sequence of valid starting points
#' 
#' @return a dataframe with starting and ending positions
#' 
#' @keywords internal
#' 
get_bait_positions <- function(valid_ranges, seed = NULL, size, n) {
	# get the bait positions
	if (!is.null(seed))
		set.seed(seed)
	recipient <- matrix(ncol = 2, nrow = n)
	recipient[, 1] <- sample(valid_ranges, size = n, replace = FALSE)
	recipient[, 2] <- recipient[, 1] + size
	colnames(recipient) <- c("Start", "Stop")
	return(recipient)
}

#' Check if any range starts/ends within an exclusion zone
#' 
#' @inheritParams region_baits
#' @param ranges a dataframe with starting and ending positions
#' 
#' @return a dataframe with starting and ending positions
#' 
#' @keywords internal
#' 
trim_ranges <- function(ranges, exclusions) {
	exclude.range <- NULL
	# for each range, check if ..
	for (i in 1:nrow(ranges)) {
		# for each exclusion area ...
		for (j in 1:nrow(exclusions)) {
			# the starting point of the range is within an exclusion
			if (ranges[i, 1] <= exclusions[j, 2] && ranges[i, 1] >= exclusions[j, 1]) {
				if (ranges[i, 2] > exclusions [j, 2]) {
					ranges[i, 1] <- exclusions[j, 2] + 1
				} else {
					exclude.range <- c(exclude.range, i)
				}
			}				
			# the ending point of the range is within an exclusion
			if (ranges[i, 2] >= exclusions[j, 1] && ranges[i, 2] <= exclusions[j, 2]) {
				if (ranges[i, 1] < exclusions [j, 1]) {
					ranges[i, 2] <- exclusions[j, 1] - 1
				} else {
					exclude.range <- c(exclude.range, i)
				}
			}
		}
	}
	if (!is.null(exclude.range)) {
		ranges <- temp.ranges[-exclude.range, ]
	}
	return(ranges)
}

#' Compare ranges' size with intended bait size
#' 
#' @inheritParams region_baits
#' @inheritParams trim_ranges
#' 
#' @return A table with valid ranges
#' 
#' @keywords internal
#' 
check_ranges <- function(ranges, size, chr) {
	ranges$range <- (ranges$stop - ranges$start) + 1
	if (any(to.exclude <- ranges$range < size + 1)) {
		# appendTo("Screen", ...)
		cat(paste0("TEMP: ", sum(to.exclude), " sub-ranges on chromosome ", chr," are too small to fit the bait size and will be excluded.\n")); flush.console()
		if (all(to.exclude))
			stop("All ranges in chromosome ", chr," are too small to fit the desired bait size. Aborting.")
		ranges <- ranges[!to.exclude, ]
	}
	return(ranges[, 1:2])
}

#' Transform the ranges dataframe into a sequence of numbers from which to sample
#' 
#' @inheritParams region_baits
#' @inheritParams trim_ranges
#' 
#' @return A vector of valid chromosome positions
#' 
#' @keywords internal
#' 
expand_ranges <- function(ranges, size) {
	output <- c()
	for (i in 1:nrow(ranges)) {
		output <- c(output, seq(from = ranges[i, 1], to = ranges[i, 2] - size, by = 1))
	}
	return(output)
}
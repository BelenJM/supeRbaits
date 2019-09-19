#' Find valid ranges within the regions
#' 
#' @param length The length of the chromosome being analysed
#' @param n number of baits to generate
#' @param size The size of each bait
#' @param tiling The minimum number of baits desired per range.
#' @param regions A subset of the regions dataframe for the target chromosome
#' @param exclusions A subset of the exclusions dataframe for the target chromosome.
#' @param chr The chromosome name, used in error messages.
#' 
#' @return a dataframe with starting and ending positions
#' 
#' @keywords internal
#' 
region_baits <- function(length, n, size, tiling = NULL, regions, exclusions = NULL, chr = NULL) {
	temp_ranges <- regions
	if (temp_ranges$stop[nrow(temp_ranges)] > length)
		temp_ranges$stop[nrow(temp_ranges)] <- length
	if (!is.null(exclusions))
		temp_ranges <- trim_ranges(ranges = ranges, exclusions = exclusions)
	valid_ranges <- check_ranges(ranges = temp_ranges, n = n, size = size, tiling = tiling, chr = chr)
	n.per.range <- check_n(ranges = valid_ranges, n = n, size = size, tiling = tiling)
	return(get_bait_positions(ranges = valid_ranges, size = size, n = n.per.range))
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
target_baits <- function(length, n, size, tiling = NULL, targets, exclusions = NULL, chr) {
	temp_ranges <- find_target_ranges(targets = targets, size = size, length = length)
	if(!is.null(exclusions))
		temp_ranges <- trim_ranges(ranges = ranges, exclusions = exclusions)
	temp_ranges <- check_ranges(ranges = temp_ranges, size = size, tiling = tiling, chr = chr)
	n.per.range <- check_n(ranges = valid_ranges, n = n, size = size, tiling = tiling)
	return(get_bait_positions(ranges = valid_ranges, size = size, n = n.per.range))
}

#' Extract n number of baits randomly from parts of the chromosome length
#' 
#' @inheritParams region_baits
#' 
#' @return a dataframe with starting and ending positions
#' 
#' @keywords internal
#' 
random_baits <- function(length, n, size, exclusions = NULL, chr) {
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
		# If there are more exclusions, break the main range appart
		if (nrow(exclusions) > 0)
			temp_ranges <- find_random_ranges(starting.point = starting.point, length = length, exclusions = exclusions)
		else
			temp_ranges <- data.frame(Start = starting.point, Stop = length)
		valid_ranges <- check_ranges(ranges = temp_ranges, n = n, size = size, tiling = 1, chr = chr)
		n.per.range <- check_n(ranges = valid_ranges, n = n, size = size, tiling = tiling)
	} else {
		valid_ranges <- data.frame(Start = 1, Stop = length)
	}
	if (length(valid_ranges) > n) {
		cat(paste0("Warning: The maximum possible number of random baits for chromosome ", chr, " is lower than the desired n.\n"))
		n <- length(valid_ranges)
	}
	return(get_bait_positions(ranges = valid_ranges, size = size, n = n))
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
check_ranges <- function(ranges, n, size, chr, tiling = 1) {
	# Range size
	ranges$range <- (ranges$stop - ranges$start) + 1
	ranges$max.baits <- ranges$range - size
	if (any(to.exclude <- ranges$max.baits < tiling)) {
		# appendTo("Screen", ...)
		cat(paste0("TEMP: ", sum(to.exclude), " sub-ranges on chromosome ", chr," are too small to fit the desired number of baits and will be excluded.\n")); flush.console()
		if (all(to.exclude))
			stop("All ranges in chromosome ", chr," are too small to fit the desired number of baits. Aborting.")
		ranges <- ranges[!to.exclude, ]
	}
	# More ranges than n/tiling
	if (n / tiling < nrow(ranges)) {
		cat(paste0("Warning: The desired n/tiling combination is not high enough to produce baits in all valid ranges in chromosome ", chr, ". Choosing a random subset of ranges.\n"))
		max.ranges <- roundDown(n / tiling, to = 1)
		ranges <- ranges[sample(1:nrow(ranges), size = max.ranges, replace = FALSE), ]
	}	
	return(ranges[, c(1:2, 4)])
}

check_n <- function(ranges, n, size, tiling = 1) {
	if (sum(ranges$max.baits) < n) {
		cat(paste0("Warning: The maximum possible number of individual targetted baits for chromosome", chr, "is lower than the desired n.\n"))
		n <- sum(ranges$max.baits)
	}
	n.per.range <- rep(roundDown(n / nrow(ranges), to = 1), nrow(ranges))
	# ensure you get the real n back by adding some here and there if needed
	while (sum(n.per.range) < n) {
		missing.n <- n - sum(n.per.range)
		expandable <- which(ranges$max.baits > n.per.range)
		max.extra.n <- min(ranges$max.baits[expandable] - n.per.range[expandable])
		# If the number of mising baits is smaller than the available ranges
		if (missing.n < length(expandable)) {
			add.here <- sample(expandable, size = missing.n, replace = FALSE)
			n.per.range[add.here] <- n.per.range[add.here] + 1
		} else {
			# If the number of maximum baits that can be allocated to the available ranges is bigger than needed
			if (length(expandable) * max.extra.n > missing.n) {
				# This will lead to a new iteration where the first IF will be triggered
				to.add <- roundDown(missing.n / length(expandable), to = 1)
				n.per.range[expandable] <- n.per.range[expandable] + to.add
			# If the number of maximum baits that can be allocated is SMALLER than the missing bait number
			} else {
				n.per.range[expandable] <- n.per.range[expandable] + max.extra.n
			}
		}
	}
	return(n.per.range)
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

#' Extract the bait start and ending points from the valid ranges
#' 
#' @inheritParams region_baits
#' @inheritParams trim_ranges
#' @param n a vector with the same length as the number of rows in ranges, containing the number of baits to extract from each range.
#' 
#' @return a dataframe with starting and ending positions
#' 
#' @keywords internal
#' 
get_bait_positions <- function(ranges, size, n) {
	# get the bait positions
	recipient <- data.frame(
		Start = integer(),
		Stop = integer()
		)
	for (i in 1:nrow(ranges)) {
		aux.start <- sample(ranges[i, 1]:(ranges[i, 2] - size), size = n[i], replace = FALSE)
		aux.df <- data.frame(
			Start = aux.start,
			Stop = aux.start + size
			)
		recipient <- rbind(recipient, aux.df)
	}
	return(recipient)
}


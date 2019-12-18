#' Find valid ranges within the regions
#' 
#' @param chr.length The length of the chromosome being analysed
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
region_baits <- function(chr.length, n, size, tiling = NULL, regions, exclusions = NULL, chr, used.baits = NULL, verbose) {
 #	cat("debug: region_baits\n"); flush.console()
	temp.ranges <- regions[, -1, drop = FALSE]
	if (temp.ranges$stop[nrow(temp.ranges)] > chr.length)
		temp.ranges$stop[nrow(temp.ranges)] <- chr.length
	if (!is.null(exclusions))
		temp.ranges <- trim_ranges(ranges = temp.ranges, exclusions = exclusions)
	valid.ranges <- check_ranges(ranges = temp.ranges, n = n, size = size, tiling = tiling, chr = chr, used.baits = used.baits, verbose = verbose)
	n.per.range <- check_n(ranges = valid.ranges, n = n, tiling = tiling, chr = chr, type = "region", verbose = verbose)
	return(get_bait_positions(ranges = valid.ranges, size = size, n = n.per.range, used.baits = used.baits))
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
target_baits <- function(chr.length, n, size, tiling = NULL, targets, exclusions = NULL, chr, used.baits = NULL, verbose) {
 #	cat("debug: target_baits\n"); flush.console()
	temp.ranges <- find_target_ranges(targets = targets, size = size, chr.length = chr.length)
	if(!is.null(exclusions))
		temp.ranges <- trim_ranges(ranges = temp.ranges, exclusions = exclusions)
	valid.ranges <- check_ranges(ranges = temp.ranges, n = n, size = size, tiling = tiling, chr = chr, used.baits = used.baits, verbose = verbose)
	n.per.range <- check_n(ranges = valid.ranges, n = n, tiling = tiling, chr = chr, type = "target", verbose = verbose)
	return(get_bait_positions(ranges = valid.ranges, size = size, n = n.per.range, used.baits = used.baits))
}

#' Extract n number of baits randomly from parts of the chromosome length
#' 
#' @inheritParams region_baits
#' 
#' @return a dataframe with starting and ending positions
#' 
#' @keywords internal
#' 
random_baits <- function(chr.length, n, size, exclusions = NULL, chr, used.baits = NULL, verbose) {
 #	cat("debug: random_baits\n"); flush.console()
	if(!is.null(exclusions)) {
		# Check if the start is excluded
		if (exclusions$start[1] == 1) {
			starting.point <- exclusions$stop[1] + 1
			exclusions <- exclusions[-1, ]
		} else {
			starting.point <- 1
		}
		# If there are more exclusions, check if length is excluded
		if (nrow(exclusions) > 0 && exclusions$stop[nrow(exclusions)] == chr.length) {
			chr.length <- exclusions$start[nrow(exclusions)] - 1
			exclusions <- exclusions[-nrow(exclusions), ] 
		}
		# If there are more exclusions, break the main range appart
		if (nrow(exclusions) > 0)
			temp.ranges <- find_random_ranges(starting.point = starting.point, chr.length = chr.length, exclusions = exclusions)
		else
			temp.ranges <- data.frame(Start = starting.point, Stop = chr.length)
	} else {
		temp.ranges <- data.frame(start = 1, stop = chr.length)
	}
	valid.ranges <- check_ranges(ranges = temp.ranges, n = n, size = size, tiling = 1, chr = chr, used.baits = used.baits, verbose = verbose)
	n.per.range <- check_n(ranges = valid.ranges, n = n, tiling = 1, chr = chr, type = "random", verbose = verbose)

	return(get_bait_positions(ranges = valid.ranges, size = size, n = n.per.range, used.baits = used.baits))
}

#' Make a table of valid ranges within the chromosome
#' 
#' @param starting.point the first valid point in the chromosome
#' @param chr.length the last valid point in the chromosome
#' @inheritParams region_baits
#' 
#' @return A table with valid ranges
#' 
#' @keywords internal
#' 
find_random_ranges <- function (starting.point, chr.length, exclusions) {
 #	cat("debug: find_random_ranges\n"); flush.console()
	temp.ranges <- data.frame(
		start = c(starting.point, rep(NA, nrow(exclusions))),
		stop = c(rep(NA, nrow(exclusions)), chr.length)
		)
	for(i in 1:nrow(exclusions)) {
		temp.ranges$stop[i] <- exclusions$start[i] - 1
		temp.ranges$start[i + 1] <- exclusions$stop[i] + 1
	}
	return(temp.ranges)
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
find_target_ranges <- function(targets, size, chr.length) {
 #	cat("debug: find_target_ranges\n"); flush.console()
	x <- data.frame(
		start = c(targets$target - size),
		stop = c(targets$target + size)
		)
	if(x$stop[nrow(x)] > chr.length)
		x$stop[nrow(x)] <- chr.length
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
 #	cat("debug: trim_ranges\n"); flush.console()
	exclude.range <- NULL
	# for each range, check if ..
	for (i in 1:nrow(ranges)) {
		# for each exclusion area ...
		for (j in 1:nrow(exclusions)) {
			# the starting point of the range is within an exclusion
			if (ranges$start[i] <= exclusions$stop[j] && ranges$start[i] >= exclusions$start[j]) {
				if (ranges$stop[i] > exclusions$stop[j]) {
					ranges$start[i] <- exclusions$stop[j] + 1
				} else {
					exclude.range <- c(exclude.range, i)
				}
			}				
			# the ending point of the range is within an exclusion
			if (ranges$stop[i] >= exclusions$start[j] && ranges$stop[i] <= exclusions$stop[j]) {
				if (ranges$start[i] < exclusions$start[j]) {
					ranges$stop[i] <- exclusions$start[j] - 1
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
check_ranges <- function(ranges, n, size, chr, tiling = 1, used.baits = NULL, verbose) {
 #	cat("debug: check_ranges\n"); flush.console()
	# Range size
	ranges$range <- (ranges$stop - ranges$start) + 1 # +1 because the first bp also counts
	ranges$used.baits <- apply(ranges, 1, function(x) sum(used.baits >= x[1] & used.baits <= x[2]))
	ranges$max.baits <- ranges$range - ranges$used.baits - size
	if (any(to.exclude <- ranges$max.baits < tiling)) {
		if (verbose)
			warning(paste0(sum(to.exclude), " sub-ranges on chromosome ", chr," are too small to fit the desired number of baits and will be excluded."), call. = FALSE, immediate. = TRUE)
		if (all(to.exclude))
			stop("All ranges in chromosome ", chr," are too small to fit the desired number of baits. Aborting.")
		ranges <- ranges[!to.exclude, ]
	}
	# More ranges than n/tiling
	if (n / tiling < nrow(ranges)) {
		if (verbose)
			warning(paste0("The desired n/tiling combination is not high enough to produce baits in all valid ranges in chromosome ", chr, ". Choosing a random subset of ranges."), call. = FALSE, immediate. = TRUE)
		max.ranges <- roundDown(n / tiling, to = 1)
		ranges <- ranges[sample(1:nrow(ranges), size = max.ranges, replace = FALSE), ]
	}	
	return(ranges[, c("start", "stop", "max.baits")])
}

check_n <- function(ranges, n, tiling = 1, chr, type = c("random", "target", "region"), verbose) {
 #	cat("debug: check_n\n"); flush.console()
	type <- match.arg(type)
	if (sum(ranges$max.baits) < n) {
		if (verbose)
			warning(paste0("The maximum possible number of unique ", type, " baits (", sum(ranges$max.baits), ") for chromosome ", chr, " is lower than the desired n (", n, ")."), call. = FALSE, immediate. = TRUE)
		n <- sum(ranges$max.baits)
	}
	# If the number of requested baits matches the available baits (most likely because the if above was triggered, finish here)
	if (sum(ranges$max.baits) == n)
		return(ranges$max.baits)
	
	# If there are too many bait spots available, distribute equally.
	if (sum(ranges$max.baits) > n) {
		# start with the minimum tiling and work from there
		n.per.range <- rep(tiling, nrow(ranges))
		# start incrementing until n is fulfilled
		while (sum(n.per.range) < n) {
			missing.n <- n - sum(n.per.range)
			expandable <- which(ranges$max.baits > n.per.range)
			max.extra.n <- min(ranges$max.baits[expandable] - n.per.range[expandable])
			# If the number of mising baits is smaller than the available ranges
			if (missing.n < length(expandable)) {
				add.here <- sample(expandable, size = missing.n, replace = FALSE)
				n.per.range[add.here] <- n.per.range[add.here] + 1
			} else {
				# If the number of maximum baits that can be allocated to the available ranges is BIGGER than the missing bait number
				if (length(expandable) * max.extra.n > missing.n) {
					# This will lead to a new iteration where the first if will be triggered
					to.add <- roundDown(missing.n / length(expandable), to = 1)
					n.per.range[expandable] <- n.per.range[expandable] + to.add
				# If the number of maximum baits that can be allocated is SMALLER OR EQUAL to the missing bait number
				} else {
					n.per.range[expandable] <- n.per.range[expandable] + max.extra.n
				}
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
 #	cat("debug: expand_ranges\n"); flush.console()
	output <- c()
	for (i in 1:nrow(ranges)) {
		output <- c(output, seq(from = ranges$start[i], to = ranges$stop[i] - size, by = 1))
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
get_bait_positions <- function(ranges, size, n, used.baits = NULL) {
 #	cat("debug: get_bait_positions\n"); flush.console()
	# get the bait positions
	recipient <- data.frame(
		Start = integer(),
		Stop = integer()
		)
	for (i in 1:nrow(ranges)) {
		# find and exclude already used baits
		the.range <- ranges$start[i]:(ranges$stop[i] - size)
		already.used <- used.baits[used.baits %in% the.range]
		if (length(already.used) == 0) {
			new.starts <- the.range
		} else {
			already.used.position <- match(already.used, the.range)
			new.starts <- the.range[-already.used.position]
		}
		# sample new baits
		aux.start <- sample(new.starts, size = n[i], replace = FALSE)
		aux.df <- data.frame(
			Start = aux.start,
			Stop = aux.start + size
			)
		recipient <- rbind(recipient, aux.df)
	}
	return(recipient)
}


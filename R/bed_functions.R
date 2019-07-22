#' Generates random baits
#'
#' Function inspired by this comment in StackOverflow: https://stackoverflow.com/questions/49149839/simulate-random-positions-from-a-list-of-intervals
#'
#' @param n number of baits to generate
#' @param size The size of each bait
#' @param lengths A file containing the chromosomes' length
#' @param exclusions A file containing regions to exclude
#' @param regions A file containing regions to target
#' @param targets A file containing points to target
#' @param seed A number to fix the randomization process, for reproducibility
#' 
#' @return A dataframe of baits
#' 
#' @export
#' 
main_function <- function(n, size, lengths, exclusions = NULL, regions = NULL, targets = NULL, seed = NULL){
	# Initial checks
	if(!is.numeric(n))
		stop("'n' must be numeric.\n")
	if(length(n) != 1)
		stop("Length of 'n' must be 1.\n")
	if(!is.numeric(size))
		stop("'size' must be numeric.\n")
	if(length(size) != 1)
		stop("Length of 'size' must be 1.\n")

	# Reduce size to include the first bp
	size <- size - 1

	# Import data
	lengths <- load_lengths(file = lengths)
	if(!is.null(exclusions))
		exclusions <- load_exclusions(file = exclusions)
	if(!is.null(regions))
		regions <- load_regions(file = regions)
	if(!is.null(targets))
		targets <- load_targets(file = targets)

	# Compatibility checks
	check_chr_names(exclusions = exclusions, regions = regions, targets = targets, lengths = lengths)
	check_chr_boundaries(exclusions = exclusions, regions = regions, targets = targets, lengths = lengths)
	if (any(size > lengths[, 2]))
		stop("'size' is larger than at least one of the chromosome lengths.\n")

	recipient <- list()
	for (i in 1:nrow(lengths)) {
		cat(paste0("debug: examining chromosome ", i,".\n")); flush.console()
		if (!is.null(exclusions)) 
			trimmed_exclusions <- subsample(input = exclusions, link = lengths[i, 1])
		else 
			trimmed_exclusions <- NULL
		if (!is.null(regions)) 
			trimmed_regions <- subsample(input = regions, link = lengths[i, 1])
		else
			trimmed_regions <- NULL
		if (!is.null(targets)) 
			trimmed_targets <- subsample(input = targets, link = lengths[i, 1])
		else
			trimmed_targets <- NULL

		if(is.null(c(trimmed_exclusions, trimmed_regions, trimmed_targets)))
			recipient[[i]] <- all_random(length = lengths[i, 2], n = n, size = size, seed = seed)

		if(!is.null(c(trimmed_regions, trimmed_targets)) & is.null(trimmed_exclusions)) {
			recipient[[i]] <- all_targetted(length = lengths[i, 2], n = n, size = size, seed = seed, 
				regions = trimmed_regions, targets = trimmed_targets)
		}

		if(is.null(c(trimmed_regions, trimmed_targets)) & !is.null(trimmed_exclusions)) {
			recipient[[i]] <- trimmed_random(length = lengths[i, 2], n = n, size = size, 
				seed = seed, exclusions = trimmed_exclusions, chr = lengths[i, 1])
		}

		if(!is.null(c(trimmed_regions, trimmed_targets)) & !is.null(trimmed_exclusions)) {
			recipient[[i]] <- trimmed_targetted(length = lengths[i, 2], n = n, size = size, seed = seed,
			  regions = trimmed_regions, targets = trimmed_targets, exclusions = trimmed_exclusions)
		}
	}
	names(recipient) <- as.character(lengths[, 1])

	if (check_python())
		processing_baitfile(input_file1 = "chrom_salmon_chunk.fasta.txt", input_file2 = "bed.CM003279.txt")
	
	# Randomly sample bed file rows, proportional to the length of each range
	# simulated.sites <- bed2[sample(.N, size = n, replace = TRUE, prob = bed2$size)]

	# Randomly sample uniformly within each chosen range
	# simulated.sites[, position := sample(start:end, size = 1), by = 1:dim(simulated.sites)[1]]

	# Remove extra columns and format as needed
	# simulated.sites[, start  := position]
	# simulated.sites[, end := position]
	# simulated.sites[, c("size", "position") := NULL]
	# out_table <- simulated.sites
	# write.table(out_table, name_file)
	# return(out_table)
	return(recipient)
}

#' Extract the values of input that match the linking element
#' 
#' @param input A dataframe whose first column will be matched to link
#' @param link a keyword to be searched for
#' 
#' @return A subset of input
#' 
#' @keywords internal
#' 
subsample <- function(input, link) {
	chr <- link
	link <- grepl(link, input[, 1])
	if (any(link)) {
		output <- input[link, ]
		output <- output[order(output[, 2]), ]
	}	else {
		output <- NULL
	}
	if (!is.null(output)) {
		if (any(table(output[, 2]) > 1))
			stop("There are duplicated starting points/targets for one of the trimming elements of chromosome ", chr, ".")
		if (ncol(output) == 3) {
			if(any(table(output[, 3]) > 1))
				stop("There are duplicated ending points for one of the trimming elements for chromosome ", chr, ".")
			if (nrow(output) > 1) {
				trigger <- NULL
				for (i in 2:nrow(output))
					trigger[i - 1] <- output[i, 2] <= output[i - 1, 3]
				if (any(trigger))
					stop("The regions to include or exclude overlap for chromosome ", chr,".")
				trigger <- NULL
				for (i in 2:nrow(output))
					trigger[i - 1] <- output[i, 2] == output[i - 1, 3] + 1
				if (any(trigger))
					stop("There are contiguous regions to include or exclude for chromosome ", chr,". Please list these as a single region.")
			}
		}
	}
	return(output)
}

#' Extract n number of baits randomly from the whole chromosome length
#' 
#' @param length The length of the chromosome being analysed
#' @param n number of baits to generate
#' @param size The size of each bait
#' @param seed A number to fix the randomization process, for reproducibility
#' 
#' @return a dataframe with starting and ending positions
#' 
#' @keywords internal
#' 
all_random <- function(length, n, size, seed = NULL) {
	if (!is.null(seed))
		set.seed(seed)
	recipient <- matrix(ncol = 2, nrow = n)
	recipient[, 1] <- sample(seq_len(length - size), size = n, replace = TRUE)
	recipient[, 2] <- recipient[, 1] + size
	colnames(recipient) <- c("Start", "Stop")
	return(as.data.frame(recipient))
}

#' Extract n number of baits from the whole chromosome length, targeting specific areas
#' 
#' @inheritParams all_random
#' @param regions A subset of the regions dataframe for the target chromosome
#' @param targets A subset of the targets dataframe for the target chromosome
#' 
#' @return a dataframe with starting and ending positions
#' 
#' @keywords internal
#' 
all_targetted <- function(length, n, size, seed = NULL, regions = NULL, targets = NULL) {
	stop("Triggered all_targetted\n")
	# MISSING: check if size fits within regions of interest, if there are any?
}

#' Extract n number of baits randomly from parts of the chromosome length
#' 
#' @inheritParams all_random
#' @param exclusions A subset of the exclusions dataframe for the target chromosome.
#' @param chr The chromosome name, to display in messages and warnings, if needed.
#' 
#' @return a dataframe with starting and ending positions
#' 
#' @keywords internal
#' 
trimmed_random <- function(length, n, size, seed = NULL, exclusions, chr) {
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
		temp_ranges <- find_ranges(starting.point = starting.point, length = length, exclusions = exclusions)
		temp_ranges <- check_ranges(ranges = temp_ranges, size = size, chr = chr)
		valid_ranges <- expand_ranges(ranges = temp_ranges, size = size)
	} else {
		valid_ranges <- seq(from = starting.point, to = length, by = 1)
	}
	# get the bait positions
	if (!is.null(seed))
		set.seed(seed)
	recipient <- matrix(ncol = 2, nrow = n)
	recipient[, 1] <- sample(valid_ranges, size = n, replace = TRUE)
	recipient[, 2] <- recipient[, 1] + size
	colnames(recipient) <- c("Start", "Stop")
	return(as.data.frame(recipient))
}

#' Extract n number of baits from parts of the chromosome length, targeting specific areas
#' 
#' @inheritParams all_random
#' @inheritParams all_targetted
#' @inheritParams trimmed_random
#' 
#' @return a dataframe with starting and ending positions
#' 
#' @keywords internal
#' 
trimmed_targetted <- function(length, n, size, seed = NULL, regions = NULL, targets = NULL, exclusions) {
	stop("Triggered trimmed_targetted\n")
}

#' Make a table of valid ranges within the chromosome
#' 
#' @param starting.point the first valid point in the chromosome
#' @param length the last valid point in the chromosome
#' @inheritParams trimmed_random
#' 
#' @return A table with valid ranges
#' 
#' @keywords internal
#' 
find_ranges <- function (starting.point, length, exclusions) {
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

#' Compare ranges' size with intended bait size
#' 
#' @inheritParams all_random
#' @inheritParams trimmed_random
#' @param ranges A table of valid ranges
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
#' @inheritParams all_random
#' @inheritParams check_ranges
#' 
#' @return A vector of chromosome positions
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
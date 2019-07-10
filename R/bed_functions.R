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
#' @output A dataframe of baits
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
	# MISSING: check if n fits within all chromosome limits
	# These two may make more sense on a chromosome to chromosome level:
	# MISSING: check if n fits within all non-excluded regions, if there are any
	# MISSING: check if n fits within targets, if there are any?

	recipient <- list()
	for (i in 1:ncol(lengths)) {
		if (!is.null(exclusions)) 
			trimmed_exclusions <- subsample(input = exclusions, link = lengths[i, 1])
		if (!is.null(regions)) 
			trimmed_regions <- subsample(input = regions, link = lengths[i, 1])
		if (!is.null(targets)) 
			trimmed_targets <- subsample(input = targets, link = lengths[i, 1])

		if(is.null(c(trimmed_exclusions, trimmed_regions, trimmed_targets)))
			recipient[[i]] <- all_random(length = lengths[i, ], n = n, size = size, seed = seed)

		if(!is.null(c(trimmed_regions, trimmed_targets)) & is.null(trimmed.exclusions)) {
			recipient[[i]] <- all_targetted(length = lengths[i, ], n = n, size = size, seed = seed, 
				regions = trimmed_regions, targets = trimmed_targets)
		}

		if(is.null(c(trimmed_regions, trimmed_targets)) & !is.null(trimmed_exclusions)) {
			recipient[[i]] <- trimmed_random(length = lengths[i, ], n = n, size = size, 
				seed = seed, exclusions = trimmed_exclusions)
		}

		if(!is.null(c(trimmed_regions, trimmed_targets)) & !is.null(trimmed.exclusions)) {
			recipient[[i]] <- trimmed_targetted(length = lengths[i, ], n = n, size = size, seed = seed,
			  regions = trimmed_regions, targets = trimmed_targets, exclusions = trimmed_exclusions)
		}

		names(recipient)[length(recipient)] <- lengths[i, 1]
	}

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
	link <- grepl(link, input[, 1])
	if (any(link))
		output <- input[link, ]
	else
		output <- NULL
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
}

#' Extract n number of baits randomly from parts of the chromosome length
#' 
#' @inheritParams all_random
#' @param exclusions A subset of the exclusions dataframe for the target chromosome
#' 
#' @return a dataframe with starting and ending positions
#' 
#' @keywords internal
#' 
trimmed_random <- function(length, n, size, seed = NULL, exclusions) {
	stop("Triggered trimmed_random\n")	
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
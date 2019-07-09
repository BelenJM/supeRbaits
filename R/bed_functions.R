#' Generates random baits
#'
#' Function inspired by this comment in StackOverflow: https://stackoverflow.com/questions/49149839/simulate-random-positions-from-a-list-of-intervals
#'
#' @param n number of baits to generate
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
main_function <- function(n, lengths, exclusions = NULL, regions = NULL, targets = NULL, seed = NULL){
	# Initial checks
	if(!is.numeric(n))
		stop("'n' must be numeric.\n")
	if(length(n) != 1)
		stop("Length of 'n' must be 1.\n")

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

		if(is.null(c(trimmed_exclusions, trimmed_regions, trimmed_targets))))
			recipient[[i]] <- all_random(lengths = lengths[i, ], n = n, seed = seed)

		if(!is.null(c(trimmed_regions, trimmed_targets)) & is.null(trimmed.exclusions))
			recipient[[i]] <- all_targetted(lengths = lengths[i, ], n = n, seed = seed, regions = trimmed_regions, targets = trimmed_targets)

		if(is.null(c(trimmed_regions, trimmed_targets)) & !is.null(trimmed_exclusions))
			recipient[[i]] <- trimmed_random(lengths = lengths[i, ], n = n, seed = seed, exclusions = trimmed_exclusions)

		if(!is.null(c(trimmed_regions, trimmed_targets)) & !is.null(trimmed.exclusions))
			recipient[[i]] <- trimmed_targetted(lengths = lengths[i, ], n = n, seed = seed, regions = trimmed_regions, targets = trimmed_targets, exclusions = trimmed_exclusions)

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

# Extract n number of baits from the whole chromosome length
all_random <- function(lengths, n, size, seed = NULL) {
	
}

all_targetted <- function(lengths, n, , size, seed = NULL, regions = NULL, targets = NULL)) {
	
}

trimmed_random <- function(lengths, n, size, seed = NULL, exclusions) {
	
}

trimmed_targetted <- function(lengths, n, size, seed = NULL, regions = NULL, targets = NULL, exclusions) {

}
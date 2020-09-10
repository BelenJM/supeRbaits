#' Check if restrict is compatible with the number of sequences present
#' 
#' @inheritParams main_function
#' @param sequences A vector with the names of the sequences in the current data.
#' 
#' @keywords internal
#' 
#' @return an updated restrict argument (numeric)
#' 
check_restrict <- function(restrict, sequences) {
	if (is.numeric(restrict)) {
		if (all(restrict > length(sequences)))
			stop("'restrict' is set to numeric values, but ALL listed values are larger than the number of available sequences.\nNumber of available sequences: ", length(sequences), "\nMin. index requested: ", min(restrict), "\nMax. index requested: ", max(restrict), "\n", call. = FALSE)
		if (max(restrict) > length(sequences)) {
			warning("'restrict' is set to numeric values, but some listed values are larger than the number of available sequences.\nDiscarded index values: ", paste(restrict[restrict > length(sequences)], collapse = " "), immediate. = TRUE, call. = FALSE)
			restrict <- restrict[restrict <= length(sequences)]
		}
	}

	if (is.character(restrict)) {
		link <- match(restrict, sequences)
		if (all(is.na(link)))
			stop("None of the sequence names listed in 'restrict' matches the available sequences.\n", call. = FALSE)
		if (any(is.na(link))) {
			warning("Some sequences listed in 'restrict' do not match the available sequences.\nMissing sequences: '", paste(restrict[is.na(link)], collapse = "', '"), "'", immediate. = TRUE, call. = FALSE)
			link <- link[!is.na(link)]
		}
		restrict <- link
	}

	return(restrict)
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

#' Check that start and end points fall within chromosome boundaries
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
				stop("Exclusion data for sequence ", the.lengths$name[link], " is off-boundaries.\n", call. = FALSE)
		})
	}
	if (!is.null(regions)) {
		capture <- lapply(unique(regions$chr), function(i) {
			link <- match(i, the.lengths$name)
			if (any(regions[regions[ ,1] == i, 2:3] > the.lengths$size[link]))
				stop("Region data for sequence ", the.lengths$name[link], " is off-boundaries.\n", call. = FALSE)
		})
	}
	if (!is.null(targets)) {
		capture <- lapply(unique(targets$chr), function(i) {
			link <- match(i, the.lengths$name)
			if (any(targets[targets$chr == i, 2] > the.lengths$size[link]))
				stop("Target data for sequence ", the.lengths$name[link], " is off-boundaries.\n", call. = FALSE)
		})
	}
}

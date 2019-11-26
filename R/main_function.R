#' Generates random baits
#'
#' Function inspired by this comment in StackOverflow: https://stackoverflow.com/questions/49149839/simulate-random-positions-from-a-list-of-intervals
#'
#' @param n number of baits to generate
#' @param size The size of each bait
#' @param database A database of chromosomes
#' @param exclusions A file containing regions to exclude
#' @param regions A file containing regions of interest.
#' @param regions.prop The proportion of baits that should overlap the regions of interest.
#' @param regions.tiling The minimum number of baits to distribute per region of interest.
#' @param targets A file containing bp's to target.
#' @param targets.prop The proportion of baits that should overlap the targets.
#' @param targets.tiling The minimum number of baits to distribute per target.
#' @param seed A number to fix the randomization process, for reproducibility
#' @param restrict A vector of chromosome names to which the analysis should be restricted.
#' 
#' @return A dataframe of baits
#' 
#' @export
#' 
main_function <- function(n, size, database, exclusions = NULL, 
	regions = NULL, regions.prop = NULL, regions.tiling = NULL,
	targets = NULL, targets.prop = NULL, targets.tiling = NULL,
	seed = NULL, restrict = NULL, debug = FALSE, gc = c(0.3, 0.5)){

  
	if (debug) {
    message("!!!--- Debug mode has been activated ---!!!")
		on.exit(save(list = ls(), file = "supeRbaits_debug.RData"), add = TRUE)
	} else {
		on.exit(unlink("temp_folder_for_supeRbaits", recursive = TRUE), add = TRUE)
	}

	if(!is.null(seed))
		set.seed(seed)

	# Initial checks
	# check_python()
	
	if (!is.numeric(n))
		stop("'n' must be numeric.\n")
	if (length(n) != 1)
		stop("Length of 'n' must be 1.\n")
	if (!is.numeric(size))
		stop("'size' must be numeric.\n")
	if (length(size) != 1)
		stop("Length of 'size' must be 1.\n")

	if (!is.null(regions.prop) && (!is.numeric(regions.prop) | (regions.prop > 1 | regions.prop < 0)))
		stop("'regions.prop' must be numeric and between 0 and 1.\n")
	if (!is.null(targets.prop) && (!is.numeric(targets.prop) | (targets.prop > 1 | targets.prop < 0)))
		stop("'targets.prop' must be numeric and between 0 and 1.\n")
	if (!is.null(regions.tiling) && !is.integer(regions.tiling))
		stop("'regions.tiling' must be an integer.\n")
	if (!is.null(targets.tiling) && !is.integer(targets.tiling))
		stop("'targets.tiling' must be an integer.\n")

	if (is.null(regions.prop) & !is.null(regions))
		stop("Please include the desired percentage of regional baits in 'regions.prop'.\n")
	if (!is.null(regions.tiling) & is.null(regions))
		stop("'regions.tiling' is set but no regions were included.\n")
	if (is.null(targets.prop) & !is.null(targets))
		stop("Please include the desired percentage of targetted baits in 'targets.prop'.\n")
	if (!is.null(targets.tiling) & is.null(targets))
		stop("'targets.tiling' is set but no targets were included.\n")

	if (sum(regions.prop, targets.prop) > 1)
		stop("The sum of 'regions.prop' and 'targets.prop' must not be greated than one.\n")
	
	if (length(gc) != 2)
		stop("Please provide two values in 'gc' (minimum and maximum percentage).\n")
	if (!is.numeric(gc))
		stop("'gc' must be numeric and contain two values between 0 and 1.\n")
	if (any(gc > 1) | any(gc < 0))
		stop("'gc' ranges must be between 0 and 1.\n")
	if (gc[1] > gc[2])
		stop("The first value of 'gc' must be smaller or equal to the second value.\n")

	if (is.null(regions.tiling))
		regions.tiling <- 1
	if (is.null(targets.tiling))
		targets.tiling <- 1

	# Reduce size to include the first bp
	size <- size - 1

	# create temp folder to dump stuff in
	if (!dir.exists("temp_folder_for_supeRbaits"))
		dir.create("temp_folder_for_supeRbaits")
		
	# extract lengths with python
	#setwd("temp_folder_for_supeRbaits")
	get_lengths(database = database, restrict = restrict)
	#setwd("..")

	# Import data
	lengths <- load_lengths(file = "temp_folder_for_supeRbaits/genome_size.txt")
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

	bait.points <- list()
	for (i in 1:nrow(lengths)) {
		if (debug)
			cat(paste0("debug: examining chromosome ", lengths[i, 1], ".\n"))
		# extract relevant parameters
		params <- trim_parameters(chr = lengths[i, 1], exclusions = exclusions, regions = regions, targets = targets)
		# region baits
		if(!is.null(regions.prop) && regions.prop > 0) {
			n.regions = n * regions.prop
			if (!is.null(params$regions)) {
				temp.regions <- region_baits(chr.length = lengths[i, 2], n = n.regions, size = size, tiling = regions.tiling,
				regions = params$regions, exclusions = params$exclusions, chr = lengths[i, 1])
				n.regions = nrow(temp.regions)
			} else {
				message("Message: No regions found for chromosome ", lengths[i, 1], ".")
				temp.regions <- NULL
				n.regions = 0
			}
		} else {
			temp.regions <- NULL
			n.regions = 0
		}
		# targetted baits
		if(!is.null(targets.prop) && targets.prop > 0) {
			n.targets = n * targets.prop
			if (!is.null(params$targets)) {
				temp.targets <- target_baits(chr.length = lengths[i, 2], n = n.regions, size = size, tiling = targets.tiling,
				targets = params$targets, exclusions = params$exclusions, chr = lengths[i, 1])
				n.targets = nrow(temp.targets)
			} else {
				message("Message: No targets found for chromosome ", lengths[i, 1], ".")
				temp.targets <- NULL
				n.targets = 0
			}
		} else {
			temp.targets <- NULL
			n.targets = 0
		}
		# random baits
		n.random <- n - (n.regions + n.targets)
		temp.random <- random_baits(chr.length = lengths[i, 2], n = n, size = size, 
				exclusions = params$exclusions, chr = lengths[i, 1])
		if (n.random > 0)
		else
			temp.random <- NULL
		# bring together the different parts
		bait.points[[i]] <- rbind(temp.regions, temp.targets, temp.random) # not sure if this works with nulls
	}
	names(bait.points) <- as.character(lengths[, 1])

	# fetch the baits' sequences (Python stuff)
	
	baits <- lapply(seq_along(bait.points), function(i) {
		write.table(bait.points[i], file = paste0("temp_folder_for_supeRbaits/", names(bait.points)[i], ".txt"), row.names = FALSE)
		retrieve_baits(chr = names(bait.points)[i], database = database)
		output <- data.table::fread(paste0("temp_folder_for_supeRbaits/", names(bait.points)[i], "_py.txt"))
		output$pGC <- output$Number_GC / (size + 1)
		return(output)
	})
	names(baits) <- names(bait.points)

	good.baits <- lapply(seq_along(baits), function(i) {
		link <- baits[[i]]$pGC > gc[1] & baits[[i]]$pGC < gc[2]
		if (all(!link)) {
			message(paste0("Message: No baits passed the GC percentage test for chromosome ", names(baits)[i], "."))
			return(NULL)
		}
		if (any(!link))
			message(paste0("Message: ", sum(!link), " baits were excluded from chr ", names(baits)[i], " due to their GC percentage."))
		return(baits[[i]][link, ])
	})
	names(good.baits) <- names(baits)

	return(list(baits = baits, good.baits = good.baits))
}

#' extract exclusions, regions, and targets relevant for the chromosome being analysed
#' 
#' @inheritParams main_function
#' @param chr The chromosome name
#' 
#' @return a list of trimmed data frames
#' 
#' @keywords internal
#' 
trim_parameters <- function(chr, exclusions = NULL, regions = NULL, targets = NULL) {
	cat("debug: trim_parameters\n"); flush.console()
	output <- list(exclusions  = NULL, regions = NULL, targets = NULL)
	if (!is.null(exclusions)) 
		output[[1]] <- subsample(input = exclusions, link = chr)
	if (!is.null(regions)) 
		output[[2]] <- subsample(input = regions, link = chr)
	if (!is.null(targets)) 
		output[[3]] <- subsample(input = targets, link = chr)
	# names(output) <- c("exclusions", "regions", "targets")
	return(output)
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
	cat("debug: subsample\n"); flush.console()
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

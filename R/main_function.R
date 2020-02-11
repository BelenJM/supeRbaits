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
#' @param restrict A vector of chromosome names OR position numbers to which the analysis should be restricted.
#' @param gc A vector of two values between 0 and 1, specifying the minimum and maximmum GC percentage allowed in the output baits.
#' @param show.times Logical: Should time spent at each step be displayed?
#' @param debug Logical: Should debug files be created?
#' @param verbose Logical: Should detailed bait processing messages be displayed per sequence?
#' 
#' @return A dataframe of baits
#' 
#' @export
#' 
main_function <- function(n, size, database, exclusions = NULL, 
	regions = NULL, regions.prop = NULL, regions.tiling = NULL,
	targets = NULL, targets.prop = NULL, targets.tiling = NULL,
	seed = NULL, restrict = NULL, show.times = FALSE, debug = FALSE, gc = c(0.3, 0.5),
	verbose = FALSE){

	if (debug) {
    message("!!!--- Debug mode has been activated ---!!!")
		on.exit(save(list = ls(), file = "supeRbaits_debug.RData"), add = TRUE)
	} 
	if (show.times)
		message("Start time: ", Sys.time())
	
	flush.console()

	if(!is.null(seed))
		set.seed(seed)

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

	# Import data
	message("M: Compiling the sequences' lengths. This process can take some seconds."); flush.console()
	getlengths.time <- system.time({
		the.lengths <- callr::r(function(path) {supeRbaits:::getChromLengths(path = path)}, args = list(path = database), spinner = TRUE)
	})

	if (show.times)
		print(getlengths.time)
	
	if (is.numeric(restrict)) {
		if (all(restrict > nrow(the.lengths)))
			stop("'restrict' is set to numeric values, but ALL listed values are larger than the number of available sequences\nNumber of available sequences: ", nrow(the.lengths), "\nMin. index requested: ", min(restrict), "\nMax. index requested: ", max(restrict), "\n", call. = FALSE)
		if (max(restrict) > nrow(the.lengths)) {
			warning("'restrict' is set to numeric values, but some listed values are larger than the number of available sequences. Running analysis on matching sequences.\nDiscarded index values: ", paste(restrict[restrict > nrow(the.lengths)], collapse = " "), immediate. = TRUE, call. = FALSE)
			restrict <- restrict[restrict <= nrow(the.lengths)]
		}
		the.lengths <- the.lengths[restrict, ]
	}

	if (is.character(restrict)) {
		link <- match(restrict, the.lengths$name)
		if (all(is.na(link)))
			stop("None of the sequence names listed in 'restrict' matches the available sequences.\n", call. = FALSE)
		if (any(is.na(link))) {
			warning("Some sequences listed in 'restrict' do not match the available sequences. Running analysis on matching sequences\nMissing sequences: '", paste(restrict[is.na(link)], collapse = "', '"), "'", immediate. = TRUE, call. = FALSE)
			link <- link[!is.na(link)]
		}
		the.lengths <- the.lengths[link, ]
		rm(link)
	}

	message("M: Loading exclusions/regions/targets (if any is present)."); flush.console()

	load.extras.time <- system.time({
		if(!is.null(exclusions))
			exclusions <- load_exclusions(file = exclusions)
		if(!is.null(regions))
			regions <- load_regions(file = regions)
		if(!is.null(targets))
			targets <- load_targets(file = targets)
	})
	if (show.times)
		print(load.extras.time)

	# Compatibility checks
	message("M: Checking exclusions/regions/targets quality (if any is present)."); flush.console()

	check.input.time <- system.time({
		if (any(size > the.lengths$size))
			stop("'size' is larger than at least one of the chromosome lengths.\n")
		recipient <- check_chr_names(exclusions = exclusions, regions = regions, targets = targets, the.lengths = the.lengths)
		exclusions <- recipient$exclusions
		regions <- recipient$regions
		targets <- recipient$targets
		rm(recipient)

		check_chr_boundaries(exclusions = exclusions, regions = regions, targets = targets, the.lengths = the.lengths)
	})
	if (show.times)
		print(check.input.time)

	message("M: Finding bait positions for each sequence."); flush.console()
	if (!verbose)
  	pb <- txtProgressBar(min = 0, max = nrow(the.lengths), initial = 0, style = 3, width = 60)

	sample.baits.time <- system.time({
		bait.points <- lapply(1:nrow(the.lengths), function(i) {
			# extract relevant parameters
			params <- trim_parameters(chr = the.lengths$name[i], exclusions = exclusions, regions = regions, targets = targets)
			# region baits
			if(!is.null(regions.prop) && regions.prop > 0) {
				n.regions = n * regions.prop
				if (!is.null(params$regions)) {
					temp.regions <- region_baits(chr.length = the.lengths$size[i], n = n.regions, size = size, tiling = regions.tiling,
					regions = params$regions, exclusions = params$exclusions, chr = the.lengths$name[i], used.baits = NULL, verbose = verbose)
					temp.regions$Type <- rep("region", nrow(temp.regions))
					n.regions = nrow(temp.regions)
					used.baits <- temp.regions$Start
				} else {
					if (verbose)
						message("M: No regions found for chromosome ", the.lengths$name[i], ".")
					temp.regions <- NULL
					n.regions = 0
					used.baits <- NULL
				}
			} else {
				temp.regions <- NULL
				n.regions = 0
				used.baits <- NULL
			}
			# targetted baits
			if(!is.null(targets.prop) && targets.prop > 0) {
				n.targets = n * targets.prop
				if (!is.null(params$targets)) {
					temp.targets <- target_baits(chr.length = the.lengths$size[i], n = n.targets, size = size, tiling = targets.tiling,
						targets = params$targets, exclusions = params$exclusions, chr = the.lengths$name[i], used.baits = used.baits, verbose = verbose)
					temp.targets$Type <- rep("target", nrow(temp.targets))
					n.targets = nrow(temp.targets)
					used.baits <- c(used.baits, temp.targets$Start)
				} else {
					if (verbose)
						message("M: No targets found for chromosome ", the.lengths$name[i], ".")
					temp.targets <- NULL
					n.targets = 0
				}
			} else {
				temp.targets <- NULL
				n.targets = 0
			}
			# random baits
			n.random <- n - (n.regions + n.targets)
			if (n.random > 0) {
				temp.random <- random_baits(chr.length = the.lengths$size[i], n = n.random, size = size, 
					exclusions = params$exclusions, chr = the.lengths$name[i], used.baits = used.baits, verbose = verbose)
				temp.random$Type <- rep("random", nrow(temp.random))
			} else {
				temp.random <- NULL
			}
			# bring together the different parts
			output <- rbind(temp.regions, temp.targets, temp.random)
	    if (!verbose) {
	    	setTxtProgressBar(pb, i) # Progress bar
	    	flush.console()
	    }
			return(output)
		})
	  if (!verbose)
			close(pb)
		names(bait.points) <- as.character(the.lengths[, 1])
	})
	if (show.times)
		print(sample.baits.time)

	message("M: Retrieving bait base pairs. This operation can take some time."); flush.console()

	message("Temp: Running getBaits for the whole content."); flush.console()
	to.fetch <- data.table::rbindlist(bait.points, use.names = TRUE, idcol = "ChromName")

	getbaits.time <- system.time({
		baits <- getBaits(db = database, df = to.fetch)
	})

	if (show.times)
		print(getbaits.time)

	message("M: Calculating GC content in the baits"); flush.console()

	calc.baits.time <- system.time({
		baits$pGC <- baits$no_GC / (size + 1)
		baits <- split(baits, baits$bait_chrom_name)
	})

	if (show.times)
		print(calc.baits.time)

	message("M: Examining GC content in the baits"); flush.console()

	assess.baits.time <- system.time({
		good.baits <- lapply(seq_along(baits), function(i) {
			link <- baits[[i]]$pGC > gc[1] & baits[[i]]$pGC < gc[2]
			if (verbose) {
				if (all(!link)) {
					message(paste0("M: No baits passed the GC percentage test for chromosome ", names(baits)[i], "."))
					return(NULL)
				}
				if (any(!link)) {
					if (sum(!link) == 1)
						message(paste0("M: ", sum(!link), " bait was excluded from chr ", names(baits)[i], " due to its GC percentage."))
					else
						message(paste0("M: ", sum(!link), " baits were excluded from chr ", names(baits)[i], " due to their GC percentage."))
				}
			}
			return(baits[[i]][link, ])
		})
		names(good.baits) <- names(baits)
	})

	if (show.times)
		print(assess.baits.time)

	message("M: Analysis completed."); flush.console()

	if (show.times)
		return(list(baits = baits, good.baits = good.baits, chr.lengths = the.lengths, exclusions = exclusions, 
			targets = targets, regions = regions, getlengths.time = getlengths.time["elapsed"], 
			sample.baits.time = sample.baits.time["elapsed"], getbaits.time = getbaits.time["elapsed"]))
	else
		return(list(baits = baits, good.baits = good.baits, chr.lengths = the.lengths, 
			exclusions = exclusions, targets = targets, regions = regions))
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
	output <- list(exclusions  = NULL, regions = NULL, targets = NULL)
	if (!is.null(exclusions)) 
		output[[1]] <- subsample(input = exclusions, link = chr)
	if (!is.null(regions)) 
		output[[2]] <- subsample(input = regions, link = chr)
	if (!is.null(targets)) 
		output[[3]] <- subsample(input = targets, link = chr)
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
	chr <- link
	link <- grepl(link, input[, 1])
	if (any(link)) {
		output <- input[link, ]
		output <- output[order(output[, 2]), ]
		rownames(output) <- 1:nrow(output)
	}	else {
		output <- NULL
	}
	if (!is.null(output)) {
		# IF we are subsampling targets, remove duplicated entries
		if (ncol(output) == 2) {
			if (any(duplicated(output[,2]))) {
				output <- output[!duplicated(output[,2]), ]
			}
		}
		if (ncol(output) == 3) {
			if (nrow(output) > 1) {
				# combine areas contained within each other or contiguous
				try.again <- TRUE
				while (try.again) {
					link <- c(FALSE, (output[-1, 2] <= output[-nrow(output), 3] | output[-1, 2] == (output[-nrow(output), 3] + 1)))
					if (any(link)) {
						rows.to.fix <- which(link)
						for(i in rows.to.fix) {
							output[i - 1, 3] <- output[i, 2]
						}
						output <- output[!link, ]
					} else {
						try.again <- FALSE
					}
				}
			}
		}
	}
	return(output)
}

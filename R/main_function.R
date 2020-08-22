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
#' @param gc A vector of two values between 0 and 1, specifying the minimum and maximum GC percentage allowed in the output baits.
#' @param verbose Logical: Should detailed bait processing messages be displayed per sequence?
#' 
#' @return A dataframe of baits
#' 
#' @export
#' 
main_function <- function(n, size, database, exclusions = NULL, 
	regions = NULL, regions.prop = 0, regions.tiling = 1,
	targets = NULL, targets.prop = 0, targets.tiling = 1,
	seed = NULL, restrict, gc = c(0, 1),
	verbose = FALSE){

	if (!is.null(options("supeRbaits_debug")[[1]]) && options("supeRbaits_debug")[[1]]) {
    message("!!!--- Debug mode has been activated ---!!!")
		on.exit(save(list = ls(), file = "supeRbaits_debug.RData"), add = TRUE)
	} 
	if (!is.null(options("supeRbaits_show_times")[[1]]) && options("supeRbaits_show_times")[[1]])
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
	size <- as.integer(size)

	if (!is.numeric(regions.prop) | (regions.prop > 1 | regions.prop < 0))
		stop("'regions.prop' must be numeric and between 0 and 1.\n")
	if (!is.null(targets.prop) && (!is.numeric(targets.prop) | (targets.prop > 1 | targets.prop < 0)))
		stop("'targets.prop' must be numeric and between 0 and 1.\n")
	if (!is.numeric(regions.tiling))
		stop("'regions.tiling' must be numeric.\n")
	if (!is.numeric(targets.tiling))
		stop("'targets.tiling' must be numeric.\n")

	if (regions.prop == 0 & !is.null(regions))
		warning("Regions were included but regions.prop = 0. No region baits will be produced.", call. = FALSE, immediate. = TRUE)
	if (targets.prop == 0 & !is.null(targets))
		warning("Regions were included but targets.prop = 0. No region baits will be produced.", call. = FALSE, immediate. = TRUE)

	if (sum(regions.prop, targets.prop) > 1)
		stop("The sum of 'regions.prop' and 'targets.prop' must not be greater than one.\n", call. = FALSE)
	
	if (length(gc) != 2)
		stop("Please provide two values in 'gc' (minimum and maximum percentage).\n", call. = FALSE)
	if (!is.numeric(gc))
		stop("'gc' must be numeric and contain two values between 0 and 1.\n")
	if (any(gc > 1) | any(gc < 0))
		stop("'gc' ranges must be between 0 and 1.\n")
	if (gc[1] > gc[2])
		stop("The first value of 'gc' must be smaller or equal to the second value.\n")

	# Import data
	message("M: Compiling the sequences' lengths. This operation can take some time."); flush.console()

	# check line ending type
	first_line <- readLines(database, n = 1)
	if (length(first_line) == 0)	
		stop("The database file appears to be empty. Are you sure this is the correct file?", call. = FALSE)

	len_first_line <- nchar(first_line)
	if (grepl("\r$", readChar(database, len_first_line + 1)))
		stop("The line endings of your database file are incompatible with supeRbaits. You can convert your database using convert_line_endings()\n", call. = FALSE)
	# --

	getlengths.time <- system.time({
		the.lengths <- callr::r(function(getChromLengths, path) 
			{
				getChromLengths(path = path)
			}, 
			args = list(getChromLengths = getChromLengths,
									path = database),
			spinner = TRUE,
			show = TRUE)
	})

	if (!is.null(options("supeRbaits_show_times")[[1]]) && options("supeRbaits_show_times")[[1]])
		print(getlengths.time)
	
	if (!missing(restrict)) {
		restrict <- check_restrict(restrict, sequences = the.lengths$name)
		the.lengths <- the.lengths[restrict, ]
	}
	
	if (any(!is.null(exclusions), !is.null(regions), !is.null(targets)))
	message("M: Loading exclusions/regions/targets."); flush.console()

	load.extras.time <- system.time({
		if(!is.null(exclusions))
			exclusions <- load_exclusions(file = exclusions)
		if(!is.null(regions))
			regions <- load_regions(file = regions)
		if(!is.null(targets))
			targets <- load_targets(file = targets)
	})
	if (!is.null(options("supeRbaits_show_times")[[1]]) && options("supeRbaits_show_times")[[1]])
		print(load.extras.time)

	# Compatibility checks
	if (any(!is.null(exclusions), !is.null(regions), !is.null(targets)))
	message("M: Checking exclusions/regions/targets quality."); flush.console()

	check.input.time <- system.time({
		if (any((size - 1) > the.lengths$size))
			stop("'size' is larger than at least one of the chromosome lengths.\n")
		recipient <- check_chr_names(exclusions = exclusions, regions = regions, targets = targets, the.lengths = the.lengths)
		exclusions <- recipient$exclusions
		regions <- recipient$regions
		targets <- recipient$targets
		rm(recipient)

		check_chr_boundaries(exclusions = exclusions, regions = regions, targets = targets, the.lengths = the.lengths)
	})
	if (!is.null(options("supeRbaits_show_times")[[1]]) && options("supeRbaits_show_times")[[1]])
		print(check.input.time)

	message("M: Finding bait positions for each sequence."); flush.console()

	sample.baits.time <- system.time({
		bait.points <- callr::r(function(sampleBaits,
																		 chrom_lens, 
																		 exclusions, 
																		 regions, 
																		 targets, 
																		 n, 
																		 size, 
																		 regions_tiling, 
																		 targets_tiling, 
																		 regions_prop, 
																		 targets_prop) 
			{
				sampleBaits(chrom_lens, 
										exclusions, 
										regions, 
										targets, 
										n, 
										size, 
										regions_tiling, 
										targets_tiling, 
										regions_prop, 
										targets_prop)
			},
			args = list(sampleBaits = sampleBaits,
									chrom_lens = the.lengths, 
									exclusions = exclusions, 
									regions = regions, 
									targets = targets, 
									n = n, 
									size = size, 
									regions_tiling = regions.tiling, 
									targets_tiling = targets.tiling, 
									regions_prop = regions.prop, 
									targets_prop = targets.prop),
			spinner = !verbose,
			show = verbose)
	})

	if (!is.null(options("supeRbaits_show_times")[[1]]) && options("supeRbaits_show_times")[[1]])
		print(sample.baits.time)

	message("M: Retrieving bait base pairs. This operation can take some time."); flush.console()

	getbaits.time <- system.time({
		baits <- callr::r(function(getBaits,
															 db, 
															 df)
			{
				getBaits(db = db, 
								 df = df)
			}, 
			args = list(getBaits = getBaits,
									db = database,
									df = bait.points), 
			spinner = TRUE,
			show = TRUE)
	})

	if (!is.null(options("supeRbaits_show_times")[[1]]) && options("supeRbaits_show_times")[[1]])
		print(getbaits.time)

	if (nrow(baits) == 0)
		stop("No baits could be generated for any of the sequences. Aborting.\n", call. = FALSE)

	message("M: Calculating GC content in the baits."); flush.console()

	calc.baits.time <- system.time({
		baits$pGC <- baits$no_GC / size
		baits <- split(baits, baits$bait_chrom_name)
	})

	if (!is.null(options("supeRbaits_show_times")[[1]]) && options("supeRbaits_show_times")[[1]])
		print(calc.baits.time)

	if (gc[1] > 0 | gc[2] < 1) {
		message("M: Examining GC content in the baits."); flush.console()

		good.baits <- list()
		bad.baits <- list()

		assess.baits.time <- system.time({
			capture <- lapply(seq_along(baits), function(i) {
				link <- baits[[i]]$pGC > gc[1] & baits[[i]]$pGC < gc[2]
				if (verbose) {
					if (all(!link)) {
						message(paste0("M: No baits passed the GC percentage test for sequence ", names(baits)[i], "."))
						return(NULL)
					}
					if (any(!link)) {
						if (sum(!link) == 1)
							message(paste0("M: ", sum(!link), " bait was excluded from sequence ", names(baits)[i], " due to its GC percentage."))
						else
							message(paste0("M: ", sum(!link), " baits were excluded from sequence ", names(baits)[i], " due to their GC percentage."))
					}
				}
				good.baits[[i]] <<- baits[[i]][link, ]
				bad.baits[[i]] <<- baits[[i]][!link, ]
			})
			names(good.baits) <- names(baits)
			names(bad.baits) <- names(baits)
		})

		remove.empty.good <- sapply(good.baits, nrow) > 0
		good.baits <- good.baits[remove.empty.good]

		remove.empty.bad <- sapply(bad.baits, nrow) > 0
		bad.baits <- bad.baits[remove.empty.bad]
		
		if (!is.null(options("supeRbaits_show_times")[[1]]) && options("supeRbaits_show_times")[[1]])
			print(assess.baits.time)
	} else {
		good.baits <- baits
		bad.baits <- NULL
	}

	message("M: Analysis completed."); flush.console()

	input.summary <- list(chr.lengths = the.lengths, 
												exclusions = exclusions, 
												targets = targets, 
												regions = regions,
												size = size)

	if (!is.null(options("supeRbaits_show_times")[[1]]) && options("supeRbaits_show_times")[[1]]) {
		times <- data.frame(
			getLengths = getlengths.time["elapsed"],
			sampleBaits = sample.baits.time["elapsed"],
			getBaits = getbaits.time["elapsed"])
		output <- list(baits = good.baits, excluded.baits = bad.baits, input.summary = input.summary, times = times)
	}	else {
		output <- list(baits = good.baits, excluded.baits = bad.baits, input.summary = input.summary)
	}

	return(output)
}

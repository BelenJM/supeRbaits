#' Generates random baits
#'
#' Function inspired by this comment in StackOverflow: https://stackoverflow.com/questions/49149839/simulate-random-positions-from-a-list-of-intervals
#'
#' @param n Number of baits to generate (distributed across the various sequences).
#' @param n.per.seq Number of baits to generate per sequence. Ignored if n is set.
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
#' @param min.per.seq Minimum number of baits per sequence. Defaults to 1.
#' @param verbose Logical: Should detailed bait processing messages be displayed per sequence?
#' 
#' @return A dataframe of baits
#' 
#' @export
#' 
main_function <- function(n, n.per.seq, size, database, exclusions = NULL, 
	regions = NULL, regions.prop = 0, regions.tiling = 1,
	targets = NULL, targets.prop = 0, targets.tiling = 1,
	seed = NULL, restrict, gc = c(0, 1), min.per.seq = 1,
	verbose = FALSE){

	if (getOption("supeRbaits_debug", default = FALSE)) {
    message("!!!--- Debug mode has been activated ---!!!")
		on.exit(save(list = ls(), file = "supeRbaits_debug.RData"), add = TRUE)
	}
	if (getOption("supeRbaits_show_times", default = FALSE))
		message("Start time: ", Sys.time())
	
	flush.console()

	if(!is.null(seed)) {
		set.seed(seed)
		on.exit(set.seed(NULL))
	}

	if (missing(n) & missing(n.per.seq))
		stop("One of 'n' or 'n.per.seq' must be set.", call. = FALSE)

	if (!missing(n) & !missing(n.per.seq))
		warning("Both 'n' and 'n.per.seq' were set. Ignoring 'n.per.seq'.", immediate. = TRUE, call. = FALSE)

	if (!missing(n) && !is.numeric(n))
		stop("'n' must be numeric.", call. = FALSE)
	if (!missing(n) && length(n) != 1)
		stop("Length of 'n' must be 1.", call. = FALSE)
	if (!missing(n) && n < min.per.seq)
		stop("'n' cannot be lower than min.per.seq.", call. = FALSE)

	if ((missing(n) & !missing(n.per.seq)) && !is.numeric(n.per.seq))
		stop("'n.per.seq' must be numeric.", call. = FALSE)
	if ((missing(n) & !missing(n.per.seq)) && length(n.per.seq) != 1)
		stop("Length of 'n.per.seq' must be 1.", call. = FALSE)
	if ((missing(n) & !missing(n.per.seq)) && n.per.seq < min.per.seq)
		stop("'n.per.seq' cannot be lower than min.per.seq.", call. = FALSE)


	if (!is.numeric(size))
		stop("'size' must be numeric.", call. = FALSE)
	if (length(size) != 1)
		stop("Length of 'size' must be 1.", call. = FALSE)
	size <- as.integer(size)

	if (!is.numeric(regions.prop) | (regions.prop > 1 | regions.prop < 0))
		stop("'regions.prop' must be numeric and between 0 and 1.", call. = FALSE)
	if (!is.null(targets.prop) && (!is.numeric(targets.prop) | (targets.prop > 1 | targets.prop < 0)))
		stop("'targets.prop' must be numeric and between 0 and 1.", call. = FALSE)
	if (!is.numeric(regions.tiling))
		stop("'regions.tiling' must be numeric.", call. = FALSE)
	if (!is.numeric(targets.tiling))
		stop("'targets.tiling' must be numeric.", call. = FALSE)

	if (regions.prop == 0 & !is.null(regions))
		warning("Regions were included but regions.prop = 0. No region baits will be produced.", call. = FALSE, immediate. = TRUE)
	if (targets.prop == 0 & !is.null(targets))
		warning("Regions were included but targets.prop = 0. No region baits will be produced.", call. = FALSE, immediate. = TRUE)

	if (sum(regions.prop, targets.prop) > 1)
		stop("The sum of 'regions.prop' and 'targets.prop' must not be greater than one.\n", call. = FALSE)
	
	if (length(gc) != 2)
		stop("Please provide two values in 'gc' (minimum and maximum percentage).\n", call. = FALSE)
	if (!is.numeric(gc))
		stop("'gc' must be numeric and contain two values between 0 and 1.", call. = FALSE)
	if (any(gc > 1) | any(gc < 0))
		stop("'gc' ranges must be between 0 and 1.", call. = FALSE)
	if (gc[1] > gc[2])
		stop("The first value of 'gc' must be smaller or equal to the second value.", call. = FALSE)

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

	# extract sequence lengths
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
	if (getOption("supeRbaits_show_times", default = FALSE))
		print(getlengths.time)
	
	# apply restrictions, if present
	if (!missing(restrict)) {
		restrict <- check_restrict(restrict, sequences = the.lengths$name)
		the.lengths <- the.lengths[restrict, ]
	}

	# Remove sequences that are too small
	link <- the.lengths$size >= (min.per.seq + size - 1)
	if (all(!link))
		stop("No sequences have enough space to fit the minimum desired number of baits.", call. = FALSE)
	if (verbose & any(link))
		warning(sum(!link), " sequence", ifelse(sum(!link) > 1, "s are", " is"), " too small to fit the minimum desired number of baits and will be discarded.", immediate. = TRUE, call. = FALSE)

	# prepare to allocate n's
	the.lengths <- the.lengths[link, ]
	the.lengths$prop <- the.lengths$size / sum(the.lengths$size)
	the.lengths$n <- 0
	
	# allocate n's
	if (!missing(n)) {
		the.lengths <- assign_n_per_seq(the.lengths, n, min.per.seq)
	} else {
		the.lengths$n <- n.per.seq
		the.lengths$n[the.lengths$size - (size - 1) < n.per.seq] <- the.lengths$size[the.lengths$size - (size - 1) < n.per.seq] - (size - 1)
		if (sum(the.lengths$n) < n.per.seq * nrow(the.lengths))
			warning("Some sequences are not long enough to allocate the desired number of baits.", immediate. = TRUE, call. = FALSE)
	}

	# remove unneeded prop column from the lengths
	the.lengths$prop <- NULL

	# load additional parameters
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
	if (getOption("supeRbaits_show_times", default = FALSE))
		print(load.extras.time)

	# Compatibility checks
	if (any(!is.null(exclusions), !is.null(regions), !is.null(targets)))
	message("M: Checking exclusions/regions/targets quality."); flush.console()

	check.input.time <- system.time({
		recipient <- check_chr_names(exclusions = exclusions, regions = regions, targets = targets, the.lengths = the.lengths)
		exclusions <- recipient$exclusions
		regions <- recipient$regions
		targets <- recipient$targets
		rm(recipient)

		check_chr_boundaries(exclusions = exclusions, regions = regions, targets = targets, the.lengths = the.lengths)
	})
	if (getOption("supeRbaits_show_times", default = FALSE))
		print(check.input.time)

	message("M: Finding bait positions for each sequence."); flush.console()

	sample.baits.time <- system.time({
		bait.points <- callr::r(function(sampleBaits,
																		 chrom_lens, 
																		 exclusions, 
																		 regions, 
																		 targets, 
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
									size = size, 
									regions_tiling = regions.tiling, 
									targets_tiling = targets.tiling, 
									regions_prop = regions.prop, 
									targets_prop = targets.prop),
			spinner = !verbose,
			show = verbose)
	})

	if (getOption("supeRbaits_show_times", default = FALSE))
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

	if (getOption("supeRbaits_show_times", default = FALSE))
		print(getbaits.time)

	if (nrow(baits) == 0)
		stop("No baits could be generated for any of the sequences. Aborting.\n", call. = FALSE)

	message("M: Calculating GC content in the baits."); flush.console()

	calc.baits.time <- system.time({
		baits$pGC <- baits$no_GC / size
		baits <- split(baits, baits$bait_chrom_name)
	})

	if (getOption("supeRbaits_show_times", default = FALSE))
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
		
		if (getOption("supeRbaits_show_times", default = FALSE))
			print(assess.baits.time)
	} else {
		good.baits <- baits
		bad.baits <- NULL
	}

	message("M: Analysis completed."); flush.console()

	input.summary <- list(chr.lengths = data.table::as.data.table(the.lengths), 
												exclusions  = data.table::as.data.table(exclusions), 
												targets     = data.table::as.data.table(targets), 
												regions     = data.table::as.data.table(regions),
												size        = size)

	if (getOption("supeRbaits_show_times", default = FALSE)) {
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

#' Distribute the overall number of baits by the multiple sequences
#' 
#' @param the.lengths a dataframe with the length of each sequence, and the respective proportion of the overall size.
#' @inheritParams main_function
#' 
#' @return An updated the.lengths dataframe, with the number of baits to be extracted per sequence.
#' 
#' @keywords internal
#' 
assign_n_per_seq <- function(the.lengths, n, min.per.seq, size) {
	# if the number of baits to distribute is lower than the number of sequences
	if (n < (nrow(the.lengths) * min.per.seq)) {
		warning("The desired n is smaller than the number of sequences times min.per.seq. Only a subset of the sequences will be used.", immediate. = TRUE, call. = FALSE)
		# decide how many sequences we can pick
		n.seqs <- floor(n / min.per.seq)
		# pick sequences
		which.seqs <- sample(1:nrow(the.lengths), n.seqs, prob = the.lengths$prop)
		# distribute min.per.seq for the chosen sequences
		the.lengths$n[which.seqs] <- min.per.seq
		# if there are baits missing
		if (sum(the.lengths$n) < n) {
			# find out which sequences still have space
			seqs.with.space <- which(the.lengths$size[which.seqs] - min.per.seq - size > 0)
			# choose which from the above will receive + 1
			add.here <- sample(which.seqs[seqs.with.space], n - sum(the.lengths$n))
			# add one bait on the above
			the.lengths$n[add.here] <- the.lengths$n[add.here] + 1
		}
	# if the number of baits to assign is bigger than the number of sequences
	} else {
		# start by assigning the minimum necessary (sequences that are too small have already been excluded)
		the.lengths$n <- min.per.seq
		# distribute the rest of the baits
		while (sum(the.lengths$n) < n) {
			n.to.distribute <- n - sum(the.lengths$n)
			# find which sequences have space
			seqs.with.space <- which(the.lengths$size - the.lengths$n - (size - 1) > 0)
			if (length(seqs.with.space) == 0) {
				warning("Ran out of space in the sequences where to allocate baits. Could not allocate ", n.to.distribute, " baits.", immediate. = TRUE, call. = FALSE)
				break()
			}
			# recalculate weights
			free.space <- the.lengths$size - the.lengths$n - (size - 1)
			the.lengths$prop <- free.space / sum(free.space)
			# draw new places where to add baits based on new weights
			add.here <- sample(seqs.with.space, n.to.distribute, replace = TRUE, prob = the.lengths$prop)
			sort.results <- table(add.here)
			# add new baits in their respective places
			the.lengths$n[as.numeric(names(sort.results))] <- the.lengths$n[as.numeric(names(sort.results))] + sort.results
			# if any of the sequences went over the limit, drop extra baits to be redistributed
			the.lengths$n <- sapply(1:nrow(the.lengths), function(i) min(the.lengths$size[i] - (size - 1), the.lengths$n[i]))
			# If the last line changed anything, then the while will be triggered again
		}
	}
	return(the.lengths)
}

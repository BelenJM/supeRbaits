#' Return bps
#' 
#' @inheritParams main_function
#' @inheritParams region_baits
#' 
#' @return logical
#' 
#' @export
#' 
retrieve_baits <- function(database) {
	# cat("debug: retrieve_baits\n"); flush.console()
	path <- paste(system.file(package = "supeRbaits"), "new_retrieveBait.py", sep="/")
	bait.file <- "temp_folder_for_supeRbaits/bait_positions.txt"
	output.file <- "temp_folder_for_supeRbaits/bait_positions_py.txt"
	aux.py <- system2("python", args = c(shQuote(path), shQuote(database), shQuote(bait.file)), stdout = TRUE)
	if(!is.null(attr(aux.py, "status")))
		stop("Python failed to retrieve the bps.")
	output <- data.table::fread(text = aux.py, sep = "\t")
	return(output)
}

#' check if Phython is installed before attempting to run python code
#' 
#' @keywords internal
#' 
#' @return logical
#' 
check_python <- function () {
 trigger <- suppressWarnings(try(reticulate::py_config(), silent = TRUE))
  if (inherits(trigger, "try-error"))
   	stop("Python modules could not be loaded. Please install python before using supeRbaits.\n", call. = FALSE)
  if (substr(trigger$version, 1, 1) == 2)
  	stop("supeRbaits requires python3, but R is currently using python2.\n", call. = FALSE)
}

#' Python handler for lengthChrom
#' 
#' Finds the length of each chromosome
#' 
#' @inheritParams main_function
#' 
#' @return a data frame with the chromosome names and their lengths
#' 
#' @keywords internal
#' 
get_lengths <- function(database, restrict = NULL) {
	cat("debug: get_lengths\n"); flush.console()
	if (file.exists("temp_folder_for_supeRbaits/genome_size.txt"))
		file.remove("temp_folder_for_supeRbaits/genome_size.txt")
	path <- paste(system.file(package = "supeRbaits"), "new_lengthChrom.py", sep="/")
	if (is.null(restrict))
		try(suppressWarnings(response <- system2("python", args = c(shQuote(path), shQuote(database)))), silent = TRUE)
	else 
		try(suppressWarnings(response <- system2("python", args = c(shQuote(path), shQuote(database)), shQquote(restrict))), silent = TRUE)
	if (!file.exists("temp_folder_for_supeRbaits/genome_size.txt"))
		stop("Python failed to retrieve the chromosome lengths.")
}

#' Return bps
#' 
#' @inheritParams main_function
#' @inheritParams region_baits
#' 
#' @return logical
#' 
#' @export
#' 
retrieve_baits <- function(chr, database) {
	cat("debug: retrieve_baits\n"); flush.console()
	path <- paste(system.file(package = "supeRbaits"), "retrieveBait.py", sep="/")
	try(suppressWarnings(response <- system2("python", args = c(path, chr, database), stdout = TRUE)), silent = T)
	if(!is.null(attr(response,"status")))
		stop("Python failed to retrieve bps from chromosome ", chr, ".")
	else 
		return(response)
}

#' Return bps new
#' 
#' @param chr the chromosome name
#' @param positions the start and stop positions of the baits to be retrieved
#' @inheritParams main_function
#' 
#' @return logical
#' 
#' @export
#' 
retrieve_baits_new <- function(chr, positions, database) {
  reticulate::source_python(paste0(system.file(package = "supeRbaits"), "/retrieveBait.py"))
	return(processing_baitfile(chr, positions))
}

#' Return bps new
#'
#' @return logical
#' 
#' @export
#' 
py_test <- function(x) {
  reticulate::source_python(paste0(system.file(package = "supeRbaits"), "/my_function.py"))
	return(my_function(x))
}

#' check if Phython is installed before attempting to run python code
#' 
#' @keywords internal
#' 
#' @return logical
#' 
check_python <- function () {
 trigger <- try(reticulate::py_config(), silent = TRUE)
  if (inherits(trigger, "try-error")) {
   	cat("Error: Python modules could not be loaded. Please install python before using baits4pop.\n")
  	return(FALSE)
  } else {
 	  return(TRUE)
  }
}

#' Return bps
#' 
#' @param chr the chromosome name
#' @param positions the start and stop positions of the baits to be retrieved
#' @inheritParams main_function
#' 
#' @return logical
#' 
#' @export
#' 
retrieve_baits <- function(chr, positions, database) {
	path <- paste(system.file(package = "baits4pop"), "retrieveBait.py", sep="/")
	command <- paste("python", path, chr, positions, sep = " ")
	try(suppressWarnings(response <- system2("python", args = c(path, chr, positions), stdout = TRUE)), silent = T)
	if(!is.null(attr(response,"status")))
		stop("Failed to retrieve bps from chromosome ", chr, ".")
	else 
		return(response)
}
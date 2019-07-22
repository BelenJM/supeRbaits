#' check if Phython is installed before attempting to run python code
#' 
#' @keywords internal
#' 
#' @return logical
#' 
check_python <- function () {
  if (Sys.which("python") == "") {
  	cat("Error: Python modules could not be loaded. Please install python before using baits4pop.\n")
  	return(FALSE)
  } else {
 	  return(TRUE)
  }
}


#' Return bps
#' 
#' @return logical
#' 
#' @export
#' 
retrieve_bps <- function(chr, positions) {
	path <- paste(system.file(package = "baits4pop"), "parse.py", sep="/")
	  
	command <- paste("python", path, chr, positions, sep = " ")
	try(suppressWarnings(response <- system(command, 
	                                        intern = T,
	                                        ignore.stderr = TRUE)), silent = T)
	if(!is.null(attr(response,"status")))
		stop("Failed to retrieve bps.")
	else 
		return(response)
}
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

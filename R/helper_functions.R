#' Convert line endings from DOS to UNIX format
#' 
#' database files must have UNIX line endings for supeRbaits to work.
#' 
#' @param input The name of the input file
#' @param output The name of the output file (the same as the input by default)
#' 
#' @return No return value. Called for side effects.
#' 
#' @export
#' 
convert_line_endings <- function(input, output = input) {
	if (input == output) {
		aux <- tempfile()
		dos2unix(input, aux)
		file.rename(aux, output)
	} else {
		dos2unix(input, output)
	}
}

#' Ensure the database nucleotide lines have the same length
#' 
#' Database files must have the same number of nucleotides per line (except for the terminal lines of each sequence).
#' 
#' @param input The name of the input file
#' @param output The name of the output file (the same as the input by default)
#' 
#' @return No return value. Called for side effects.
#' 
#' @export
#' 
standardize_lines <- function(input, output = input) {
	if (input == output) {
		aux <- tempfile()
		fixLineLengths(input, aux)
		file.rename(aux, output)
	} else {
		fixLineLengths(input, output)
	}
}

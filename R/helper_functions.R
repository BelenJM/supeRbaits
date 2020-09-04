#' Convert line endings from DOS to UNIX format
#' 
#' Database files must have UNIX line endings for supeRbaits to work.
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

#' Extract specific nucleotide sequences
#' 
#' @param database A database of sequences
#' @param seq The name of the sequence
#' @param from the first nucleotide position to be extracted
#' @param to the last nucleotide position to be extracted
#' 
#' @return the nucleotide sequence
#' 
#' @export
#' 
extract_nucleotides <- function(database, seq, from, to) {
	if (length(database) > 1)
		stop("Please select only one database.")
	
	if (length(seq) > 1)
		stop("Please select only one sequence.")

	if (from <= 0)
		stop("'from' must be positive.")
	
	if (to <= 0)
		stop("'to' must be positive.")

	if (length(from) > 1)
		stop("Please select only one starting point.")

	if (length(from) > 1)
		stop("Please select only one ending point.")

	if (to < from)
		stop("'to' should not be lower than 'from'.")
	
	extractNucleotides(database, seq, from, to)
}

#' Forcefully round a number down
#'
#' Forces the rounding of the input to the next lower rounded value.
#' 
#' @param input The value to be rounded
#' @param to The level of rounding to be applied (i.e. to=10 will round 14.8 to 10; to=1 will round i to 14)
#' 
#' @return The rounded value
#' 
#' @export
#' 
roundDown <- function(input, to = 10) {
  to * (input%/%to)
}

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
	dos2unix(input, output)
}
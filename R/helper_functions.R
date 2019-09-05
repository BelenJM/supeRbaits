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

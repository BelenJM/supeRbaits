.onAttach <- function(libname, pkgname){
  packageStartupMessage("You need python to run this!")
}

.onLoad <- function(libname, pkgname) {
  SeqIO <<- reticulate::import("SeqIO", delay_load = TRUE)
}

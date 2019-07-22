.onAttach <- function(libname, pkgname){
  packageStartupMessage("You need python to run this!")
}

.onLoad <- function(libname, pkgname) {
  SeqIO <<- reticulate::import("SeqIO", delay_load = TRUE)
  reticulate::source_python(paste0(system.file(package = "baits4pop"), "/baits4pop_extract_sequence_genome2.py"))
}

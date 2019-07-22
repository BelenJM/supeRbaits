.onAttach <- function(libname, pkgname){
  trigger <- try(reticulate::py_config(), silent = TRUE)
  if (inherits(trigger, "try-error")) {
  	packageStartupMessage("Error: Python modules could not be loaded. Please install python before using baits4pop.")
  } else {
 	  packageStartupMessage("Python modules successfully loaded!")
  }
}

.onLoad <- function(libname, pkgname) {
  SeqIO <<- reticulate::import("SeqIO", delay_load = TRUE)
  trigger <- try(reticulate::py_config(), silent = TRUE)
  if (!inherits(trigger, "try-error"))
  	reticulate::source_python(paste0(system.file(package = "baits4pop"), "/baits4pop_extract_sequence_genome2.py"))
}

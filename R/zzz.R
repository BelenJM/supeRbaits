.onAttach <- function(libname, pkgname){
  if (reticulate::py_available()) {
  	reticulate::source_python(paste0(system.file(package = "baits4pop"), "/baits4pop_extract_sequence_genome2.py"))
 	  packageStartupMessage("Python modules successfully loaded!")
 	  cat("test")
  } else {
  	packageStartupMessage("Error: Python modules could not be loaded. Please install python before using baits4pop.")
  }
}

.onLoad <- function(libname, pkgname) {
  SeqIO <<- reticulate::import("SeqIO", delay_load = TRUE)
}

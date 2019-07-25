.onAttach <- function(libname, pkgname){
 trigger <- try(reticulate::py_config(), silent = TRUE)
  if (inherits(trigger, "try-error")) {
   	packageStartupMessage("Error: Could not find Python modules. Please install python before using baits4pop.")
  } else {
  	reticulate::source_python(paste0(system.file(package = "baits4pop"), "/baits4pop_extract_sequence_genome2.py"))
 	  packageStartupMessage("Python modules successfully loaded!")
 	  cat("test")
  }
  packageStartupMessage("Welcome to baits4pop!")
}

# .onLoad <- function(libname, pkgname) {
#   SeqIO <<- reticulate::import("SeqIO", delay_load = TRUE)
# }

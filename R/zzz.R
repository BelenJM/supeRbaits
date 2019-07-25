.onAttach <- function(libname, pkgname){
 # trigger <- try(reticulate::py_config(), silent = TRUE)
 #  if (inherits(trigger, "try-error")) {
 #   	packageStartupMessage("Error: Could not find Python modules. Please install python before using baits4pop.")
 #  } else {
 # 	  packageStartupMessage("Python modules successfully loaded!")
 # 	  cat("test")
 #  }
  packageStartupMessage("Welcome to baits4pop!")

}

# .onLoad <- function(libname, pkgname) {
#   SeqIO <<- reticulate::import("SeqIO", delay_load = TRUE)
#   pyth <<-  reticulate::source_python(paste0(system.file(package = "baits4pop"), "/retrieveBait.py"))
# }

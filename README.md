# supeRbaits
## An R-package to design baits for capture sequencing experiments 

## What can supeRbaits help you with?
You are a researcher/lab manager wanting to carry out a capture sequencing experiment, where you want to design your own set of baits from a genome sequence of your species of interest. The R-package supeRbaits, written in C++ and implemented in R, can help you with this task.


<img src="vignettes/mb_arrays.svg" alt="drawing" width="870"/>

### Main functions:

**1. explore()**

 explore() allows you to quickly get a summary of your data. You can use explore() to get a general feel for the study results, and check if the input files are behaving as expected. It is also a good candidate if you just want to validate your detections for later use in other analyses.
 
**2. migration()**
```
## install the package
setwd("C:/Users/bmen/Documents/GitHub/supeRbaits/")
require("devtools")
install()

setwd("C:\\Users\\bmen\\Documents\\GitHub\\supeRbaits\\testdata")
# devtools::load_all()
```

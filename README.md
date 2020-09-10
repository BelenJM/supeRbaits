# supeRbaits<img src="vignettes/mb_arrays.svg" alt="drawing" width="870"/>
## An R-package to design baits for capture sequencing experiments 

## What can supeRbaits help you with?
You are a researcher/lab manager wanting to carry out a capture sequencing experiment, where you want to design your own set of baits from a genome sequence of your species of interest. The R-package supeRbaits, written in C++ and implemented in R, can help you with this task.

### How to start:
**Installation:
```
## install the package
setwd("C:/Users/bmen/Documents/GitHub/supeRbaits/")
require("devtools")
install()

setwd("C:\\Users\\bmen\\Documents\\GitHub\\supeRbaits\\testdata")
# devtools::load_all()
```

### Main function and its options:
main() is the main function of the package. The different options are:
**1.n**
 With the n option, you can specify the number of baits you would like to generate, in total.
 
**2.size**
 Size of each of the baits to generate. An example: 120 base pairs.
 
**3.database**

**4.exclusions**
**5.regions**
**6.targets**
**7.seed**
**8.restrict**
**9.gc**
**10.debug**
**11.verbose**

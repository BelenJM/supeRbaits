# supeRbaits<img src="supeRbaits.png" align="left" width="100" />

## An R-package to design baits for capture sequencing experiments 

## What can supeRbaits help you with?
You are a researcher/lab manager wanting to carry out a capture sequencing experiment, where you want to design your own set of baits from a genome sequence of your species of interest. The R-package supeRbaits, written in C++ and implemented in R, can help you with this task.

### How to start:
**Installation:**
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

**1.Number of total baits (n)**
 
 With the n option, you can specify the number of total baits you would like to generate.
 
**2.Size of bait (size)**

Desired size of each of the baits to generate. An example: 120 base pairs.
 
**3.Genomic dataset to use as a reference (database)**

The genomic information available from the species of interest that you want the baits to be designed from. The database can be of any format, and it can be at a chromosome or contig level. 

**4.Areas of the database to exclude (exclusions)**

Use this option if you want to exclude certain areas from your genomic database and not generate baits from those. An example would be to exclude areas with a specific GC content, or repeated regions, depending on the type of analysis you are interested in performing. 

**5.Areas of the database that you want to specifically include (regions)**

This option allows you to specify regions of the genomic database that you are very interested in including within your baits. A region file could be a set of genes, where for each gene you have an interval of base pairs where you are interested in having baits designed from.

**6.targets**

You can include another type of regions of interest where you want some baits to be designed from. Targets would typically consist on Single Nucleotide Polymorphisms (SNP) where you know the position at the genomic database where they are located, and you want to design a bait from such area.

**7.seed**

You can specify a seed number in order to reproduce your bait design. 

**8.restrict**


**9.gc**

Specific range of GC content where you want the baits to fall.

**10.debug**

If set up as TRUE, you have the option to debug the script. Default set to FALSE.

**11.verbose**

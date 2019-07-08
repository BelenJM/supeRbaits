# Implementation BEDtools in R

# bedTools.random<-function(functionstring="random",l,n,seed,g,out,opt.string="")
# {
#   #create temp files
#   n = tempfile()
#   l = tempfile()
#   seed = tempfile()
#   g = tempfile()
#   out = tempfile()
#   options(scipen =99) # not to use scientific notation when writing out
#
#   #write bed formatted dataframes to tempfile
#   write.table(g,file=g,quote=F,sep="\t",col.names=F,row.names=F)
#
#   # create the command string and call the command using system()
#   command=paste(functionstring,"-n",n,"-l",l,"-seed", seed,
#                 "-g", g, ">", out, sep=" ")
#   cat(command,"\n")
#   try(system(command))
#
#   res=read.table(out,header=F)
#   unlink(l);unlink(n);unlink(seed);unlink(g);unlink(out)
#   return(res)
# }

#### FUNCTION 1 ####
# Function that generates random sequences in chromosome positions
# Input:
# - nsites: number of sites to generate
# - chr_file: file provided by the user with three columns: chromosome name, start (bp), end (bp)
# - name of the output file
# Function inspired by this comment in StackOverflow: https://stackoverflow.com/questions/49149839/simulate-random-positions-from-a-list-of-intervals


# Function
random_areas <- function(nsites, chr_file,name_file){
nSites <- nsites
# read the bed file
bed <- read.table(file=chr_file, header=T, sep='\t')
# convert it to a data.table
bed2 <- data.table::as.data.table(bed, keep.colnames=TRUE)

# calculate size of range
bed2[, size := 1 + end-start]

# Randomly sample bed file rows, proportional to the length of each range
simulated.sites <- bed2[sample(.N, size=nSites, replace=TRUE, prob=bed2$size)]

# Randomly sample uniformly within each chosen range
simulated.sites[, position := sample(start:end, size=1), by=1:dim(simulated.sites)[1]]

# Remove extra columns and format as needed
simulated.sites[, start  := position]
simulated.sites[, end := position]
simulated.sites[, c("size", "position") := NULL]
out_table <- simulated.sites
write.table(out_table, name_file)
return(out_table)
}


#### FUNCTION 2 ####
# bedtools_intersect("random_areas", out1=TRUE, bedtools=12)
#
#
# bedTools.random("random",100,10,3569,"chrSize_example.txt","test")
#
# install.packages("HelloRanges")
# library(HelloRanges)
# library(HelloRangesData)
# code <- bedtools_intersect("-a cpg.bed -b exons.bed")
# code{genome <- Seqinfo(genome = NA_character_)gr_a <- import("cpg.bed", genome = genome)gr_b <- import("exons.bed", genome = genome)pairs <- findOverlapPairs(gr_a, gr_b, ignore.strand = TRUE)ans <- pintersect(pairs, ignore.strand = TRUE)ans}

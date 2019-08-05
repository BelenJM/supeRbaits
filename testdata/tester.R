setwd("C:\\Users\\hdmfla\\Documents\\GitHub\\baits4pop\\testdata")
# devtools::load_all()
library(supeRbaits)
# all_random test
main_function(n = 10, size = 20, database = "chrom_salmon_chunk.fasta.txt")
# all_targetted test
main_function(n = 10, size = 20, lengths = "chrSize_example.txt", 
	regions = "repeatMasker_example.txt")
# trimmed_random test
main_function(n = 10, size = 20, lengths = "chrSize_example.txt", 
	exclusions = "repeatMasker_example.txt")
# trimmed_targetted test
main_function(n = 10, size = 20, lengths = "chrSize_example.txt", 
	regions = "repeatMasker_example.txt", exclusions = "repeatMasker_example.txt")

# tests to the subsample function's fail-safes
input <- data.frame(name = c("A", "A"), Start = c(2, 4))
subsample(input, "A") # all good
#   name Start
# 1    A     2
# 2    A     4
input <- data.frame(name = c("A", "A"), Start = c(2, 2))
subsample(input, "A") # error - duplicated start
# Error in subsample(input, "A") : 
#   There are duplicated starting points/targets for one of the trimming elements of chromosome A.
input <- data.frame(name = c("A", "A"), Start = c(3, 2))
subsample(input, "A") # Should reorder rows
#   name Start
# 2    A     2
# 1    A     3
input <- data.frame(name = c("A", "A"), Start = c(1, 4), stop = c(2, 5))
subsample(input, "A") # all good
#   name Start stop
# 1    A     1    2
# 2    A     4    5
input <- data.frame(name = c("A", "A"), Start = c(1, 4), stop = c(5, 5))
subsample(input, "A") # error - duplicated end
# Error in subsample(input, "A") : 
#  There are duplicated ending points for one of the trimming elements for chromosome A.
input <- data.frame(name = c("A", "A"), Start = c(1, 4), stop = c(3, 5))
subsample(input, "A") # error - contiguous
# Error in subsample(input, "A") : 
#   There are contiguous regions to include or exclude for chromosome A. Please list these as a single region.
input <- data.frame(name = c("A", "A"), Start = c(1, 3), stop = c(3, 5))
subsample(input, "A") # error - overlapping
# Error in subsample(input, "A") : 
#   The regions to include or exclude overlap for chromosome A.

main_function(n = 10, size = 20, lengths = "chrSize_example_simple.txt",
	exclusions = "repeatMasker_example.txt")


size <- 20
lengths <- baits4pop:::load_lengths("chrSize_example.txt")
exclusions <- baits4pop:::load_exclusions("repeatMasker_example.txt")

# for chr 1
chr <- "ssa01"
length <- lengths[1, 2]
exclusions <- exclusions[1:3, ]

# for chr 2
chr <- "ssa02"
length <- lengths[2, 2]
exclusions <- exclusions[4, ]


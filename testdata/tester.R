# Start by moving to the main repository folder (i.e. GitHub\supeRbaits)

library("supeRbaits")

cd(".."); Rcpp::compileAttributes(); devtools::document(); cd("testdata/")

# all random test
	x <- main_function(n = 10, size = 20, gc = c(0.2, 0.8), database = "chrom_salmon_chunk.fasta.txt", regions.prop = 0, targets.prop = 0)

# exclusions tests
	x <- main_function(n = 10, size = 20, gc = c(0.2, 0.8), database = "chrom_salmon_chunk.fasta.txt", exclusions = "exclusion_example-bad_name.txt")
	# should fail due to wrong names in the exclusions file

	x <- main_function(n = 10, size = 20, gc = c(0.2, 0.8), database = "chrom_salmon_chunk.fasta.txt", exclusions = "exclusion_example-bad_length.txt")
	# should fail due to wrong lengths in the exclusions file

	x <- main_function(n = 10, size = 20, gc = c(0.2, 0.8), database = "chrom_salmon_chunk.fasta.txt", exclusions = "exclusion_example.txt")
	test_exclusions <- read.table("exclusion_example.txt")
	any(apply(test_exclusions, 1, function(e) x[[1]][[1]]$Start_bp >= e[2] & x[[1]][[1]]$Start_bp <= e[3])) # Should return false
	any(apply(test_exclusions, 1, function(e)   x[[1]][[1]]$End_bp >= e[2] & x[[1]][[1]]$End_bp   <= e[3])) # Should return false
	# Works nicely.

	# Try larger n
	x <- main_function(n = 100, size = 20, gc = c(0.2, 0.8), database = "chrom_salmon_chunk.fasta.txt", exclusions = "exclusion_example.txt")
	any(apply(test_exclusions, 1, function(e) x[[1]][[1]]$Start_bp >= e[2] & x[[1]][[1]]$Start_bp <= e[3])) # Should return false
	any(apply(test_exclusions, 1, function(e)   x[[1]][[1]]$End_bp >= e[2] & x[[1]][[1]]$End_bp   <= e[3])) # Should return false
	# Works nicely as well.

# regions tests
	x <- main_function(n = 10, size = 20, gc = c(0.2, 0.8), database = "chrom_salmon_chunk.fasta.txt", regions = "region_example.txt")
	# Returns error, as expected
	# Please include the desired percentage of regional baits in 'regions.prop'.
	x <- main_function(n = 10, size = 20, gc = c(0.2, 0.8), database = "chrom_salmon_chunk.fasta.txt", 
		regions = "region_example.txt", regions.prop = 1)
	test_regions <- read.table("region_example.txt")
	any(apply(test_regions, 1, function(r) x[[1]][[1]]$Start_bp > r[3] & x[[1]][[1]]$Start_bp < r[2])) # Should return false
	any(apply(test_regions, 1, function(r)   x[[1]][[1]]$End_bp > r[3] & x[[1]][[1]]$End_bp   < r[2])) # Should return false
	# Works nicely.

	# Try larger n
	x <- main_function(n = 100, size = 20, gc = c(0.2, 0.8), database = "chrom_salmon_chunk.fasta.txt", 
		regions = "region_example.txt", regions.prop = 1)
	# The first 32 baits should be region, the rest should be random
	any(apply(test_regions, 1, function(r) x[[1]][[1]]$Start_bp[1:32] > r[3] & x[[1]][[1]]$Start_bp[1:32] < r[2])) # Should return false
	any(apply(test_regions, 1, function(r)   x[[1]][[1]]$End_bp[1:32] > r[3] & x[[1]][[1]]$End_bp[1:32]   < r[2])) # Should return false
	# Works nicely as well.

# region + exclusions test
	x <- main_function(n = 100, size = 20, gc = c(0.2, 0.8), database = "chrom_salmon_chunk.fasta.txt", 
		exclusions = "exclusion_example.txt", regions = "region_example.txt", regions.prop = 1)
	# The first 12 baits should be region, the rest should be random
	any(apply(test_regions, 1, function(r) x[[1]][[1]]$Start_bp[1:12] > r[3] & x[[1]][[1]]$Start_bp[1:12] < r[2])) # Should return false
	any(apply(test_regions, 1, function(r)   x[[1]][[1]]$End_bp[1:12] > r[3] & x[[1]][[1]]$End_bp[1:12]   < r[2])) # Should return false
	# No baits should be within the exclusion zones
	any(apply(test_exclusions, 1, function(e) x[[1]][[1]]$Start_bp >= e[2] & x[[1]][[1]]$Start_bp <= e[3])) # Should return false
	any(apply(test_exclusions, 1, function(e)   x[[1]][[1]]$End_bp >= e[2] & x[[1]][[1]]$End_bp   <= e[3])) # Should return false
	# all good.

	# random baits should not be duplicates of the region baits
	any(duplicated(x[[1]][[1]]$Start_bp)) # should return false

# targets test
	x <- main_function(n = 10, size = 20, gc = c(0.2, 0.8), database = "chrom_salmon_chunk.fasta.txt", 
		targets = "targets_example.txt", targets.prop = 1)
	test_targets <- read.table("targets_example.txt")
	!all(sapply(x[[1]][[1]]$Start_bp, function(x_i) {
		any(sapply(test_targets[, 2], function(t) {
			x_i[1] >= (t - 19) & x_i[1] <= t
		}))
	})) # Should return false

# targets + exclusions test
	x <- main_function(n = 10, size = 20, gc = c(0.2, 0.8), database = "chrom_salmon_chunk.fasta.txt", 
		exclusions = "exclusion_example.txt", targets = "targets_example.txt", targets.prop = 1)
	!all(sapply(x[[1]][[1]]$Start_bp, function(x_i) {
		any(sapply(test_targets[, 2], function(t) {
			x_i[1] >= (t - 19) & x_i[1] <= t
		}))
	})) # Should return false
	any(apply(test_exclusions, 1, function(e) x[[1]][[1]]$Start_bp >= e[2] & x[[1]][[1]]$Start_bp <= e[3])) # Should return false
	any(apply(test_exclusions, 1, function(e)   x[[1]][[1]]$End_bp >= e[2] & x[[1]][[1]]$End_bp   <= e[3])) # Should return false
	
	#try larger n
	x <- main_function(n = 100, size = 20, gc = c(0.2, 0.8), database = "chrom_salmon_chunk.fasta.txt", 
		exclusions = "exclusion_example.txt", targets = "targets_example.txt", targets.prop = 1)
	# only the first 24 are targetted
	!all(sapply(x[[1]][[1]]$Start_bp[1:24], function(x_i) {
		any(sapply(test_targets[, 2], function(t) {
			x_i[1] >= (t - 19) & x_i[1] <= t
		}))
	})) # Should return false
	any(apply(test_exclusions, 1, function(e) x[[1]][[1]]$Start_bp >= e[2] & x[[1]][[1]]$Start_bp <= e[3])) # Should return false
	any(apply(test_exclusions, 1, function(e)   x[[1]][[1]]$End_bp >= e[2] & x[[1]][[1]]$End_bp   <= e[3])) # Should return false
	# all looking good

# regions + targets test
	x <- main_function(n = 10, size = 20, gc = c(0.2, 0.8), database = "chrom_salmon_chunk.fasta.txt", 
		targets = "targets_example.txt", targets.prop = 1, regions = "region_example.txt", regions.prop = 1)
	# Returns error: The sum of 'regions.prop' and 'targets.prop' must not be greated than one.

	x <- main_function(n = 10, size = 20, gc = c(0.2, 0.8), database = "chrom_salmon_chunk.fasta.txt", 
		targets = "targets_example.txt", targets.prop = 0.5, regions = "region_example.txt", regions.prop = 0.5)
	# the first 5 should be regions
	any(apply(test_regions, 1, function(r) x[[1]][[1]]$Start_bp[1:5] > r[3] & x[[1]][[1]]$Start_bp[1:5] < r[2])) # Should return false
	any(apply(test_regions, 1, function(r)   x[[1]][[1]]$End_bp[1:5] > r[3] & x[[1]][[1]]$End_bp[1:5]   < r[2])) # Should return false
	# the last 5 should be targets
	!all(sapply(x[[1]][[1]]$Start_bp[6:10], function(x_i) {
		any(sapply(test_targets[, 2], function(t) {
			x_i[1] >= (t - 19) & x_i[1] <= t
		}))
	})) # Should return false

	# try with greater n
	x <- main_function(n = 100, size = 20, gc = c(0.2, 0.8), database = "chrom_salmon_chunk.fasta.txt", 
		targets = "targets_example.txt", targets.prop = 0.5, regions = "region_example.txt", regions.prop = 0.5)
	# the first 32 should be regions
	any(apply(test_regions, 1, function(r) x[[1]][[1]]$Start_bp[1:32] > r[3] & x[[1]][[1]]$Start_bp[1:32] < r[2])) # Should return false
	any(apply(test_regions, 1, function(r)   x[[1]][[1]]$End_bp[1:32] > r[3] & x[[1]][[1]]$End_bp[1:32]   < r[2])) # Should return false
	# the following 40 should be targets
	!all(sapply(x[[1]][[1]]$Start_bp[33:72], function(x_i) {
		any(sapply(test_targets[, 2], function(t) {
			x_i[1] >= (t - 19) & x_i[1] <= t
		}))
	})) # Should return false
	# the rest should be random.

# regions + targets + exclusions test
	x <- main_function(n = 10, size = 20, gc = c(0.2, 0.8), database = "chrom_salmon_chunk.fasta.txt", exclusions = "exclusion_example.txt", 
		targets = "targets_example.txt", targets.prop = 0.5, regions = "region_example.txt", regions.prop = 0.5)
	# the first 5 should be regions
	any(apply(test_regions, 1, function(r) x[[1]][[1]]$Start_bp[1:5] > r[3] & x[[1]][[1]]$Start_bp[1:5] < r[2])) # Should return false
	any(apply(test_regions, 1, function(r)   x[[1]][[1]]$End_bp[1:5] > r[3] & x[[1]][[1]]$End_bp[1:5]   < r[2])) # Should return false
	# the last 5 should be targets
	!all(sapply(x[[1]][[1]]$Start_bp[6:10], function(x_i) {
		any(sapply(test_targets[, 2], function(t) {
			x_i[1] >= (t - 19) & x_i[1] <= t
		}))
	})) # Should return false
	# and none should be in the exclusions
	any(apply(test_exclusions, 1, function(e) x[[1]][[1]]$Start_bp >= e[2] & x[[1]][[1]]$Start_bp <= e[3])) # Should return false
	any(apply(test_exclusions, 1, function(e)   x[[1]][[1]]$End_bp >= e[2] & x[[1]][[1]]$End_bp   <= e[3])) # Should return false
	# looking good

	# try with larger n
	x <- main_function(n = 100, size = 20, gc = c(0.2, 0.8), database = "chrom_salmon_chunk.fasta.txt", exclusions = "exclusion_example.txt", 
		targets = "targets_example.txt", targets.prop = 0.5, regions = "region_example.txt", regions.prop = 0.5)
	# the first 12 should be regions
	any(apply(test_regions, 1, function(r) x[[1]][[1]]$Start_bp[1:12] > r[3] & x[[1]][[1]]$Start_bp[1:12] < r[2])) # Should return false
	any(apply(test_regions, 1, function(r)   x[[1]][[1]]$End_bp[1:12] > r[3] & x[[1]][[1]]$End_bp[1:12]   < r[2])) # Should return false
	# the following 24 should be targets
	!all(sapply(x[[1]][[1]]$Start_bp[13:36], function(x_i) {
		any(sapply(test_targets[, 2], function(t) {
			x_i[1] >= (t - 19) & x_i[1] <= t
		}))
	})) # Should return false
	# the rest should be random.
	# None should be in the exclusions
	any(apply(test_exclusions, 1, function(e) x[[1]][[1]]$Start_bp >= e[2] & x[[1]][[1]]$Start_bp <= e[3])) # Should return false
	any(apply(test_exclusions, 1, function(e)   x[[1]][[1]]$End_bp >= e[2] & x[[1]][[1]]$End_bp   <= e[3])) # Should return false
	# and finally none should be duplicated
	any(duplicated(x[[1]][[1]]$Start_bp)) # should return false
	
# All seems to be up and running.



# Test print functions

	gc_table(baits = x$good.baits)
	gc_table(baits = x$good.baits, combine = TRUE)

	# only works working on debug environment
	print_coverage(chr.lengths = the.lengths, baits = good.baits, exclusions = exclusions, targets = targets)
	coverage(chr.lengths = the.lengths, baits = good.baits, exclusions = exclusions)


# Test temporary variable useR

	x <- main_function(n = 100, size = 20, gc = c(0.2, 0.8), database = "chrom_salmon_chunk.fasta.txt", exclusions = "exclusion_example.txt", 
		targets = "targets_example.txt", targets.prop = 0.5, regions = "region_example.txt", regions.prop = 0.5, useR = FALSE)


# ----------------	

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

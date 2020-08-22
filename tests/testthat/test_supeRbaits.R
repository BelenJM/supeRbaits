path.aux <- system.file(package = "supeRbaits")[1]

# set up test directory
home.wd <- getwd()
setwd(tempdir())
dir.create("test_supeRbaits")
setwd("test_supeRbaits")
file.copy(paste0(path.aux, "/", list.files(path.aux)), list.files(path.aux), overwrite = TRUE)

convert_line_endings(paste0(path.aux, "/", "sequences.txt"), "sequences.txt")
# ---

# all random test
test_that("random extration is working", {
	x <- main_function(n = 10, size = 20, database = "sequences.txt")
	expect_equal(names(x), c("baits", "excluded.baits", "input.summary"))
	expect_equal(x$baits$CMF$bait_seq, c("GCATATCCCAAAATTTCFFF", "CATATCCCAAAATTTCFFFF", "ATATCCCAAAATTTCFFFFF"))
})

test_that("saturated results are equal to reference", {
	x <- main_function(n = 200, size = 10, database = "sequences.txt")
	## RUN THESE LINES ONLY TO RESET THE REFERENCE!
	# saturated_simple_results <- x
	# save(saturated_simple_results, file = paste0(home.wd, "/saturated_simple_results.RData"))
	load(paste0(home.wd, "/saturated_simple_results.RData"))
	expect_equal(x, saturated_simple_results)
})

# exclusions tests
test_that("extraction with exclusions is working", {
	expect_warning(main_function(n = 10, size = 20, gc = c(0.2, 0.8), database = "sequences.txt", exclusions = paste0(home.wd, "/exclusions_bad_name.txt")),
			"Not all of the sequences' names in the exclusions match the names listed in the database. Removing orphan exclusions.", fixed = TRUE)

	expect_error(main_function(n = 10, size = 20, gc = c(0.2, 0.8), database = "sequences.txt", exclusions = paste0(home.wd, "/exclusions_bad_length.txt")),
		"Exclusion data for sequence CM003279 is off-boundaries.", fixed = TRUE)

	x <- main_function(n = 10, size = 20, gc = c(0.2, 0.8), database = "sequences.txt", exclusions = "exclusions.txt")
	test_exclusions <- read.table("exclusions.txt")
	expect_false(any(apply(test_exclusions, 1, function(e) x$baits$CM003279$Start_bp >= e[2] & x$baits$CM003279$Start_bp <= e[3])))
	expect_false(any(apply(test_exclusions, 1, function(e)   x$baits$CM003279$End_bp >= e[2] & x$baits$CM003279$End_bp   <= e[3])))

	# Try larger n
	x <- main_function(n = 100, size = 20, gc = c(0.2, 0.8), database = "sequences.txt", exclusions = "exclusions.txt")
	expect_false(any(apply(test_exclusions, 1, function(e) x$baits$CM003279$Start_bp >= e[2] & x$baits$CM003279$Start_bp <= e[3])))
	expect_false(any(apply(test_exclusions, 1, function(e)   x$baits$CM003279$End_bp >= e[2] & x$baits$CM003279$End_bp   <= e[3])))
})

# regions tests
test_that("extraction with regions is working", {
	expect_warning(main_function(n = 10, size = 20, gc = c(0.2, 0.8), database = "sequences.txt", regions = "regions.txt"),
	"Regions were included but regions.prop = 0. No region baits will be produced.", fixed = TRUE)

	x <- main_function(n = 10, size = 20, gc = c(0.2, 0.8), database = "sequences.txt", regions = "regions.txt", regions.prop = 1)
	test_regions <- read.table("regions.txt")
	expect_false(any(apply(test_regions, 1, function(r) x$baits$CM003279$Start_bp > r[3] & x$baits$CM003279$Start_bp < r[2])))
	expect_false(any(apply(test_regions, 1, function(r)   x$baits$CM003279$End_bp > r[3] & x$baits$CM003279$End_bp   < r[2])))

	# Try larger n
	x <- main_function(n = 100, size = 20, gc = c(0.2, 0.8), database = "sequences.txt", 
		regions = "regions.txt", regions.prop = 1)
	# The first 32 baits should be region, the rest should be random
	aux <- table(x$baits$CM003279$bait_type)
	expect_equal(names(aux), c("random", "region"))
	expect_equal(as.vector(aux), c(68, 32))
	expect_false(any(apply(test_regions, 1, function(r) x$baits$CM003279$Start_bp[1:32] > r[3] & x$baits$CM003279$Start_bp[1:32] < r[2])))
	expect_false(any(apply(test_regions, 1, function(r)   x$baits$CM003279$End_bp[1:32] > r[3] & x$baits$CM003279$End_bp[1:32]   < r[2])))
})

# region + exclusions test
test_that("extraction with exclusions and regions is working", {
	x <- main_function(n = 100, size = 20, gc = c(0.2, 0.8), database = "sequences.txt", 
		exclusions = "exclusions.txt", regions = "regions.txt", regions.prop = 1)
	# The first 12 baits should be region, the rest should be random
	aux <- table(x$baits$CM003279$bait_type)
	expect_equal(names(aux), c("random", "region"))
	expect_equal(as.vector(aux), c(68, 32))
	expect_false(any(apply(test_regions, 1, function(r) x$baits$CM003279$Start_bp[1:12] > r[3] & x$baits$CM003279$Start_bp[1:12] < r[2])))
	expect_false(any(apply(test_regions, 1, function(r)   x$baits$CM003279$End_bp[1:12] > r[3] & x$baits$CM003279$End_bp[1:12]   < r[2])))
	# No baits should be within the exclusion zones
	expect_false(any(apply(test_exclusions, 1, function(e) x$baits$CM003279$Start_bp >= e[2] & x$baits$CM003279$Start_bp <= e[3])))
	expect_false(any(apply(test_exclusions, 1, function(e)   x$baits$CM003279$End_bp >= e[2] & x$baits$CM003279$End_bp   <= e[3])))

	# random baits should not be duplicates of the region baits
	expect_false(any(duplicated(x$baits$CM003279$Start_bp)))
})

# targets test
test_that("extraction with targets is working", {
	x <- main_function(n = 10, size = 20, gc = c(0.2, 0.8), database = "sequences.txt", 
		targets = "targets.txt", targets.prop = 1)
	test_targets <- read.table("targets.txt")
	expect_true( # all baits should be within target range
		all(sapply(x$baits$CM003279$Start_bp, function(x_i) {
			any(sapply(test_targets[, 2], function(t) {
				x_i[1] >= (t - 19) & x_i[1] <= t
			}))
		}))
	)
})

# targets + exclusions test
test_that("extraction with exclusions and targets is working", {
	x <- main_function(n = 10, size = 20, gc = c(0.2, 0.8), database = "sequences.txt", 
		exclusions = "exclusions.txt", targets = "targets.txt", targets.prop = 1)
	expect_true( # all baits should be within target range
		all(sapply(x$baits$CM003279$Start_bp, function(x_i) {
			any(sapply(test_targets[, 2], function(t) {
				x_i[1] >= (t - 19) & x_i[1] <= t
			}))
		}))
	)
	# no baits should be in the exclusion zones
	expect_false(any(apply(test_exclusions, 1, function(e) x$baits$CM003279$Start_bp >= e[2] & x$baits$CM003279$Start_bp <= e[3])))
	expect_false(any(apply(test_exclusions, 1, function(e)   x$baits$CM003279$End_bp >= e[2] & x$baits$CM003279$End_bp   <= e[3])))
	
	#try larger n
	x <- main_function(n = 100, size = 20, gc = c(0.2, 0.8), database = "sequences.txt", 
		exclusions = "exclusions.txt", targets = "targets.txt", targets.prop = 1)
	aux <- table(x$baits$CM003279$bait_type)
	expect_equal(names(aux), c("random", "target"))
	expect_equal(as.vector(aux), c(80, 20))

	# only the first 20 are targetted
	expect_true( # all target baits should contain the target
		all(sapply(x$baits$CM003279$Start_bp[1:20], function(x_i) {
			any(sapply(test_targets[, 2], function(t) {
				x_i[1] >= (t - 19) & x_i[1] <= t
			}))
		}))
	)
	# no baits should be in the exclusion zones
	expect_false(any(apply(test_exclusions, 1, function(e) x$baits$CM003279$Start_bp >= e[2] & x$baits$CM003279$Start_bp <= e[3])))
	expect_false(any(apply(test_exclusions, 1, function(e)   x$baits$CM003279$End_bp >= e[2] & x$baits$CM003279$End_bp   <= e[3])))
})

# regions + targets test
test_that("extraction with regions and targets is working", {
	expect_error(main_function(n = 10, size = 20, gc = c(0.2, 0.8), database = "sequences.txt", 
		targets = "targets.txt", targets.prop = 1, regions = "regions.txt", regions.prop = 1),
	"The sum of 'regions.prop' and 'targets.prop' must not be greater than one.", fixed = TRUE)

	x <- main_function(n = 10, size = 20, gc = c(0.2, 0.8), database = "sequences.txt", 
		targets = "targets.txt", targets.prop = 0.5, regions = "regions.txt", regions.prop = 0.5)
	aux <- table(x$baits$CM003279$bait_type)
	expect_equal(names(aux), c("region", "target"))
	expect_equal(as.vector(aux), c(5, 5))

	# the first 5 should be regions
	expect_false(any(apply(test_regions, 1, function(r) x$baits$CM003279$Start_bp[1:5] > r[3] & x$baits$CM003279$Start_bp[1:5] < r[2])))
	expect_false(any(apply(test_regions, 1, function(r)   x$baits$CM003279$End_bp[1:5] > r[3] & x$baits$CM003279$End_bp[1:5]   < r[2])))

	# the last 5 should be targets
	expect_true( # all target baits should contain the target
		all(sapply(x$baits$CM003279$Start_bp[6:10], function(x_i) {
			any(sapply(test_targets[, 2], function(t) {
				x_i[1] >= (t - 19) & x_i[1] <= t
			}))
		}))
	)

	# try with greater n
	x <- main_function(n = 100, size = 20, gc = c(0.2, 0.8), database = "sequences.txt", 
		targets = "targets.txt", targets.prop = 0.5, regions = "regions.txt", regions.prop = 0.5)

	aux <- table(x$baits$CM003279$bait_type)
	expect_equal(names(aux), c("random", "region", "target"))
	expect_equal(as.vector(aux), c(40, 28, 32))

	# the first 32 should be regions
	expect_false(any(apply(test_regions, 1, function(r) x$baits$CM003279$Start_bp[1:32] > r[3] & x$baits$CM003279$Start_bp[1:32] < r[2])))
	expect_false(any(apply(test_regions, 1, function(r)   x$baits$CM003279$End_bp[1:32] > r[3] & x$baits$CM003279$End_bp[1:32]   < r[2])))
	# the following 40 should be targets
	expect_true(
		all(sapply(x$baits$CM003279$Start_bp[33:72], function(x_i) {
			any(sapply(test_targets[, 2], function(t) {
				x_i[1] >= (t - 19) & x_i[1] <= t
			}))
		}))
	)
	# the rest should be random.
})

# regions + targets + exclusions test
test_that("extraction with exclusions, regions and targets is working", {
	x <- main_function(n = 10, size = 20, gc = c(0.2, 0.8), database = "sequences.txt", exclusions = "exclusions.txt", 
		targets = "targets.txt", targets.prop = 0.5, regions = "regions.txt", regions.prop = 0.5)

	aux <- table(x$baits$CM003279$bait_type)
	expect_equal(names(aux), c("region", "target"))
	expect_equal(as.vector(aux), c(5, 5))

	# the first 5 should be regions
	expect_false(any(apply(test_regions, 1, function(r) x$baits$CM003279$Start_bp[1:5] > r[3] & x$baits$CM003279$Start_bp[1:5] < r[2])))
	expect_false(any(apply(test_regions, 1, function(r)   x$baits$CM003279$End_bp[1:5] > r[3] & x$baits$CM003279$End_bp[1:5]   < r[2])))
	# the last 5 should be targets
	expect_true(
		all(sapply(x$baits$CM003279$Start_bp[6:10], function(x_i) {
			any(sapply(test_targets[, 2], function(t) {
				x_i[1] >= (t - 19) & x_i[1] <= t
			}))
		}))
	)
	# and none should be in the exclusions
	expect_false(any(apply(test_exclusions, 1, function(e) x$baits$CM003279$Start_bp >= e[2] & x$baits$CM003279$Start_bp <= e[3])))
	expect_false(any(apply(test_exclusions, 1, function(e)   x$baits$CM003279$End_bp >= e[2] & x$baits$CM003279$End_bp   <= e[3])))
	
	# try with larger n
	x <- main_function(n = 100, size = 20, gc = c(0.2, 0.8), database = "sequences.txt", exclusions = "exclusions.txt", 
		targets = "targets.txt", targets.prop = 0.5, regions = "regions.txt", regions.prop = 0.5)

	aux <- table(x$baits$CM003279$bait_type)
	expect_equal(names(aux), c("random", "region", "target"))
	expect_equal(as.vector(aux), c(64, 12, 24))

	# the first 12 should be regions
	expect_false(any(apply(test_regions, 1, function(r) x$baits$CM003279$Start_bp[1:12] > r[3] & x$baits$CM003279$Start_bp[1:12] < r[2])))
	expect_false(any(apply(test_regions, 1, function(r)   x$baits$CM003279$End_bp[1:12] > r[3] & x$baits$CM003279$End_bp[1:12]   < r[2])))
	
	# the following 24 should be targets
	expect_true(
		all(sapply(x$baits$CM003279$Start_bp[13:36], function(x_i) {
			any(sapply(test_targets[, 2], function(t) {
				x_i[1] >= (t - 19) & x_i[1] <= t
			}))
		}))
	)
	# the rest should be random.

	# None should be in the exclusions
	expect_fale(any(apply(test_exclusions, 1, function(e) x$baits$CM003279$Start_bp >= e[2] & x$baits$CM003279$Start_bp <= e[3])))
	expect_false(any(apply(test_exclusions, 1, function(e)   x$baits$CM003279$End_bp >= e[2] & x$baits$CM003279$End_bp   <= e[3])))

	# and finally none should be duplicated
	expect_false(any(duplicated(x$baits$CM003279$Start_bp)))
})


results <- main_function(n = 100, size = 20, gc = c(0.2, 0.8), database = "sequences.txt", exclusions = "exclusions.txt", 
		targets = "targets.txt", targets.prop = 0.5, regions = "regions.txt", regions.prop = 0.5)

# Test print functions

test_that("gc_table is working", {
	x <- gc_table(results)
	expect_true(is.list(x))
	expect_equal(names(x), c("CM003279", "CMF"))
	expect_equal(x$CM003279$n, as.vector(table(results$baits$CM003279$bait_type)))

	x <- gc_table(results, combine = TRUE)
	expect_true(is.data.frame(x))
	# need to expand this once the large n bug has been fixed
})

test_that("coverage is working", {
	x <- coverage(results)
	expect_true(is.list(x))
	expect_equal(names(x), c("CM003279", "CMF"))
	expect_equal(x$CM003279$bait_type, "All")

	x <- coverage(results, combine = FALSE)
	expect_true(is.list(x))
	expect_equal(names(x), c("CM003279", "CMF"))
	# expand once the large n bug has been fixed
})


	# print_coverage(results, "CMF")
	# print_coverage(results, "CM003279")

setwd(home.wd)
rm(list = ls())

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/do_baits.R
\name{do_baits}
\alias{do_baits}
\title{Generates random baits}
\usage{
do_baits(
  n,
  n.per.seq,
  size,
  database,
  lengths,
  exclusions = NULL,
  regions = NULL,
  regions.prop = 0,
  regions.tiling = 1,
  targets = NULL,
  targets.prop = 0,
  targets.tiling = 1,
  seed = NULL,
  restrict,
  gc = c(0, 1),
  min.per.seq = 1,
  verbose = FALSE,
  force = FALSE
)
}
\arguments{
\item{n}{Number of baits to generate (distributed across the various sequences).}

\item{n.per.seq}{Number of baits to generate per sequence. Ignored if n is set.}

\item{size}{The size of each bait}

\item{database}{A database of sequences}

\item{lengths}{Optional: the lengths of the sequences, compiled through \code{\link{load_lengths}}. If missing, lengths will be compiled from the database on-the-fly.}

\item{exclusions}{A file containing regions to exclude}

\item{regions}{A file containing regions of interest.}

\item{regions.prop}{The proportion of baits that should overlap the regions of interest.}

\item{regions.tiling}{The minimum number of baits to distribute per region of interest.}

\item{targets}{A file containing bp's to target.}

\item{targets.prop}{The proportion of baits that should overlap the targets.}

\item{targets.tiling}{The minimum number of baits to distribute per target.}

\item{seed}{A number to fix the randomization process, for reproducibility}

\item{restrict}{A vector of chromosome names OR position numbers to which the analysis should be restricted.}

\item{gc}{A vector of two values between 0 and 1, specifying the minimum and maximum GC percentage allowed in the output baits.}

\item{min.per.seq}{Minimum number of baits per sequence. Defaults to 1.}

\item{verbose}{Logical: Should detailed bait processing messages be displayed per sequence?}

\item{force}{Logical: Proceed even if the number of baits requested is very large?}
}
\value{
A dataframe of baits
}
\description{
Function inspired by this comment in StackOverflow: https://stackoverflow.com/questions/49149839/simulate-random-positions-from-a-list-of-intervals
}

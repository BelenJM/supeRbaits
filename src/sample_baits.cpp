// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <limits>
#include <vector>
#include <fstream>
#include <iostream>
#include <unordered_map>

/*
  test_df_targets <- data.frame("ChromName" = c("CM003279", "CMF", "CMF", "CMF", "CMF"),
  "Target" = c(72, 75, 111, 75, 11))

  test_df_ranges <- data.frame("ChromName" = c("CM003279", "CM003279"),
  "Start" = c(10, 70),
  "Stop" = c(50, 300))
*/

const int DF_NAME_INDEX   = 0,
          DF_LEN_INDEX    = 0,
          DF_START_INDEX  = 1,
          DF_STOP_INDEX   = 2,
          DF_TARGET_INDEX = 1;

std::vector<size_t> filter_targets(std::vector<size_t> targets, size_t len) {
  std::vector<size_t> uniq_targets;
  std::sort(targets.begin(), targets.end());

  if (targets[0] > len) {
    return uniq_targets;
  }
  size_t prev = targets[0];
  uniq_targets.push_back(prev);
  for (size_t i = 1; i < targets.size(); i++) {
    if (targets[i] > len) {
      break;
    }
    if (targets[i] != prev) {
      prev = targets[i];
      uniq_targets.push_back(prev);
    }
  }

  return uniq_targets;
}

std::vector<std::pair<size_t, size_t>> filter_ranges(std::vector<std::pair<size_t, size_t>> ranges, size_t len) {
  std::vector<std::pair<size_t, size_t>> uniq_ranges;
  std::sort(ranges.begin(), ranges.end());
  ranges.push_back(std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));

  if (ranges[0].first > len) {
    return uniq_ranges;
  }
  std::pair<size_t, size_t> prev = std::make_pair(ranges[0].first, std::min(ranges[0].second, len));
  for (size_t i = 1; i < ranges.size(); i++) {
    if (ranges[i].first <std::numeric_limits<size_t>::max() && ranges[i].first> len) {
      uniq_ranges.push_back(prev);
      break;
    }
    if (prev.second + 1 >= ranges[i].first) {
      prev.second = std::min(ranges[i].second, len);
    } else {
      uniq_ranges.push_back(prev);
      prev = std::make_pair(ranges[i].first, std::min(ranges[i].second, len));
    }
  }

  return uniq_ranges;
}

std::vector<size_t> subsample_targets(std::string name, size_t len, Rcpp::DataFrame df) {
  Rcpp::StringVector df_names = df[DF_NAME_INDEX];

  std::vector<size_t> targets;
  Rcpp::NumericVector df_targets = df[DF_TARGET_INDEX];
    
  for (size_t i = 0; i < df_names.size(); i++) {
    if (df_names[i] == name) {
      targets.push_back(df_targets[i]);
    }
  }
  return filter_targets(targets, len);
}

std::vector<std::pair<size_t, size_t>> subsample_ranges(std::string name, size_t len, Rcpp::DataFrame df) {
  Rcpp::StringVector df_names = df[DF_NAME_INDEX];
  
  std::vector<std::pair<size_t, size_t>> ranges;
  Rcpp::NumericVector df_starts = df[DF_START_INDEX], df_stops = df[DF_STOP_INDEX];
    
  for (size_t i = 0; i < df_names.size(); i++) {
    if (df_names[i] == name) {
      ranges.push_back(std::make_pair(df_starts[i], df_stops[i]));
    }
  }
  return filter_ranges(ranges, len);
}

// [[Rcpp::export]]
void sample_baits(Rcpp::DataFrame chrom_lens, Rcpp::DataFrame exclusions, Rcpp::DataFrame regions, Rcpp::DataFrame targets) {
  Rcpp::StringVector df_names = chrom_lens[DF_NAME_INDEX];
  Rcpp::NumericVector df_lens = chrom_lens[DF_LEN_INDEX];
  for (size_t i = 0; i < df_names.size(); i++) { // for each chrom/len in chrom_lens
    std::vector<size_t> target_subsamples = subsample_targets(Rcpp::as<std::string>(df_names[i]), df_lens[i], targets);
    std::vector<std::pair<size_t, size_t>> exclusion_subsamples = subsample_ranges(Rcpp::as<std::string>(df_names[i]), df_lens[i], exclusions);
    std::vector<std::pair<size_t, size_t>> regions_subsamples = subsample_ranges(Rcpp::as<std::string>(df_names[i]), df_lens[i], regions);
    
    // TODO: calculate baits for each type
  }
  return;
}

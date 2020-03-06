// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <limits>
#include <vector>
#include <fstream>
#include <iostream>
#include <unordered_map>

/*
  test_df_chrom_lens <- data.frame("ChromName" = c("CM003279", "CMF"),
  "Length" = c(250, 150))

  test_df_targets <- data.frame("ChromName" = c("CM003279", "CMF", "CMF", "CMF", "CMF"),
  "Target" = c(72, 75, 111, 75, 11))

  test_df_regions <- data.frame("ChromName" = c("CM003279", "CM003279"),
  "Start" = c(10, 70),
  "Stop" = c(50, 300))

  test_df_exclusions <- data.frame("ChromName" = c("CM003279", "CM003279"),
  "Start" = c(40, 80),
  "Stop" = c(70, 200))

  sampleBaits(test_df_chrom_lens, test_df_exclusions, test_df_regions, test_df_targets, 2, 2, 1, 1, 1, 0)
*/

typedef std::vector<size_t> vec;
typedef std::vector<std::pair<size_t, size_t>> vec_pair;

const int DF_NAME_INDEX   = 0,
          DF_LEN_INDEX    = 1,
          DF_START_INDEX  = 1,
          DF_STOP_INDEX   = 2,
          DF_TARGET_INDEX = 1;

struct RegionBait {
  size_t start;
  size_t stop;
  std::string type;
};

vec_pair trim_ranges(vec_pair ranges, vec_pair exclusions) {
  vec_pair trimmed_ranges;
  while(ranges.size() > 0) {
    std::pair<size_t, size_t> cur_range = ranges.back(); ranges.pop_back();
    for (auto cur_exclusion : exclusions) {
      if (cur_exclusion.first > cur_range.first) {
	if (cur_exclusion.first > cur_range.second) {
	    continue; // no intersection
	}
	if (cur_exclusion.second > cur_range.second) {
	  cur_range.second = cur_exclusion.first-1;
	} else {
	  ranges.push_back(std::make_pair(cur_exclusion.second+1, cur_range.second));
	  cur_range.second = cur_exclusion.first-1;
	}
      } else {
	if (cur_exclusion.second < cur_range.first) {
	  continue; // no intersection
	}
	if (cur_exclusion.second < cur_range.second) {
	  cur_range.first = cur_exclusion.second+1;
	} else {
	  cur_range.first = -1; break; // complete overlap
	} 
      }
    }
    if (cur_range.first != -1) {
      trimmed_ranges.push_back(cur_range);
    }
  }
  sort(trimmed_ranges.begin(), trimmed_ranges.end());
  return trimmed_ranges;
}

// std::vector<RegionBait> region_baits(std::string chrom_name,
// 				     size_t n,
// 				     size_t size,
// 				     size_t tiling,
// 				     vec_pair ranges,
// 				     vec_pair exclusions) {
//   if (exclusions.size() > 0) {
//     vec_pair trimmed_regions = trim_ranges(ranges, exclusions);
//   }
// }


vec filter_targets(vec targets, size_t len) {
  vec uniq_targets;
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

vec_pair filter_ranges(vec_pair ranges, size_t len) {
  vec_pair uniq_ranges;
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

vec subsample_targets(std::string name, size_t len, Rcpp::DataFrame df) {
  Rcpp::StringVector df_names = df[DF_NAME_INDEX];

  vec targets;
  Rcpp::NumericVector df_targets = df[DF_TARGET_INDEX];
    
  for (size_t i = 0; i < df_names.size(); i++) {
    if (df_names[i] == name) {
      targets.push_back(df_targets[i]);
    }
  }
  return filter_targets(targets, len);
}

vec_pair subsample_ranges(std::string name, size_t len, Rcpp::DataFrame df) {
  Rcpp::StringVector df_names = df[DF_NAME_INDEX];
  
  vec_pair ranges;
  Rcpp::NumericVector df_starts = df[DF_START_INDEX], df_stops = df[DF_STOP_INDEX];
    
  for (size_t i = 0; i < df_names.size(); i++) {
    if (df_names[i] == name) {
      ranges.push_back(std::make_pair(df_starts[i], df_stops[i]));
    }
  }
  return filter_ranges(ranges, len);
}

// [[Rcpp::export]]
void sampleBaits(Rcpp::DataFrame chrom_lens,
		 Rcpp::DataFrame exclusions,
		 Rcpp::DataFrame regions,
		 Rcpp::DataFrame targets,
		 size_t n,
		 size_t size,
		 size_t regions_tiling,
		 size_t targets_tiling,
		 double regions_prop,
		 double targets_prop) {
  Rcpp::StringVector df_names = chrom_lens[DF_NAME_INDEX];
  Rcpp::NumericVector df_lens = chrom_lens[DF_LEN_INDEX];
  for (size_t i = 0; i < df_names.size(); i++) { // for each chrom/len in chrom_lens
    vec targets_subsample = subsample_targets(Rcpp::as<std::string>(df_names[i]), df_lens[i], targets);
    vec_pair exclusions_subsample = subsample_ranges(Rcpp::as<std::string>(df_names[i]), df_lens[i], exclusions);
    vec_pair regions_subsample = subsample_ranges(Rcpp::as<std::string>(df_names[i]), df_lens[i], regions);

    // region baits
    if (regions_prop > 0) {
      size_t n_regions = n * regions_prop;
      if (regions_subsample.size() > 0) {
	Rcpp::Rcout << "Trimmed regions:" << std::endl;
	for (auto t : trim_ranges(regions_subsample, exclusions_subsample)) {
	  Rcpp::Rcout << "(" << t.first << ", " << t.second << ")" << std::endl;
	};
      } else {
    	Rcpp::Rcout << "M: No regions found for chromosome " << df_names[i] << "." << std::endl;
      }

    }
    // target baits
    if (targets_prop > 0) {
      size_t n_targets = n * targets_prop;
      if (targets_subsample.size() > 0) {
    	// TODO
      } else {
    	Rcpp::Rcout << "M: No targets found for chromosome " << df_names[i] << "." << std::endl;
      }
    }
  }
  return;
}

// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <limits>
#include <vector>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <stdlib.h>
#include <math.h>

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
typedef std::vector<std::pair<std::pair<size_t, size_t>, size_t>> vec_pair_value;
typedef std::vector<std::pair<std::pair<size_t, size_t>, std::string>> vec_pair_string;

const int DF_NAME_INDEX   = 0,
          DF_LEN_INDEX    = 1,
          DF_START_INDEX  = 1,
          DF_STOP_INDEX   = 2,
          DF_TARGET_INDEX = 1;

struct SampleBait {
  std::string chrom;
  std::string type;
  size_t start;
  size_t stop;
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
    
    if (cur_range.first != (size_t) -1) {
      trimmed_ranges.push_back(cur_range);
    }
  }
  sort(trimmed_ranges.begin(), trimmed_ranges.end());
  
  return trimmed_ranges;
}

vec filter_targets(vec targets, size_t len) {
  if (!targets.size())
    return vec(0);

  vec uniq_targets;
  std::sort(targets.begin(), targets.end());

  if (targets[0] > len)
    return uniq_targets;

  size_t prev = targets[0];
  uniq_targets.push_back(prev);
  
  for (size_t i = 1; i < targets.size(); i++) {
    if (targets[i] > len) break;
    if (targets[i] != prev) {
      prev = targets[i];
      uniq_targets.push_back(prev);
    }
  }

  return uniq_targets;
}

vec_pair filter_ranges(vec_pair ranges, size_t len) {
  if (!ranges.size())
    return vec_pair(0);
  
  vec_pair uniq_ranges;
  std::sort(ranges.begin(), ranges.end());

  ranges.push_back(std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));

  if (ranges[0].first > len)
    return uniq_ranges;

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
  if (!df.nrows())
    return vec(0);
  
  Rcpp::StringVector df_names = df[DF_NAME_INDEX];

  vec targets;
  Rcpp::NumericVector df_targets = df[DF_TARGET_INDEX];

  for (size_t i = 0; i < (size_t) df_names.size(); i++)
    if (df_names[i] == name)
      targets.push_back(df_targets[i]);

  return filter_targets(targets, len);
}

vec_pair subsample_ranges(std::string name, size_t len, Rcpp::DataFrame df) {
  if (!df.nrows())
      return vec_pair(0);

  Rcpp::StringVector df_names = df[DF_NAME_INDEX];
  
  vec_pair ranges;
  Rcpp::NumericVector df_starts = df[DF_START_INDEX], df_stops = df[DF_STOP_INDEX];
    
  for (size_t i = 0; i < (size_t) df_names.size(); i++)
    if (df_names[i] == name)
      ranges.push_back(std::make_pair(df_starts[i], df_stops[i]));

  return filter_ranges(ranges, len);
}

vec_pair_value check_ranges(vec_pair ranges, size_t n, size_t size, std::string chrom, size_t tiling, std::unordered_set<size_t> used_baits, std::string type) {
  vec_pair_value output_ranges;
  
  for (auto r : ranges) {
    size_t range = r.second - r.first + 1;
    size_t count_used = 0;
    
    for (auto u : used_baits) // check if used_bait start is in range
      if (u >= r.first && u <= (r.second-(size-1)))
	count_used += 1;

    size_t max_baits = 0;
    if (range > count_used+(size-1)) {
      max_baits = range-count_used-(size-1);
    }
    
    if (max_baits >= tiling)
      output_ranges.push_back(std::make_pair(std::make_pair(r.first, r.second), max_baits));
  }

  if (output_ranges.size()) { // only perform checks if there are valid ranges
    if (output_ranges.size() < ranges.size())
      Rcpp::Rcout << "Warning: " << (ranges.size() - output_ranges.size()) << " sub-ranges on sequence " << chrom << " are too small to fit the desired number of baits and will be excluded." << std::endl;

    if (n / tiling < output_ranges.size()) { // if there are too many ranges, randomly select needed
      Rcpp::Rcout << "Warning: The desired n/tiling combination is not high enough to produce baits in all valid ranges in sequence " << chrom << ". Choosing a random subset of ranges." << std::endl;
    
      size_t max_ranges = std::floor((double)n / tiling);
      std::random_shuffle(output_ranges.begin(), output_ranges.end());
    
      output_ranges.resize(max_ranges);
      sort(output_ranges.begin(), output_ranges.end());
    }
  } else {
   Rcpp::Rcout << "Warning: All " << type << " ranges in sequence " << chrom << " are too small to fit the desired number of baits per range. Skipping..." << std::endl;
  }
  
  return output_ranges;
}

vec check_n(vec_pair_value ranges, size_t n, size_t tiling, std::string chrom, std::string type) {
  vec n_per_range;
  size_t sum_max_baits = 0;
  
  for (auto r : ranges)
    sum_max_baits += r.second;
  
  if (sum_max_baits <= n) {
    if (sum_max_baits < n) // not enough space for all baits
      Rcpp::Rcout << "Warning: The maximum possible number of unique " << type << " baits (" << sum_max_baits << ") for sequence " << chrom << " is lower than the desired n (" <<  n << ")." << std::endl;

    for (auto r : ranges)
      n_per_range.push_back(r.second);
  } else { // there are many bait spots available, distribute proportionally
    size_t n_after_tiling = n;
    
    for (size_t i = 0; i < ranges.size(); i++) {
      n_per_range.push_back(tiling);
      n_after_tiling -= tiling;
      ranges[i].second -= tiling;
    }
    size_t remaining_n = n_after_tiling;

    sum_max_baits = 0;
    for (auto r : ranges)
      sum_max_baits += r.second;
    
    for (size_t i = 0; i < ranges.size(); i++) {
      size_t n_baits_r = std::floor(((double) ranges[i].second/sum_max_baits)*n_after_tiling);
      remaining_n -= n_baits_r;
      n_per_range[i] += n_baits_r;
    }
    
    if (remaining_n > 0) {
      vec seq;
      for (size_t i = 0; i < n_per_range.size(); i++)
	if ((ranges[i].second+tiling) > n_per_range[i])
	  seq.push_back(i);

      std::random_shuffle(seq.begin(), seq.end());
      for (size_t i = 0; i < remaining_n; i++)
	n_per_range[seq[i]] += 1;
    }

    size_t check_n_per_range = 0;
    for (auto r : n_per_range)
      check_n_per_range += r;
   
    if (check_n_per_range != n)
      Rcpp::stop("Unhadled error: contact developers :)");
  }
  
  return n_per_range;
}

vec_pair_string get_bait_positions(vec_pair_value ranges, size_t size, vec n_per_range, std::unordered_set<size_t> used_baits, std::string type) {
  vec_pair bait_positions;
  for (size_t i = 0; i < ranges.size(); i++) {
    std::set<size_t> rands;
    while (rands.size() < n_per_range[i]) {
      size_t max = ranges[i].first.second - (size - 1);
      size_t min = ranges[i].first.first;
      size_t rand_num = std::ceil(R::runif(min, max));
      if (!used_baits.count(rand_num))
	rands.insert(rand_num);
    }

    for (size_t myrand : rands)
      bait_positions.push_back(std::make_pair(myrand, myrand + (size - 1)));
  }

  vec_pair_string bait_pos_type;
  for (auto bait_pos : bait_positions)
    bait_pos_type.push_back(std::make_pair(bait_pos, type));

  return bait_pos_type;
}

Rcpp::DataFrame buildSamplesRdf(std::vector<SampleBait> &baits) {
  std::vector<std::string> names, types;
  std::vector<size_t> starts, stops;

  for (SampleBait b : baits) {
    names.push_back(b.chrom);
    types.push_back(b.type);
    starts.push_back(b.start);
    stops.push_back(b.stop);
  }

  Rcpp::DataFrame df = Rcpp::DataFrame::create(Rcpp::_["ChromName"] = names,
					       Rcpp::_["Type"]      = types,
					       Rcpp::_["Start"]     = starts,
					       Rcpp::_["Stop"]      = stops);
  return df;
}

// [[Rcpp::export]]
Rcpp::DataFrame sampleBaits(Rcpp::DataFrame chrom_lens,
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

  std::vector<SampleBait> all_baits;
  for (size_t i = 0; i < (size_t) df_names.size(); i++) { // for each chrom/len in chrom_lens
    vec targets_subsample = subsample_targets(Rcpp::as<std::string>(df_names[i]), df_lens[i], targets);
    vec_pair exclusions_subsample = subsample_ranges(Rcpp::as<std::string>(df_names[i]), df_lens[i], exclusions);
    vec_pair regions_subsample = subsample_ranges(Rcpp::as<std::string>(df_names[i]), df_lens[i], regions);

    std::unordered_set<size_t> used_baits;
    size_t n_regions = 0, n_targets = 0;
    vec_pair regions, targets, random;
    vec_pair_string region_bait_positions, target_bait_positions, random_bait_positions;
    
    Rcpp::Rcout << "M: Sampling sequence " << df_names[i] << "." << std::endl;

    // region baits
    if (regions_prop > 0) {
      n_regions = std::floor(n * regions_prop);
      if (regions_subsample.size() > 0) {
	regions = trim_ranges(regions_subsample, exclusions_subsample);

	vec_pair_value valid_regions = check_ranges(regions, n_regions, size, (std::string) df_names[i], regions_tiling, used_baits, "region");

	if (valid_regions.size()) {
	  vec n_per_range = check_n(valid_regions, n_regions, regions_tiling, (std::string) df_names[i], "region");
	
	  region_bait_positions = get_bait_positions(valid_regions, size, n_per_range, used_baits, "region");

	  for (auto t : region_bait_positions) {
	    all_baits.push_back(SampleBait {(std::string) df_names[i], "region", t.first.first, t.first.second});
	    used_baits.insert(t.first.first);
	  }
	  n_regions = region_bait_positions.size();
	} else {
	  n_regions = 0;
	}
      } else {
	n_regions = 0;
    	Rcpp::Rcout << "M: No regions found for sequence " << df_names[i] << "." << std::endl;
      }
    }

    // target baits
    if (targets_prop > 0) {
      n_targets = std::floor(n * targets_prop);
      if (targets_subsample.size() > 0) {
	vec_pair target_ranges;
	for (auto t : targets_subsample) // create ranges from targets
	  target_ranges.push_back(std::make_pair(std::max(t-(size-1), (size_t)1), std::min(t+size-1, (size_t)df_lens[i])));

	targets = trim_ranges(target_ranges, exclusions_subsample);

	vec_pair_value valid_targets = check_ranges(targets, n_targets, size, (std::string) df_names[i], targets_tiling, used_baits, "target");

	if (valid_targets.size()) {
	  vec n_per_range = check_n(valid_targets, n_targets, targets_tiling, (std::string) df_names[i], "target");

	  target_bait_positions = get_bait_positions(valid_targets, size, n_per_range, used_baits, "target");

	  for (auto t : target_bait_positions) {
	    all_baits.push_back(SampleBait {(std::string) df_names[i], "target", t.first.first, t.first.second});
	    used_baits.insert(t.first.first);
	  }

	  n_targets = target_bait_positions.size();
	} else {
	  n_targets = 0;
	}
      } else {
	n_targets = 0;
    	Rcpp::Rcout << "M: No targets found for sequence " << df_names[i] << "." << std::endl;
      }
    }

    // random baits
    if (n_regions + n_targets < n) {
      size_t n_random = n - (n_regions + n_targets);

      vec_pair random_ranges; random_ranges.push_back(std::make_pair(1, df_lens[i]));

      for (auto t : trim_ranges(random_ranges, exclusions_subsample))
	random.push_back(std::make_pair(t.first, t.second));
      
      vec_pair_value valid_random = check_ranges(random, n_random, size, (std::string) df_names[i], 1, used_baits, "random");

      if (valid_random.size()) {
	vec n_per_range = check_n(valid_random, n_random, 1, (std::string) df_names[i], "random");
	
	random_bait_positions = get_bait_positions(valid_random, size, n_per_range, used_baits, "random");

	for (auto t : random_bait_positions)
	  all_baits.push_back(SampleBait {(std::string) df_names[i], "random", t.first.first, t.first.second});
      }
    }
  }
  
  return buildSamplesRdf(all_baits);
}

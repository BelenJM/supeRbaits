// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <unordered_map>

/* 
Change <unordered_map> to <map> if keeping ordering of keys in output is necessary.
Accordingly, change unorderd_map<size_t, size_t> to map<size_t, size_t>.
*/

// [[Rcpp::export]]
Rcpp::DataFrame importBlastResults(std::string path) {
  std::ifstream input(path.c_str());
  
  if(!input.good()) {
    Rcpp::stop("Error opening file '%s'. Exiting...", path);
  }

  std::unordered_map<size_t, size_t> map;
  std::string line;
  while(getline(input, line).good()) {
    size_t bait_no = 0;
    for (size_t i = 0; line[i] != ',' && i < (size_t) line.size(); i++) {
      bait_no = 10*bait_no + (line[i] - '0');
    }

    auto elem = map.find(bait_no);
    if (elem == map.end()) { // key not in map
      map.insert(std::make_pair(bait_no, 1));
    } else {
      map[elem->first] = elem->second + 1;
    }    
  }

  std::vector<size_t> bait_no;
  std::vector<size_t> no_matches;
  for (auto it = map.begin(); it != map.end(); it++) {
    bait_no.push_back(it->first);
    no_matches.push_back(it->second);
  }
  
  Rcpp::DataFrame df =  Rcpp::DataFrame::create(Rcpp::_["bait_no"] = bait_no,
						Rcpp::_["no_matches"] = no_matches);

  return df;
}

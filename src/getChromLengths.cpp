// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <vector>
#include <fstream>
#include <iostream>

// [[Rcpp::export]]
Rcpp::DataFrame getChromLengths(std::string path) {
  std::ifstream input(path.c_str());
  
  if(!input.good()) {
    Rcpp::stop("Error opening file '%s'. Exiting...", path);
  }

  size_t size;
  std::string line, name;
  std::vector<std::string> names;
  std::vector<size_t> sizes;
  while (getline(input, line).good()) {
    if (line.empty() || line[0] == '>') {
      if(!name.empty()) {
	names.push_back(name);
	sizes.push_back(size);
	name.clear();
      }
      if(!line.empty()) {
	for (size_t i=1; line[i] != ' ' && i < line.size(); i++) {
	  name += line[i];
	}
      }
      size = 0;
    } else if(!name.empty()) {
      if(line.find(' ') != std::string::npos){
	name.clear();
	size = 0;
      } else {
	size += line.size();
      }
    }
  }
  if(!name.empty()) {
    names.push_back(name);
    sizes.push_back(size);
  }

  Rcpp::DataFrame df = Rcpp::DataFrame::create(Rcpp::_["name"] = names, Rcpp::_["size"] = sizes);
  
  return df;
}

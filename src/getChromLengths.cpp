#include <Rcpp.h>
#include <vector>
#include <fstream>
#include <iostream>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
DataFrame getChromLengths(string path) {
  ifstream input(path);
  
  if(!input.good()) {
    stop("Error opening file '%s'. Exiting...", path);
  }

  int size;
  string line, name;
  vector<string> names;
  vector<int> sizes;

  while (getline(input, line).good()) {
    if (line.empty() || line[0] == '>') { // identifier marker
      if(!name.empty()) {
	names.push_back(name);
	sizes.push_back(size);
	name.clear();
      }
      if(!line.empty()) {
	for (int i=1; line[i] != ' ' && i < line.length(); i++) {
	  name += line[i];
	}
      }
      size = 0;
    } else if(!name.empty()) {
      if(line.find(' ') != string::npos){ // invalid: sequence with spaces
	name.clear();
	size = 0;
      } else {
	size += line.length();
      }
    }
  }
  if( !name.empty() ) {
    names.push_back(name);
    sizes.push_back(size);
  }

  DataFrame df = DataFrame::create(_["name"] = names, _["size"] = sizes);
  
  return df;
}

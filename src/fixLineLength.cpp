// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <string>
#include <fstream>
#include <iostream>

const int DEFAULT_LENGTH = 80;

// [[Rcpp::export]]
void fixLineLengths(std::string fin_path, std::string fout_path) {
  std::ifstream fin(fin_path.c_str());
  std::ofstream fout(fout_path.c_str());

  if(!fin.good()) {
    Rcpp::stop("Error opening file '%s'. Exiting...", fin_path);
  }

  if(!fout.good()) {
    Rcpp::stop("Error opening file '%s'. Exiting...", fout_path);
  }

  std::string line = "", seq = "";
  while(getline(fin, line).good()) {
    if (line[0] == '>') {
      if (seq.size() > 0) {
	for (size_t i = 1; i < seq.size()+1; i++) {
	  if (!(i % DEFAULT_LENGTH)) {
	    fout << seq[i-1] << "\n";
	  } else {
	    fout << seq[i-1];
	  }
	}
	if (seq.size() % DEFAULT_LENGTH) {
	  fout << std::endl;
	}
      }
      seq = "";
      fout << line << "\n";
    }
    else {
      seq += line;
    }
  }
  
  if (seq.size() > 0) {
    for (size_t i = 1; i < seq.size()+1; i++) {
      if (!(i % DEFAULT_LENGTH)) {
	fout << seq[i-1] << "\n";
      } else {
	fout << seq[i-1];
      }
    }
    if (seq.size() % DEFAULT_LENGTH) {
      fout << std::endl;
    }
  }

  fin.close();
  fout.close();
}

// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

// [[Rcpp::export]]
void removeNLs(std::string fin_path, std::string fout_path) {
  std::ifstream fin(fin_path.c_str());
  std::ofstream fout(fout_path.c_str());

  if(!fin.good()) {
    Rcpp::stop("Error opening file '%s'. Exiting...", fin_path);
  }

  if(!fout.good()) {
    Rcpp::stop("Error opening file '%s'. Exiting...", fout_path);
  }

  bool first_line = true;
  std::string line;
  while(getline(fin, line).good()) {
    if (line[0] == '>') {
      if (first_line) {
	first_line = false;
      } else {
	fout << "\n";
      }
      fout << line << "\n";
    } else {
      fout << line;
    }
  }
  fout << "\n";
    
  fin.close();
  fout.close();
}

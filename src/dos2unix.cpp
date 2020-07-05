// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

// [[Rcpp::export]]
void dos2unix(std::string fin_path, std::string fout_path) {
  // only for converting CRLF to LF
  ifstream fin(fin_path.c_str());
  ofstream fout(fout_path.c_str());

  if(!fin.good()) {
    Rcpp::stop("Error opening file '%s'. Exiting...", fin_path);
  }

  if(!fout.good()) {
    Rcpp::stop("Error opening file '%s'. Exiting...", fout_path);
  }

  char c;
  while (fin.get(c))
    if (c != '\r')
      fout << c;

  fin.close();
  fout.close();
}

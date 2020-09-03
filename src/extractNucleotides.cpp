// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <string>
#include <fstream>
#include <iostream>

// [[Rcpp::export]]
std::string extractNucleotides(std::string db_path, std::string chrom_name, size_t start, size_t stop) {
  std::ifstream db(db_path.c_str());
  
  if (!db.good()) {
    Rcpp::stop("Error opening file '%s'. Exiting...", db_path);
  }

  std::string line = "", name = "", seq = "";
  bool seq_found = false;
  while(getline(db, line).good()) {
    if (line[0] == '>') {
      name = "";
      if (seq_found) {
	break;
      }
      for (size_t i=1; line[i] != ' ' && i < line.size(); i++) {
	name += line[i];
      }
      if (name == chrom_name) {
	seq_found = true;
      }
    }
    else {
      if (seq_found) {
	seq += line;
      }
    }
  }

  std::string subseq = "";
  if (seq_found) {
    try {
      subseq = seq.substr(start-1, stop-start+1);
    } catch (const std::out_of_range& oor) {
      Rcpp::stop("The selected nucleotide subsequence (%d, %d) starts after the selected chromosome ends.", start, stop);
    }
  } else {
    Rcpp::stop("Chromosome %s was not found.", chrom_name);
  }

  if (subseq.size() < stop-start+1) {
    Rcpp::warning("The selected nucleotide subsequence (%d, %d) ends after the selected chromosome ends. Ending subsequence at the end of the selected chromosome...", start, stop);
  }

  return subseq;
}

// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <vector>
#include <fstream>
#include <iostream>

const int DF_BAIT_NO = 0,
          DF_BAIT_SEQ = 1;

// [[Rcpp::export]]
void writeFasta(Rcpp::DataFrame baits, std::string fasta_path) {
  Rcpp::NumericVector df_bait_no = baits[DF_BAIT_NO];
  Rcpp::StringVector  df_bait_seq = baits[DF_BAIT_SEQ];

  std::ofstream fout(fasta_path.c_str(), std::ios::binary);

  if(!fout.good()) {
    Rcpp::stop("Error opening file '%s'. Exiting...", fasta_path);
  }

  for (size_t i = 0; i < (size_t) df_bait_no.size(); i++) {
    fout << ">" << df_bait_no[i] << std::endl;
    fout << df_bait_seq[i] << std::endl;
  }

  fout.close();
}

// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <unordered_map>

/*
test_df <- data.frame("ChromName" = c("CM003279", "CM003279", "CMF", "CMF"),
 		      "Start" = c(53, 56, 2, 1),
 		      "Stop" = c(72, 75, 21, 20),
 		      "Type" = c("region", "region", "random", "random"))
*/

const int DF_NO_COLS = 4,
          DF_NAME_INDEX = 0,
          DF_START_INDEX = 1,
          DF_STOP_INDEX = 2,
          DF_TYPE_INDEX = 3;

struct Bait {
  size_t no;
  std::string name;
  std::string type;
  std::string seq;
  size_t start;
  size_t stop;
  size_t no_A;
  size_t no_T;
  size_t no_G;
  size_t no_C;
  size_t no_UNK;
};
  

Bait getBait(std::ifstream &db,
	     std::unordered_map<std::string,
	     size_t> &map,
	     size_t no,
	     std::string name,
	     std::string type,
	     size_t start,
	     size_t stop) {
  size_t fpos = map[name];
  if (fpos == 0) {
    Rcpp::stop("Error: '%s' not found. Exiting...", name);
  }

  db.clear(); db.seekg(fpos+(start-1), std::ios::beg);
  
  char c;
  std::string seq = "";
  size_t no_A = 0, no_T = 0, no_C = 0, no_G = 0, no_UNK = 0;
  for (size_t i = 0; i <= stop-start && db.get(c); i++) {
    if (c == '>') {
      Rcpp::stop("Error: sequence stop overflow for '%s'. Exiting...", name);
    }
    if (c == '\n') {
      i--;
    } else {
      switch(toupper(c)) {
      case 'A': no_A++; break;
      case 'T': no_T++; break;
      case 'C': no_C++; break;
      case 'G': no_G++; break;
      default: no_UNK++; break;
      }
      seq += c;
    }
  }
  if (db.eof()) {
    Rcpp::stop("Error: sequence stop overflow for '%s'. Exiting...", name);
  }

  Bait b {no, name, type, seq, start, stop, no_A, no_T, no_G, no_C, no_UNK};
  return b;
}

std::unordered_map<std::string, size_t> preProcDB(std::ifstream &db) {
  std::unordered_map<std::string, size_t> map;
  
  std::string line, name;
  while(getline(db, line).good()) {
    if(line[0] == '>') {
      for (size_t i=1; line[i] != ' ' && i < line.size(); i++) {
	name += line[i];
      }
      map[name] = db.tellg();
      name.clear();
    }
  }

  return map;
}

Rcpp::DataFrame buildRdf(std::vector<Bait> &baits) {
  std::vector<std::string> names, types, seqs;
  std::vector<size_t> nos, starts, stops, sizes,
    no_As, no_Ts, no_Gs, no_Cs, no_UNKs, no_ATs, no_GCs;

  for (Bait b : baits) {
    nos.push_back(b.no);
    names.push_back(b.name);
    types.push_back(b.type);
    seqs.push_back(b.seq);
    starts.push_back(b.start);
    stops.push_back(b.stop);
    sizes.push_back(b.stop-b.start+1);
    no_As.push_back(b.no_A);
    no_Ts.push_back(b.no_T);
    no_Gs.push_back(b.no_G);
    no_Cs.push_back(b.no_C);
    no_UNKs.push_back(b.no_UNK);
    no_ATs.push_back(b.no_A + b.no_T);
    no_GCs.push_back(b.no_G + b.no_C);
  }

  Rcpp::DataFrame df = Rcpp::DataFrame::create(Rcpp::_["bait_no"]         = nos,
					       Rcpp::_["bait_chrom_name"] = names,
					       Rcpp::_["bait_type"]       = types,
					       Rcpp::_["bait_seq"]        = seqs,
					       Rcpp::_["bait_start"]      = starts,
					       Rcpp::_["bait_stop"]       = stops,
					       Rcpp::_["bait_seq_size"]   = sizes,
					       Rcpp::_["no_A"]            = no_As,
					       Rcpp::_["no_T"]            = no_Ts,
					       Rcpp::_["no_G"]            = no_Gs,
					       Rcpp::_["no_C"]            = no_Cs,
					       Rcpp::_["no_UNK"]          = no_UNKs,
					       Rcpp::_["no_AT"]           = no_ATs,
					       Rcpp::_["no_GC"]           = no_GCs);
  return df;
}
			 
// [[Rcpp::export]]
Rcpp::DataFrame getBaits(std::string db_path, Rcpp::DataFrame df) {
  std::vector<Bait> baits;
  std::ifstream db(db_path.c_str());
  
  if(!db.good()) {
    Rcpp::stop("Error opening file '%s'. Exiting...", db_path);
  }
  
  std::unordered_map<std::string, size_t> map = preProcDB(db);

  if (df.size() != DF_NO_COLS) {
    Rcpp::stop("Error: data frame has the wrong number of columns. Exiting...");
  }

  Rcpp::StringVector df_names = df[DF_NAME_INDEX], df_types = df[DF_TYPE_INDEX];
  Rcpp::NumericVector df_starts = df[DF_START_INDEX], df_stops  = df[DF_STOP_INDEX];

  for (size_t i = 0; i < df_names.size(); i++) {
    baits.push_back(getBait(db,
			    map,
			    i,
			    Rcpp::as<std::string>(df_names[i]),
			    Rcpp::as<std::string>(df_types[i]),
			    df_starts[i],
			    df_stops[i]));
  }

  return buildRdf(baits);
}

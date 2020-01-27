// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <unordered_map>

using namespace Rcpp;
using namespace std;

/*
test_df <- data.frame("ChromName" = c("CM003279", "CM003279", "CMF", "CMF"),
 		      "Start" = c(53, 56, 2, 1),
 		      "Stop" = c(72, 75, 21, 20),
 		      "Type" = c("region", "region", "random", "random"))
*/

const int NO_COLUMNS = 4;

enum {BAIT_CHROM_NAME = 0,
      BAIT_START      = 1,
      BAIT_STOP       = 2,
      BAIT_TYPE       = 3};

// [[Rcpp::export]]
DataFrame getBaits(std::string gen_path, DataFrame bait_df) {
  ifstream gen_input(gen_path.c_str());
  if(!gen_input.good()) {
    stop("Error opening file '%s'. Exiting...", gen_path);
  }
  unordered_map<string, string> memo; // chromosome cache
  vector<string> bait_chrom_names, bait_types, bait_seqs;
  vector<long long> bait_nos, bait_starts, bait_stops, bait_seq_sizes,
                    no_As, no_Ts, no_Gs, no_Cs, no_UNKs, no_ATs, no_GCs;

  if (bait_df.size() != NO_COLUMNS) {
      stop("Error parsing baits file. Wrong number of columns...");
  }
  StringVector bait_df_chrom_names = bait_df[BAIT_CHROM_NAME];
  NumericVector bait_df_starts = bait_df[BAIT_START];
  NumericVector bait_df_stops = bait_df[BAIT_STOP];
  StringVector bait_df_types = bait_df[BAIT_TYPE];

  long long index = 0;
  while (index < bait_df_chrom_names.size()) {
    string bait_chrom_name = as<string>(bait_df_chrom_names[index]);
    long long bait_start   = bait_df_starts[index];
    long long bait_stop    = bait_df_stops[index];
    string bait_type       = as<string>(bait_df_types[index]);
      
    string gen_chrom_name, gen_chrom;
    unordered_map<string, string>::iterator it = memo.find(bait_chrom_name); 
    if (it != memo.end()) { // cache hit
      gen_chrom_name = it->first;
      gen_chrom      = it->second;
    } else { // cache miss
      string gen_line;
      while (getline(gen_input, gen_line).good()) {
	if (gen_line.empty() || gen_line[0] == '>') { // identifier marker
	  if (!gen_chrom_name.empty() && gen_chrom_name == bait_chrom_name) {
	    memo[gen_chrom_name] = gen_chrom;
	    break;
	  }
	  if(!gen_line.empty()) {
	    gen_chrom_name.clear();
	    gen_chrom.clear();
	    for (unsigned long long i=1; gen_line[i] != ' ' && i < gen_line.length(); i++) {
	      gen_chrom_name += gen_line[i];
	    }
	  }
	} else if(!gen_chrom_name.empty()) {
	  if(gen_line.find(' ') != string::npos){ // invalid: sequence with spaces
	    gen_chrom_name.clear();
	    gen_chrom.clear();
	  } else {
	    gen_chrom += gen_line;
	  }
	}
      }
      if (!gen_chrom_name.empty() && gen_chrom_name == bait_chrom_name) {
	memo[gen_chrom_name] = gen_chrom;
      }
      gen_input.clear(); gen_input.seekg(0, ios::beg);
    }

    if (!gen_chrom_name.empty() && gen_chrom_name == bait_chrom_name) { // chromosome found
      string bait_seq = gen_chrom.substr(bait_start-1, bait_stop-(bait_start-1));
      long long no_A = 0, no_T = 0, no_G = 0, no_C = 0, no_UNK = 0;
      for (char c : bait_seq) {
	c = toupper(c);
	switch(c) {
	case 'A':
	  no_A++; break;
	case 'T':
	  no_T++; break;
	case 'G':
	  no_G++; break;
	case 'C':
	  no_C++; break;
	default:
	  no_UNK++;
	  break;
	}
      }
      bait_nos.push_back(index);
      bait_chrom_names.push_back(bait_chrom_name);
      bait_types.push_back(bait_type);
      bait_starts.push_back(bait_start);
      bait_stops.push_back(bait_stop);
      bait_seqs.push_back(bait_seq);
      bait_seq_sizes.push_back(bait_seq.size());
      no_As.push_back(no_A);
      no_Ts.push_back(no_T);
      no_Gs.push_back(no_G);
      no_Cs.push_back(no_C);
      no_UNKs.push_back(no_UNK);
      no_ATs.push_back(no_A + no_T);
      no_GCs.push_back(no_G + no_C);
    } else { // chromosome not found
      Rcout << "Could not find bait chromosome " << bait_chrom_name << " in genome file " << gen_path << ". Skipping..." << endl;
    }
    index++;
  }

  DataFrame df = DataFrame::create(_["bait_no"]         = bait_nos,
				   _["bait_chrom_name"] = bait_chrom_names,
				   _["bait_type"]       = bait_types,
				   _["bait_start"]      = bait_starts,
				   _["bait_stop"]       = bait_stops,
				   _["bait_seq"]        = bait_seqs,
				   _["bait_seq_size"]   = bait_seq_sizes,
				   _["no_A"]            = no_As,
				   _["no_T"]            = no_Ts,
				   _["no_G"]            = no_Gs,
				   _["no_C"]            = no_Cs,
				   _["no_UNK"]          = no_UNKs,
				   _["no_AT"]           = no_ATs,
				   _["no_GC"]           = no_GCs);
				   
  return df;
}

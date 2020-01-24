// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <unordered_map>

using namespace Rcpp;
using namespace std;

const int NO_COLUMNS = 4;

enum {BAIT_CHROM_NAME = 0,
      BAIT_START      = 1,
      BAIT_STOP       = 2,
      BAIT_TYPE       = 3};

// [[Rcpp::export]]
DataFrame getBaits(std::string gen_path, std::string bait_path) {
  ifstream gen_input(gen_path.c_str()),
           bait_input(bait_path.c_str());
  if(!gen_input.good()) {
    stop("Error opening file '%s'. Exiting...", gen_path);
  }
  if(!bait_input.good()) {
    stop("Error opening file '%s'. Exiting...", bait_path);
  }
  unordered_map<string, string> memo; // chromosome cache
  vector<string> bait_chrom_names, bait_types, bait_seqs;
  vector<long long> bait_nos, bait_starts, bait_stops, bait_seq_sizes,
                    no_As, no_Ts, no_Gs, no_Cs, no_UNKs, no_ATs, no_GCs;

  string bait_line;
  long long index = 0;
  while (getline(bait_input, bait_line).good()) {
    if (index == 0) { // skip header
      index++;
      continue;
    }
    istringstream iss(bait_line);
    std::vector<std::string> columns{std::istream_iterator<std::string>(iss), 
				     std::istream_iterator<std::string>()};
    if (columns.size() != NO_COLUMNS) {
      stop("Error parsing baits file. Wrong number of columns...");
    }
    string bait_chrom_name = columns[BAIT_CHROM_NAME];
    long long bait_start   = stoll(columns[BAIT_START]);
    long long bait_stop    = stoll(columns[BAIT_STOP]);
    string bait_type       = columns[BAIT_TYPE];

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

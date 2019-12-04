#include <iostream>
#include <fstream>
#include <vector>


//////////////////////////////////////////////////////
// TODO: remove this comment			    //
// 						    //
// Compile this program like this:		    //
// g++ getChromLengths.cpp -o getChromLengths	    //
// this create an executable called getChromLengths //
// 						    //
// Run the executable like this:		    //
// ./getChromLengths "FASTA file"		    //
//////////////////////////////////////////////////////

using namespace std;

vector<pair<string, int>> getChromLengths(string pathToFasta) {
  ifstream input(pathToFasta);
  
  if(!input.good()) {
    cerr << "Error opening '" << pathToFasta << "'. Exiting..." << endl;
    exit(-1);
  }

  int size;
  string line, name;
  vector<pair<string, int>> buff;
  
  while (getline( input, line ).good()) {
    // Identifier marker
    if (line.empty() || line[0] == '>') { 
      if(!name.empty()) {
	buff.push_back(make_pair(name, size));
	name.clear();
      }
      if(!line.empty()) {
	for (int i=1; line[i] != ' ' && i < line.length(); i++) {
	  name += line[i];
	}
      }
      size = 0;
    } else if(!name.empty()) {
      // Invalid sequence--no spaces allowed
      if(line.find(' ') != string::npos){ 
	name.clear();
	size = 0;
      } else {
	size += line.length();
      }
    }
  }
  if( !name.empty() ) {
    buff.push_back(make_pair(name, size));
  }

  return buff;
}

int main(int argc, char **argv) {
  ios_base::sync_with_stdio(false);
  
  if( argc <= 1 ) {
    cerr << "Usage: " << argv[0] << " FASTA file " << endl;
    exit(-1);
  }

  auto lens = getChromLengths(argv[1]);
  
  for (auto e : lens) {
    cout << e.first << " " << e.second << endl;
  }

  exit(0);
}



#include <iostream>
#include <iostream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdarg>
#include <sstream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <vector>

using namespace std;

unsigned long int getCurrentLine(std::istream& is);

int main(int argc, char**argv) {
  istream& in = cin;
  string s;
  while(in >> s) {
   	if (s=="bad") {
	  cout << "Current position: " << in.fail() << ", " <<  in.tellg() << endl;
      cout << "Error! Bad on line " << getCurrentLine(in) << endl;
 	  break;
   	}
  }
}

unsigned long int getCurrentLine(std::istream& is) {
    unsigned long int line = 1;
    is.clear(); // need to clear error bits otherwise tellg returns -1.
    std::streampos originalPos = is.tellg();
    if (originalPos < 0) return 0;
    is.seekg(0);
    char c;
	// TODO: Update to using sentry and rdbuf to improve efficiency.
	// TODO: Check for old-style Mac line endings ('\r')
    while ((is.tellg() < originalPos) && is.get(c))
        if (c == '\n') ++line;
    
    // ok but if we are AT a newline, subtract one.
    if (c == '\n') --line;

    return line;
}

/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Authors: Chris Connett and Andrew Yohn, 2006
 * Contributor: Josh Keegan, 2006
 */

#include "configfile.h"

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

using namespace std;

ConfigFile::ConfigFile() {
}

ConfigFile::ConfigFile(const string &filename)
  : valid(true) {
  
  ifstream config(filename.c_str());

  if(!config) {
    
    // Could not open config file.
    valid = false;
    cerr << "Could not open config file named '" << filename << "'." << endl;

  } else {

    // Config file was successfully opened.
    int lineCounter = 0;
    string line;
    while(getline(config, line)) {

      lineCounter++;

      // Skip blank lines and comments.
      istringstream temp(line);
      string token;
      if (temp >> token && token[0] != '#') {
      
        // Try to find a '=' character.
        string::size_type equalSignIndex = line.find('=');
        if(equalSignIndex == string::npos) {

          // No '=' character found, syntax error.
          valid = false;
          cerr << "Skipping line " << lineCounter << " in config file due to "
               << "syntax error: expected '='." << endl;
        
        } else {

          // Get the key and value from the line (separated by '=').
          string unparsedkey = line.substr(0, equalSignIndex);
          string unparsedvalue = line.substr(equalSignIndex + 1);
          istringstream keystream(unparsedkey);
          istringstream valuestream(unparsedvalue);
          string key;
          string value;
          if(keystream >> key && valuestream >> value &&
             !(keystream >> token) && !(valuestream >> token)) {

            // Got the key/value pair
	          string lowerKey = key;
	          int length = lowerKey.size();
	          for( int i = 1; i <= length; i++ ) {
		          lowerKey[i-1] = tolower( lowerKey[i-1] );
	          }
	          options[lowerKey] = value;

          } else {
            // Invalid key/value pair
            valid = false;
            cerr << "Skipping line " << lineCounter 
                 << " in config file due to syntax error: invalid key/value pair." << endl 
                 << "'" << line << "'" << endl;
            
          }
        }
      }
    }
    config.close();
  }
}

bool ConfigFile::isValid() const {
  return valid;
}

bool ConfigFile::contains(const string &option) const {
	string lowerOption = option;
	int length = lowerOption.size();
	for( int i = 1; i <= length; i++ ) {
		lowerOption[i-1] = tolower( lowerOption[i-1] );
	}

  return options.find(lowerOption) != options.end();
}

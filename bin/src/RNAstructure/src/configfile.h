/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Authors: Chris Connett and Andrew Yohn, 2006
 * Contributor: Josh Keegan, 2006
 */

#ifndef CONFIGFILE_H
#define CONFIGFILE_H

#include <ctype.h>
#include <sstream>
#include <map>
#include <string>

#ifndef NDEBUG
#include <iostream>
#endif

class ConfigFile {
private:
  std::map<std::string, std::string> options;
  bool valid;
  ConfigFile(); // Prevent construction with no filename

public:
  ConfigFile(const std::string &filename);

  bool isValid() const;
  bool contains(const std::string &option) const;
  
  template<typename T>
  T getOption(const std::string &option) const;

  template<typename T>
  bool getOptionByRef(const std::string &option, T &setting) const;
};

//! Get the value of the setting. 
//! If the setting was NOT found, the return value is the default value of the specified type. (e.g. 0 for ints, "" for strings etc).
template<typename T>
T ConfigFile::getOption(const std::string &option) const {
  T setting;
  getOptionByRef<T>(option, setting);
#ifndef NDEBUG
    std::cout << "Set " << option  << " = " << setting << std::endl;
#endif
  return setting;
}

//! Set the value of setting (passed in by reference) but ONLY if the setting was found. 
//! otherwise the setting retains its original value.
//! /return True if the setting was found. False otherwise.
template<typename T>
bool ConfigFile::getOptionByRef(const std::string &option, T &setting) const {
  std::string lowerOption = option;
  int length = lowerOption.size();
  for( int i = 1; i <= length; i++ )
	  lowerOption[i-1] = tolower( lowerOption[i-1] );
  std::map<std::string, std::string>::const_iterator i = options.find(lowerOption);
  if (i != options.end()) {
    std::istringstream temp((*i).second);
    temp >> setting;
    return true;
  }
  return false;
}

#endif

/*
 * An implementation file for a program that parses Unix-like command lines to
 * determine the data being given for input.
 *
 * (c) 2008  Mathews Lab, University of Rochester Medical Center
 * Redone in 2012.
 * Written by Jessica S. Reuter
 */

#include "ParseCommandLine.h"
#include "common_utils.h"
#include "version.h" // Defines the version and build date for the RNAstructure package.
#include <ctype.h>

#define PROG_ERR(MSG) const char * const error_message = "Programming error: " MSG ; cerr << error_message; throw error_message;

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
ParseCommandLine::ParseCommandLine( string name ) {
	// Note: the 'version' field has been replaced by a pre-processor definition in version.h
	
	// Initialize the error flag and the specialized usage flags.
	error = false;
	helpFlag = false;
	specializedUsage = false;
	specializedUsageSet = false;

	// Initialize the interface name and usage string.
	interfaceName = name;
	usageString = "USAGE: " + name + " ";

	// Add the help flag, which is found in all text interfaces.
	standardHelpFlags.push_back( "-h" );
	standardHelpFlags.push_back( "--help" );
	addOptionFlagsNoParameters( standardHelpFlags, "Display the usage details message." );

	// Add the version flag, which is found in all text interfaces.
	standardVersionFlags.push_back( "-v" );
	standardVersionFlags.push_back( "--version" );
	addOptionFlagsNoParameters( standardVersionFlags, "Display version and copyright information for this interface." );
}

///////////////////////////////////////////////////////////////////////////////
// Add a group of option flags that don't have parameters.
///////////////////////////////////////////////////////////////////////////////
void ParseCommandLine::addOptionFlagsNoParameters(const vector<string> &list, const string &text) {
	if (list.size()==0) { PROG_ERR("addOptionFlagsNoParameters was called with an empty list of parameters."); }

	// Build the vector list into a string.
	string flags;
	for( unsigned int i = 1; i <= list.size(); i++ )
		flags += list[i-1] + " ";

	// Add the full option description.
	descriptionsOfOptionsNoFlags[flags] = text;

	// Transform flags to lower case.
	for( unsigned int i = 0; i < list.size(); i++ ) {
		string copy = list[i];
		lowerOptionsNoFlags.insert( toLower(copy) );
	}
}

///////////////////////////////////////////////////////////////////////////////
// Add a group of option flags that have parameters.
///////////////////////////////////////////////////////////////////////////////
void ParseCommandLine::addOptionFlagsWithParameters(const vector<string> &list, const string &text) {
	if (list.size()==0) { PROG_ERR("addOptionFlagsWithParameters was called with an empty list of parameters."); }

	// Build the vector list into a string.
	string flags;
	for( unsigned int i = 1; i <= list.size(); i++ )
		flags += list[i-1] + " ";

	// Add the full option description.
	descriptionsOfOptionsWithFlags[flags] = text;

	// Transform flags to lower case.
	for( unsigned int i = 0; i < list.size(); i++ ) {
		string copy = list[i];
		lowerOptionsWithFlags.insert( toLower(copy) );
	}
}

vector<string> ParseCommandLine::addFlag(const bool hasParameters, const string& flags, const string& text) {
    vector<string> tokens = split(flags, " ", false/*include-blanks*/); // split defined in src/common_utils.h
    if (hasParameters)
    	addOptionFlagsWithParameters(tokens, text);
    else
    	addOptionFlagsNoParameters(tokens, text);
    return tokens;
}

///////////////////////////////////////////////////////////////////////////////
// Add a parameter description.
///////////////////////////////////////////////////////////////////////////////
void ParseCommandLine::addParameterDescription(
	const string& id, const string &description ) {

	// Create the parameter string.
	string paramString = "<" + id + ">";

	// Update the usage string.
	usageString += ( paramString + " " );

	// Add the full parameter description.
	descriptionsOfParameters.push_back( make_pair( paramString, description ) );
}

///////////////////////////////////////////////////////////////////////////////
// Check whether the parser contains an option in a list.
///////////////////////////////////////////////////////////////////////////////
bool ParseCommandLine::contains( const vector<string> &list ) {
	// Go through the possible flags that denote the option.
	// If one is found, return true.
	const vector<string>::const_iterator end=list.end();
	for(vector<string>::const_iterator it=list.begin(); it!=end;it++) {
		// convert all characters in the options list to lower case.
		if( parsedData.find(toLower(*it)) != parsedData.end() )
			return true; 
	}
	// Return false if no option was found.
	return false;
}

///////////////////////////////////////////////////////////////////////////////
// Get an optional value as a string.
///////////////////////////////////////////////////////////////////////////////
string ParseCommandLine::getOptionString( const vector<string> &list, bool verifyFileExists  /* default: true, for historic reasons */) {
	string data;
	setOptionString(list, data, verifyFileExists);
	return data;
}

///////////////////////////////////////////////////////////////////////////////
// Get an optional value as a string.
///////////////////////////////////////////////////////////////////////////////
bool ParseCommandLine::setOptionString( const vector<string> &list, string& defaultValue, const bool verifyFileExists /* default: false */) {

	// Go through the possible flags that denote the option.
	if (error) return false;
	for( int i = 1; i <= list.size(); i++ ) {
		// Get the next option if it exists.
		string option = toLower(list[i-1]);
		if( parsedData.find( option ) != parsedData.end() ) {

			// If string is meant to be an existing file name, but the file
			// doesn't exist, set an error and return an empty string.
			string &data = parsedData[option];
			if( verifyFileExists && !fileExists(data.c_str()) ) {
				//cerr << "File required for option " << option << " does not exist. (Path: " << data.c_str() << ")" << endl;
				setErrorSpecialized(sfmt("File required for option %s does not exist. (Path: %s)", option.c_str(), data.c_str()));
				return "";
			}

			// Return the file.
			defaultValue = data;
		}
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////
// Get a parameter.
///////////////////////////////////////////////////////////////////////////////
string ParseCommandLine::getParameter( const int number, const bool verifyFileExists /* default: false */) {
	if (number > descriptionsOfParameters.size()) {
		cerr << "Programming error: invalid index for required parameter in ParseCommandLine::getParameter. 1-based index is " << number << ", but number of parameters is " << descriptionsOfParameters.size() << "." << endl;
		setError();
	}
	if (error) return "";

	// Create a string that would be this parameter's key in the data map.
	stringstream stream( stringstream::in | stringstream::out );
	stream << "param" << number;

	// Return the parameter's value.
	string data = parsedData[stream.str()];
	if( verifyFileExists && !fileExists(data.c_str()) ) {
		cerr << "File required for parameter " << number << " " << descriptionsOfParameters[number-1].first << " does not exist. (Path: " << data.c_str() << ")" << endl;
		setError();
		return "";
	}

	return data;
}

///////////////////////////////////////////////////////////////////////////////
// Get whether the parser encountered an error.
///////////////////////////////////////////////////////////////////////////////
bool ParseCommandLine::isError() {

	return error;
}

///////////////////////////////////////////////////////////////////////////////
// Get whether the parser encountered the -h (--help) flag.
///////////////////////////////////////////////////////////////////////////////
bool ParseCommandLine::isHelp() {
	return helpFlag;
}

///////////////////////////////////////////////////////////////////////////////
// Get whether the class doesn't need to handle usage.
///////////////////////////////////////////////////////////////////////////////
bool ParseCommandLine::isSpecializedUsage() {

	return specializedUsageSet;
}

// returns true if the argument looks like a flag.
// Examples: -d and --dna are flags while -6 and '-' are not.
bool isFlagLike(const string& test) {
	return test.size() > 1 && test[0] == '-' && !::isdigit(test[1]);
}
///////////////////////////////////////////////////////////////////////////////
// Parse the command line into its pieces.
///////////////////////////////////////////////////////////////////////////////
void ParseCommandLine::parseLine( int argc, char* argv[] ) {

	// Add the options piece to the usage string, because there always are at
	// least two types of options (help and version information).
	usageString += "[options]";

	// Convert the command line into strings.
	// During this process, transform all flags into lower case.
	vector<string> cmdLine;
	for( int i = 1; i < argc; i++ ) {
		string next( argv[i] );
		if (isFlagLike(next)) // it is NOT an option if it is "", "-" or "-#..." where # is a number.
			toLower(next);
		cmdLine.push_back( next );
	}
	//unsigned int length = cmdLine.size();

	// Check if a help flag exists.
	// If one does, print out a usage message, set an error, and return.
	// This isn't a real error, it's just to stop parsing and running later.
	vector<string>::iterator helpIt;
	for( unsigned int i = 1; i <= standardHelpFlags.size(); i++ ) {
		string next = standardHelpFlags[i-1];
		helpIt = find( cmdLine.begin(), cmdLine.end(), next );
		if( helpIt != cmdLine.end() ) {
			helpFlag = true;
			usage();
			setError();
			return;
		}
	}

	// Check if a version flag exists.
	// If one does, print out a version message, set an error, and return.
	// This isn't a real error, it's just to stop parsing and running later.
	vector<string>::iterator versionIt;
	for( unsigned int i = 1; i <= standardVersionFlags.size(); i++ ) {
		string next = standardVersionFlags[i-1];
		versionIt = find( cmdLine.begin(), cmdLine.end(), next );
		if( versionIt != cmdLine.end() ) {
			cout << interfaceName << ": Version " << RNASTR_BUILD_VERSION << " (" << RNASTR_BUILD_DATE << ")." << endl
			     << "Copyright Mathews Lab, University of Rochester." << endl;
			setError();
			return;
		}
	}

	// Populate the parsed data map.
	unsigned int numParams = 0;
	for( unsigned int i = 1; i <= cmdLine.size(); i++ ) {

		// Get the next command line argument.
		string next = cmdLine[i-1];

		// If the argument doesn't start with "-", call it a parameter.
		if( !isFlagLike(next) ) {
			numParams++;
			stringstream paramStream( stringstream::in | stringstream::out );
			paramStream << "param" << numParams;
			parsedData[paramStream.str()] = next;
		}

		// Otherwise, process the flag.
		else {
			// Create a variable to track whether the flag exists.
			bool exists = false;

			// Check the flag against the possible flags with no parameters.
			// Set the flag's value if found.
			set<string>::iterator it1;
			set<string> s1 = lowerOptionsNoFlags;
			for( it1 = s1.begin(); it1 != s1.end(); it1++ ) {
				string flag = *it1;
				if( next == flag ) {
					parsedData[next] = "";
					exists = true;
				}
			}

			// Check the flag against the possible options with parameters.
			// Set the flag's value if found.
			// Only do this if a flag wasn't previously found.
			if( exists == false ) {
				set<string>::iterator it2;
				set<string> s2 = lowerOptionsWithFlags;
				for( it2 = s2.begin(); it2 != s2.end(); it2++ ) {
					string flag = *it2;
					if( next == flag ) {
						bool noParam = false;
						if ( i == cmdLine.size() )
							noParam = true;
						else if ( cmdLine[i][0] == '-' ) {  //could be another flag OR a negative number: e.g. " -s -3 "
							double dummy;
							stringstream testStream( cmdLine[i] );
							noParam = !( testStream >> dummy );
						}
						
						if( noParam ) {
							cerr << "Option missing for flag: " << next << endl;
							setError();
							return;
						} else {
							parsedData[next] = cmdLine[i];
							exists = true;
							i++;
						}
					}
				}
			}

			// If the flag doesn't exist, show an error and return.
			if( exists == false ) {
				cerr << "Flag " << next << " does not exist." << endl;
				setError();
				return;
			}
		}
	}

	// If the number of parameters in the command line isn't what it should
	// be, show an error and return.
	if( numParams != descriptionsOfParameters.size() ) {
		cerr << "Incorrect number of required parameters given." 
		     << " (Found " << numParams << " but expected " << descriptionsOfParameters.size() << ".)" << endl
		     << usageString << endl 
		     << "Use any of the following options to get a help message: ";
		
		for( unsigned int i = 1; i <= standardHelpFlags.size(); i++ ) {
			cerr << standardHelpFlags[i-1] << " ";
		}

		cerr << endl;
		setError();
		return;
	}
}

///////////////////////////////////////////////////////////////////////////////
// Set that an error occurred in parsing.
///////////////////////////////////////////////////////////////////////////////
void ParseCommandLine::setError() {

	error = true;
}

///////////////////////////////////////////////////////////////////////////////
// Set that an error occurred in parsing.
///////////////////////////////////////////////////////////////////////////////
void ParseCommandLine::setError( const string &type ) {

	setError();
	cerr << "Invalid " << type << " given." << endl;
}

///////////////////////////////////////////////////////////////////////////////
// Set that an error occurred in parsing.
///////////////////////////////////////////////////////////////////////////////
void ParseCommandLine::setErrorSpecialized( const string& errString ) {

	setError();
	cerr << errString << endl;
}

///////////////////////////////////////////////////////////////////////////////
// Set an optional value as a double.
// Returns true if the option was found and was parseable, and false otherwise.
///////////////////////////////////////////////////////////////////////////////
bool ParseCommandLine::setOptionDouble(
	const vector<string> &list, double& defaultValue ) {

	// Go through the possible flags that denote the option.
	// If one is found, attempt to convert it to a double and return that value.
	if (error) return false;
	for( int i = 1; i <= list.size(); i++ ) {
		string option = toLower(list[i-1]);
		if( parsedData.find( option ) != parsedData.end() ) {
			// Try to read the given value as a double.
			// If that can't be done, show an error.
			if (parseVal(parsedData[option], defaultValue) ) return true;
			setErrorSpecialized(sfmt("Non-numeric input given for flag \"%s\". Expected a floating-point number, but found \"%s\".", option.c_str(), parsedData[option].c_str()));
			break;
		}
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////
// Set an optional value as a float.
// Returns true if the option was found and was parseable, and false otherwise.
///////////////////////////////////////////////////////////////////////////////
bool ParseCommandLine::setOptionFloat(
	const vector<string> &list, float& defaultValue ) {
	if (error) return false;
	// Go through the possible flags that denote the option.
	// If one is found, attempt to convert it to a float and return that value.
	for( int i = 1; i <= list.size(); i++ ) {
		string option = toLower(list[i-1]);
		if( parsedData.find( option ) != parsedData.end() ) {
			// Try to read the given value as a double.
			// If that can't be done, show an error.
			if (parseVal(parsedData[option], defaultValue) ) return true;
			setErrorSpecialized(sfmt("Non-numeric input given for flag \"%s\". Expected a floating-point number, but found \"%s\".", option.c_str(), parsedData[option].c_str()));
			break;
		}
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////
// Set an optional value as an integer.
// Returns true if the option was found and was parseable, and false otherwise.
///////////////////////////////////////////////////////////////////////////////
bool ParseCommandLine::setOptionInteger(
	const vector<string> &list, int& defaultValue ) {
	if (error) return false;
	// Go through the possible flags that denote the option.
	// If one is found, attempt to convert it to an int and return that value.
	for( int i = 1; i <= list.size(); i++ ) {
		string option = toLower(list[i-1]);
		if( parsedData.find( option ) != parsedData.end() ) {
			// Try to read the given value as an int.
			// If that can't be done, show an error.
			if (parseVal(parsedData[option], defaultValue) ) return true;
			setErrorSpecialized(sfmt("Non-numeric input given for flag \"%s\". Expected an integer, but found \"%s\".", option.c_str(), parsedData[option].c_str()));
			break;
		}
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////
// Set an optional value as a long integer.
// Returns true if the option was found and was parseable, and false otherwise.
///////////////////////////////////////////////////////////////////////////////
bool ParseCommandLine::setOptionLong(
	const vector<string> &list, long& defaultValue ) {
	if (error) return false;
	// Go through the possible flags that denote the option.
	// If one is found, attempt to convert it to an int and return that value.
	for( int i = 1; i <= list.size(); i++ ) {
		string option = toLower(list[i-1]);
		if( parsedData.find( option ) != parsedData.end() ) {
			// Try to read the given value as an int.
			// If that can't be done, show an error.
			if (parseVal(parsedData[option], defaultValue) ) return true;
			setErrorSpecialized(sfmt("Non-numeric input given for flag \"%s\". Expected an integer, but found \"%s\".", option.c_str(), parsedData[option].c_str()));
			break;
		}
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////
// Tell the class that it doesn't need to handle usage.
///////////////////////////////////////////////////////////////////////////////
void ParseCommandLine::setSpecializedUsage() {

	specializedUsage = true;
}

///////////////////////////////////////////////////////////////////////////////
// Print out a detailed usage message.
///////////////////////////////////////////////////////////////////////////////
void ParseCommandLine::usage() {

	// If the standard usage message should not be printed out, return.
	if( specializedUsage == true ) {
		specializedUsageSet = true;
		return;
	}

	// Print out a short usage string and general flag information.
	cout << usageString << endl
	     << "All flags are case-insensitive, "
	     << "and grouping of flags is not allowed." << endl << endl;

	// Print out the required parameters.
	cout << "=============================" << endl
	     << "==== Required Parameters ====" << endl
	     << "=============================" << endl;
	unsigned int numParams = descriptionsOfParameters.size();
	for( unsigned int i = 1; i <= numParams; i++ ) {
		pair<string, string> nextPair = descriptionsOfParameters[i-1];
		cout << nextPair.first << endl;
		printDescription( nextPair.second );
	}

	// Print out option flags without parameters, if they exist.
	unsigned int numOptionsNoParams = descriptionsOfOptionsNoFlags.size();
	if( numOptionsNoParams >= 1 ) {
		cout << "=========================================" << endl
		     << "==== Option Flags Without Parameters ====" << endl
		     << "=========================================" << endl;
		map<string, string, compareOptions>::iterator it2;
		map<string, string, compareOptions> m2 = descriptionsOfOptionsNoFlags;
		for( it2 = m2.begin(); it2 != m2.end(); it2++ ) {
			cout << it2->first << endl;
			printDescription( it2->second );
		}
	}

	// Print out option flags with parameters, if they exist.
	unsigned int numOptionsWithParams = descriptionsOfOptionsWithFlags.size();
	if( numOptionsWithParams >= 1 ) {
		cout << "======================================" << endl
		     << "==== Option Flags With Parameters ====" << endl
		     << "======================================" << endl
		     << "All parameters must follow their associated flag directly."
		     << endl
		     << "Failure to do so may result in aberrant program behavior."
		     << endl << endl;
		map<string, string, compareOptions>::iterator it1;
		map<string, string, compareOptions> m1 = descriptionsOfOptionsWithFlags;
		for( it1 = m1.begin(); it1 != m1.end(); it1++ ) {
			cout << it1->first << endl;
			printDescription( it1->second );
		}
	}
}

void replaceAll(std::string& str, const std::string& find, const std::string& replaceWith) {
    if(find.empty()) return;
    size_t pos = 0;
    while((pos=str.find(find, pos)) != std::string::npos) {
        str.replace(pos, find.length(), replaceWith);
        pos += replaceWith.length();
    }
}

void addLines(vector<string> &lines, string line, const int maxlen) {
		// if the line is too long, split it into smaller lines.
		// Try to break at spaces, but if there are no spaces, add a hypen.

		// Do special treatment of lines that begin with a tab. Interpret these as a list of options.
		// Remove the tab character, but indent the line and also indent any wrapped segments.
		string indent;
		if (!line.empty()&&line[0]=='\t') {
			indent = "  ";
			line.erase(0,1);
		}
		while(line.length() > maxlen-indent.size()) {
			size_t end=line.find_last_of(" \t", maxlen+1); // use maxlen+1 because the found character will not actually be part of the line.
			if (end==string::npos) {
				end = maxlen-1; // subtract 1 because we add with a hyphen and the total length should be maxlen
				lines.push_back(indent+line.substr(0, end)+"-"); 
			} else {
				//Specify how SHAPE data should be handled when multiple values exist for the
				lines.push_back(indent+line.substr(0, end)); // use end (instead of end+1) because we do not include the space or tab in this line.
			}
			// if the character was a space, skip it. If it was a tab or other character (in the case of a hyphenated word), put it at the start of the next line.
			line = line.substr(line[end]==' '?end+1:end);
			if (indent.size()==2) indent="    "; // indent subsequent lines further.
		}
		if (line.size() != 0)
			lines.push_back(indent+line);
}

///////////////////////////////////////////////////////////////////////////////
// Word wrap a string.
///////////////////////////////////////////////////////////////////////////////
void ParseCommandLine::printDescription( const string &text ) {
	// The descriptions are wrapped at a 78-character limit.
	// Previously, this function wrapped text after periods, but that isn't necessary and 
	//     wrapped lines at wrong places (e.g. in numbers "3.5" or abbreviations "i.e." etc)  
	// Now lines are wrapped irrelevant of periods, except in two cases (discussed below).
	// A program can still "force" text onto a new line using any of the following techniques:
	//   1) Add a newline character (\n) in the text
	//   2) Follow a period with two spaces. e.g. "Here is a line.  This will be on another line."
	//   3) For legacy reasons, a sentence starting with the text " Default" is also pushed onto a new line. 
	//      e.g. "Here is some parameter. Default is 3."

	// Determine the maximum width of the string to wrap.
	const size_t bufferSize = 78;
	const string spacer("    ");
	const size_t lineSize = bufferSize - spacer.length(); // maximum length of a line of content (after spacer)
	vector<string> lines;
	size_t pos;
	size_t prev = 0;
    while(prev<text.size()) {

		pos=min(text.find("\n", prev), text.find(".  ", prev));
		pos=min(pos, text.find(". Default", prev));
		if (pos==string::npos) {
			// add all the remaining text (after the last . or \n )
			addLines(lines, text.substr(prev), lineSize);
			break;
		} else {
			// if the char is a dot, keep it, otherwise use the part of the string BEFORE the \n, thus discarding it.
			addLines(lines, text.substr(prev, text[pos]=='\n'?pos-prev:pos-prev+1), lineSize);
			if (text[pos]=='.') {
				// remove up to two leading spaces from the next line following a period.
				if (pos+1<text.size()&&text[pos+1]==' ') pos++;
				if (pos+1<text.size()&&text[pos+1]==' ') pos++;
			}
		}
		prev=pos+1; // outer loop searching for "." or "\n"
    }
	for(vector<string>::iterator it=lines.begin(); it!=lines.end();it++)
		cout << spacer << *it << endl;

	cout << endl;
}

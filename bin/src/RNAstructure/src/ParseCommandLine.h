/*
 * A header file for a program that parses Unix-like command lines to determine
 * the data being given for input.
 *
 * (c) 2008  Mathews Lab, University of Rochester Medical Center
 * Redone in 2012.
 * Written by Jessica S. Reuter
 */

#ifndef PARSE_COMMAND_LINE_H
#define PARSE_COMMAND_LINE_H

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

class ParseCommandLine {
 public:
	// Public methods.

	/*
	 * Name:        Constructor
	 * Description: Initializes the parsing object.
	 * Arguments:
	 *     1. name
	 *        The interface name.
	 */
	ParseCommandLine( string name );

	/*
	 * Name:        addOptionFlagsNoParameters
	 * Description: Add a group of option flags that don't have parameters to
	 *              the parser.
	 *              These options always set a boolean, an on-off switch.
	 * Arguments:
	 *     1. list
	 *        The list of flags.
	 *     2. text
	 *        The string describing these parameters.
	 */
	void addOptionFlagsNoParameters( const vector<string> &list, const string &text );

	/*
	 * Name:        addOptionFlagsWithParameters
	 * Description: Add a group of option flags that don't have parameters to
	 *              the parser.
	 * Arguments:
	 *     1. list
	 *        The list of flags.
	 *     2. text
	 *        The string describing these parameters.
	 */
	void addOptionFlagsWithParameters( const vector<string> &list, const string &text );

	/*
	 * Name:        addFlags
	 * Description: A more convenient form of the addOptionFlags... functions.
	 *              This function splits a space-delimited string of flags, 
	 *              then calls either addOptionFlagsNoParameters or addOptionFlagsWithParameters
	 *              and returns the vector of parsed flags.
	 * Arguments:
	 *     1. Whether or not these flags have parameters.
	 *     2. A string containing the flag and any alias flags. Flags must be separated by spaces.
	 *     3. text
	 *        The string describing these parameters.
	 */
	vector<string> addFlag(const bool hasParameters, const string &flags, const string &text);

	/*
	 * Name:        addFlags
	 * Description: A more convenient form of the addOptionFlags... functions.
	 *              This function splits a space-delimited string of flags, 
	 *              then calls either addOptionFlagsNoParameters or addOptionFlagsWithParameters
	 *              and returns the vector of parsed flags.
	 * Arguments:
	 *     1. Whether or not these flags have parameters.
	 *     2. A string containing the flag and any alias flags. Flags must be separated by spaces.
	 *     3. text
	 *        The string describing these parameters.
	 *     4. defaultValue
	 *        The default value of the parameter. This is only for display purposes. 
	 *        It adds the text "(default: ...)" to the description of the flag.
	 */
	template<typename T>
	vector<string> addFlag(const bool hasParameters, const string& flags, const string& text, T defaultValue) {
		std::stringstream ss;
		ss << text << " (default: " << defaultValue << ").";
		return addFlag(hasParameters, flags, ss.str());
	}

	/*
	 * Name:        addParameterDescription
	 * Description: Add a description of the next required parameter.
	 * Arguments:
	 *     1. id
	 *        The next parameter's id (shows up in brackets)
	 *     2. description
	 *        The next parameter's description.
	 */
	void addParameterDescription( const string &id, const string &description );

	/*
	 * Name:        contains
	 * Description: Determine if the parser contains a particular option.
	 * Arguments:
	 *     1. list
	 *        The list of option flags to check.
	 * Returns:
	 *     True if the command line contains the option, false if not.
	 */
	bool contains( const vector<string> &list );

	/*
	 * Name:        getOptionString
	 * Description: Get a particular option as a string.
	 * Arguments:
	 *     1. list
	 *        The list of option flags to check.
	 *     2. exists
	 *        True if the string is a file name that should exist, false if not.
	 *        Default is true.
	 * Returns:
	 *     The option set by these flags, as a string.
	 */
	string getOptionString( const vector<string> &list, const bool verifyFileExists = true);

	/*
	 * Name:        setOptionString
	 * Description: Get a particular option as a string.
	 * Arguments:
	 *     1. list
	 *        The list of option flags to check.
	 *     2. exists
	 *        True if the string is a file name that should exist, false if not.
	 *        Default is true.
	 *     3. defaultValue
	 *        A reference to the string that should receive the option value.
	 *        If the option is not found, the value will be left as-is, so it can be pre-set to a default value before this call.
	 * Returns:
	 *     The option set by these flags, as a string.
	 */
	bool setOptionString( const vector<string> &list, string &defaultValue, const bool verifyFileExists = false);

	/*
	 * Name:        getParameter
	 * Description: Get a required parameter.
	 * Arguments:
	 *     1. number
	 *        The parameter number, one-indexed.
	 */
	string getParameter( const int number, const bool verifyFileExists = false );

	/*
	 * Name:        isError
	 * Description: Determine if the parser encountered an error.
	 * Returns:
	 *     True if the parser encountered an error, false if not.
	 */
	bool isError();

	/*
	* Name:        isHelp
	* Description: Determine if the parser encountered the help flag (-h or --help).
	* Useful for distinguishing this state from other errors 
	* (because the presence of the help flag sets the error state to true).
	*/
	bool isHelp();

	/**
	 * Name:        getSpecializedUsage
	 * Description: Get whether the class doesn't handle its own usage messages.
	 * Returns:
	 *     True if the class doesn't handle its own messages, false if it does.
	 */
	bool isSpecializedUsage();

	/*
	 * Name:        parseLine
	 * Description: Parse the command line for parameter and option data.
	 * Arguments:
	 *     1. argc
	 *        The number of command line arguments.
	 *     2. argv
	 *        The command line arguments themselves.
	 */
	void parseLine( int argc, char* argv[] );

	/*
	 * Name:        setError
	 * Description: Set that an error occurred in parsing, but don't show any
	 *              specifics about it.
	 */
	void setError();

	/*
	 * Name:        setError
	 * Description: Set that an error occurred in parsing, and show it.
	 *              The error is printed out as "Invalid <type> given."
	 * Arguments:
	 *     1. type
	 *        The type of error that occurred.
	 */
	void setError( const string &type );

	/**
	 * Name:        setErrorSpecialized
	 * Description: Set that an error occurred, and print out a specialized
	 *              error given by the user.
	 * Arguments:
	 *     1. errString
	 *        The error string to print out.
	 */
	void setErrorSpecialized( const string &errString );

	/*
	 * Name:        setOptionDouble
	 * Description: Set a particular option as a double.
	 * Arguments:
	 *     1. list
	 *        The list of option flags to check.
	 *     2. defaultValue
	 *        The default value this option has.
	 * Returns:
	 *		True if the option was found and was parseable, and false otherwise.
	 */
	bool setOptionDouble( const vector<string> &list, double& defaultValue );

	/*
	 * Name:        setOptionFloat
	 * Description: Set a particular option as a float.
	 * Arguments:
	 *     1. list
	 *        The list of option flags to check.
	 *     2. defaultValue
	 *        The default value this option has.
	 * Returns:
	 *		True if the option was found and was parseable, and false otherwise.
	 */
	bool setOptionFloat( const vector<string> &list, float& defaultValue );

	/*
	 * Name:        setOptionInteger
	 * Description: Set a particular option as an integer.
	 * Arguments:
	 *     1. list
	 *        The list of option flags to check.
	 *     2. defaultValue
	 *        The default value this option has.
	 * Returns:
	 *		True if the option was found and was parseable, and false otherwise.
	 */
	bool setOptionInteger( const vector<string> &list, int& defaultValue );

	/*
	 * Name:        setOptionLong
	 * Description: Set a particular option as a long integer.
	 * Arguments:
	 *     1. list
	 *        The list of option flags to check.
	 *     2. defaultValue
	 *        The default value this option has.
	 * Returns:
	 *		True if the option was found and was parseable, and false otherwise.
	 */
	bool setOptionLong( const vector<string> &list, long& defaultValue );

	/**
	 * Name:        setSpecializedUsage
	 * Description: Sets a flag in the ParseCommandLine class that says an
	 *              external class, not this one, will handle usage messages.
	 */
	void setSpecializedUsage();

	/**
	 * Name:        printDescription
	 * Description: Output a string wrapped so it's no wider than 79 characters.
	 */
	void printDescription( const string &text );

 private:
	// Private usage method.

	/*
	 * Name:        usage
	 * Description: Print out a long usage message with detailed descriptions.
	 */
	void usage();

 private:
	// A private struct that holds a comparator for option flag maps.
	// This makes flags ignore case and number of dashes when sorting.
	struct compareOptions {
		bool operator()( string a, string b ) const {
			string one = a.substr( a.find_first_not_of( '-' ) );
			transform( one.begin(), one.end(), one.begin(), ::tolower );

			string two = b.substr( b.find_first_not_of( '-' ) );
			transform( two.begin(), two.end(), two.begin(), ::tolower );

			return ( one < two );
		}
	};

 private:
	// Private variables.

	// Boolean that tracks whether the parser encountered an error.
	bool error;

	// Boolean that tracks whether the parser encountered the help flag (to distinguish this state from other errors)
	bool helpFlag;

	// The interface name.
	string interfaceName;

	// The map of descriptions for options without flags.
	map<string, string, compareOptions> descriptionsOfOptionsNoFlags;

	// The map of descriptions for options with flags.
	map<string, string, compareOptions> descriptionsOfOptionsWithFlags;

	// The vector of descriptions for parameters.
	vector< pair<string, string> > descriptionsOfParameters;

	// The set of options without parameters transformed to lower case.
	set<string> lowerOptionsNoFlags;

	// The set of options with parameters transformed to lower case.
	set<string> lowerOptionsWithFlags;

	// The map of parsed data.
	map<string, string> parsedData;

	// Boolean flag, true if a specialized usage message should be implemented
	// outside of the class, false if not.
	bool specializedUsage;

	// Boolean flag, true if a specialized usage message has been asked for,
	// false if not.
	bool specializedUsageSet;

	// The vector of standardized help flags.
	vector<string> standardHelpFlags;

	// The vector of standardized version flags.
	vector<string> standardVersionFlags;

	// The simple usage string.
	string usageString;
};

#endif /* PARSE_COMMAND_LINE_H */

#ifndef COMMON_UTILS_H
#define COMMON_UTILS_H

#include "defines.h"
#include <iostream>
#include <string>
#include <iterator>
#include <vector>
#include <sstream>
#include <cstring>

using namespace std;

#define MAX_FILENAME_LEN 255 // maximum filename length (in Windows) used by createSafeFilename

//###################### FILE Utilities ##################################
//! Determine if the path <directory>/<file> exists (and is not itself a directory)
bool fileExists(const char* const directory, const char* const file);
//! Determine if the file specified by filePath exists and is not itself a directory.
bool fileExists(const char* const filePath, const bool verifyReadable = false);
//! Determine if the file specified by filePath exists and is not itself a directory.
inline bool fileExists(const string& filePath, const bool verifyReadable = false) { return fileExists(filePath.c_str(), verifyReadable); }

//! Determine if the specified directory exists.
bool dirExists(const char* const fullDirPath);
//! Determine if the file specified by filePath exists and is not itself a directory.
inline bool dirExists(const string& fullDirPath) { return dirExists(fullDirPath.c_str()); }

//! Determine whether the specified file path signifies that STDIO should
//! be used instead of a regular file. Currently this means the file path must be
//! identically "-" (as commonly used in POSIX programs for the same purpose).
bool isStdIoFile(const char* const fullPath);

//! Gets the portion of a path that represents the filename.
//! I.e. removes the directory portion of a path. 
//! If removeExtension is true, the file extension is also removed.
string getFileName(const char * const path, bool removeExtension = false);

//! Gets the portion of a path that represents the file extension 
//! i.e. the part of the file name after the last dot (.)
string getFileExt(const string& filePath);

//! Gets the directory portion of a path. 
//! It is assumed that the path represents a file.
string getDirName(const char * const path);
//###################### STRING Utilities ##################################

//! Create a valid file name from a string label (which could contain invalid characters)
//! and an optional extension.
//! If the label is too long for a filename, the name is truncated.
string createSafeFilename(const string& label, const string& extension="", bool replaceSpaceChar=false);

// Utility function that creates a string and fills it by sprintf-formatting the arguments. 
string sfmt(const char* const format, ...);


// The following string functions operate in-place on their subjects
// and then return a reference to that same subject (for coding fluency).

//! Converts control characters in subject into escape sequences.
//! (useful for debugging strings that might contain control characters.)
//! e.g. "Hello\n" becomes "Hello\\n" 
//! i.e. The linefeed is converted into a literal slash and the letter 'n'
//! (operates in-place on subject and returns it).
string& escapeChars(string &subject);

//! Finds invalid file name characters in a string and replaces them with the 
//!   specified character (usually underscore '_').
//! Operates in-place on the subject string.
string& replaceInvalidFileNameChars(string &subject, char replaceWith='_', bool replaceSpace=false);

//! Trim leading whitespace from a string (operates in-place on subject and returns it).
string& trimLeft(string &subject);

//! Trim trailing whitespace from a string (operates in-place on subject and returns it).
string& trimRight(string &subject);

//! Trim leading and trailing whitespace from a string (operates in-place on subject and returns it).
string& trim(string &subject);

//! Converts a string to lower-case (operates in-place on subject and returns it).
string& toLower(string &subject);


// The following string functions operate on const strings and 
// return a copy of the string with the desired modifiations.

//! Converts a string to upper-case (operates in-place on subject and returns it).
string& toUpper(string &subject);

//! Returns a copy of the subject in which leading whitespace has been removed.
string trimLeft(const string &s);

//! Returns a copy of the subject in which trailing whitespace has been removed.
string trimRight(const string &s);

//! Returns a copy of the subject in which leading and trailing whitespace has been removed.
string trim(const string &s);

//! Returns a copy of the subject in which all characters have been converted to lower case.
string toLower(const string &subject);

//! Returns a copy of the subject in which all characters have been converted to upper case.
string toUpper(const string &subject);

//! Find a character in subject (a c-string)
//! Returns the 0-based index where the char is found or string::npos if it is not found.
std::size_t findchr(const char* const subject, const char find);

//! Returns true if a c-string is NULL or empty ("")
//! Equivalent to (cstr==NULL||strlen(cstr)==0)
inline bool is_blank(const char*const cstr) { return cstr==NULL||*cstr=='\0'; } // if the first character is the null-char, '\0', the string length is 0.

//! Make a copy of a c-string (i.e. a char* or char[] etc) 
//! There's a lot of code duplicated around RNAstructure that does this (and should use
//! this function instead).
//! \returns A pointer to a newly allocated buffer that is a copy of the source buffer.
//!          The returned pointer must be deleted later on.
char* copy_cstr(const char* src);

//! Converts a cstring (char*) into a string, making sure not to
//! dereference a NULL pointer.
//! \return A string copy of cstr or an empty string if cstr is NULL.
inline string as_str(const char*const cstr) { return cstr==NULL?"":cstr; }

// Find the next non-whitespace character and return a pointer to it.
// Note:  '\0' is not considered whitespace, so every valid 
// (null-terminated) c-string has at least one non-whitespace character --
// namely the terminal one.
inline const char* nextNonWhitespaceChar(const char* input) { while(std::isspace(*input)) input++; return input; }
// Find the next non-whitespace character and return a pointer to it.
inline char* nextNonWhitespaceChar(char* input) { while(std::isspace(*input)) input++; return input; } // non-const version

//! Attempts to convert a (null-terminated) c-string into an int (e.g. " 12 " becomes 12)
//! Returns true if the conversion succeeded or false if it failed
//!  (e.g. the value was NOT a number or the number was out of range).
//! Whitespace before or after the number is ignored, but extra text after 
//! the number will cause the conversion to fail if readToEnd is true (the default).
bool parseInt(const char* const input, int& intResult, const bool readToEnd=true);


//! Attempts to convert a (null-terminated) c-string into an double (e.g. "6.02E23" becomes (double)6.02E23)
//! Returns true if the conversion succeeded or false if it failed
//!  (e.g. the value was NOT a number or the number was out of range).
//! Whitespace before or after the number is ignored, but extra text after 
//! the number will cause the conversion to fail if readToEnd is true (the default).
bool parseDbl(const char* const input, double& dblResult, const bool readToEnd=true);


// Parse a c-string to get an integer.
// If conversion succeeds, the parsed value is returned. Otherwise defaultValue is returned.
inline int tryParseInt(const char* const input, const int defaultValue) { int value=defaultValue; parseInt(input, value); return value; }

//! Parse a c-string to get a double. 
//! If conversion succeeds, the parsed value is returned. Otherwise defaultValue is returned.
inline bool tryParseDbl(const char* input, const double defaultValue){ double value=defaultValue; parseDbl(input, value); return value; }

//! Attempts to convert a string into an int (e.g. " 12 " becomes 12)
//! Returns true if the conversion succeeded or false if it failed
//!  (e.g. the value was NOT a number or the number was out of range).
//! Whitespace before or after the number is ignored, but extra text after 
//! the number will cause the conversion to fail.
inline bool parseInt(const string& s, int& intResult) { return parseInt(s.c_str(), intResult); /*stringstream ss(s); return !(ss >> intVal).fail();*/ }

//! Attempts to convert a string into a double (e.g. "6.02E23" becomes (double)6.02E23)
//! Returns true if the conversion succeeded or false if it failed 
//! (e.g. the value was NOT a number or the number was out of range).
//! Whitespace before or after the number is ignored, but extra text after 
//! the number will cause the conversion to fail.
inline bool parseDbl(const string& s, double& dblResult) { return parseDbl(s.c_str(), dblResult); /*stringstream ss(s); return !(ss >> dblVal).fail();*/ }

//! Attempts to convert a string into another type (e.g. short, char, long int, etc)
//! returns true if the conversion succeded or false if it failed (e.g. the value was NOT convertible, using a stringstream).
template<typename T>
bool parseVal(const string &subject, T& result) {
	stringstream ss(subject); return !(ss >> result).fail();
}

//! write vector contents to an output stream. used by join and the 'ofstream<<vector' operator.
template<class T>
void join(ostream &out, const vector<T> &v, const char* const delim=",") {
	if(v.size() > 1) copy(v.begin(),v.end()-1, ostream_iterator<T>(out, delim));
	if(!v.empty()) out << v.back();
}
//! join any vector into a string with elements separated by the specified delimiter.
//! (Useful for cout debugging.)
template<class T>
string join(const vector<T> &v, const char* const delim=",") {
	ostringstream oss; join(oss,v,delim); return oss.str();
}

//! Split a string using a delimiter.
vector<string> split(const string& str, const string& delim, const bool includeEmptyValues=true);

//! Define the ostream << vector operator. (useful for debugging)
//! e.g. cout << "Vector Contents: " << v << endl;
template <typename T>
std::ostream& operator<<(ostream &out, const vector<T> &v) {
  out << '[';
  join(out,v,", ");
  return out << "]";
}

//! Class used to implement a no-op output stream buffer.
class NullBuffer : public std::streambuf { 	public: int overflow(int c) { return c; } };
//! Class used to implement a no-op output stream.
class NullStream : public std::ostream {
	public: NullStream() : std::ostream(&m_sb) {}
	static NullStream Default;
	private: NullBuffer m_sb;
};

//###################### Pointer and Memory Utilities ##################################

//! The auto_delete class provides a very limited implementation of auto_ptr (precursor to unique_ptr etc)
//! It only does one thing -- deletes a pointer when the auto_delete variable goes out of scope.
//! It is not meant to be copied or passed by value etc.
//! One use is to make it easier to return from a function without having to worry about calling delete for 
//! dynamically allocated objects at the end of the function.
//! Another use is to prevent memory analysis tools (e.g. valgrind) from complaining about 
//! file-level variables (aka global variables) that are dynamically allocated singletons that 
//! live permanently until the end of program execution.
//! With the advent of C++11 there are better ways to handle this and when RNAstructure
//! can routinely be built with c++11 (or later), this class can be removed/replaced.
template<typename T,bool IS_ARRAY=false> 
struct auto_delete {
	auto_delete(T *p): _p(p){}
	~auto_delete() { if (_p!=NULL) { if (IS_ARRAY) delete[] _p; else delete _p;} }
	T& operator*() const { return *_p; }
	operator T*() const { return _p; }
	T* get() const { return _p; }
private:
	T*_p;
	auto_delete(const auto_delete &copy){}  // disable the copy constructor
	auto_delete& operator=(const auto_delete &copy){} // disable the assignment operator
};

//! Cross-Platform `getline` 
//! Accepts DOS ("\r\n"), Old-Mac ("\r"), and Unix ("\n") line endings.
//! This allows e.g. DOS files to be read on Linux etc. 
//! (Note that even as of 2018, some programs on Mac OSX still 
//! use the old-style carriage-return-only line endings.)
std::istream& getlineXP(std::istream& in, std::string& line);

//###################### Operating System Compatibility ##################################
//! Windows's printf function incorrectly outputs 3-digit exponents when only 2 digits 
//! are required, violating the C99+ spec. e.g. Windows: "1e-002" vs Linux: "1e-02"
//! This causes issues with regression tests due to different output on each OS.
//!
//! Calling this function (before relevant output) fixes the problem.
void FixWindowsExponentFormat();

#endif // COMMON_UTILS_H

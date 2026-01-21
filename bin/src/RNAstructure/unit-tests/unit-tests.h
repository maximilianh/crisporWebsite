// very basic unit testing of the RNAstructure library

#ifndef UNIT_TESTS_H
#define UNIT_TESTS_H

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <typeinfo>
#include "../RNA_class/RNA.h"

#define TEST_EQUAL(ACTUAL,EXPECTED)  test_result(ACTUAL, EXPECTED, __LINE__)
using namespace std;

void post_result(bool passed);

// typeName is a function that outputs the dynamic type of a variable.
template<typename T>
const string typeName(const T &value) { return typeid(value).name(); }
template<>
const string typeName<string>(const string &value) {  return "string"; }

template<typename T>
const string typeName(const vector<T> &value) { 
	const string& info = value.size()!=0 ? typeName(value[0]) : string(typeid(T).name());
	return "vector<" + info + ">";
}

// FVal is a class to customize ostream (e.g. cout) output of unit-test results.
template<class T> struct FVal { const T value; FVal(const T &value) : value(value) { } };
// The generic implementation of cout << FVal<T>(x) is to output x followed by its typeid in brackets.
template<class T> ostream& operator<<(ostream& out, const FVal<T>& fval) { return out << fval.value << flush; } // << " [" << typeName(fval.value) << "]";  }
// A convenience function to convert a value (x) into an FVal<X>(x). This simply aleviates the need to specify <T> directly.
template<class T> FVal<T> fmt(const T &value) { return FVal<T>(value); }

// Explicitly create `fmt(const char *value)`, because otherwise the  
// fmt<T> template will be instantiated for each different size of char[N]
//FVal<const char *> fmt(const char *value) { cerr << endl << "USING FMT CHAR*USING FMT CHAR*USING FMT CHAR*" << endl; return FVal<const char *>(value); }

// a function to facilitate the display of strings 
// used by both char* and string ostreams<< operators
string dispstr(const string& str) {
	string copy = str;
	escapeChars(copy);
	if (copy.size() > 90) { copy.resize(86); copy.append("<..>"); }
	return copy;
}

// Customize the output of a const char*
ostream& operator<<(ostream& out, const FVal<const char*>& fval) { 
	if (fval.value==NULL)
		return out << "<NULL> [CHAR*]";
	else
		return out << "\"" << dispstr(fval.value) << "\" [" << strlen(fval.value) << " char*]";
}

// Customize the output of a string
ostream& operator<<(ostream& out, const FVal<string>& fval) { 
	return out << "\"" << dispstr(fval.value) << "\" [" << fval.value.size() << " str]";
}

// Customize the output of a vector
template<class T>
ostream& operator<<(ostream& out, const FVal<vector<T> >& fval) { 
	const vector<T> &v = fval.value;
	const size_t size=v.size();
	out << typeName(v) << "[" << size << "] { ";
	for(size_t i=0; i<size;i++)
		out << (i==0 ? "" : ",")
	        << i << ":" << fmt(v[i]);
	return out << "}";
}

template<class T>
bool compare(const T& actual, const T& expected) { return actual==expected; }
bool compare(const char* actual, const char* expected) { return actual==expected || (actual!=NULL&&expected!=NULL&&strcmp(actual,expected)==0); }
template<class T>
bool compare(const vector<T> &actual, const vector<T> &expected) { 
	if (actual.size()!=expected.size())
		return false;
	const size_t size = actual.size();
	for(size_t i=0; i<size; i++)
		if (!compare(actual[i], expected[i])) return false;
	return true;
}

template<class T>
void test_result(const T& actual, const T& expected, const int line) {
	bool pass = compare(actual,expected);
	post_result(pass);
	// if (pass)
	// 	std::cout << "PASSED TEST[" << line << "] Expected==Actual==" << actual << std::endl; 
	// else
	// 	std::cerr << "!! FAILED TEST[" << line << "] Expected: " << expected << " Actual: " << actual << std::endl; 
	if (pass)
		std::cout << "PASSED TEST[" << line << "] Expected==Actual==" << fmt(actual) << std::endl; 
	else
		std::cerr << "!! FAILED TEST[" << line << "] Expected: " << fmt(expected) << " Actual: " << fmt(actual) << std::endl; 
}

// Deal with mismatched template arguments
// Note that all string/char* cases can be handled with a single
// function: `test_result(string,string)` 
// however this would ignore cases where the actual
// result is (char*)NULL.
// A smaller subset of the following functions are also 
// sufficient for well-authored test cases, but for 
// highest robustness, we include all of them:

// #1 required for (char[N], char[M]) etc.
void test_result(const char* actual, const char* expected, const int line) {
	test_result<const char*>(actual, expected, line);
}
// #2 required for (char[N], string) etc.  -- not necessary 
//    if all tests use char* expected values when the actual value is char*
void test_result(const char* actual, const string& expected, const int line) {
	test_result(actual, expected.c_str(), line); // uses #1
}
// #3 required for (string, char[N]) etc.
void test_result(const string& actual, const char* expected, const int line) {
  test_result(actual, (string)expected, line); // uses test_result<string> from template
}


#endif //UNIT_TESTS_H
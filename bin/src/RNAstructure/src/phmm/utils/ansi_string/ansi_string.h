#ifndef _ANSI_STRING_
#define _ANSI_STRING_

//#include <iostream>
#include <vector>
using namespace::std;

// This is the initial buffer size for string with no initializer.
#define IBS (5000)

class t_string;

typedef vector<t_string*> t_string_tokens;

/*
The ansi string library:
Portable string representation for extendible, memory managed string with big capabilities,
This is necessary for not having to implement a new version of this code every time a new application is
implemented. tokenizer and sprintf are two vital functions to implement.
*/

class t_string
{
public:
	t_string(const char* string);
	t_string(t_string* string);
	t_string();
	~t_string();

	char* obj_string; // null terminated string that contains the string information.
	int obj_str_mem_size; // memory size of the string, not the actual size. This number is greater than or equal to length of the string.

	// Tokenize the string of this string object and return a
	// vector of string objects.
	t_string_tokens* tokenize_by_chars(const char* delimiter_list);
	t_string_tokens* tokenize_by_str(const char* delimiter_string);
	char* substring(int i, int j);

	void copy(const char* string);
	void copy(t_string* string);

	// Static string library functions.
	static void copy(char* dest_string, const char* src_string);
	static int string_length(const char* string);
	static int string_length(t_string* string);
	static t_string* num2str(int num, int base);
	static int str2num(char* num_str, int base);
	static int str2num(t_string* num_str, int base);
	static bool compare_strings(t_string* str1, t_string* str2);
    static bool compare_strings(char* str1, char* str2);
    static bool compare_strings_ci(t_string* str1, t_string* str2);
    static bool compare_strings_ci(char* str1, char* str2);
	static void to_upper(char* string);
	static void clean_tokens(t_string_tokens* tokens);

	static bool is_balanced(char* str, char* left_pars, char* right_pars);

	static void replace_avoid_list(char* str, char* avoided_char_list, char char_to_replace);

	// Parse the consercutive numbers in the string and return them in a vector.
	vector<int>* get_integers_in_string();

	void remove_beginning_spaces();

	bool is_balanced(char* left_pars, char* right_pars);

	// Concatenation functions. These are the basis of other functions because
	// these do the memory scaling. No other functions should need the memory scaling for the string.
	void concat_char(char _char);
	void concat_string(char* string);
	void concat_string(t_string* string);
	void concat_int(int i_num);
	void concat_float(double f_num);

	// sprintf declaration.
	void sprintf(char* fmt_string, ...);

	bool compare(char* string);
	bool compare(t_string* string);
	bool compare_ci(t_string* string);
	bool compare_ci(char* string);
	bool starts_with(char* string);
	bool starts_with(t_string* string);

	char& x(int i);
	char* str();
	int length();
	void empty();
	void revert();	
	void to_upper();
};

#endif // _ANSI_STRING_




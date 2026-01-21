#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h> // For va_... functions.
#include "ansi_string.h"
#include "ctype.h"

// DO NOT INCLUDE string.h!

t_string::t_string(const char* string)
{
        // Allocate the memory.
        this->obj_string = (char*)malloc(string_length(string) + 3);
        this->obj_str_mem_size = string_length(string) + 1;

        this->copy(string);
}

t_string::t_string(t_string* string)
{
        // Allocate the memory.
        this->obj_string = (char*)malloc(string_length(string) + 3);
        this->obj_str_mem_size = string_length(string) + 1;

        this->copy(string);
}

t_string::t_string()
{
        this->obj_string = (char*)malloc(IBS + 3);
        this->obj_str_mem_size = IBS;

        // Initialize the string.
        this->empty();
}

t_string::~t_string()
{
        //printf("Freeing a string object.\n");
        free(this->obj_string);
}

vector<int>* t_string::get_integers_in_string()
{
        vector<int>* int_vector = new vector<int>();

        // The number string.
        t_string* cur_int_str = new t_string();

        for(int char_i = 0; char_i < this->length(); char_i++)
        {
                if(this->x(char_i) >= '0' && this->x(char_i) <= '9')
                {
                        // Add this character.
                        cur_int_str->concat_char(this->x(char_i));
                }
                else
                {
                        if(cur_int_str->length() > 0)
                        {
                                int new_num = str2num(cur_int_str, 10);

                                // Add the new number
                                int_vector->push_back(new_num); 

                                cur_int_str->copy("");
                        }
                }
        } // char index.

        // Process the last string if it is a number.
        if(cur_int_str->length() > 0)
        {
                int new_num = str2num(cur_int_str, 10);

                // Add the new number
                int_vector->push_back(new_num); 

                cur_int_str->copy("");
        }

        // Free integer storing string.
        delete(cur_int_str);

        return(int_vector);
}

void t_string::copy(char* dest_string, const char* src_string)
{
        int l_src_string = string_length(src_string);
        
        // Do the actual the copying, this loop copies the null character!
        for(int i = 0; i <= l_src_string; i++)
        {
                dest_string[i] = src_string[i];
        }
}

/*
Copy does not require an initialized string in the object's buffer.
*/
void t_string::copy(const char* string)
{
        int l_src_string = string_length(string);
        if(this->obj_str_mem_size < l_src_string + 3)
        {
                free(this->obj_string);
                this->obj_string = (char*)malloc(l_src_string + 3);
                this->obj_str_mem_size = l_src_string + 3;
        }
        
        // Do the actual the copying, this loop copies the null character!
        for(int i = 0; i <= l_src_string; i++)
        {
                this->x(i) = string[i];
        }
}

void t_string::copy(t_string* string)
{
        this->copy(string->obj_string);
}

// Static string library functions: These are thread safe functions.
int t_string::string_length(const char* string)
{
        int str_length = 0;
        while(string[str_length] != 0)
        {
                str_length++;
        }

        return(str_length);
}

// Static string library functions: These are thread safe functions.
int t_string::string_length(t_string* string)
{
        return(string->length());
}

int t_string::length()
{
        return(string_length(this->obj_string));
}

// empty the string.
void t_string::empty()
{
        this->x(0) = 0;
}

// Getter function.
char& t_string::x(int i)
{
        return(this->obj_string[i]);
}

char* t_string::str()
{
        return(this->obj_string);
}

// Return the substring at substring starting at i and ending at j, inclusive.
char* t_string::substring(int i, int j)
{
        //t_string* sub_str = new t_string();
        char* sub_str;

        // Validity checks on the arguments.
        if(i > j ||
                i > this->length() || 
                j > this->length())
        {
                return(NULL);
        }
        else
        {
                sub_str = (char*)malloc(sizeof(char) * (j-i+1+2));
                sub_str[0] = 0; // Make this an empty string.
        }

        // Copy the substring.
        int ip = i;
        int i_sub_str = 0;
        while(ip <= j)
        {
                sub_str[i_sub_str] = this->x(ip);

                i_sub_str++;
                ip++;
        }

        // End the string.
        sub_str[i_sub_str] = 0;

        return(sub_str);
}

void t_string::replace_avoid_list(char* str, char* avoided_char_list, char char_to_replace)
{
        for(int i_str = 0; i_str < string_length(str); i_str++)
        {
                // Look for avoided chars. if found, replace.
                bool found_avoided_char = false;
                for(int i_avoided = 0; i_avoided < string_length(avoided_char_list); i_avoided++)
                {
                        if(str[i_str] == avoided_char_list[i_avoided])
                        {
                                found_avoided_char = true;
                        }
                } // i_avoided loop

                if(found_avoided_char)
                {
                        str[i_str] = char_to_replace;
                }
        } // i_str loop
}

// Tokenizer: Tokenizes with respect to the characters in the delimiter list, which is a null terminated list.
t_string_tokens* t_string::tokenize_by_chars(const char* delimiter_list)
{
        //printf("String: %sFIN\n", this->str());
        t_string_tokens* token_vector = new vector<t_string*>();

        // Current token is nothing at the beginning.
        t_string* current_token = new t_string();

        // Go over the obj_string.
        int str_len = string_length(this->str());
        for(int i_str = 0; i_str < str_len; i_str++)
        {
                //printf("Processing: %c\n", this->x(i_str));
                // This is set to true if i_str'th character is a delimiter.
                bool is_delimiter = false;

                // Check if current character in the string is a delimiter character.
                for(int i_del = 0; i_del < string_length(delimiter_list); i_del++)
                {
                        if(this->x(i_str) == delimiter_list[i_del])
                        {
                                is_delimiter = true;

                                // Push current token to vector (if it is appropriate for pushing.) and start a new token.
                                // Push the previous token. Cannot push the next token since there is no way of knowing that
                                // it is a valid token beforehand.
                                if(current_token != NULL && current_token->length() != 0)
                                {
                                        // Push this token to token vector.
                                        //printf("pushing %s\n", current_token->obj_string);
                                        token_vector->push_back(current_token);

                                        // Do not free current token now, allocate a new one.
                                        current_token = new t_string();
                                }
                                else // The token is not a good token (i.e., an empty token), free its memory and 
                                {
                                        current_token->empty();
                                }

                                break; // Break from the delimiter char loop.
                        }
                } // delimiter char loop.

                // If this character is not a delimiter, concatenate the character to the current_token.
                if(!is_delimiter)
                {
                        current_token->concat_char(this->x(i_str));
                }

        } // string char loop.

        // The last token is added to the tokens vector after whole loop is finished.
        if(current_token != NULL && current_token->length() != 0)
        {
                // Push this token to token vector.
                token_vector->push_back(current_token);
                //printf("pushing %s\n", current_token->obj_string);
        }
        else    
        {
                delete(current_token);
        }

        return(token_vector);
}

bool t_string::is_balanced(char* left_pars, char* right_pars)
{
        return(t_string::is_balanced(this->str(), left_pars, right_pars));
}

bool t_string::is_balanced(char* str, char* left_pars, char* right_pars)
{
        int n_pars = t_string::string_length(left_pars);

        // This is the list of left parenthesis characters in the string.
        char* left_par_chars = (char*)malloc(sizeof(char) * (t_string::string_length(str) + 1));

        // Count the left and right parenthesis characters, they should be equaal in number, watch the structure along the way.
        int l_str = t_string::string_length(str);
        int l_left_par_chars = 0;
        left_par_chars[l_left_par_chars] = 0; // Initialize.
        for(int i = 0; i < l_str; i++)
        {
                // Check all the parenthesis characters to check if this is a parenthesis character.
                for(int i_par = 0; i_par < n_pars; i_par++)
                {
                        // Check first for pairs that close a frame.
                        if(right_pars[i_par] == str[i])
                        {
                                // Make sure that this does close the last left parenthesis character.
                                if((l_left_par_chars < 1) ||
                                        (l_left_par_chars >= 1 && left_par_chars[l_left_par_chars-1] != left_pars[i_par]))
                                {
                                        // This character may actually be a left parenthesis, whose left and right characters are the same.
                                        if(left_pars[i_par] == right_pars[i_par])
                                        {
                                                // Interpret this as the left parenthesis characters.
                                                left_par_chars[l_left_par_chars] = str[i];
                                                l_left_par_chars++;
                                        }
                                        else
                                        {
                                                printf("Knotted parenthesis: Right parenthesis character %c is closing left parenthesis character %c.\n", right_pars[i_par], left_par_chars[l_left_par_chars-1]);
                                                return(false);
                                        }
                                }
                                else
                                {
                                        // Note that for the pairs of parenthesis characters whose left and right elements are the same, nesting is not possible.

                                        // Remove this left parenthesis character since its matching right parenthesis is found.
                                        left_par_chars[l_left_par_chars-1] = 0;
                                        l_left_par_chars--;
                                }
                        }
                        else if(left_pars[i_par] == str[i])
                        {
                                // Add this left parenthesis character to the end of list of left parenthesis characters encountered so far.
                                left_par_chars[l_left_par_chars] = str[i];
                                l_left_par_chars++;
                        }
                } // i_par loop.

                // Update par_cnt if this is a parenthesis character, depending on it being left or right parenthesis character.

                // If the last parenthesis character was a left one, the new character cannot be a right one and vice versa. (No knotted pairs.)
        } // i loop

        if(l_left_par_chars == 0)
        {
                return(true);
        }
        else
        {
                return(false);
        }
}

void t_string::clean_tokens(t_string_tokens* tokens)
{
        for(int i = 0; i < (int)tokens->size(); i++)
        {
                delete(tokens->at(i));
        }

        tokens->clear();
        delete(tokens);
}

void t_string::remove_beginning_spaces()
{
        int i_last_space = 0;

        for(int i = 0; i < this->length(); i++)
        {
                if(this->x(i) != ' ')
                {
                        i_last_space = i;
                        break;
                }
        }

        char* reloc_copy = (char*)malloc(sizeof(char*) * (this->length() + 2));

        copy(reloc_copy, this->str());

        // Relocate.
        copy(this->str(), &reloc_copy[i_last_space]);
        free(reloc_copy);
}

// Tokenizer: Tokenizes with respect to the characters in the delimiter list, which is a null terminated list.
t_string_tokens* t_string::tokenize_by_str(const char* delimiter_string)
{
        t_string_tokens* token_vector = new vector<t_string*>();

        // Current token is nothing at the beginning.
        t_string* current_token = new t_string();

        // Go over the obj_string.
        int str_len = this->length();
        for(int i_str = 0; i_str < str_len; i_str++)
        {
                //printf("Processing: %c\n", this->x(i_str));
                // This is set to true if i_str'th character is a delimiter.
                bool is_delimiter = false;

                int i_str_search = i_str;
                // Check if current character in the string is a delimiter character.
                for(int i_del = 0; i_del < string_length(delimiter_string); i_del++)
                {
                        if(i_str_search == this->length() || this->x(i_str_search) != delimiter_string[i_del])
                        {
                                is_delimiter = false;
                                break;
                        }                               

                        // If the current substring matches the delimiter string, set delimiter to true.
                        if(i_del == string_length(delimiter_string) - 1)
                        {
                                // Push current token to vector (if it is appropriate for pushing.) and start a new token.
                                // Push the previous token. Cannot push the next token since there is no way of knowing that
                                // it is a valid token beforehand.
                                if(current_token != NULL && current_token->length() != 0)
                                {
                                        // Push this token to token vector.
                                        //printf("pushing %s\n", current_token->obj_string);
                                        token_vector->push_back(current_token);

                                        // Do not free current token now, allocate a new one.
                                        current_token = new t_string();
                                }
                                else // The token is not a good token (i.e., an empty token), free its memory and 
                                {
                                        current_token->empty();
                                }

                                is_delimiter = true;

                                // Note that the last char that is at i_str_search is the last char
                                // of delimiter string but that will be jumped over by for loop.
                                // This is the same in tokenize_by_chars function.
                                i_str = i_str_search;
                                break;
                        }

                        // Increment substring counter.
                        i_str_search++;

                } // delimiter char loop.

                // If this character is not a delimiter, concatenate the character to the current_token.
                if(!is_delimiter)
                {
                        current_token->concat_char(this->x(i_str));
                }

        } // string char loop.

        // The last token is added to the tokens vector after whole loop is finished.
        if(current_token != NULL && current_token->length() != 0)
        {
                // Push this token to token vector.
                token_vector->push_back(current_token);
                //printf("pushing %s\n", current_token->obj_string);
        }
        else
        {
                delete(current_token);
        }

        return(token_vector);
}

// Variable argumented sprintf: Print the string of this object using fmt_string as formatter string.
void t_string::sprintf(char* fmt_string, ...)
{
        va_list pars;
        va_start(pars, fmt_string);

        // Initialize string as having nothing.
        this->copy("");

        // Go over the string: Every time a formatting string is found by %.., it is processed.
        for(int i_str = 0; i_str < string_length(fmt_string); i_str++)
        {
                if(fmt_string[i_str] == '%')
                {
                        // Check the formatter character.
                        i_str++;

                        if(fmt_string[i_str] == 'd')
                        {
                                int _num = va_arg(pars,int);
                                this->concat_int(_num);
                        }
                        else if(fmt_string[i_str] == 'c')
                        {
                                int _char = va_arg(pars,int);
                                this->concat_char(_char);
                        }
                        else if(fmt_string[i_str] == 's')
                        {
                                char* _str = (char*)va_arg(pars, void*);
                                this->concat_string(_str);
                        }
                        else if(fmt_string[i_str] == '%')
                        {
                                this->concat_char(fmt_string[i_str]);
                        }

                        // The for loop jumps over formatter character.
                        //i_str++;
                }
                else // This is a regular character, just copy it.
                {
                        this->concat_char(fmt_string[i_str]);
                }
        }

        va_end(pars);
}

void t_string::concat_char(char _char)
{
        // A memory check.
        int str_len = this->length();
        while(this->obj_str_mem_size <= (str_len + 10))
        {
                //printf("Doubling buffer size.\n");
                //this->obj_str_mem_size = (this->length() + 10);
                this->obj_str_mem_size *= 2;
                char* temp_buf = this->obj_string;
                this->obj_string = (char*)malloc(this->obj_str_mem_size);
                this->copy(temp_buf);
                free(temp_buf);
        }

        int unconcat_length = str_len;
        this->obj_string[unconcat_length] = _char;

        // This is vert important. Finish the string!
        this->obj_string[unconcat_length + 1] = 0;
}

void t_string::concat_string(char* string)
{
        // Add all the chars in the string.
        int str_len = string_length(string);
        for(int i = 0; i <= str_len; i++)
        {
                this->concat_char(string[i]);
        }
}

void t_string::concat_string(t_string* string)
{
        this->concat_string(string->obj_string);
}

void t_string::concat_int(int i_num)
{
        t_string* num_str = num2str(i_num, 10);
        this->concat_string(num_str);
        delete(num_str);
}

void t_string::concat_float(double f_num)
{
}

// Take the reverse of the string in this object's string buffer.
void t_string::revert()
{
        t_string* temp_buf = new t_string(this->str());
        int str_len = temp_buf->length();
        for(int i = 0; i < str_len; i++)
        {
                this->x(i) = temp_buf->x(str_len - i - 1);
        }

        delete(temp_buf);
}

bool t_string::compare(t_string* string)
{
        return(this->compare(string->str()));
}

bool t_string::compare_strings(t_string* str1, t_string* str2)
{
        return(compare_strings(str1->str(), str2->str())); 
}

bool t_string::compare_strings(char* str1, char* str2)
{
        int str_len1 = string_length(str1);
        int str_len2 = string_length(str2);
        if(str_len1 != str_len2)
        {
                return(false);
        }

        for(int i = 0; i < str_len1; i++)
        {
                if(str1[i] != str2[i])
                {
                        return(false);
                }
        }

        return(true);
}

bool t_string::compare_strings_ci(t_string* str1, t_string* str2)
{
        return(compare_strings_ci(str1->str(), str2->str()));
}

bool t_string::compare_strings_ci(char* str1, char* str2)
{
        int str_len1 = string_length(str1);
        int str_len2 = string_length(str2);

        if(str_len1 != str_len2)
        {
                return(false);
        }

        for(int i = 0; i < string_length(str1); i++)
        {
                if(toupper(str1[i]) != toupper(str2[i]))
                {
                        return(false);
                }
        }

        return(true);
}

bool t_string::compare_ci(char* string)
{
        return(compare_strings_ci(this->str(), string));
}

bool t_string::compare_ci(t_string* string)
{
        return(compare_strings_ci(this->str(), string->str()));
}


bool t_string::compare(char* string)
{
        return(compare_strings(this->str(), string));
}

bool t_string::starts_with(char* string)
{
        if(this->length() < string_length(string))
        {
                return(false);
        }

        for(int i_str = 0; i_str < string_length(string); i_str++)
        {
                if(this->x(i_str) != string[i_str])
                {
                        return(false);
                }
        }       

        return(true);
}

bool t_string::starts_with(t_string* string)
{
        return(this->starts_with(string->str()));
}


t_string* t_string::num2str(int num, int base)
{
        t_string* num_str = new t_string();
        int divident = num;
        int residual = divident % base;

        do
        {
                num_str->concat_char((char)(residual + 48)); // Convert to ascii value.         
                divident = divident / base; // Shift right.
                residual = divident % base; // Compute the last digit.

                if(residual > 9)
                {
                        printf("The residual greater than 9!\n");
                }
        }
        while(divident != 0);

        // Reverse the num string.
        num_str->revert();

        return(num_str);
}

int t_string::str2num(char* num_str, int base)
{
        int num = 0;
        int extended_base = 1;
        for(int i = (int)string_length(num_str)-1; i >= 0; i--)
        {
                int cur_digit = (int)num_str[i];
                if(cur_digit >= '0' && cur_digit <= '9')
                        num += extended_base * (cur_digit - (int)'0');
                else if(cur_digit >= 'A' && cur_digit <= 'F')
                         num += extended_base * (cur_digit - (int)'A' + 10);
                else if(cur_digit >= 'a' && cur_digit <= 'f')
                         num += extended_base * (cur_digit - (int)'a' + 10);
                else
                {
                        printf("Could not resolve character as number in %s for base %d\n", num_str, base);
                        exit(0);
                }

                //printf("%d = %d * (%d - %d)\n", num, extended_base, cur_digit, 48);
                extended_base *= base;
        }

        return(num);
}

void t_string::to_upper()
{
        t_string::to_upper(this->str());
}

void t_string::to_upper(char* string)
{
        printf("%s->", string);
        //int upper_diff = (int)'A' - (int)'a';

        int str_len = string_length(string); 
        for(int i = 0; i < str_len; i++)
        {
                if(!(string[i] <= 'Z' &&  string[i] >= 'A') && !(string[i] <= 'z' && string[i] >= 'a'))
                {
                }
                else if(string[i] <= 'z' && string[i] >= 'a')
                {
                        string[i] = string[i] + ((int)'A' - (int)'a');
                }
                else
                {
                }
        }
        printf("%s\n", string);

}       
int t_string::str2num(t_string* num_str, int base)
{
        return(str2num(num_str->str(), base));
}



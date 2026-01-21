#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <cstdlib>
#include "io.h"
#include <iterator>
using std::string;
using std::vector;
using std::cout;
using std::endl;

const options op;

//some functions for internal use
std::vector<std::string> split_on(std::string tobesplit, const char* delims);
std::vector<int> _read_ct(std::string name, const int sz);
std::vector<std::vector<std::string> >tokenize_file(std::ifstream& infile);
//read the whole file at infile_location into a string
std::string slurp(std::ifstream& infile);

//some handy specializations of split_on_delim
std::vector<std::string> split_on_newline(std::string input);
std::vector<std::string> split_on_tab(std::string tobesplit);

constraints IO::read_constraints(string name, const int sz){
	std::ifstream infile(name,std::ifstream::in);
	if(!infile.is_open()) std::cout <<"file "<<name<<" not found!"<<std::endl;
	//string data = slurp(infile);
    vector<vector<string>> tok = tokenize_file(infile);
    if(sz==0) return constraints();
    vector<int> cst(sz);
    for(int k=0;k<sz;k++){
        cst[k] = k;
    }
    for(size_t k=0;k<tok.size();k++){
        const int i = atoi(tok[k][0].c_str());
        const int j = atoi(tok[k][1].c_str());
        cst[j] = i;
        cst[i] = j;
    }
    return constraints(cst);
}

constraints from_db_string(const string& db){
  vector<int> cst(db.size());
  for(size_t k=0;k<db.size();k++){
    cst[k] = k;
  }
  vector<int> stack;
  for(size_t i=0;i<db.size();i++){
    if (db[i] == '(' ){
      stack.push_back(i);
    }
    else if (db[i] == ')' ){
      if(stack.empty()) throw string("unbalanced parens in dot-bracket string");
      int partner = stack.back();
      stack.pop_back();
      cst[i] = partner;
      cst[partner] = i;
    }
    else if (db[i] == '.' ){
      continue;
    }
    else {
      char error[50];
      sprintf(error, "unrecognized character in dot-bracket string: %c", db[i]);
      throw string(error);
    }
  }
  if (!stack.empty()){
    throw string("unbalanced parens in dot-bracket string");
  }
  return constraints(cst);
}

constraints IO::read_unpairing_constraints(string name, const int sz){
	std::ifstream infile(name,std::ifstream::in);
	if(!infile.is_open()) std::cout <<"file "<<name<<" not found!"<<std::endl;
	//string data = slurp(infile);
    vector<vector<string>> tok = tokenize_file(infile);
    if(sz==0) return constraints();
    vector<bool> cst(sz,false);
    for(size_t k=0;k<tok.size();k++){
        const int i = atoi(tok[k][0].c_str());
        cst[i] = true;
    }
    return constraints(cst);
}

vector<int> _read_ct(string name, const int sz){
	std::ifstream infile(name,std::ifstream::in);
	if(!infile.is_open()){
        std::cout <<"file "<<name<<" not found!"<<std::endl;
        throw "file not found";
    }
    vector<vector<string>> tok = tokenize_file(infile);
    vector<int> cst(sz);
    for(int k=1;k<=sz;k++){
        const int i = atoi(tok[k][0].c_str());
        const int j = atoi(tok[k][4].c_str());
        if(j==0){
            cst[i-1] = i-1;
        }
        else {
            cst[j-1] = i-1;
            cst[i-1] = j-1;
        }
    }
	return cst;
}

constraints IO::read_ct(string name, const int sz){
    return constraints(_read_ct(name,sz));
}

constraints IO::read_unpairing_ct(string name, const int sz){
    vector<int> pairs = _read_ct(name,sz);
	vector<bool> ret;
	for(int i=0;i<(int)pairs.size();i++){
		ret.push_back(i!=pairs[i]);
	}
	return constraints(ret);
}


vector <vector <string> > IO::read_datafile(string name, const options& op)
{
	string path = op.parameter_files + name;
	std::ifstream infile (path,std::ifstream::in);
	if(!infile.is_open()){
		string errormessage = string("data file ")+path+string(" not found!\nset CYCLEFOLD_DATAPATH environment variable to /path/to/RNAstructure/CycleFold/datafiles");
		throw errormessage;
	}
	return tokenize_file(infile);
}

bool contains(string text, string pattern){
  return text.find(pattern) != std::string::npos;
}

bool match(char c, const string& cmp){
	for(string::const_iterator it=cmp.begin();it!=cmp.end();++it){
		if(c==*it) return true;
	}
	return false;
}

vector<string> split_on_delim(const string& text, const string& delims){
	//first, identify the position of every non-delimiter
	vector<size_t> positions;
	if(!text.empty() && match(text.at(0), delims)) positions.push_back(0);
	for(size_t i=0;i<text.size();i++){
		if(!match(text.at(i),delims)){
			positions.push_back(i);
		}
	}
	if(positions.size() == 0) return vector<string>();

	//then find stretches of adjacent non-delimiters
	vector<string> ret;
	vector<char> cs = {text.at(*positions.begin())};
	for(auto it=positions.begin()+1;it!=positions.end();++it){
		size_t i = *(it-1);
		size_t j = (*it);
		if(i+1!=j){
			ret.push_back(string(cs.begin(), cs.end()));
			cs.clear();
		}
		cs.push_back(text.at(j));
	}
	if(cs.size()>0) ret.push_back(string(cs.begin(),cs.end()));
	return ret;
}

string trim(const string& text, const string& whitespace){
	size_t lpos = text.find_first_not_of(whitespace);
	size_t rpos = text.find_last_not_of(whitespace);
	if(lpos==string::npos)
		return string("");
	return text.substr(lpos, rpos-lpos+1);
}


vector<sequence> IO::read_fasta(std::string name)
{
	std::ifstream infile(name,std::ifstream::in);
	if(!infile.is_open()) std::cout <<"file "<<name<<" not found!"<<std::endl;
//check if open!!!
	string data = slurp(infile);
	//std::cout<<data<<std::endl;
	vector<sequence> seqs;
	vector<string> temp = split_on_delim(data,string(">"));//split on > to get each sequence
	//first entry in temp is empty string, so we start at temp.begin()+1
	for (std::vector<string>::iterator it = temp.begin(); it != temp.end(); ++it)
	{
		vector<string> lines = split_on_newline(*it);
		lines[0] = trim(lines[0]);
		string tag = string(lines[0].begin()+1, lines[0].end());
		string seq="";
		for (std::vector<string>::iterator it2 = lines.begin()+1; it2 != lines.end(); ++it2)
		{
			*it2 = trim(*it2);
      if(contains(*it2, string(".")) || contains(*it2, string(")")) || contains(*it2, string("("))){
        continue;
      }
			seq+=*it2;
		}
		seqs.push_back(sequence(tag,seq));
	}
	return seqs;
}

vector<sequence> IO::read_dotseq(std::string name)
{
	std::ifstream infile(name,std::ifstream::in);
	if(!infile.is_open()) std::cout <<"file "<<name<<" not found!"<<std::endl;
	string data = slurp(infile);
	vector<string> lines = split_on_delim(data,string("\n"));//split on > to get each sequence
	//remove whitespace
	for(auto it = lines.begin(); it != lines.end(); ++it){
		*it = trim(*it);
	}
	auto it = lines.begin();
	while((*it).at(0) == ';'){
		++it;
	}
	if(it==lines.end()) throw "could not read seq file\n";
	string tag = *it;
	string seq = "";
	for(++it; it !=lines.end(); ++it){
		seq += *it;
	}
	if(seq.at(seq.size()-1) != '1') throw "couldn't interpret seq file: no '1' character at the end of the sequence\n";
	seq.erase(seq.end()-1);
	return { sequence(tag,seq) };
}

vector<constraints> IO::read_fasta_constraints(std::string name)
{
	std::ifstream infile(name,std::ifstream::in);
	if(!infile.is_open()) std::cout <<"file "<<name<<" not found!"<<std::endl;
  //check if open!!!
	string data = slurp(infile);
	//std::cout<<data<<std::endl;
	vector<sequence> seqs;
	vector<constraints> structures;
	vector<string> temp = split_on_delim(data, string(">"));//split on > to get each sequence
	//first entry in temp is empty string, so we start at temp.begin()+1
	for (std::vector<string>::iterator it = temp.begin(); it != temp.end(); ++it)
    {
      vector<string> lines = split_on_newline(*it);
	  lines[0] = trim(lines[0]);
      string tag = lines[0];
      string seq="";
      string structure="";
      for (std::vector<string>::iterator it2 = lines.begin()+1; it2 != lines.end(); ++it2)
        {
		  *it2 = trim(*it2);
          if(contains(*it2, string(".")) || contains(*it2, string(")")) || contains(*it2, string("(")))
            structure+=*it2;
          else {
            seq+=*it2;
          }
        }
      seqs.push_back(sequence(tag,seq));
	  constraints c = from_db_string(structure);
//	  cout<<"db string "<<structure;
//	  c.show();
      structures.push_back(c);
    }
	return structures;
}

vector<vector<string>> tokenize_file(std::ifstream& infile)
{
    vector<vector<string>> ret;
    string line;
    std::istringstream buffer;
    while(std::getline(infile, line)){
        buffer.clear();
        buffer.str(line);
        vector<string> tok{std::istream_iterator<string>(buffer),
                           std::istream_iterator<string>()};
        ret.push_back(tok);
    }
    return ret;
}

string slurp(std::ifstream& infile) {
    std::stringstream sstr;
    sstr << infile.rdbuf();
    return sstr.str();
}

vector<string> split_on_newline(string input)
{
	return split_on_delim(input, string("\n"));
}

vector<string> split_on_tab(string tobesplit)
{
	return split_on_delim(tobesplit, string("\n"));
}

vector<string> split_on(string tobesplit, const char* delims)
{
	return split_on_delim(tobesplit, string(delims));
}

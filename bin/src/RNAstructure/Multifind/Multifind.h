#include "../RNA_class/Multifind_object.h"
#include "../src/ParseCommandLine.h"

#include <vector>
#include <string>
#include <sstream>

class Multifind_Interface{
 public:
  

  Multifind_Interface();
  
  bool parse(int argc,char** argv);
  void read_fasta();
  void run();
  string compact(string& sequence);
 private:
  int processors;
  int structure_output;
  string fasta_file;
  string multifind_file;
  int number_seq;
  vector<string> ct_files;
  vector<string> input_alignment;
  vector<string> input_sequences;
};


#ifndef CMD_OPTIONS_H
#define CMD_OPTIONS_H
#include <vector>
#include <string>

class options{
 public:
  options();
  std::vector<std::string> input_files;//input file paths
  std::string output_file;
  std::string parameter_files;//directory with parameters
  std::string constraint_ct;
  std::string constraint_file;
	bool partition;
	bool maxexpect;
  bool turbo;
  bool set;//indicates valid options set and calculation can proceed
  int iter;
  double gamma;
  bool bigloops;
  bool unpair;
  bool fasta_constraints;
  bool dotseq_format;

};

options parse(int ac, char* av[]);

namespace parser{
	options parseCommandLine(int ac,char* av[]);
}

#endif //CMD_OPTIONS_H

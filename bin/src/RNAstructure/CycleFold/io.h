#ifndef __FILEREADER_INCLUDE__
#define __FILEREADER_INCLUDE__

#include "options.h"
#include "sequence.h"
#include "constraints.h"
#include <string>
#include <vector>
#include <fstream>

namespace IO {
  std::vector<std::vector<std::string> > read_datafile(std::string infile_name, const options& op);
  std::vector<sequence> read_fasta(std::string infile_location);
  std::vector<sequence> read_dotseq(std::string infile_location);
  std::vector<constraints> read_fasta_constraints(std::string infile_location);
  constraints read_constraints(std::string infile_location,const int sz);
  constraints read_unpairing_constraints(std::string infile_location,const int sz);
  constraints read_ct(std::string infile_location,const int sz);
  constraints read_unpairing_ct(std::string infile_location,const int sz);
}
std::vector<std::string> split_on_delim(const std::string& text, const std::string& delims);
std::string trim(const std::string& text, const std::string& whitespace=std::string(" \r\r\n"));

#endif //__FILEREADER_INCLUDE__

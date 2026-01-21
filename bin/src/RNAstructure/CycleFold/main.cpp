#ifndef TEST
#include "mainloop.h"
#include "options.h"
#include "arrays.h"
#include "io.h"
#include "NCM_parameters.h"
#include "constraints.h"
#include "maxexpect.h"
#include "turbo_calculation.h"
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <vector>

using std::vector;
using std::cout;
using std::pair;

void prettyprint(std::vector< std::pair<int,int> > input, const int   n){
  std::vector<char> output(n,'.');
  for(auto it : input){
    output[it.first] = '(';
    output[it.second] = ')';
    printf("%d %d \n", it.first, it.second);
  }
  for (auto i: output)
    std::cout << i;
  std::cout<<"\n";
}

int main(int ac,char* av[]){
	try{
		options op = parse(ac,av);
		bool mfe = !(op.turbo || op.partition || op.maxexpect);
		// mfe = false;
		bool turbo = op.turbo;
	  constraints rst = constraints();

	  //Store the sequence information in class sequence.  Use a vector to read multiple sequences.
      vector<sequence> seqs = op.dotseq_format? IO::read_dotseq(op.input_files[0]) : IO::read_fasta(op.input_files[0]);

	  //Store a single sequence, the first sequence, for use in a deafult calculation (not TurboFold)
      sequence seq =  seqs[0];
	  vector<constraints> c;
	  if (op.fasta_constraints){
		c = IO::read_fasta_constraints(op.input_files[0]);
	  }

	  if (op.constraint_ct != ""){
			if(!op.unpair)
				rst = IO::read_ct(op.constraint_ct.c_str(), seq.getLength());
			else
				rst = IO::read_unpairing_ct(op.constraint_ct.c_str(), seq.getLength());
	  }
	  if (op.constraint_file != ""){
			if(!op.unpair)
				rst = IO::read_constraints(op.constraint_file.c_str(), seq.getLength());
			else
				rst = IO::read_unpairing_constraints(op.constraint_file.c_str(), seq.getLength());
	  }
	  //rst.show();
	  rst.build_allowedpairs(seq.getLength());
		if(mfe) {
			parameters<int> par = parameters<int>(op);
		for(sequence& ss : seqs)
		  calculate_mfe(ss, par, rst, op);
			return 0;
		}
		else if (turbo) {
			parameters<real_t> par = parameters<real_t>(op);

			//create a vector of table_t, pair probabilities for a single sequence.
			vector<table_t> turbo = turbo_calculation(seqs, par, op.iter, op.gamma, c);
			for(int i=0;i<(int)seqs.size();i++){
				//cout<<seqs[i].getTag()<<std::endl;
				if(op.partition){
					show_pair_probs_arrayview(turbo[i]);
				}
				else {
				  show_ct(maxexpect(turbo[i]), seqs[i]);
			  }
			}
		//show_ct(maxexpect(turbo[0]), seqs[0]);
		}
		else {
			parameters<real_t> par = parameters<real_t>(op);
			table_t pair_probs = calculate_pairing_probabilities(seq,par,rst);
			if(op.partition) {
				show_pair_probs_arrayview(pair_probs);
			}
			if(op.maxexpect) {
				vector<pair<int,int>> pairs = maxexpect(pair_probs);
				show_ct(pairs,seq);
			}
			return 0;
		}
	}
	catch(const char* oops){
		cout<<oops<<"\n";
		return 1;
	}
	catch(std::string oops){
		cout<<oops<<"\n";
		return 1;
	}
	catch (int code) { // parse will throw 0 for the -h flag.
		return code;
	}
}

#else //compile for testing
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "test.h"
#endif //TEST

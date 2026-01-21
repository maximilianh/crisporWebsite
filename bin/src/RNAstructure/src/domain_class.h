#include "../RNA_class/RNA.h"
#include "common_utils.h"

#include "defines.h"
#include <iostream>

struct pair_data{
	int i,j;
	double prob;
};

class domain_calc {
 	public:

		domain_calc(bool circularize, double outside_weight, int leng);
		~domain_calc();

		int calc_domains();

		/*
		Function Name;	fill_domain_matrices

		Function:
			Fills the dynamic progranning arrays used to calculate the inner and total base pairing probabilities

		Arguments:
			None
		*/
		void fill_domain_matrices();

		/*
		Function Name;	writeconstraints

		Function:
			Writes a file with the optimum domain segmentation.

		Arguments:
			std::string filename:	The path to the constraint file to be written
		*/
		void writeconstraints(std::string filename);

		/*
		Function Name;	read_pp_file

		Function:
			Reads a pair probability file that is output from ProbabilityPlot or partition-cuda.  It needs a sequence file to provide the nucleotide sequence.

		Arguments:
			std::string filename:	The path to the pair probability file
		*/
		void read_pp_file(std::string filename);

		/*
		Function Name;	read_pfs_file

		Function:
			Imports pair probability information when reading information from a pfs file.

		Arguments:
			RNA *strand:	An RNA object.  This strand object has already been filled with 
							information from the pfs file.
		*/
		void read_pfs_file(RNA *strand);

		/*
		Function Name;	complete_domain_assignment

		Function:
			Finds optimum domain segmentation, including circular domains

		Arguments:
			int min_domain_length:	The minimum allowed domain size.
		*/
		void complete_domain_assignment(int min_domain_length);

		/*
		Function Name;	calc_domain_scores

		Function:
			Calculates all possible domain scores. The results are used to find the optimum domain
			segmentation.

		Arguments:
		*/
		void calc_domain_scores();

		/*
		Function Name;	copy_constraints

		Function:
			Copies the optimum domain segmentation as folding constraints in the structure object.
			Currently, circular domains are not added, but are implictly considered.

		Arguments:
			structure *ct:	Pointer to the structure object that the domains will be copied to.
		*/
		void copy_constraints(structure *ct);

		/*
		Function Name;	score

		Function:
			Returns the score for the domain between i and j nucleotides.  If i > j, then the score
			for the domain that spans the 5' and 3' end of the sequence (the circular domain) is 
			returned.

		Arguments:
			i:  Integer index of the 5' most nucleotide in the domain
			j:	Integer index of the 3' most nucleotide in the domain
		*/
		double score(int i, int j);
	
	private:

 		//base_pairs lists all possible base pairs, with their pairing probability, listed by their index
		vector<pair_data> base_pairs;

		//The sequence of the RNA to be segmented
		const char *seq;

		//The length of the RNA, in nucleotides
		int length;

		int domain_phase;

		//A number of dynamic programming arrays used by the domain findng algorithm
		DynProgArrayU<double> *prob_matrix, *row_sum, *col_sum, *inner, *total, *scores, *scores2;
		
		//constraints contains a number of i,j pairs, where (i,j) defines a base pair closing a folding domain
		vector<int> constraints5;
		vector<int> constraints3;

		bool circ;

		double w2;
};


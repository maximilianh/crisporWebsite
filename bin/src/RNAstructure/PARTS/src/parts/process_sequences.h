#ifndef _PROC_SEQS_
#define _PROC_SEQS_

#define MAX_SEQ_PATH (500)

// Do sequence data validation, 
// Do sequence data querying,
// Couple raw sequence data and internal recursion, allocation, etc. code.
// This part is the interface between sequence inputs and other parts of the code. Whenever a part of code needs to access
// sequence data (e.g. query for length of sequences) it refers to this class, sequences are only processed here.
// For instance, 
// 	do randomization of unknown nucleotides here.
//	do correction for indexing here.
// Encapsulates 2 structure structs.
class t_structure;
class t_ppf_cli;

// Sequence manager class.
class t_seq_man
{
	public:
	// Constructor: Do input validation and creation of structures.
	t_seq_man(t_ppf_cli* _ppf_cli); // Construct
	~t_seq_man(); // Destruct!!!

	// File paths of 2 sequences.
	char seq1_fp[MAX_SEQ_PATH];
	char seq2_fp[MAX_SEQ_PATH];

	t_ppf_cli* ppf_cli;

	// Encapsulated sequences.
	t_structure* seq1;
	t_structure* seq2;

	// Get lengths of sequences.
	int get_l_seq1();
	int get_l_seq2();

	// Get nucleotides from sequences, note that those are after 
	// Ambiguous nucleotide corrections (by replacing with random nucleotides).
	// Returns 0, 1, 2, 3 corresponding to A, C, G, U.
	// index argument start from 1 to N1 and N2 (including N1 and N2).
	char get_nuc_seq1(int index);
	char get_nuc_seq2(int index);
};

#endif 

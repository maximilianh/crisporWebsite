#ifndef _TEMPLATE_PF_ARRAY_
#define _TEMPLATE_PF_ARRAY_

class t_seq_man; // Sequence manager to interface with sequence data.
class t_spf_array; // Single partition function array.
class t_aln_priors;
class t_ppf_cli;

class t_template_pf_array
{
	public:
	t_template_pf_array(t_seq_man* _seq_man, bool mallocate); // Constructor.

	// Constructor, seq1_fold_env and seq2_fold_env are coming from coming from spf array, should not be 
	// played with by template array class.
	// Dont touch, only watch principle.
	t_template_pf_array(t_seq_man* seq_man, 
						short** _seq1_ptr_rel_map,
						short** _seq2_ptr_rel_map,
						bool mallocate);

	~t_template_pf_array(); // Destructor.

	t_seq_man* seq_man;

	int N1;
	int N2;

	static int n_arrays;

	double n_alloced_bytes;

	double zero;

	// This is the essential array which is calculated.
	double**** pf_array;

	// fold envelopes, if there are any.
	// Use those for allocation only, the pairibility is handled by spf classes,
	// not by these classes. Actually only seq1_fold_env is enough for allocation.
	bool** seq1_fold_env; // template_pf_array only uses fold envelope of first sequence for allocation of arrays.
	bool** seq2_fold_env; // template_pf_array only uses fold envelope of first sequence for allocation of arrays.
	bool** seq1_str_coinc_env; // template_pf_array only uses fold envelope of first sequence for allocation of arrays.
	bool** seq2_str_coinc_env; // template_pf_array only uses fold envelope of first sequence for allocation of arrays.

	short** seq1_ptr_rel_map;
	short** seq2_ptr_rel_map;

	// Return number of bytes allocated.
	void alloc_init_pf_array(t_seq_man* seq_man, bool mallocate);	

	// Access V pf array using sequence indices.
	double& x(int i, int j, int k, int l);

	// Following are for managing loop limits for V coming from alignment module.
	static int* low_limits; // v_low_limit[i] = lowest sequence index for k. those are sequence indices in 2nd sequence.
	static int* high_limits; // v_high_limit[i] = highest sequence index for k. those are sequence indices in 2nd sequence.
	static void alloc_init_loop_limits(t_seq_man* seq_man); // Sets loop limits as if no alignment envelope exists.
	static void alloc_init_loop_limits(t_seq_man* seq_man, int* _low_limits, int* _high_limits);
	static void update_loop_limits(t_seq_man* seq_man, int* _low_limits, int* _high_limits);
	static void update_loop_limits(t_seq_man* seq_man, char* lls_fp);
	bool check_boundary(int i1, int i2); // Checks if i2, seq2 sequence index, is in loop limit of i1, seqwence index of 1st sequence.

	// Check accessibility of indices for general cases.
	bool check_4D_ll(int i, int j, int k, int l);
};

#endif // _TEMPLATE_PF_ARRAY_

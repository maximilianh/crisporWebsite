#ifndef _W_
#define _W_

/*
W Class: Handles W array for pairwise partition function with pair scoring.
Array is the partition function for subsequences with common secondary structures
where i-j is a pair and k-l is a pair and i~k and j~l where - corresponds to pairing 
and ~ corresponds to alignment.

Definition this object from implementation point of view (i.e. what a w object can do and what it cannot do):

Loop limits for w array are handled from that class.
*/

// Forward declare the objects that W depends on calculation.
class t_seq_man; // Sequence manager to interface with sequence data.
class t_ppf_WL;
class t_ppf_WR;
class t_ppf_V_mhe;
class t_ppf_V;
class t_aln_priors;
class t_spf_array;
class t_template_pf_array;
class t_ppf_tb_stack;
class t_MAP_alignment;
class t_stoch_sampled_str_aln_sample_set;
class t_stoch_sampling_math;
class t_ppf_loops;

class t_ppf_W
{
public:
	t_ppf_W(t_ppf_loops* _ppf_loops); // Constructor.
	~t_ppf_W(); // Destructor.

	int N1;
	int N2;

	t_ppf_loops* ppf_loops;

	t_aln_priors* aln_priors;
	t_spf_array* seq1_spf;
	t_spf_array* seq2_spf;

	t_template_pf_array* pf_array;
	t_template_pf_array* ext_pf_array;

	void alloc_init_ppf_W_pf_array(t_seq_man* seq_man);

	// Main calculation unit for W, calculate one cell of 4D pf array.
	void calculate_W(int i, int j, int k, int l, bool bt = false); // Arrays to use.
	void calculate_ext_dependencies(int i, int j, int k, int l); // Arrays to use.

	// Access W pf array using sequence indices.
	double& x_setter(int i, int j, int k, int l);
	double& x_ext(int i, int j, int k, int l);
	double x(int i, int j, int k, int l);

	// Following are for managing loop limits for W coming from alignment module.
	bool check_boundary(int i1, int i2); // Checks if i2, seq2 sequence index, is in loop limit of i1, seqwence index of 1st sequence.
};

#endif // _W_

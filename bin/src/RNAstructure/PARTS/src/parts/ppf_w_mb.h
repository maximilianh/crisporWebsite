#ifndef _WMB_
#define _WMB_

/*
W Class: Handles W array for pairwise partition function with pair scoring.
Array is the partition function for subsequences with common secondary structures
where i-j is a pair and k-l is a pair and i~k and j~l where - corresponds to pairing 
and ~ corresponds to alignment.

Definition this object from implementation point of view (i.e. what a w object can do and what it cannot do):

Loop limits for w array are handled from that class.
*/

// Forward declare the objects that WMB depends on calculation.
class t_seq_man; // Sequence manager to interface with sequence data.
class t_ppf_WMBL; 
class t_aln_priors;
class t_spf_array;
class t_template_pf_array;
class t_ppf_tb_stack;
class t_MAP_alignment;
class t_stoch_sampled_str_aln_sample_set;
class t_stoch_sampling_math;
class t_ppf_loops;

class t_ppf_WMB
{
	public:
	t_ppf_WMB(t_ppf_loops* _ppf_loops); // Alignment prior scores.

	~t_ppf_WMB(); // Destructor.

	t_ppf_loops* ppf_loops;

	int N1;
	int N2;

	t_aln_priors* aln_priors;
	t_spf_array* seq1_spf;
	t_spf_array* seq2_spf;

	t_template_pf_array* pf_array;
	t_template_pf_array* ext_pf_array;

	void alloc_init_ppf_WMB_pf_array(t_seq_man* seq_man);

	// Main calculation unit for WMB, calculate one cell of 4D pf array.
	void calculate_WMB(int i, int j, int k, int l, bool bt = false); // Arrays to use.
	void calculate_ext_dependencies(int i, int j, int k, int l); // Arrays to use.

	// Access WMB pf array using sequence indices.
	double& x_setter(int i, int j, int k, int l);
	double x(int i, int j, int k, int l);
	double& x_ext(int i, int j, int k, int l);

	// Following are for managing loop limits for WMB coming from alignment module.
	bool check_boundary(int i1, int i2); // Checks if i2, seq2 sequence index, is in loop limit of i1, seqwence index of 1st sequence.
};

#endif // _WMB_

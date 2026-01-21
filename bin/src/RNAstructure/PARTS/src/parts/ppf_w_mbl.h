#ifndef _WMBL_
#define _WMBL_

/*
class t_ppf_WMBL: Enumerates scores of structures that have 5' unpaired nucleotides only, 
that is, 3' ends are always paired to an inner nucleotide in subsequence. This is 
formed by doing WMB-V or W-V concatenations.
*/

// Forward declare the objects that WMBL depends on calculation.
class t_seq_man; // Sequence manager to interface with sequence data.
class t_ppf_V; // Single partition function array.
class t_ppf_W;
class t_ppf_WL;
class t_ppf_WMB;
class t_spf_array;
class t_aln_priors;
class t_template_pf_array;
class t_ppf_tb_stack;
class t_stoch_sampled_str_aln_sample_set;
class t_stoch_sampling_math;

class t_ppf_WMBL
{
	public:
	t_ppf_WMBL(t_ppf_loops* _ppf_loops); // Constructor.

	t_ppf_loops* ppf_loops;

	~t_ppf_WMBL(); // Destructor.

	int N1;
	int N2;

	t_aln_priors* aln_priors;
	t_spf_array* seq1_spf;
	t_spf_array* seq2_spf;

	t_template_pf_array* pf_array;
	t_template_pf_array* ext_pf_array;

	// Main calculation unit for WMB, calculate one cell of 4D pf array.
	void calculate_WMBL(int i, int j, int k, int l, bool bt = false); // Arrays to use.
	void calculate_ext_dependencies(int i, int j, int k, int l); // Arrays to use.

	// Access WMB pf array using sequence indices.
	double& x_setter(int i, int j, int k, int l);
	double x(int i, int j, int k, int l);
	double& x_ext(int i, int j, int k, int l);

	// Following are for managing loop limits for WMB coming from alignment module.
	bool check_boundary(int i1, int i2); // Checks if i2, seq2 sequence index, is in loop limit of i1, seqwence index of 1st sequence.
};

#endif // _WMBL_

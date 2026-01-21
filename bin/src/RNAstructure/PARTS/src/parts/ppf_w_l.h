#ifndef _WL_
#define _WL_

/*
WL Class:  A closed structure with only possible unpaired bases at 5' side.
*/

// Forward declare the objects that WL depends on calculation.
class t_seq_man; // Sequence manager to interface with sequence data.
class t_ppf_V;
class t_ppf_V_bpi;
class t_ppf_V_bp_aln_up;
class t_aln_priors;
class t_spf_array;
class t_template_pf_array;
class t_ppf_W_bpi;
class t_ppf_tb_stack;
class t_ppf_V_mhe;
class t_MAP_alignment;
class t_stoch_sampled_str_aln_sample_set;
class t_stoch_sampling_math;
class t_ppf_loops;

class t_ppf_WL
{
public:
	t_ppf_WL(t_ppf_loops* _ppf_loops); // Constructor.
	~t_ppf_WL(); // Destructor.

	t_ppf_loops* ppf_loops;

	int N1;
	int N2;

	t_aln_priors* aln_priors;
	t_spf_array* seq1_spf;
	t_spf_array* seq2_spf;

	t_template_pf_array* pf_array;
	t_template_pf_array* ext_pf_array;

	// Main calculation unit for WL, calculate one cell of 4D pf array.
	void calculate_WL(int i, int j, int k, int l, bool bt = false); // Arrays to use.
	void calculate_ext_dependencies(int i, int j, int k, int l); // Arrays to use.

	// Access WL pf array using sequence indices.
	double x(int i, int j, int k, int l);
	double& x_ext(int i, int j, int k, int l);
	double& x_setter(int i, int j, int k, int l);

	// Following are for managing loop limits for WL coming from alignment module.
	bool check_boundary(int i1, int i2); // Checks if i2, seq2 sequence index, is in loop limit of i1, seqwence index of 1st sequence.
};

#endif // _WL_

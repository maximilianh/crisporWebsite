#ifndef _V_mhe_
#define _V_mhe_

// Forward declare the objects that V_mhe depends on calculation.
class t_seq_man; // Sequence manager to interface with sequence data.
class t_ppf_V;
class t_aln_priors;
class t_template_pf_array;
class t_spf_array;
class t_ppf_W;
class t_ppf_WMB;
class t_SS_str;
class t_ppf_tb_stack;
class t_MAP_structures;
class t_ppf_V_bpi;
class t_ppf_W_bpi;
class t_ppf_W_bp_aln_up;
class t_ppf_V_bp_aln_up;
class t_ppf_W_mhi;
class t_stoch_sampled_str_aln_sample_set;
class t_stoch_sampling_math;
class t_ppf_loops;

class t_ppf_V_mhe
{
public:
	t_ppf_V_mhe(t_ppf_loops* _ppf_loops); // Constructor.

	~t_ppf_V_mhe(); // Destructor.

	t_ppf_loops* ppf_loops;

	int N1;
	int N2;

	t_aln_priors* aln_priors;
	t_spf_array* seq1_spf;
	t_spf_array* seq2_spf;

	t_template_pf_array* pf_array;
	t_template_pf_array* ext_pf_array;

	// Main calculation unit for W_bi, calculate one cell of 4D pf array.
	void calculate_V_mhe(int i, int j, int k, int l, bool bt = false); // Arrays to use.
	void calculate_ext_dependencies(int i, int j, int k, int l); // Arrays to use.
	//void map_tb_V_mhe(int i, int j, int k, int l); // Arrays to use.
	//void stoch_tb_V_mhe(int i, int j, int k, int l);

	// Access W_bi pf array using sequence indices.
	double& x_setter(int i, int j, int k, int l);
	double x(int i, int j, int k, int l);
	double& x_ext(int i, int j, int k, int l);

	// Following are for managing loop limits for W coming from alignment module.
	bool check_boundary(int i1, int i2); // Checks if i2, seq2 sequence index, is in loop limit of i1, seqwence index of 1st sequence.

	double calc_ij_mhe_bp_score(int i, int j, int k, int l); // Arrays to use.
	double calc_kl_mhe_bp_score(int i, int j, int k, int l); // Arrays to use.
};

#endif // _V_mhe_

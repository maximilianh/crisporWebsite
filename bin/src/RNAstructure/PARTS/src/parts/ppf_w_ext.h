#ifndef _WEXT_
#define _WEXT_

// Forward declare the classes that W_ext depends on calculation.
class t_ppf_W;
class t_ppf_WMB;
class t_seq_man; // Sequence manager to interface with sequence data.
class t_spf_array; // Single partition function array.
class t_aln_priors;
class t_template_pf_array;
class t_frag_aln_enum_array;
class t_ppf_V_bpi;
class t_SS_str;
class t_ppf_W_bpi;
class t_ppf_tb_stack;
class t_MAP_structures;
class t_ppf_W_bp_aln_up;
class t_ppf_W_mhi;
class t_MAP_alignment;
class t_stoch_sampled_str_aln_sample_set;
class t_stoch_sampling_math;

class t_ppf_WEXT
{
public:
	t_ppf_WEXT(t_ppf_loops* _ppf_loops);
	~t_ppf_WEXT();

	t_ppf_loops* ppf_loops;

	// 2D pointer.
	double** pf_array;
	double** ext_pf_array;

	void alloc_init_w_ext_array();
	void alloc_init_ext_w_ext_array();

	double& x(int i, int k);
	double& x_ext(int i, int k);
	bool check_boundary(int i1, int i2);

	int N1;
	int N2;

	t_seq_man* seq_man;
	t_aln_priors* aln_priors;
	t_spf_array* seq1_spf;
	t_spf_array* seq2_spf;

	// This is a 2D array with no constraints on the maximum separation of structurally coincidable nucleotides.
	//t_template_pf_array* w_ext_pf_array;

	void calculate_W_ext();
	void calculate_ext_W_ext();
	void calculate_ext_dependencies(int i, int j, int k, int l); // Arrays to use.
	void backtrack(int i, int k);

	void rescale_internal(bool up_scale);
	void rescale_external(bool up_scale);
};

#endif // _WEXT_




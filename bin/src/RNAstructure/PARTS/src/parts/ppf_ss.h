#ifndef _PPF_SS_
#define _PPF_SS_

// Single stranded sequence containing structure handler.
class t_ppf_V;
class t_ppf_W;
class t_ppf_WL;
class t_ppf_WR;
class t_ppf_WMB;
class t_ppf_WMBL;
class t_ppf_WMBR;
class t_frag_aln_enum_array;
class t_seq_man;
class t_aln_priors;
class t_spf_array;
class t_MAP_alignment;
class t_stoch_sampled_str_aln_sample_set;
class t_stoch_sampling_math;
class t_ppf_loops;

// single stranded sequence containing structures, those are all exterior structures where j > N1 and l > N2.
class t_ppf_SS 
{
public:
	t_ppf_SS(t_ppf_loops* _ppf_loops); // Constructor.
	~t_ppf_SS(); // Destructor.

	double n_alloced_bytes;

	t_ppf_loops* ppf_loops;

	double**** pf_array;
	double**** ext_pf_array;

	int N1;
	int N2;

	t_seq_man* seq_man;
	t_spf_array* seq1_spf;
	t_spf_array* seq2_spf;
	t_aln_priors* aln_priors;

	short** seq1_ptr_reloc_map;
	short** seq2_ptr_reloc_map;

	// Compute all the alignments.
	void compute();
	void backtrack(int i, int j, int k, int l);

	void alloc_init_pf_array();
	void alloc_init_ext_pf_array();

	void rescale(bool up_scale);

	bool check_boundary(int i1, int i2);
	bool check_str_coinc_ll(int i, int j, int k, int l);

	double x(int i, int j, int k, int l);
	double& x_ext(int i, int j, int k, int l);
	double& x_setter(int i, int j, int k, int l);
};

#endif // _PPF_SS_

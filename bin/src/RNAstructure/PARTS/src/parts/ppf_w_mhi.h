#ifndef _W_mhi_
#define _W_mhi_

// Forward declare the objects that v_bpi depends on calculation.
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
class t_ppf_W_bpi;
class t_ppf_W_bp_aln_up;
class t_ppf_V_bp_aln_up;
class t_ppf_V_mhe;

class t_stoch_sampled_str_aln_sample_set;
class t_stoch_sampling_math;

class t_ppf_loops;

class t_ppf_W_mhi
{
public:
	t_ppf_W_mhi(t_ppf_loops* _ppf_loops); // Constructor.

	~t_ppf_W_mhi(); // Destructor.

	int N1;
	int N2;

	t_ppf_loops* ppf_loops;

	t_spf_array* seq1_spf;
	t_spf_array* seq2_spf;
	t_aln_priors* aln_priors;

	t_template_pf_array* ij_pf_array;
	t_template_pf_array* kl_pf_array;

	t_template_pf_array* ext_ij_pf_array;
	t_template_pf_array* ext_kl_pf_array;

	bool check_boundary(int i1, int i2);

	// Main calculation unit for W_bi, calculate one cell of 4D pf array.
	void calculate_ij_W_mhi(int i, int j, int k, int l, bool bt = false); // Arrays to use. 
	void calculate_ij_ext_dependencies(int i, int j, int k, int l); // Arrays to use.
	void calculate_kl_W_mhi(int i, int j, int k, int l, bool bt = false); // Arrays to use. 
	void calculate_kl_ext_dependencies(int i, int j, int k, int l); // Arrays to use.

	void map_tb_W_mhi(int i, int j, int k, int l); // Arrays to use.
	void map_tb_W_ij_mhi(int i, int j, int k, int l); // Arrays to use.
	void map_tb_W_kl_mhi(int i, int j, int k, int l); // Arrays to use.
	void stoch_tb_W_mhi(int i, int j, int k, int l);
	void stoch_tb_W_ij_mhi(int i, int j, int k, int l);
	void stoch_tb_W_kl_mhi(int i, int j, int k, int l);

	// Note that following are still used in map calculations.
	double& x_ij_setter(int i, int j, int k, int l);
	double x_ij(int i, int j, int k, int l);

	double& x_kl_setter(int i, int j, int k, int l);
	double x_kl(int i, int j, int k, int l);

	double& x_ext_ij(int i, int j, int k, int l);
	double& x_ext_kl(int i, int j, int k, int l);
};

#endif // _W_mhi_

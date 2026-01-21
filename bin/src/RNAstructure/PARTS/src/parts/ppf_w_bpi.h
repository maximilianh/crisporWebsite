#ifndef _w_bpi_
#define _w_bpi_

/*
	w_bpi: Base pair insertions on open structures: Insert base pairs on open hairpin ss's, open mbl's and W's.

	Note that those structures include a closed structure of one sequence and open structure of other sequence.
	This alignment does not go into any of the array definitions so that w_bpi array is needed. 

	-> w_bpi array turns into V_bpi with a base pair insertion which closes the open structure of sequence.
		This insertion can be handled in V_bpi. V_bpi must also recuse on w_bpi.

	-> w_bpi array turns into V by an aligned base pair addition in V calculation. V must also recurse on w_bpi.
		
	-> Other than those base pair insertions can be carried on in V_bpi. However W must recurse directly on V_bpi.
		That ensures insertion of base pairs at opening of internal loops for interior calculations and closing
		of internal loops for exterior calculations.
	
	-> Exclusion of isolated base pair insertion runs can be acCOMPARElished by making sure that base pairs are inserted 
		next to an aligned base pair. Those are introduced by recursing W on w_bpi's, which should not happen. 

	-> w_bpi recurses on ss_str, W, WMB.
		* ss_str is is hairpins for internal calculations and single stranded structures for external calculations.
		os_ij_bpi(i,j,k,l) refers to bpi at i-j and and open structure at k-l, which corresponds to 
		(W+WMB+ss_str)(i+1, j-1, k, l). 

*/

// Forward declare the objects that w_bpi depends on calculation.
class t_ppf_V;
class t_ppf_W;
class t_ppf_WMB;
class t_seq_man; // Sequence manager to interface with sequence data.
class t_spf_array; // Single partition function array.
class t_aln_priors;
class t_template_pf_array;
class t_frag_aln_enum_array;
class t_ppf_V_bpi;
class t_ppf_V_bp_aln_up;
class t_SS_str;
class t_ppf_tb_stack;
class t_MAP_structures;
class t_ppf_W_bp_aln_up;
class t_ppf_V_mhe;
class t_ppf_W_mhi;
class t_MAP_alignment;
class t_stoch_sampled_str_aln_sample_set;
class t_stoch_sampling_math;
class t_ppf_loops;

class t_ppf_W_bpi
{
public:
	t_ppf_W_bpi(t_ppf_loops* _ppf_loops);
	~t_ppf_W_bpi(); // Destructor.

	int N1;
	int N2;

	t_ppf_loops* ppf_loops;

	t_aln_priors* aln_priors;
	t_spf_array* seq1_spf;
	t_spf_array* seq2_spf;

	t_template_pf_array* os_ij_bpi_pf_array;
	t_template_pf_array* os_kl_bpi_pf_array;

	// Main calculation unit for W_bi, calculate one cell of 4D pf array.
	double calculate_os_ij_bpi(int i, int j, int k, int l); // Arrays to use.
	double calculate_os_kl_bpi(int i, int j, int k, int l); // Arrays to use.
	void map_tb_os_ij_bpi(int i, int j, int k, int l); // Arrays to use.
	void map_tb_os_kl_bpi(int i, int j, int k, int l); // Arrays to use.
	void stoch_tb_os_ij_bpi(int i, int j, int k, int l); // Arrays to use.
	void stoch_tb_os_kl_bpi(int i, int j, int k, int l); // Arrays to use.

	// Access W_bi pf array using sequence indices.
	double x_bpi(int i, int j, int k, int l);
	double x_map_bpi(int i, int j, int k, int l);

	double& x_ij_bpi(int i, int j, int k, int l);
	double& x_kl_bpi(int i, int j, int k, int l);

	// Following are for managing loop limits for W coming from alignment module.
	bool check_boundary(int i1, int i2); // Checks if i2, seq2 sequence index, is in loop limit of i1, seqwence index of 1st sequence.

	// Calculate scores of base pairing with base pair insertions on both sides.
	double calc_ij_bpi_bp_score(int i, int j, int k, int l); // Arrays to use.
	double calc_kl_bpi_bp_score(int i, int j, int k, int l); // Arrays to use.
};

#endif // _w_bpi_

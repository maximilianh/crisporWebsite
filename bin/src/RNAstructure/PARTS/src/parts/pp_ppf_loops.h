#ifndef _PPF_LOOPS_
#define _PPF_LOOPS_

class t_seq_man;
class t_spf_array;
class t_aln_priors;

class t_ppf_V;
class t_ppf_WL;
class t_ppf_WR;
class t_ppf_W;
class t_ppf_WMBL;
class t_ppf_WMBR;
class t_ppf_WMB;
class t_ppf_W_mhi;
class t_ppf_V_mhe;
class t_ppf_V_bp_aln_up;
class t_ppf_W_bp_aln_up;
class t_ppf_W_bpi;
class t_ppf_V_bpi;
class t_SS_str;
class t_ppf_WEXT;

class t_ppf_cli;

class t_pp_results;

// Calculate pairwise partition function based on base pairing 
// probabilities from individual single partition functions.
// Call this function as main interface to pairwise partition function.
t_pp_results* calculate_bp_ppf(t_ppf_cli* ppf_cli, // Need access to command line interface input.
					  t_seq_man* seq_man, 
					  t_spf_array* seq1_spf, 
					  t_spf_array* seq2_spf, 
					  t_aln_priors* aln_priors);

void stochastic_traceback(t_ppf_cli* ppf_cli, t_seq_man* seq_man, 
						  t_ppf_V* V, 
						  t_ppf_V_bp_aln_up* v_bp_aln_up, 
						  t_ppf_W_bp_aln_up* w_bp_aln_up, 
						  t_ppf_V_mhe* V_mhe, t_ppf_W_mhi* w_mhi, 
						  t_ppf_W* W, t_ppf_WL* WL, t_ppf_WR* WR, 
						  t_ppf_WMB* WMB, t_ppf_WMBL* WMBL, t_ppf_WMBR* WMBR, 
						  t_ppf_V_bpi* V_bpi, t_ppf_W_bpi* w_bpi, 
						  t_SS_str* ss_scorer,
						  t_ppf_WEXT* W_ext);

void rescale_all_ppf_arrays_by_indices(int i, int j, int k, int l,
								t_ppf_V* V,
								t_ppf_WL* WL,
								t_ppf_WR* WR,
								t_ppf_W* W,
								t_ppf_WMBL* WMBL,
								t_ppf_WMBR* WMBR,
								t_ppf_WMB* WMB,
								t_ppf_W_mhi* w_mhi,
								t_ppf_V_mhe* V_mhe,
								bool up_scale);

void scale_pp_value(double& value, double factor_de_scale, bool up_scale);


#endif

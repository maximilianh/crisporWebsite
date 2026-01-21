#ifndef _POSTERIOR_PAIRING_PROB_LOOPS_
#define _POSTERIOR_PAIRING_PROB_LOOPS_

class t_ppf_V;
class t_ppf_WL;
class t_ppf_WR;
class t_ppf_W;
class t_ppf_WMBL; 
class t_ppf_WMBR;
class t_ppf_WMB;
class t_seq_man;
class t_SS_str;
class t_ppf_V_bpi;
class t_ppf_W_bpi;
class t_spf_array;
class t_aln_priors;
class t_ppf_W_bp_aln_up;
class t_ppf_V_bp_aln_up;
class t_ppf_W_mhi;
class t_ppf_V_mhe;

class t_ppf_cli;

class t_pp_results;

t_pp_results* calculate_bp_scores_mixed_bp_inserts(t_ppf_cli* ppf_cli, t_ppf_V* V, t_ppf_W* W, t_ppf_V_mhe* v_mhe, t_ppf_W_mhi* w_mhi, t_ppf_WMB* WMB, t_SS_str* ss_scorer);

#endif // _POSTERIOR_PAIRING_PROB_LOOPS_


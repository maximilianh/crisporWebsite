#ifndef _PFE_SAMPLED_STRUCTURES_
#define _PFE_SAMPLED_STRUCTURES_

class t_seq_man;
class t_ppf_cli;

class t_stoch_sampled_structures
{
public:
	t_stoch_sampled_structures(t_seq_man* seq_man, t_ppf_cli* ppf_cli);
	~t_stoch_sampled_structures();

	int* seq1_sampled_ct_bps;
	int* seq2_sampled_ct_bps;
	int N1;
	int N2;
	t_seq_man* seq_man;
	t_ppf_cli* ppf_cli;

	void add_bp_sampled_ct1(int i, int j);
	void add_bp_sampled_ct2(int k, int l);

	void dump_sampled_cts(int _sample_id);

	void reset_bps();
};

#endif // _PFE_SAMPLED_STRUCTURES_



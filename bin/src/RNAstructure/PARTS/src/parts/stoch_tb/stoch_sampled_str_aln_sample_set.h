#ifndef _STRUCTURAL_ALIGNMENT_SAMPLE_SET_
#define _STRUCTURAL_ALIGNMENT_SAMPLE_SET_

// This class encapsulates the structural alignment sample set.
// Handles both sampled alignments and also structures.
class t_stoch_sampled_structures;
class t_stoch_sampled_alignment;
class t_ppf_cli;
class t_seq_man;
class t_ppf_loops;

class t_spf_array;
class t_aln_priors;

class t_stoch_sampled_str_aln_sample_set
{
public:
	//t_stoch_sampled_str_aln_sample_set(t_ppf_cli* _ppf_cli, t_seq_man* _seq_man, int _n_samples);
	t_stoch_sampled_str_aln_sample_set(t_ppf_loops* _ppf_loops);

	~t_stoch_sampled_str_aln_sample_set();

	t_ppf_loops* ppf_loops;

	char sampled_ops_dir[1000];
	t_ppf_cli* ppf_cli;
	t_seq_man* seq_man;
	t_spf_array* seq1_spf;
	t_spf_array* seq2_spf;
	t_aln_priors* aln_priors;
	int n_samples;

	// Indices of of min and max prob. str alignments.
	int min_pfe_str_aln_index;
	int max_pfe_str_aln_index;

	int i_current_sample; // Index of current sample being processed.

	t_stoch_sampled_structures* current_sampled_structures();
	t_stoch_sampled_alignment* current_sampled_alignment();

	t_stoch_sampled_structures* sampled_structures;
	t_stoch_sampled_alignment* sampled_alignment;

	// Inferences of base pairing and alignment probabilities.
	void dump_strs_aln(int _sample_id);
	void reset_strs_aln();

	void marginalize_alignment_probabilities();
	void marginalize_pairing_probabilities();

	double** seq1_pairing_probs;
	double** seq2_pairing_probs;
	double** state_aln_probs;
	double** state_ins1_probs;
	double** state_ins2_probs;
};

#endif // _STRUCTURAL_ALIGNMENT_SAMPLE_SET_


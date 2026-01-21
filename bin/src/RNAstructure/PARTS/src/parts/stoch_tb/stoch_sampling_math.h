#ifndef _SAMPLING_MATH_
#define _SAMPLING_MATH_

#include "../parts_compilation_directives.h"

class CRandomMother;
class t_ppf_cli;
class t_seq_man;
class t_spf_array;
class t_aln_priors;
class t_rng;
class t_ppf_loops;

class t_stoch_sampling_math
{
public:
	t_stoch_sampling_math(t_ppf_loops* _ppf_loops);

	t_ppf_loops* ppf_loops;

	/*
	t_stoch_sampling_math(t_ppf_cli* _ppf_cli, 
				t_seq_man* _seq_man, 
				t_spf_array* _seq1_spf,
				t_spf_array* _seq2_spf,
				t_aln_priors* _aln_priors,
				int _prn_gen_seed);
				*/

	~t_stoch_sampling_math();

	t_ppf_cli* ppf_cli;
	t_seq_man* seq_man;
	t_spf_array* seq1_spf;
	t_spf_array* seq2_spf;
	t_aln_priors* aln_priors;

	t_rng* prn_gen; // ran3 random number generator based on numerical recipes book.

	int prn_gen_seed; // Pseudo random number generator seed.
	double* global_score_array; // Global scoring array which is used to store all the scores. Allocate this once and do not reallocate.

	int sample_size;

	double random_double_1_0();
	int sample_score_array(double* scores, int n_scores); // Sample score array.
};

#endif // _SAMPLING_MATH_




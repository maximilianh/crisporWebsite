#ifndef _PHMM_
#define _PHMM_

#include <algorithm>
#define MIN(x,y) std::min(x,y)
#define MAX(x,y) std::max(x,y)

/*
t_phmm class encapsulates the hidden markov model parameters which define the current model used.
It is initialized from a hmm parameter file or a set of parameters.
All the hidden markov model specific stuff goes here.
This class does not specify a singleton because the probability computation uses more than 1 set of parameters.
*/

#define N_STATES (3)
#define N_OUTPUTS (27)
#define N_BINZ (10)

enum {STATE_INS1, STATE_INS2, STATE_ALN};
#define STATE_START STATE_ALN
#define STATE_END STATE_ALN

#define GAP_SYM '.'

// This file defines the hmm model for alignment hmm model.
// This file is only for alignment hmm model, the hmm model 
// is a pairwise alignment model where 3 states can output 
// 2 sequence symbols (nucleotides) or gaps. 
// Following are the parameters trained with 10 full families.
// Following parameters are used for ML calculation. They should be loaded.
static double ML_emit_probs[N_OUTPUTS][N_STATES] = 
{
	{0.000000, 0.000000, 0.134009}, // AA
	{0.000000, 0.000000, 0.027164}, // AC
	{0.000000, 0.000000, 0.049659}, // AG
	{0.000000, 0.000000, 0.028825}, // AU
	{0.211509, 0.000000, 0.000000}, // A.
	{0.000000, 0.000000, 0.027164}, // CA
	{0.000000, 0.000000, 0.140242}, // CC
	{0.000000, 0.000000, 0.037862}, // CG
	{0.000000, 0.000000, 0.047735}, // CU
	{0.257349, 0.000000, 0.000000}, // C.
	{0.000000, 0.000000, 0.049659}, // GA
	{0.000000, 0.000000, 0.037862}, // GC
	{0.000000, 0.000000, 0.178863}, // GG
	{0.000000, 0.000000, 0.032351}, // GU
	{0.271398, 0.000000, 0.000000}, // G.
	{0.000000, 0.000000, 0.028825}, // UA
	{0.000000, 0.000000, 0.047735}, // UC
	{0.000000, 0.000000, 0.032351}, // UG
	{0.000000, 0.000000, 0.099694}, // UU
	{0.259744, 0.000000, 0.000000}, // U.
	{0.000000, 0.211509, 0.000000}, // .A
	{0.000000, 0.257349, 0.000000}, // .C
	{0.000000, 0.271398, 0.000000}, // .G
	{0.000000, 0.259744, 0.000000}, // .U
	{0.000000, 0.000000, 0.000000}, // ..
	{0.000000, 0.000000, 1.000000}, // START
	{0.000000, 0.000000, 1.000000}  // END
};

static double ML_trans_probs[N_STATES][N_STATES] = 
{
	{0.666439, 0.041319, 0.292242}, // INS1
	{0.041319, 0.666439, 0.292242}, // INS2
	{0.022666, 0.022666, 0.954668}  // ALIGN
};

static char pw_state_names[][40] = {"PW_STATE_INS1", "PW_STATE_INS2", "PW_STATE_ALN"};

class t_phmm
{
public:
	// Replace the emission and transition probabilities.
	t_phmm(double new_emit_probs[N_OUTPUTS][N_STATES], double new_trans_probs[N_STATES][N_STATES]);
	t_phmm(const char* phmm_pars_file);

	// Set the pointers to parameters using the pointers in the arguments. 
	// This constructor is aimed at saving memory.
	t_phmm(double** _emit_probs, double** _trans_probs, double* _fam_thresholds);
	~t_phmm();

	void set_parameters_by_sim(double similarity);
	int get_bin_index(double similarity, int n_bins);
	double get_fam_threshold(double similarity);

	double get_emit_prob(int sym_index, int state);
	double get_trans_prob(int prev, int next);

	void alloc_init_params();
	void free_params();

	//double emission_probs[N_OUTPUTS][N_STATES];
	//double trans_probs[N_STATES][N_STATES];
	double** emission_probs;
	double** trans_probs;

	//double fam_hmm_pars[N_BINZ * (N_STATES + N_OUTPUTS) * N_STATES];
	//double fam_thresholds[N_BINZ];
	double* fam_hmm_pars;
	double* fam_thresholds;

	static char state_names[N_STATES][100];

	void dump_parameters();
};

#endif // _PHMM_


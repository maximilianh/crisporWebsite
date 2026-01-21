#ifndef _HMM_ALN_
#define _HMM_ALN_
#include <string>
#include <iostream>
#include <vector>
using namespace::std;

#define DEFAULT_ALN_ENV_BAND_CONSTRAINT (15)
#define DEFAULT_PHMM_BAND_CONSTRAINT (0x1ffff)


enum{BANDED_ALN_ENV,
FULL_ALN_ENV,
PROB_ALN_ENV,
MANUAL_ALN_ENV};

class t_structure;
class t_phmm;
class t_p_alignment;
class t_phmm_array;
class t_matrix;

struct t_pp_result
{
	double** aln_probs;
	double** ins1_probs;
	double** ins2_probs;
	double similarity;
	double op_prob;
	double ml_similarity;
	double fam_threshold;
	double phmm_band_constraint_size;
};

struct t_ML_result
{
	double ml_prob;
	double ml_similarity;
	t_p_alignment* ml_alignment;
	vector<char>* seq1_aln_line;
	vector<char>* seq2_aln_line;
};

// Results of an alignment envelope computation.
struct t_aln_env_result
{
	//t_pp_result* pp_result; // Results of the pp computation.
	int* low_limits;
	int* high_limits;
};

// t_phmm_aln object encapsulates a hidden markov model alignment.
// This can be a posterior probability computation or an ML alignment computation.
class t_phmm_aln
{
public:
	//t_phmm_aln(char* seq1_fp, char* seq2_fp, int _band_size);
	//t_phmm_aln(t_structure* seq1_fp, t_structure* seq2_fp, int _band_size);
	t_phmm_aln(char* seq1_fp, char* seq2_fp);
	t_phmm_aln(t_structure* seq1_fp, t_structure* seq2_fp);
	//t_phmm_aln(t_structure* seq1_fp, t_structure* seq2_fp, double** _aln_coinc_priors);
	//t_phmm_aln(t_structure* seq1_fp, t_structure* seq2_fp, t_matrix* _aln_coinc_priors);
	~t_phmm_aln();

	// Setters, for decreasing number of constructors.
	void set_aln_coinc_priors(t_matrix* _aln_priors);
	void set_aln_coinc_priors(double** _aln_priors);
	void set_match_priors(t_matrix* _match_priors);
	void set_match_priors(double** _match_priors);
	void set_phmm_band_size(int _phmm_band_size);
	void set_aln_constraints(int* _aligned_positions);
    void set_match_score_priors(t_matrix* _match_score_priors);
	void set_match_score_priors(double** _match_score_priors);

	void check_set_seqs();

	// The priors.
	t_matrix* aln_coinc_priors;
	double get_coinc_prior(int i, int k);
	t_matrix* match_priors;
        double get_match_prior(int i, int k, int n1, int n2);
        
        void set_indel_priors(vector<double>* _indel_priors);
        vector<double>* indel_priors;
        double get_indel_prior(int i, int k, int n1,int n2);

    t_matrix* match_score_priors;
        double get_match_score_prior(int i, int k, int n1, int n2);

	// The sequences.
	t_structure* seq1;
	t_structure* seq2;

	int l1();
	int l2();

	// The model that this alignment uses.
	t_phmm* phmm; 

	// Applies to all computations done via this phmm_aln object.
	int phmm_band_constraint_size;

	// This is the alignment constraint: Constraint on aligned positions.
	int* seq1_alignment_constraints;
	int* seq2_alignment_constraints;

	// For replacing unknown nucleotides.
	char generate_random_nuc();
	double get_aln_similarity(char* aln_line1, char* aln_line2);

	void init_forward_array(t_phmm_array* fore_array);
	void init_backward_array(t_phmm_array* back_array);
	void init_ML_array(t_phmm_array* ml_array);

	// ML and pp functions.
	void compute_forward_array(t_phmm_array* fore_array);
	void compute_backward_array(t_phmm_array* back_array);
	t_phmm_array* compute_ML_array(t_ML_result* ml_res);
	void traceback_ml_array(t_phmm_array* ml_array, t_ML_result* ml_res);

	// aln_env_utils.
	bool check_backward_connection(int i, int k, bool** pruned_aln_env);
	bool check_forward_connection(int i, int k, bool** pruned_aln_env);
	bool** prune_aln_env(bool** aln_env);

	bool check_connection(bool** aln_env);
	bool check_connection(t_aln_env_result* aln_env_res);

	void load_map_limits_from_map(const char* aln_map_fn, int* low_limits, int* high_limits);

	// Check the loop limits for making sure that 
	void check_ins1_ins2(t_aln_env_result* aln_env_result);

	double get_trans_emit_prob(int prev_state, int current_state, int i, int k);
	int nuc2num(char nuc);

	// These are the main services that this object offers.
	//t_pp_result* compute_posterior_probs(int* aligned_positions, double** coinc_priors);
	t_pp_result* compute_posterior_probs();
	void* free_pp_result(t_pp_result* pp_result);

	//t_ML_result* compute_ML_alignment(int* aligned_positions, double** coinc_priors);
	t_ML_result* compute_ML_alignment();
	t_aln_env_result* compute_alignment_envelope(int aln_env_type, int par);
	t_aln_env_result* compute_alignment_envelope(int aln_env_type, 
															 t_pp_result* pp_result, 
															 double log_threshold, 
															 int par);

	void* free_aln_env_result(t_aln_env_result* aln_env_result);
	void* free_ML_result(t_ML_result* ml_result);

	void dump_data_files(t_pp_result* pp_results);
	void dump_data_files(t_ML_result* ml_results);

	int* get_seq2_aln_const(int* seq1_aln_const);
	void get_aln_permissions(bool& forbid_STATE_ALN, 
								bool& forbid_STATE_INS1, 
								bool& forbid_STATE_INS2, 
								int i, 
								int k);
};

t_phmm_aln* create_phmm_aln(const vector<char>& seq1_nucs, const vector<char>& seq2_nucs);
void delete_phmm_aln(t_phmm_aln* a);
void write_probability_array(t_pp_result* pp_result, const char* outputfile,int l1,int l2,bool logprobs);
void write_ML_alignment(t_ML_result* ML_result, const char* outputfile, int l1, int l2, const char* seq1, const char* seq2);


#endif // _HMM_ALN_


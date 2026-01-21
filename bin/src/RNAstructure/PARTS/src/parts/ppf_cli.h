#ifndef _PPF_CLI_
#define _PPF_CLI_

/*
t_ppf_cli: 
ppf command line interface class to parse command line arguments that are input to ppf program. This class basically controls
how program will run based on command line arguments.
Format is:

pairwise_part_func [sequence 1 path] [sequence 2 path] -map/-pp -ap [aln prior file path] -fp [folding prior] -ct1 [correct ct 1 path] -ct2 [correct ct 2 path]

All parts after -map/-pp are optional. Correct ct files are needed to dump to plots and calculate structural 
prediction accuracies.

-map/-pp flag is used to indicate which calculation to make, that is, maximum likelihood alignment or posterior probability 
calculation to get pairing probabilities. One of those is necessary.
*/

#define PARTS_RUN_MODE_INIT 0
#define PARTS_RUN_MODE_MAP 1
#define PARTS_RUN_MODE_PP 2
#define PARTS_RUN_MODE_STOCH_SAMPLE 16

// Error codes for t_ppf_loops class:
enum{NO_CLI_ERROR,
ERR_CONF_FILE_NOT_EXISTENT,
ERR_SEQ1_NOT_EXISTENT,
ERR_SEQ2_NOT_EXISTENT,
ERR_MODE_NOT_VALID,
N_CLI_ERRS};

static char ppf_cli_error_msgs[N_CLI_ERRS][100] = {"No error.",
"Could not open configuration file.",
"Could not open sequence 1.",
"Could not open sequence 2.",
"Invalid mode parameter."};

class t_ppf_cli
{
	public:

	t_ppf_cli(int argc, char* argv[]); // ctor: parse command line arguments.
	t_ppf_cli(t_ppf_cli* _ppf_cli_to_copy);
	t_ppf_cli(char* seq1_fp, 
			 char* seq2_fp, 
			 int mode);

	t_ppf_cli(char* conf_fp);

	~t_ppf_cli();
	
	//!  Set new output file prefixes for sequence 1 and sequence 2 even after a call to the constructor.
        //!  Typically, this function will be used to set new output file names after calling the:
        //!  t_ppf_cli(char* seq1_fp, char* seq2_fp, int mode)
        //!  version of the constructor.  It takes two null terminated character strings, deletes the contents of
        //!  seq1_op_file_prefix and seq2_op_file_prefix, reallocates each, and then fills them with the contents of the
        //!  two supplied strings.
	void set_output_prefixes(const char* seq1_prefix, const char* seq2_prefix);

	// Input sequence paths.
	char* seq1_path;
	char* seq2_path;

	bool use_array_files;
	bool save_array_files;

	void init_paths_vars();

	int phmm_band_constraint_size;

	char* seq1_op_file_prefix;
	char* seq2_op_file_prefix;

	char* seq1_map_ct_op;
	char* seq2_map_ct_op;
	char* map_aln_op;
	char* seq1_pp_op;
	char* seq2_pp_op;
	char* seq1_sample_ct_op;
	char* seq2_sample_ct_op;
	char* sample_aln_op;

	char* seq1_SHAPE_path;
	char* seq2_SHAPE_path;

	// alignment and sequence ids for determining 
	// an identifier for the alignment.
	char* seq1_id;
	char* seq2_id;
	char* aln_id;

	 // Is this a map calculation?
	//bool map;
	static int mode; 

	// This is the value of alignment weight if there is any, this is going to
	// be the log weight of alignment in pseudo free energy model.
	double log_aln_weight_per_aln_pos;

	int stoch_sampling_seed;	

	// This is the overriding rescaling factor per nucleotide.
	// Input by -rpfn flag.
	double cmd_rescaling_increment_factor_per_nucleotide;

	// aln prior file path, if existing.
	char* aln_prior_path;
	char* ins1_prior_path;
	char* ins2_prior_path;

	double array_mem_limit_in_megs;

	// This is the propagation constant of priors,
	// it affects how the priors before affect the posterior
	// probabilities in the pfe sampling. Assuming that these probabilities
	// are going to be used as priors for next algorithm, this is 
	// necessary to keep continuity and smooth out transitions
	// This parameter determines what fraction of the sample size
	// is going to be used as pseuodocount in posteriors.
	double prior_to_posterior_prop_rate;

	// fold prior 1 file path, if existing.
	char* fold_prior1_path;
	char* str1_coinc_info_path;

	// fold prior 2 file path, if existing.
	char* fold_prior2_path;

	// unpairing paths.
	char* i_unpairing_prior1_path;
	char* i_unpairing_prior2_path;
	char* j_unpairing_prior1_path;
	char* j_unpairing_prior2_path;

	// correct ct1 path.
	char* correct_ct1_path;

	// correct ct2 path.
	char* correct_ct2_path;

	bool check_file(char* file_path);
	void print_usage_msg();
	bool have_correct_ct1();
	bool have_correct_ct2();
	bool have_fold_prior1();
	bool have_fold_prior2();
	bool have_aln_prior();

	void get_seq_aln_ids();

	// Sampling mode; number of structural alignment to sample.
	int n_samples;

	// The path of array dump file.
	char* array_dump_fp;

	int max_n_separation_between_nucs;

	double fold_env_prob_treshold;
	double str_coinc_env_prob_treshold;

	char* loop_limits_fn;

	void load_conf_file(char* conf_file);

	// Constructor error message. 0 signals no errors.
	int GetErrorCode();

	// Return the error message for error_code.
	char* GetErrorMessage(const int error_code);

private:
	int last_error_code;
};

#endif // _PPF_CLI_

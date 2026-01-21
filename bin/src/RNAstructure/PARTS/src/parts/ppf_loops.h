#ifndef _PPF_LOOPS_
#define _PPF_LOOPS_

#include <stdio.h>


class t_ppf_cli;
class t_seq_man;
class t_spf_array;
class t_aln_priors;
class t_ppf_W;
class t_ppf_WL;
class t_ppf_WMB;
class t_ppf_WMBL;
class t_ppf_V_mhe;
class t_ppf_W_mhi;
class t_ppf_tb_stack;
class t_ppf_SS;
class t_ppf_WEXT;
class t_seq_man; // Sequence manager to interface with sequence data.
class t_map_results;
class t_pp_results;
class t_ppf_cli;
class t_ppf_tb_stack;
class t_stoch_sampling_math;
class t_stoch_sampled_str_aln_sample_set;
class t_MAP_structures;
class t_MAP_alignment;
class t_ppf_scale;

class t_map_results;

// Error codes for t_ppf_loops class:
enum{NO_PPF_ERROR,
ERR_CLI_ERROR,
ERR_MAP_COMPUTATIONS_NOT_PERFORMED,
ERR_PPF_COMPUTATIONS_NOT_PERFORMED,
ERR_INVALID_ARGUMENT,
N_ERRS};



static char ppf_loops_error_msgs[N_ERRS][100] = {"No error.",
"Cli error.",
"MAP computations not performed, yet.",
"PPF computations not performed, yet.",
"Invalid argument is supplied to the function."};

class t_ppf_loops
{
public:
	// Use configuration file.
	t_ppf_loops(char* conf_file);

	// Use sequence paths and mode to do computations, do not dump.
	t_ppf_loops(char* seq1_fp,
				char* seq2_fp,
				int mode);

	t_ppf_loops(char* seq1_fp,
				char* seq2_fp);

	t_ppf_loops(int argc, char* argv[]);

	void init_loops();
	
	//!  Set new output file prefixes for sequence 1 and sequence 2 even after a call to the constructor.
        //!  Typically, this function will be used to set new output file names after calling the:
        //!  t_ppf_loops(char* seq1_fp, char* seq2_fp)
        //!  version of the constructor.  It takes two null terminated character strings, and passes them to t_ppf_cli::set_output_prefixes.  By extention it then deletes the contents of
        //!  t_ppf_cli::seq1_op_file_prefix and t_ppf_cli::seq2_op_file_prefix, reallocates each, and then fills them with the contents of the
        //!  two supplied strings.
        
	void set_output_prefixes(const char* seq1_prefix, const char* seq2_prefix);

	~t_ppf_loops();

	t_ppf_cli* ppf_cli;
	t_seq_man* seq_man; 
	t_spf_array* seq1_spf;
	t_spf_array* seq2_spf; 
	t_aln_priors* aln_priors;

	// Arrays:
	t_ppf_V_mhe* V_mhe;
	t_ppf_W* W;
	t_ppf_WL* WL;
	t_ppf_WMB* WMB;
	t_ppf_WMBL* WMBL;
	t_ppf_W_mhi* W_mhi;
	t_ppf_SS* SS;
	t_ppf_WEXT* W_ext;

	// Loop computation.
	//int compute_pseudo_free_energy_loops();
	int compute_pseudo_free_energy_loops();
	void compute_internal_pseudo_free_energy_loops();
	void compute_external_pseudo_free_energy_loops();

	int mallocate_arrays();
	void delete_arrays();

	void compute_posterior_pp();
	double** seq1_pp;
	double** seq2_pp;

	// Loop backtracking.
	void* backtrack_pseudo_free_energy_loops();

	// Stack that manages MAP/Stochastic backtracking.
	t_ppf_tb_stack* tb_stack;

	// Stoch sampled strs/alns.
	t_stoch_sampling_math* stoch_sampler;
	t_stoch_sampled_str_aln_sample_set* stoch_sampled_set;

	// MAP structures/alignment
	t_MAP_structures* map_strs;
	t_MAP_alignment* map_alignment;

	t_ppf_scale* ppf_scaler;

	void rescale_all_ppf_arrays_by_indices(int i_lim, int j_lim, int k_lim, int l_lim, bool up_scale);
	void rescale_all_external_ppf_arrays_by_indices(int i_lim, int j_lim, int k_lim, int l_lim, bool up_scale);
	void scale_pp_value(double& value, double factor_de_scale, bool up_scale);

	// Changes depending on the type of computation requested.
	double (*MAX_SUM)(double,double);
	bool (*BT_OP)(double,double);

	void set_aln(int i, int k);
	void set_seq1_ins(int i, int k);
	void set_seq2_ins(int i, int k);

	void add_ct1_bp(int i, int j);
	void add_ct2_bp(int k, int l);

	// Query individual probability values from a PP computation.	
	// Return <0 if there is an error.
	double GetPairProbabilitySeq1(const int i, const int j);
	double GetPairProbabilitySeq2(const int k, const int l);

	// Constructor error message. 0 signals no errors.
	int GetErrorCode();

	// Return the error message for error_code.
	char* GetErrorMessage(const int error_code);

	// Return the length of the alignment. Return -1 on error.
	signed int Get_MAP_Alignment_Length();
	//signed int Get_Seq1_Alignment(const int i);
	//signed int Get_Seq2_Alignment(const int k);

	// 
	signed int GetAlignmentSeq1(int alignmentcolumn);
	signed int GetAlignmentSeq2(int alignmentcolumn);

	// Return the 
	signed int GetAlignedNucleotideSeq1(int nucleotide1);
	signed int GetAlignedNucleotideSeq2(int nucleotide2);

	// Return sequence lengths.
	signed int GetLengthSequence1();
	signed int GetLengthSequence2();

	void dump_internal_arrays();
	void load_internal_arrays(char* array_dir);
	bool read_array_indices_value(FILE* f_array, int& i, int& j, int& k, int& l, double& val);

	// Perform MAP and PPF computations, return 0 on success.
	int compute_MAP();
	int compute_PPF();

private:
	int last_error_code;
};

#endif

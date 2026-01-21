#ifndef _SINGLE_PF_ARRAY_
#define _SINGLE_PF_ARRAY_

#include <stdio.h>
#include <stdlib.h>

// Single pf array class: Interface between single partition function calculation and double pf calculation.
// Read resulting probability arrays from single partition function into this class or take arrays directly from that calculation.
// The array is used for getting base pairing probabilities from single pf calculation. So the class is just encapsulation
// of one array of pairing probabilities.

class t_ppf_cli;
class t_folding_constraints;

class t_spf_array
{
public:

	int N; // Needs those in order to check if j < N1 is accessed.

	// The first nucleotide of each sequence is indexed by 1, not 0.
	double** pairing_array;
	double* ind_unpairing_array; // This is the probability that a nucleotide is not paired, ind_ prefix refers to the fact that array is 1D.
	double* ind_pairing_array; // This is the marginalized pairing probability, that is, probability that a nucleotide is paired.

	double** i_unpairing_array;
	double** j_unpairing_array;

	double n_bytes_alloced;

	double* ppf_scales; // Scales for linear calculations this is stored in order to make accession operation fast.

	// This always exists to at least get rid of allocations where
	// base pairing prior is 0, it can be set alternatively to 
	// filter base pairs with pairing prior lower than a threshold.
	bool** fold_env;
	bool** str_coinc_env;

	t_ppf_cli* ppf_cli;

	t_spf_array(int seq_length, bool mallocate); // Default constructor, inits all scores to CONVERT_FROM_LIN(1.0).
	t_spf_array(int seq_length, 
				 char* seq_path, 
				 t_ppf_cli* _ppf_cli, 
				 char* pairing_probs_file,
				 bool mallocate);

	~t_spf_array();

	// This is a 2D array where pointer_relocation_map[i_seq][j_seq] points to [i_seq][j_mem].
	// This helps to generate a nonlinear allocation. The non-linear allocation enables decrease in
	// allocated memory by avoiding allocation of non-structurally co-incident locations.
	//short** coinc_pointer_relocation_map; // Use to allocate unpaired (only co-incident locations)
	//short** paired_pointer_relocation_map; // Use to allocate paired arrays.
	t_folding_constraints* folding_constraints;

	// Calculate unpairing probabilities after spf is loaded.
	void calculate_unpairing_probs();

	// Dump single pairing probabilities (score) in a matlab loadable format.
	//void dump_spf_plane();
	//void dump_fold_env();

	// Do pfunction and pfunction_probs calls in class.
	void calculate_spf(char* seq_path);

	// Set ppf scaler.
	//void set_scaler(double* ppf_scaler);
	void set_scaler(double* _ppf_scales);

	// Check fold envelope.
	bool is_pairable(int i, int j);

	double px(int i1, int i2);
	double ux_5p(int i1, int i2);
	double ux_3p(int i1, int i2);
};

#endif

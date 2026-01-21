#include <string.h>
#include <limits.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "../../../RNA_class/RNA.h"
#include "parts_compilation_directives.h"
#include "single_pf_array.h"
#include "ppf_math.h"
#include "process_sequences.h"
#include "parts_paths.h"

#include "ppf_scale.h"
#include "ppf_cli.h"
#include "../../../src/phmm/structure/folding_constraints.h"
#include "../../../src/phmm/utils/file/utils.h"

#include <math.h>

#include <iostream>

using namespace std;

#include <algorithm>
#define MIN(x,y) min(x,y)
#define MAX(x,y) max(x,y)

bool _DUMP_SPF_MESSAGES_ = false;
bool _DUMP_SPF_PLANES_ = false;

// SPF array deals with log probabilities!!!
t_spf_array::t_spf_array(int seq_length, bool mallocate)
{
	// Copy sequence length.
	this->N = seq_length;
	this->n_bytes_alloced = 0.0f;

	if(mallocate)
	{
		this->pairing_array = (double**)malloc(sizeof(double*) * (N + 1));
		this->ind_unpairing_array = (double*)malloc(sizeof(double) * (N + 1));
		this->ind_pairing_array = (double*)malloc(sizeof(double) * (N + 1));
		this->fold_env = (bool**)malloc(sizeof(bool*) * (N + 1));
		this->str_coinc_env = (bool**)malloc(sizeof(bool*) * (N + 1));
	}
	else
	{
		this->pairing_array = NULL;
		this->ind_unpairing_array = NULL;
		this->ind_pairing_array = NULL;
		this->fold_env = NULL;
	}

	this->n_bytes_alloced += (sizeof(double*) * (N + 1) + 
								sizeof(double) * (N + 1) + 
								sizeof(double) * (N + 1) +
								sizeof(bool*) * (N + 1) +
								sizeof(bool*) * (N + 1));

	// SPF array is indexed starting from 1 to N1 and to N2.
	// utilize the constraints.
	for(int i1 = 1; i1 <= N; i1++)
	{
		if(mallocate)
		{
			this->fold_env[i1] = (bool*)malloc(sizeof(bool) * (N - i1 + 1));
			this->fold_env[i1] -= i1; // Do pointer shift for fold envelope.
		}
		this->n_bytes_alloced += (sizeof(bool) * (N - i1 + 1));

		if(mallocate)
		{
			this->pairing_array[i1] = (double*)malloc(sizeof(double) * (N - i1 + 1));
			this->pairing_array[i1] -= i1; // Do pointer shift.
		}
		this->n_bytes_alloced += (sizeof(double) * (N - i1 + 1));

		if(mallocate)
		{
			this->str_coinc_env[i1] = (bool*)malloc(sizeof(bool) * (N - i1 + 1));
			this->str_coinc_env[i1] -= i1; // Do pointer shift for fold envelope.
		}
		this->n_bytes_alloced += (sizeof(bool) * (N - i1 + 1));

		if(mallocate)
		{
			this->ind_unpairing_array[i1] = CONVERT_FROM_LIN(1.0);
			this->ind_pairing_array[i1] = CONVERT_FROM_LIN(1.0);

			for(int i2 = i1 + 1; i2 <= N; i2++)
			{
				this->pairing_array[i1][i2] = CONVERT_FROM_LIN(1.0); // Initialize the probabilities to 1. For uniform priors of folding.
				this->fold_env[i1][i2] = true; // Set all possible pairs.
				this->str_coinc_env[i1][i2] = true;
			}
		}
	} // i1 loop

if(_DUMP_SPF_MESSAGES_)
	printf("t_spf_array allocated %lf bytes\n", this->n_bytes_alloced);
}

// Destructor: free all arrays.
t_spf_array::~t_spf_array()
{
	if(this->ind_unpairing_array == NULL)
	{
		return;
	}

	// Allocate pairing and unpairing spf arrays.
	for(int i1 = 0; i1 <= N; i1++)
	{
		// Include the max_separation criterion in the allocation function.
		int min_i2 = i1;
		int max_i2 = MIN(i1 + ppf_cli->max_n_separation_between_nucs, N);

		this->fold_env[i1] += i1; // Do pointer shift for fold envelope.
		free(this->fold_env[i1]);
		
		this->pairing_array[i1] += i1; // Do pointer shift to access the array using sequence indices.
		free(this->pairing_array[i1]);
		
		this->str_coinc_env[i1] += i1; // Do pointer shift for fold envelope.
		free(this->str_coinc_env[i1]);
	} // i1 loop

	free(this->pairing_array);
	free(this->ind_unpairing_array);
	free(this->ind_pairing_array);
	free(this->fold_env);
	free(this->str_coinc_env);

	free(this->ppf_scales);

	delete(this->folding_constraints);
}

// Run spf calculation on sequence and load folding priors.
// THE PRIORS COMING FROM SPF PROGRAM ARE IN LINEAR DOMAIN!!!
t_spf_array::t_spf_array(int seq_length, 
						 char* seq_path, 
						 t_ppf_cli* _ppf_cli, 
						 char* pairing_probs_file, 
						 bool mallocate)
{
	this->ppf_cli = _ppf_cli;

	this->n_bytes_alloced = 0.0f;

	// Copy sequence length.
	this->N = seq_length;

	if(mallocate)
	{
		this->pairing_array = (double**)malloc(sizeof(double*) * (N + 1));
		this->ind_unpairing_array = (double*)malloc(sizeof(double) * (N + 1));
		this->ind_pairing_array = (double*)malloc(sizeof(double) * (N + 1));
		this->fold_env = (bool**)malloc(sizeof(bool*) * (N + 1));
		this->str_coinc_env = (bool**)malloc(sizeof(bool*) * (N + 1));

		this->n_bytes_alloced += ((sizeof(double*) * (N + 1)) + 
									(sizeof(double) * (N + 1)) + 
									(sizeof(double) * (N + 1)) +
									(sizeof(bool*) * (N + 1)) + 
									(sizeof(bool*) * (N + 1)) + 
									(sizeof(short*) * (N + 1)));
	}
	else
	{
		this->pairing_array = NULL;
		this->ind_unpairing_array = NULL;
		this->ind_pairing_array = NULL;
		this->fold_env = NULL;
		this->str_coinc_env = NULL;
	}

	this->n_bytes_alloced += ((sizeof(double*) * (N + 1)) + 
								(sizeof(double) * (N + 1)) + 
								(sizeof(double) * (N + 1)) +
								(sizeof(bool*) * (N + 1)) + 
								(sizeof(bool*) * (N + 1)) + 
								(sizeof(short*) * (N + 1)));

	// Allocate pairing and unpairing spf arrays.
	for(int i1 = 0; i1 <= N; i1++)
	{
		// Include the max_separation criterion in the allocation function.
		int min_i2 = i1;
		int max_i2 = MIN(i1 + ppf_cli->max_n_separation_between_nucs, N);

		if(mallocate)
		{
			this->fold_env[i1] = (bool*)malloc(sizeof(bool) * (max_i2 - min_i2 + 2));
			this->fold_env[i1] -= i1; // Do pointer shift for fold envelope.
		}

		this->n_bytes_alloced += (sizeof(bool) * (max_i2 - min_i2 + 2));

		if(mallocate)
		{
			this->pairing_array[i1] = (double*)malloc(sizeof(double) * (max_i2 - min_i2 + 2)); // Allocate pairing prob.
			this->pairing_array[i1] -= i1; // Do pointer shift to access the array using sequence indices.
		}

		this->n_bytes_alloced += (sizeof(double) * (max_i2 - min_i2 + 2));

		if(mallocate)
		{
			this->str_coinc_env[i1] = (bool*)malloc(sizeof(bool) * (max_i2 - min_i2 + 2));
			this->str_coinc_env[i1] -= i1; // Do pointer shift for fold envelope.
		}

		this->n_bytes_alloced += (sizeof(double) * (max_i2 - min_i2 + 2));

		if(mallocate)
		{
			this->ind_pairing_array[i1] = ZERO;
			this->ind_unpairing_array[i1] = ZERO;

			for(int i2 = min_i2; i2 <= max_i2; i2++)
			{
				this->pairing_array[i1][i2] = CONVERT_FROM_LIN(0.0); // Initialize the probabilities to 0.
				this->fold_env[i1][i2] = false; // Set all possible pairs.
				this->str_coinc_env[i1][i2] = false;
			}
		}
	} // i1 loop

	if(!mallocate)
	{
		return;
	}

	// Now arrays are allocated, do single partition function calculation for that sequence
	if(pairing_probs_file == NULL)
	{
        	RNA* rna = new RNA(seq_path, FILE_SEQ);
	        rna->PartitionFunction();

		// Load pairing array.
		for(int i = 1; i <= this->N; i++)
		{
			int min_j = i+1;
			int max_j = MIN(i + ppf_cli->max_n_separation_between_nucs, N);

			for(int j = min_j; j <= max_j; j++)
			{
				this->pairing_array[i][j] = rna->GetPairProbability(i, j);
			}
		}
	}
	else
	{
		// Read spf file:
		/*
		Read spf array file, the format is as following:
		1 2 0.000000000000000000000000000000
		1 3 0.000000000000000000000000000000
		1 4 0.000000000000000000000000000000
		1 5 0.000000000000000000000000000000
		1 6 0.000000000000000000000000000000
		1 7 0.000000000000000000000000000000
		1 8 0.000000000000000000000000000000
		1 9 0.000000000000000000000000000000
		1 10 0.000000000000000000000000000000
		1 11 0.000002568332606195572231287628
		...

		where each line consists of
		[index 1] [index 2] [pairing probability of two nucleotides]
		*/
		char spf_array_fn[1000];
		strcpy(spf_array_fn, pairing_probs_file);

		// Read file, read all lines, # of lines read must be equal to # of nucleotides in sequence.
		FILE* spf_file = open_f(spf_array_fn, "rb");

		if(spf_file == NULL)
		{
			printf("Could not open single partition function %s @ %s(%d)\n", spf_array_fn, __FILE__, __LINE__);
			exit(0);
		}

		// SPF file do not contain all the (i1, i2) pairs, it rather includes 
		// i1, i2 pairs where i1 < i2. However in ppf calculations, 
		int i1 = 0;
		int i2 = 0;
		double current_lin_prob = ZERO;
		int n_samples = 0;
		int n_curr_pp_cnt = 0;

		if(fread(&n_samples, sizeof(int), 1, spf_file) != 1)
		{
			printf("Could not read number of samples from %s\n", spf_array_fn);
			exit(0);
		}
		else
		{
			printf("%d samples are processed to estimate base pairing probabilities.\n", n_samples);
		}

		// 1 11 0.000002568332606195572231287628
		while(true)
		{
			if(fread(&i1, sizeof(int), 1, spf_file) != 1)
			{
				break;
			}
		
			if(fread(&i2, sizeof(int), 1, spf_file) != 1)
			{
				printf("Could not read i2 for i1=%d in %s\n", i1, spf_array_fn);
				exit(0);
			}

			if(fread(&n_curr_pp_cnt, sizeof(int), 1, spf_file) != 1)
			{
				printf("Could not read i2 for i1=%d in %s\n", i1, spf_array_fn);
				exit(0);
			}

			current_lin_prob = (double)n_curr_pp_cnt / (double)n_samples;

			// Check max_separation criterion.
			if(i2 > i1 && (i2 - i1) <= ppf_cli->max_n_separation_between_nucs)
			{
				// It should be noted that probabilities in spf file are linear, might need to change them.
				this->pairing_array[i1][i2] = CONVERT_FROM_LIN(current_lin_prob);

	if(_DUMP_SPF_MESSAGES_)
				printf("pp(%d, %d) = %.25f\n", i1, i2, current_lin_prob);

				// If the pairing probability is smaller than fold_env_prob_treshold, set fold envelope for this to 0.
				if(this->pairing_array[i1][i2] >= CONVERT_FROM_LIN(ppf_cli->fold_env_prob_treshold))
				{
					this->fold_env[i1][i2] = true;
				}
			}
			else
			{
				// Out of bounds, do not set the value here since it is not allocated.
			}

	if(_DUMP_SPF_MESSAGES_)
			printf("P_pair(%d, %d) = %.25f\n", i1, i2, this->pairing_array[i1][i2]);

			//fscanf(spf_file, "%d %d %lf", &i1, &i2, &current_lin_prob);
			//printf("read %d %d\n", i1, i2);
		}

		fclose(spf_file);
	} // read the pairing probabilities from external file.

	// Compute the pairing and coincidence ptr relocation maps with base pairing enforced for pairs that have 0.999 or higher probability of pairing.
	this->folding_constraints = new t_folding_constraints(seq_path, this->pairing_array, 0.999f);

	// Weigh all pairing probabilities with a factor in log domain to
	// decrease affect of positive feedback.
	this->calculate_unpairing_probs();

if(_DUMP_SPF_MESSAGES_)
	printf("t_spf_array allocated %lf bytes\n", this->n_bytes_alloced);

//	// Dump spf plane if desired.
//if(_DUMP_SPF_PLANES_)
//{
//	this->dump_spf_plane();
//	this->dump_fold_env();
//}
}

// Calculate unpairing probabilities, call this function after loading spf file.
// This shouldnt be called after uniform prior score setting but it will give an error 
// and exit anyway, i.e., must be called after reading PROBABILITIES as scores from spf calculations,
// because this code assumes that ind_unpairing_prob[i] + ind_pairing_prob[i] = 1.0
void t_spf_array::calculate_unpairing_probs()
{
	// Calculate individual pairing probability of all bases.
	for(int i1 = 1; i1 <= this->N; i1++)
	{	
		// Note that individual pairing and unpairing probabilities are already initied to CONVERT_FROM_LIN(0.0) in allocation.
		// for pairing_array[i1][i2], this has two pairing probabilites: i1-i2 pair and i2-i1 pair. Since two
		// loops do not account for cases where i1 > i2, that is, only accounts half of i1xi2 plane, must account
		// for both pairs in one case.
		for(int i2 = i1 + 1; i2 <= this->N; i2++)
		{
			// It is important not to use px for accessing pairing probabilities since px returns
			// scaled probabilities in linear domain.
			this->ind_pairing_array[i1] = SUM(this->ind_pairing_array[i1], this->pairing_array[i1][i2]);
			this->ind_pairing_array[i2] = SUM(this->ind_pairing_array[i2], this->pairing_array[i1][i2]);
		}

if(_DUMP_SPF_MESSAGES_)
		printf("ind pairing prob of %d: %.25f\n", i1, this->ind_pairing_array[i1]);

		// individual unpairing probability is 1-ind_pairing_probability.
		this->ind_unpairing_array[i1] = SUB(CONVERT_FROM_LIN(1.0), this->ind_pairing_array[i1]);

if(_DUMP_SPF_MESSAGES_)
		printf("Unpairing probability of %d. nucleotide is %.25f (%.5f)\n", i1, this->ind_unpairing_array[i1], (this->ind_unpairing_array[i1]));
	}
}

// Dump single pairing probabilities (score) in a matlab loadable format.
//void t_spf_array::dump_spf_plane()
//{
//	return;
//}

//void t_spf_array::dump_fold_env()
//{
//	char seq_fold_env_fp[500];
//	if(this->id == 0)
//	{
//		//sprintf(seq_fold_env_fp, "seq1_fold_env");
//		//dump_entropy(this->N, this->pairing_array, "seq1_spf_entropy.txt");
//		concat_path_fn(seq_fold_env_fp, DUMPING_DIR, "seq1_fold_env");
//	}
//	else if(this->id == 1)
//	{
//		//sprintf(seq_fold_env_fp, "seq2_fold_env");
//		//dump_entropy(this->N, this->pairing_array, "seq2_spf_entropy.txt");
//		concat_path_fn(seq_fold_env_fp, DUMPING_DIR, "seq2_fold_env");
//	}
//	else
//	{
//		printf("Illegal number of spf's instantiated, %d of'em...", n_spfs);
//		exit(0);
//	}
//
//	FILE* fold_env_dump_file = open_f(seq_fold_env_fp, "w");
//
//	// Dump pairing_array.
//	for(int i1 = 1; i1 <= this->N; i1++)
//	{
//		for(int i2 = 1; i2 <= this->N; i2++)
//		{
//			if(i2 > i1)
//			{
//				fprintf(fold_env_dump_file, "%d ", this->fold_env[i1][i2]);
//			}
//			else if(i1 > i2)
//			{
//				fprintf(fold_env_dump_file, "%d ", this->fold_env[i2][i1]);
//			}
//			else // i1 == i2
//			{
//				fprintf(fold_env_dump_file, "0 ");
//			}
//		}
//
//		fprintf(fold_env_dump_file, "\n");
//	}
//
//	fclose(fold_env_dump_file);
//}

void t_spf_array::set_scaler(double* _ppf_scales)
{
	this->ppf_scales = (double*)malloc(sizeof(double) * (this->N + 2));

	for(int cnt = 0; cnt <= this->N; cnt++)
	{
		this->ppf_scales[cnt] = 1.0;
	}

	if(_ppf_scales != NULL)
	{
		for(int cnt = 0; cnt <= this->N; cnt++)
		{
			this->ppf_scales[cnt] = _ppf_scales[cnt];
		}
	}
}

// Return pairability of nucs at indices i1 and i2.
bool t_spf_array::is_pairable(int i1, int i2)
{
	if(i2 > this->N)
	{
		return(this->fold_env[i2-N][i1]);
	}
	else
	{
		return(this->fold_env[i1][i2]);
	}
}

	// Note that acessing function does not have to return a reference since these 
	// values are coming from single partition funtion calculations.

	// px: Returns pairing score for an alignment position (i1, i2).
#ifdef _LINEAR_COMPUTATIONS_
double t_spf_array::px(int i1, int i2) 
{
	if(i1 >= i2)
	{
		double* p = 0;
		*p = 0;
		return(0.0);
		//return(0.0);			
	}

	if(i1 > N)
	{
		double* p = 0;
		*p = 0;
		return(0);
	}

	if(i2 <= N)
	{
		return(this->pairing_array[i1][i2] * this->ppf_scales[i1] * this->ppf_scales[i2]);
	}
	else
	{
		double* p = 0;
		*p = 0;
	}
}

double t_spf_array::ux_5p(int i1, int i2) 
{
	if(i1 <= N)
	{
		return(this->ind_unpairing_array[i1] * this->ppf_scales[i1]);
	}
	else
	{
		return(this->ind_unpairing_array[i1 - this->N] * this->ppf_scales[i1 - this->N]);
	}
}

// ux_j: Returns unpairing score for a nucleotide position.
double t_spf_array::ux_3p(int i1, int i2) 
{
	if(i2 <= N)
	{
		return(this->ind_unpairing_array[i2] * this->ppf_scales[i2]);
	}
	else
	{
		return(this->ind_unpairing_array[i2 - this->N] * this->ppf_scales[i2 - this->N]);
	}
}

#endif // _LINEAR_COMPUTATIONS_

#ifdef _LOG_COMPUTATIONS_
double t_spf_array::px(int i1, int i2)
{
	if(i1 >= i2)
	{
		double* p = 0;
		*p = 0;
	}

	if(i1 > N)
	{
		double* p = 0;
		*p = 0;
	}

	if((i2 - i1) > this->ppf_cli->max_n_separation_between_nucs)
	{
		return(ZERO);
	}

	if(i2 <= N)
	{
		return(this->pairing_array[i1][i2]);
	}
	else if(i2 > N)
	{
		double* p = 0;
		*p = 0;
	}
}

// ux_i: Returns unpairing score for a nucleotide position.
double t_spf_array::ux_5p(int i1, int i2) 
{
	if(i2 <= N)
	{
		if((i2 - i1) > this->ppf_cli->max_n_separation_between_nucs)
		{
			return(ZERO);
		}

		return(this->ind_unpairing_array[i1]);
	}
	else
	{
		double* p = 0;
		*p = 0;
	}
}

// ux_j: Returns unpairing score for a nucleotide position.
double t_spf_array::ux_3p(int i1, int i2) 
{
	if(i2 <= N)
	{
		if((i2 - i1) > this->ppf_cli->max_n_separation_between_nucs)
		{
			return(ZERO);
		}

		return(this->ind_unpairing_array[i2]);
	}
	else
	{
		double* p = 0;
		*p = 0;
	}
}
#endif // _LOG_COMPUTATIONS_




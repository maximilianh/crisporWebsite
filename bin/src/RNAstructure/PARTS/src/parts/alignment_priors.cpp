#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parts_compilation_directives.h"
#include "alignment_priors.h"
#include "ppf_math.h"
#include "ppf_cli.h"
#include "phmm_parameters.h"
#include "template_pf_array.h"
#include "process_sequences.h"
#include "parts_paths.h"
#include "math.h"

#include <string.h>
#include <limits.h>

#include "../../../src/phmm/utils/file/utils.h"

#include <iostream>

using namespace std;

/*
class t_aln_priors
{
public:
	int N1;
	int N2;
	double*** priors; // Accession is the same as forward-backward arrays: [seq1 index][seq2 index][state].

	// Constructor can allocate and initialize prior array.
	t_aln_priors(int N1, int N2); // Allocate empty priors.
	t_aln_priors(char* priors_file); // Load priors from a file, does not need lengths since it determines lengths from 

	inline double x(int i1, int i2, int state)
	{
		return(priors[i1][i2][state]);
	}
};
*/

bool _DUMP_PPF_ALN_PRIOR_MESSAGES_ = false;

t_aln_priors::t_aln_priors(t_seq_man* seq_man, bool mallocate)
{
	// Get sequence length.
	this->N1 = seq_man->get_l_seq1();
	this->N2 = seq_man->get_l_seq2();

if(_DUMP_PPF_ALN_PRIOR_MESSAGES_)
	printf("Allocating alignment priors with %d X %d array with uniform priors on alignments.\n", N1, N2);

	this->alloc_array(mallocate);
}

t_aln_priors::~t_aln_priors()
{
	this->free_array();
}

// Loads the priors from files.
t_aln_priors::t_aln_priors(t_seq_man* seq_man, t_ppf_cli* _ppf_cli, bool mallocate)
{
	this->ppf_cli = _ppf_cli; 

	// Get sequence length.
	this->N1 = seq_man->get_l_seq1();
	this->N2 = seq_man->get_l_seq2();

if(_DUMP_PPF_ALN_PRIOR_MESSAGES_)
	printf("Allocating alignment priors with %d X %d array with phmm priors on alignments.\n", N1, N2);

	this->alloc_array(mallocate);

	char aln_prob_matrix_fn[__MAX_PATH];
	char ins1_prob_matrix_fn[__MAX_PATH];
	char ins2_prob_matrix_fn[__MAX_PATH];

	// Set up file names to read alignment prioris.
	strcpy(aln_prob_matrix_fn, seq_man->ppf_cli->aln_prior_path);
	strcpy(ins1_prob_matrix_fn, seq_man->ppf_cli->ins1_prior_path);
	strcpy(ins2_prob_matrix_fn, seq_man->ppf_cli->ins2_prior_path);
	
	// Open priori files for reading.
	FILE* f_aln_probs = open_f(aln_prob_matrix_fn, "rb");
	FILE* f_ins1_probs = open_f(ins1_prob_matrix_fn, "rb");
	FILE* f_ins2_probs = open_f(ins2_prob_matrix_fn, "rb");

	while(1)
	{
		int i;
		int k;
		double current_prob;
		if(fread(&i, sizeof(int), 1, f_aln_probs) != 1)
		{
			break;
		}

		if(fread(&k, sizeof(int), 1, f_aln_probs) != 1)
		{
			printf("Could not read k for i = %d in %s\n", i, aln_prob_matrix_fn);
			exit(0);
		}

		if(fread(&current_prob, sizeof(double), 1, f_aln_probs) != 1)
		{
			printf("Could not read prob for (i,k) = (%d, %d) in %s\n", i, k, aln_prob_matrix_fn);
			exit(0);
		}

		printf("ALN(%d, %d) = %.15f\n", i,k, current_prob);
		int min_k = t_template_pf_array::low_limits[i];
		int max_k = t_template_pf_array::high_limits[i];

		if(k <= max_k && k >= min_k)
		{
			this->priors[i][k][STATE_ALN] = CONVERT_FROM_LOG(current_prob * this->ppf_cli->log_aln_weight_per_aln_pos);
		}
	}

	while(1)
	{
		int i;
		int k;
		double current_prob;
		if(fread(&i, sizeof(int), 1, f_ins1_probs) != 1)
		{
			break;
		}

		if(fread(&k, sizeof(int), 1, f_ins1_probs) != 1)
		{
			printf("Could not read k for i = %d in %s\n", i, aln_prob_matrix_fn);
			exit(0);
		}

		if(fread(&current_prob, sizeof(double), 1, f_ins1_probs) != 1)
		{
			printf("Could not read prob for (i,k) = (%d, %d) in %s\n", i, k, aln_prob_matrix_fn);
			exit(0);
		}

		printf("INS1(%d, %d) = %.15f\n", i,k, current_prob);
		int min_k = t_template_pf_array::low_limits[i];
		int max_k = t_template_pf_array::high_limits[i];

		if(k <= max_k && k >= min_k)
		{
			this->priors[i][k][STATE_INS1] = CONVERT_FROM_LOG(current_prob * this->ppf_cli->log_aln_weight_per_aln_pos);
		}
	}

	while(1)
	{
		int i;
		int k;
		double current_prob;
		if(fread(&i, sizeof(int), 1, f_ins2_probs) != 1)
		{
			break;
		}

		if(fread(&k, sizeof(int), 1, f_ins2_probs) != 1)
		{
			printf("Could not read k for i = %d in %s\n", i, aln_prob_matrix_fn);
			exit(0);
		}

		if(fread(&current_prob, sizeof(double), 1, f_ins2_probs) != 1)
		{
			printf("Could not read prob for (i,k) = (%d, %d) in %s\n", i, k, aln_prob_matrix_fn);
			exit(0);
		}

		printf("INS2(%d, %d) = %.15f\n", i,k, current_prob);
		int min_k = t_template_pf_array::low_limits[i];
		int max_k = t_template_pf_array::high_limits[i];

		if(k <= max_k && k >= min_k)
		{
			this->priors[i][k][STATE_INS2] = CONVERT_FROM_LOG(current_prob * this->ppf_cli->log_aln_weight_per_aln_pos);
		}
	}

	fclose(f_aln_probs);
	fclose(f_ins1_probs);
	fclose(f_ins2_probs);

	// Open priori files for reading.
	f_aln_probs = open_f("aln_dump.txt", "w");
	f_ins1_probs = open_f("ins1_dump.txt", "w");
	f_ins2_probs = open_f("ins2_dump.txt", "w");

	for(int i = this->N1 - 100; i <= this->N1; i++)
	{
		for(int k = this->N2 - 100; k <= this->N2; k++)
		{
			int min_k = t_template_pf_array::low_limits[i];
			int max_k = t_template_pf_array::high_limits[i];

			if(k <= max_k && k >= min_k)		
			{
				fprintf(f_aln_probs, "%.5f", this->priors[i][k][STATE_ALN]);
				fprintf(f_ins1_probs, "%.5f", this->priors[i][k][STATE_INS1]);
				fprintf(f_ins2_probs, "%.5f", this->priors[i][k][STATE_INS2]);
			}
		}

		fprintf(f_aln_probs, "\n");
		fprintf(f_ins1_probs, "\n");
		fprintf(f_ins2_probs, "\n");
	}

	fclose(f_aln_probs);
	fclose(f_ins1_probs);
	fclose(f_ins2_probs);
}

t_aln_priors::t_aln_priors(t_seq_man* seq_man, 
						   t_ppf_cli* _ppf_cli, 
						   double** aln_probs, 
						   double** ins1_probs, 
						   double** ins2_probs, 
						   bool mallocate)
{
	// Get sequence length.
	this->N1 = seq_man->get_l_seq1();
	this->N2 = seq_man->get_l_seq2();

	this->ppf_cli = _ppf_cli;

	this->alloc_array(mallocate);

	if(!mallocate)
	{
		return;
	}

	for(int i1 = 0; i1 <= N1; i1++)
	{
		//this->priors[i1] = (double**)malloc(sizeof(double*) * (N2 + 3));
		int min_i2 = t_template_pf_array::low_limits[i1];
		int max_i2 = t_template_pf_array::high_limits[i1];

		for(int i2 = min_i2; i2 <= max_i2; i2++)
		{
			// Apply log_aln_weight_per_aln_pos.
			// Note that this weight is multiplied with alignment position log-score of each alignment plane position.
			// So since the scores that are read from hmm output are already logarithmic, the logarithmic weight can be applied
			// to score of each alignment position directly with following. Note that multiplication is not MUL but 
			// normal linear multiplication by definition of log_aln_weight_per_aln_pos.
			this->priors[i1][i2][STATE_ALN] = aln_probs[i1][i2] * this->ppf_cli->log_aln_weight_per_aln_pos;
			this->priors[i1][i2][STATE_INS1] = ins1_probs[i1][i2] * this->ppf_cli->log_aln_weight_per_aln_pos;
			this->priors[i1][i2][STATE_INS2] = ins2_probs[i1][i2] * this->ppf_cli->log_aln_weight_per_aln_pos;

			// Convert read values, which are in log domain into ppf calculation domain.
			// This is also required for spf priors.
			this->priors[i1][i2][STATE_ALN] = CONVERT_FROM_LOG(this->priors[i1][i2][STATE_ALN]);
			this->priors[i1][i2][STATE_INS1] = CONVERT_FROM_LOG(this->priors[i1][i2][STATE_INS1]);
			this->priors[i1][i2][STATE_INS2] = CONVERT_FROM_LOG(this->priors[i1][i2][STATE_INS2]);

			//printf("%d %d %.15f\n", i1, i2, aln_probs[i1][i2]);
			//printf("%d %d %.15f\n", i1, i2, ins1_probs[i1][i2]);
			//printf("%d %d %.15f\n", i1, i2, ins2_probs[i1][i2]);

		} // i2 loop.
	} // i1 loop.
}

// Exits on error, utilizes the alignment constraints from the phmm envelope computation.
void t_aln_priors::alloc_array(bool mallocate)
{
	this->priors = (double***)malloc(sizeof(double**) * (this->N1 + 3));

	this->n_bytes_alloced = 0.0f;

	if(this->priors == NULL)
	{
		printf("Could not allocate alignment prior array for %d nucleotides, exiting @ %s(%d)\n", N1, __FILE__, __LINE__);
		getc(stdin);
		exit(0);
	}

	for(int i1 = 0; i1 <= N1; i1++)
	{
		int min_i2 = t_template_pf_array::low_limits[i1];
		int max_i2 = t_template_pf_array::high_limits[i1];

		if(mallocate)
		{
			this->priors[i1] = (double**)malloc(sizeof(double*) * (max_i2 - min_i2 + 1));

			if(this->priors[i1] == NULL)
			{
				printf("Could not allocate alignment prior array, exiting @ %s(%d)\n", __FILE__, __LINE__);
				exit(0);
			}
		}

		this->n_bytes_alloced += (sizeof(double*) * (max_i2 - min_i2 + 1));

		this->priors[i1] -= min_i2;

		for(int i2 = min_i2; i2 <= max_i2; i2++)
		{
			if(mallocate)
			{
				this->priors[i1][i2] = (double*)malloc(sizeof(double) * N_STATES);

				if(this->priors[i1][i2] == NULL)
				{
					printf("Could not allocate alignment prior array, exiting @ %s(%d)\n", __FILE__, __LINE__);
					exit(0);
				}

				for(int state_cnt = 0; state_cnt < N_STATES; state_cnt++)
				{				
					this->priors[i1][i2][state_cnt] = CONVERT_FROM_LIN(0.0); // LOG SCORES!!!
				}
			} // mallocate = true

			this->n_bytes_alloced += (sizeof(double) * N_STATES);
		} // i2 loop 
	} // i1 loop

if(_DUMP_PPF_ALN_PRIOR_MESSAGES_)
	printf("t_aln_priors allocated %lf bytes\n", this->n_bytes_alloced);
}

void t_aln_priors::free_array()
{
	//this->priors = (double***)malloc(sizeof(double**) * (this->N1 + 3));

	////this->n_bytes_alloced = 0.0f;

	//if(this->priors == NULL)
	//{
	//	printf("Could not allocate alignment prior array for %d nucleotides, exiting @ %s(%d)\n", N1, __FILE__, __LINE__);
	//	getc(stdin);
	//	exit(0);
	//}

	for(int i1 = 0; i1 <= N1; i1++)
	{
		int min_i2 = t_template_pf_array::low_limits[i1];
		int max_i2 = t_template_pf_array::high_limits[i1];

		//if(mallocate)
		//{
		//	this->priors[i1] = (double**)malloc(sizeof(double*) * (max_i2 - min_i2 + 1));

		//	if(this->priors[i1] == NULL)
		//	{
		//		printf("Could not allocate alignment prior array, exiting @ %s(%d)\n", __FILE__, __LINE__);
		//		exit(0);
		//	}
		//}

		//this->priors[i1] -= min_i2;

		for(int i2 = min_i2; i2 <= max_i2; i2++)
		{
			//if(mallocate)
			//{
				//this->priors[i1][i2] = (double*)malloc(sizeof(double) * N_STATES);
				free(this->priors[i1][i2]);

				//if(this->priors[i1][i2] == NULL)
				//{
				//	printf("Could not allocate alignment prior array, exiting @ %s(%d)\n", __FILE__, __LINE__);
				//	exit(0);
				//}

				//for(int state_cnt = 0; state_cnt < N_STATES; state_cnt++)
				//{				
				//	this->priors[i1][i2][state_cnt] = CONVERT_FROM_LIN(0.0); // LOG SCORES!!!
				//}
			//} // mallocate = true
		} // i2 loop 

		this->priors[i1] += min_i2;
		free(this->priors[i1]);
	} // i1 loop

	free(this->priors);
}

double t_aln_priors::x(int i1, int i2, int state)
{
	// Give an error and exit if requested alignment prior is aligning an external nuc and an internal nuc.
	// This means that there is a problem with the loop limits.
	if((i1 > N1 && i2 < N2) || (i1 < N1 && i2 > N2))
	{
		printf("Algorithm is aligning an internal base (%d) and an external base (%d) @ %s(%d).\n", i1, i2, __FILE__, __LINE__);
		double* p = NULL;
		*p = 0;
		while(1){};
		//exit(0);
	}
/*
	if(i1 > len1 && i2 > len2)
	{
		printf("Aligning an external to an external or what??? %d <-> %d\n", i1, i2);
	}
*/
	if(i1 > N1)
	{
		i1 -= N1;
	}

	if(i2 > N2)
	{
		i2 -= N2;
	}

	return(priors[i1][i2][state]);
}

#include <stdio.h>
#include <stdlib.h>
#include "phmm.h"
#include <string.h>

#include "structure/structure_object.h"
#include "utils/xmath/log/xlog_math.h"
#include "utils/file/utils.h"

char t_phmm::state_names[N_STATES][100] = {"STATE_INS1", "STATE_INS2", "STATE_ALN"};

t_phmm::t_phmm(double new_emission_probs[N_OUTPUTS][N_STATES], double new_trans_probs[N_STATES][N_STATES])
{
	this->alloc_init_params();

	// Copy transition matrix.
	for(int cnt1 = 0; cnt1 < N_STATES; cnt1++)
	{
		for(int cnt2 = 0; cnt2 < N_STATES; cnt2++)
		{
			this->trans_probs[cnt1][cnt2] = xlog(new_trans_probs[cnt1][cnt2]);
		} // cnt2 loop 
	} // cnt1 loop 

	// Copy emission probabilities.
	for(int cnt1 = 0; cnt1 < N_OUTPUTS; cnt1++)
	{
		for(int cnt2 = 0; cnt2 < N_STATES; cnt2++)
		{
			this->emission_probs[cnt1][cnt2] = xlog(new_emission_probs[cnt1][cnt2]);
		} // cnt2 loop 
	} // cnt1 loop 
}

t_phmm::t_phmm(const char* phmm_pars_file)
{
	this->alloc_init_params();

	// Read the parameters file. This parameters file contains 10 different parameter sets and thresholds.
	FILE* fam_par_file = open_f(phmm_pars_file, "r");

	if(fam_par_file == NULL)
	{
		double* p=0;
		*p = 0;
		printf("Cannot find phmm parameters file, exiting @ %s(%d).\n", __FILE__, __LINE__);
		exit(0);
	}

	// Load all parameters from file.
	for(int cnt = 0; cnt < N_BINZ * (N_STATES + N_OUTPUTS) * N_STATES; cnt++)
	{
		fscanf(fam_par_file, "%lf", &fam_hmm_pars[cnt]);
		//printf("%d: %.3f\n", cnt, fam_hmm_pars[cnt]);
	}

	//printf("\n\n\nReading thresholds!!!\n");
	// Read thresholds.
	for(int cnt = 0; cnt < N_BINZ; cnt++)
	{
		fscanf(fam_par_file, "%lf", &fam_thresholds[cnt]);
		//printf("%d: %f\n", cnt, fam_thresholds[cnt]);
	}

	fclose(fam_par_file);
}

void t_phmm::alloc_init_params()
{
	// Copy transition matrix.
	this->trans_probs = (double**)malloc(sizeof(double*) * (N_STATES + 2));
	for(int cnt1 = 0; cnt1 < N_STATES; cnt1++)
	{
		this->trans_probs[cnt1] = (double*)malloc(sizeof(double) * (N_STATES + 2));
		for(int cnt2 = 0; cnt2 < N_STATES; cnt2++)
		{
			trans_probs[cnt1][cnt2] = xlog(0.0f);
		} // cnt2 loop 
	} // cnt1 loop

	// Copy emission probabilities.
	this->emission_probs = (double**)malloc(sizeof(double*) * (N_OUTPUTS + 2));
	for(int cnt1 = 0; cnt1 < N_OUTPUTS; cnt1++)
	{
		this->emission_probs[cnt1] = (double*)malloc(sizeof(double) * (N_STATES + 2));
		for(int cnt2 = 0; cnt2 < N_STATES; cnt2++)
		{
			emission_probs[cnt1][cnt2] = xlog(0.0f);
		} // cnt2 loop
	} // cnt1 loop

	this->fam_hmm_pars = (double*)malloc(sizeof(double) * (N_BINZ * (N_STATES + N_OUTPUTS) * N_STATES + 2));
	this->fam_thresholds = (double*)malloc(sizeof(double) * (N_BINZ + 2));
}

void t_phmm::free_params()
{
	// Free transition matrix.
	for(int cnt1 = 0; cnt1 < N_STATES; cnt1++)
	{
		free(this->trans_probs[cnt1]);
	} // cnt1 loop
	free(this->trans_probs);

	// Free emission probabilities.
	for(int cnt1 = 0; cnt1 < N_OUTPUTS; cnt1++)
	{
		free(this->emission_probs[cnt1]);
	} // cnt1 loop
	free(this->emission_probs);

	free(this->fam_hmm_pars);

	free(this->fam_thresholds);
}

t_phmm::~t_phmm()
{
	this->free_params();
}

void t_phmm::dump_parameters()
{
	// Dump emission probabilities.
	for(int cnt1 = 0; cnt1 < N_OUTPUTS; cnt1++)
	{
		for(int cnt2 = 0; cnt2 < N_STATES; cnt2++)
		{
			printf("%.3f ", xexp(emission_probs[cnt1][cnt2]));
		}

		printf("\n");
	}

	// Dump transition probabilities.
	printf("\n");
	for(int cnt1 = 0; cnt1 < N_STATES; cnt1++)
	{
		for(int cnt2 = 0; cnt2 < N_STATES; cnt2++)
		{
			printf("%.3f ", xexp(trans_probs[cnt1][cnt2]));
		}

		printf("\n");
	}
}

void t_phmm::set_parameters_by_sim(double similarity)
{
	//int fam_par_set_index = (int)(similarity * (double)N_BINZ);
	int fam_par_set_index = get_bin_index(similarity, N_BINZ);

	// Load emission probabilities.
	// Each parameter set is (N_STATES + N_OUTPUTS) * N_STATES doubles long.
	int start_linear_index = (N_STATES + N_OUTPUTS) * N_STATES * get_bin_index(similarity, N_BINZ);
	double* par_ptr = fam_hmm_pars + start_linear_index;

	for(int cnt1 = 0; cnt1 < N_OUTPUTS; cnt1++)
	{
		for(int cnt2 = 0; cnt2 < N_STATES; cnt2++)
		{
			//emission_probs[cnt1][cnt2] = *(par_ptr + cnt1 * N_STATES + cnt2);
			emission_probs[cnt1][cnt2] = xlog(par_ptr[cnt1 * N_STATES + cnt2]);
		}
	}

	start_linear_index = (N_STATES + N_OUTPUTS) * N_STATES * fam_par_set_index + N_STATES * N_OUTPUTS;
	par_ptr = fam_hmm_pars + start_linear_index;

	// Load trans probabilities.
	for(int cnt1 = 0; cnt1 < N_STATES; cnt1++)
	{
		for(int cnt2 = 0; cnt2 < N_STATES; cnt2++)
		{
			//trans_probs[cnt1][cnt2] = *(par_ptr + cnt1 * N_STATES + cnt2);
			trans_probs[cnt1][cnt2] = xlog(par_ptr[cnt1 * N_STATES + cnt2]);
		}
	}		
}

// Get index of bin of parameters for a sequence alignment.
int t_phmm::get_bin_index(double similarity, int n_bins)
{
	if(similarity == 1.0)
	{
		return(n_bins - 1);
	}
	else
	{
		return((int)(n_bins * similarity));
	}
}

double t_phmm::get_fam_threshold(double similarity)
{
	int bin_index = get_bin_index(similarity, N_BINZ);

	return(fam_thresholds[bin_index]);
}

double t_phmm::get_trans_prob(int prev, int next)
{
	return(this->trans_probs[prev][next]);
}

double t_phmm::get_emit_prob(int sym_index, int state)
{
	return(this->emission_probs[sym_index][state]);
}


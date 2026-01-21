#include <string.h>
#include <limits.h>
#include "stoch_sampling_math.h"
#include "../parts_compilation_directives.h" // Usign same computation configuration as rest of the algorithm.
#include <stdio.h>
#include <stdlib.h>
#include "../ppf_cli.h"
#include "../phmm_parameters.h"
#include "../single_pf_array.h"
#include "../alignment_priors.h"
#include "../process_sequences.h"
#include <time.h> // For generating sed for pseudo random number generator.
#include <math.h>
#include "../ppf_math.h"
#include "../ppf_timer.h"
#include "../../../../src/phmm/utils/rng/rng.h"
#include "../../../../src/phmm/utils/rng/seed_manager.h"
#include "../ppf_loops.h"

bool _DUMP_STOCH_SAMPLER_MESSAGES_ = false;

t_stoch_sampling_math::t_stoch_sampling_math(t_ppf_loops* _ppf_loops)
{
	this->ppf_loops = _ppf_loops;

	// Copy all pointers.
	this->ppf_cli = ppf_loops->ppf_cli;
	this->seq_man = ppf_loops->seq_man;
	this->seq1_spf = ppf_loops->seq1_spf;
	this->seq2_spf = ppf_loops->seq2_spf;
	this->aln_priors = ppf_loops->aln_priors;

	// Set up sample size.
	this->sample_size = ppf_cli->n_samples;
	if(this->sample_size == 0)
		this->sample_size = 1;

	// Set the seed from time(0).
	//this->prn_gen_seed = t_seed_manager::seed_me();
	if(this->ppf_loops->ppf_cli->stoch_sampling_seed != 0xfffff)
	{
		this->prn_gen_seed = this->ppf_loops->ppf_cli->stoch_sampling_seed;
	}	
	else
	{
		this->prn_gen_seed = t_seed_manager::seed_me();
	}	
	//char seed_dump_command[100];
	//sprintf(seed_dump_command, "echo %d >> seeds.txt", this->prn_gen_seed);
	//system(seed_dump_command);

	// Create pseudo random number generator.
	this->prn_gen = new t_rng(this->prn_gen_seed);

if(_DUMP_STOCH_SAMPLER_MESSAGES_)
	printf("Allocated global_score_array with length %d\n", (this->seq_man->get_l_seq1() * this->seq_man->get_l_seq2() + 3));
}

t_stoch_sampling_math::~t_stoch_sampling_math()
{
	delete(this->prn_gen);
}

double t_stoch_sampling_math::random_double_1_0()
{
	return(CONVERT_FROM_LIN(this->prn_gen->random_double_ran3()));
}


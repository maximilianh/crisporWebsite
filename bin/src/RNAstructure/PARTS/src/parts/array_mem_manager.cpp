#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "array_mem_manager.h"
#include "process_sequences.h"
#include "template_pf_array.h"
#include "../../../src/phmm/phmm.h"
#include "../../../src/phmm/phmm_aln.h"
#include "../../../src/phmm/phmm_array.h"
#include "ppf_math.h"
#include "../../../src/phmm/utils/xmath/log/xlog_math.h"
#include "ppf_cli.h"
#include "alignment_priors.h"
#include "single_pf_array.h"
#include "../../../src/phmm/structure/structure_object.h"
#include "../../../src/phmm/structure/folding_constraints.h"

t_array_mem_manager::t_array_mem_manager(t_seq_man* _seq_man, t_ppf_cli* _ppf_cli)
{
	this->seq_man = _seq_man;
	this->ppf_cli = _ppf_cli;
}

t_aln_env_result* t_array_mem_manager::get_env_threshold_per_mem_limit(double array_mem_limit)
{	
	t_structure* seq1 = new t_structure(this->seq_man->seq1_fp);
	t_structure* seq2 = new t_structure(this->seq_man->seq2_fp);

	t_spf_array* spf1 = new t_spf_array(seq1->numofbases, this->seq_man->seq1_fp, this->ppf_cli, NULL, false);
	t_spf_array* spf2 = new t_spf_array(seq2->numofbases, this->seq_man->seq2_fp, this->ppf_cli, NULL, false);

	// Compute 3 alignment planes, and the loop limits.
	//t_phmm_aln* phmm_aln = new t_phmm_aln(this->seq_man->seq1_fp, this->seq_man->seq2_fp, this->seq_man->ppf_cli->phmm_band_constraint_size);
	t_phmm_aln* phmm_aln = new t_phmm_aln(this->seq_man->seq1_fp, this->seq_man->seq2_fp);

	// No alignment constraints.
	//t_aln_env_result* aln_env_result = phmm_aln->compute_alignment_envelope(PROB_ALN_ENV, 0);
        t_pp_result* pp_result = phmm_aln->compute_posterior_probs();
        double log_threshold = pp_result->fam_threshold;
        t_aln_env_result* aln_env_result = phmm_aln->compute_alignment_envelope(PROB_ALN_ENV, pp_result, log_threshold, 0);

	if(aln_env_result == NULL)
	{
		printf("Alignment envelope computation failed.\n");
		exit(0);
	}
		
	// Have to initialize loop limits before allocating and initing pf array.
	t_template_pf_array::alloc_init_loop_limits(seq_man, aln_env_result->low_limits, aln_env_result->high_limits);

	double per_array_mem = this->get_mem_per_array();

	t_aln_priors* aln_priors = new t_aln_priors(this->seq_man, false);

	double current_log_threshold = log_threshold; 

	while(((per_array_mem * N_ARRAYS) + spf1->n_bytes_alloced + spf2->n_bytes_alloced + aln_priors->n_bytes_alloced) > array_mem_limit)
	{
		printf("Current total usage is %lf bytes (%lf x %d + %lf + %lf + %lf) for allowed size of %lf.\n", 
			((per_array_mem * N_ARRAYS) + spf1->n_bytes_alloced + spf2->n_bytes_alloced + aln_priors->n_bytes_alloced),
			per_array_mem, N_ARRAYS,
			spf1->n_bytes_alloced,
			spf2->n_bytes_alloced,
			aln_priors->n_bytes_alloced,
			array_mem_limit);

		// Free current alignment envelope.
		//phmm_aln->free_aln_env_result(aln_env_result);
		free(aln_env_result->high_limits);
		free(aln_env_result->low_limits);

		// Increase threshold by value of threshold update.
		current_log_threshold = xlog_sum(current_log_threshold, xlog(LOG_THRESHOLD_UPDATE));

		// Recompute alignment envelope utilizing the probability planes computed in previous step, do not recompute the 
		// probability plane.
		aln_env_result = phmm_aln->compute_alignment_envelope(PROB_ALN_ENV, pp_result, current_log_threshold, 0);		

		// Note that alignment envelope computation returns NULL in the case the alignment envelope is not connected.
		if(aln_env_result == NULL)
		{
			// Decrement log threshold to make alignment envelope connected.
			current_log_threshold = xlog_sub(current_log_threshold, xlog(LOG_THRESHOLD_UPDATE));

			// Update alignment envelope result.
			aln_env_result = phmm_aln->compute_alignment_envelope(PROB_ALN_ENV, pp_result, current_log_threshold, 0);

			printf("Detected disconnected alignment envelope, returning with total usage of %lf bytes (%lf x %d + %lf + %lf + %lf) for allowed size of %lf.\n", 
					((per_array_mem * N_ARRAYS) + spf1->n_bytes_alloced + spf2->n_bytes_alloced + aln_priors->n_bytes_alloced),
					per_array_mem, N_ARRAYS,
					spf1->n_bytes_alloced,
					spf2->n_bytes_alloced,
					aln_priors->n_bytes_alloced,
					array_mem_limit);

			delete(phmm_aln);

			return(aln_env_result);
		}

		// Update loop limits.
		t_template_pf_array::update_loop_limits(seq_man, aln_env_result->low_limits, aln_env_result->high_limits);

		// Allocate arrays and check memory.
		per_array_mem = this->get_mem_per_array();
		delete(aln_priors);
		aln_priors = new t_aln_priors(this->seq_man, false);
	}

	printf("Total usage of %lf bytes (%lf x %d + %lf + %lf + %lf) for allowed size of %lf.\n", 
			((per_array_mem * N_ARRAYS) + spf1->n_bytes_alloced + spf2->n_bytes_alloced + aln_priors->n_bytes_alloced),
			per_array_mem, N_ARRAYS,
			spf1->n_bytes_alloced,
			spf2->n_bytes_alloced,
			aln_priors->n_bytes_alloced,
			array_mem_limit);

	//for(int i = 0; i <= this->seq_man->get_l_seq1(); i++)
	//{
	//	printf("%d: %d -> %d\n", i, aln_env_result->low_limits[i], aln_env_result->high_limits[i]);
	//}

	delete(seq1);
	delete(seq2);

	// Free current alignment envelope.
	delete(phmm_aln);
	return(aln_env_result);
}

double t_array_mem_manager::get_total_mem_size()
{
	t_structure* seq1 = new t_structure(this->seq_man->seq1_fp);
	t_structure* seq2 = new t_structure(this->seq_man->seq2_fp);

	printf("%s: %d nucleotides\n%s: %d nucleotides.\n", seq1->ctlabel, seq1->numofbases, seq2->ctlabel, seq2->numofbases);

	t_spf_array* spf1 = NULL;
	if(this->seq_man->ppf_cli->fold_prior1_path != NULL)
	{
		spf1 = new t_spf_array(seq1->numofbases, this->seq_man->seq1_fp, this->ppf_cli, this->seq_man->ppf_cli->fold_prior1_path, false);
		printf("1 %s: %lf megabytes allocated in spf1.\n", this->seq_man->seq1_fp, spf1->n_bytes_alloced / MEG_BYTES);
	}
	else
	{
		spf1 = new t_spf_array(seq1->numofbases, this->seq_man->seq1_fp, this->ppf_cli, NULL, false);
		printf("1 %s: %lf megabytes allocated in spf1.\n", this->seq_man->seq1_fp, spf1->n_bytes_alloced / MEG_BYTES);
	}

	t_spf_array* spf2 = NULL;
	if(this->seq_man->ppf_cli->fold_prior2_path != NULL)
	{	
		spf2 = new t_spf_array(seq2->numofbases, this->seq_man->seq2_fp, this->ppf_cli, this->seq_man->ppf_cli->fold_prior2_path, false);
		printf("2 %s: %lf megabytes allocated in spf2.\n", this->seq_man->seq2_fp, spf2->n_bytes_alloced / MEG_BYTES);
	}
	else
	{
		spf2 = new t_spf_array(seq2->numofbases, this->seq_man->seq2_fp, this->ppf_cli, NULL, false);
		printf("2 %s: %lf megabytes allocated in spf2.\n", this->seq_man->seq2_fp, spf2->n_bytes_alloced / MEG_BYTES);
	}

	if(this->ppf_cli->loop_limits_fn == NULL)
	{
		// Compute 3 alignment planes, and the loop limits.
		t_phmm_aln* phmm_aln = new t_phmm_aln(this->seq_man->seq1_fp, this->seq_man->seq2_fp);

		// No alignment constraints.
		t_aln_env_result* aln_env_result = phmm_aln->compute_alignment_envelope(PROB_ALN_ENV, 0);

		if(aln_env_result == NULL)
		{
			printf("Alignment envelope computation failed @ %s(%d).\n", __FILE__, __LINE__);
			return(-1);
		}

		// Update loop limits.
		t_template_pf_array::update_loop_limits(seq_man, aln_env_result->low_limits, aln_env_result->high_limits);
	}
	else
	{
		t_template_pf_array::update_loop_limits(seq_man, this->ppf_cli->loop_limits_fn);
	}

	t_aln_priors* aln_priors = new t_aln_priors(this->seq_man, false);

	double per_array_mem_size = this->get_mem_per_array();
	double total_mem_size = ((per_array_mem_size * N_ARRAYS) + spf1->n_bytes_alloced + spf2->n_bytes_alloced + aln_priors->n_bytes_alloced);

	if(total_mem_size > GIG_BYTES)
	{
		printf("Total usage of %lf gigabytes (%lf x %d + %lf + %lf + %lf).\n", 
				total_mem_size / GIG_BYTES,
				per_array_mem_size, N_ARRAYS,
				spf1->n_bytes_alloced,
				spf2->n_bytes_alloced,
				aln_priors->n_bytes_alloced);
	}
	else if(total_mem_size > MEG_BYTES)
	{
		printf("Total usage of %lf megabytes (%lf x %d + %lf + %lf + %lf).\n", 
				total_mem_size / MEG_BYTES,
				per_array_mem_size, N_ARRAYS,
				spf1->n_bytes_alloced,
				spf2->n_bytes_alloced,
				aln_priors->n_bytes_alloced);
	}
	else
	{
		printf("Total usage of %lf bytes (%lf x %d + %lf + %lf + %lf).\n", 
				total_mem_size,
				per_array_mem_size, N_ARRAYS,
				spf1->n_bytes_alloced,
				spf2->n_bytes_alloced,
				aln_priors->n_bytes_alloced);
	}

	delete(aln_priors);
	delete(spf1);
	delete(spf2);
	return(total_mem_size);
}

// Assuming that the loop limits are set.
double t_array_mem_manager::get_mem_per_array()
{
	t_structure* seq1 = new t_structure(this->seq_man->seq1_fp);
	t_structure* seq2 = new t_structure(this->seq_man->seq2_fp);

	// Following does not utilize the folding envelope for first sequence.
	t_spf_array* seq1_spf = NULL;
	if(this->seq_man->ppf_cli->fold_prior1_path != NULL)
	{
		seq1_spf = new t_spf_array(seq1->numofbases, this->seq_man->seq1_fp, this->ppf_cli, this->seq_man->ppf_cli->fold_prior1_path, true);
		printf("%s: %lf megabytes allocated in spf1.\n", this->seq_man->seq1_fp, seq1_spf->n_bytes_alloced / MEG_BYTES);
	}
	else
	{
		seq1_spf = new t_spf_array(seq1->numofbases, this->seq_man->seq1_fp, this->ppf_cli, NULL, true);
		printf("%s: %lf megabytes allocated in spf1.\n", this->seq_man->seq1_fp, seq1_spf->n_bytes_alloced / MEG_BYTES);
	}

	t_spf_array* seq2_spf = NULL;
	if(this->seq_man->ppf_cli->fold_prior2_path != NULL)
	{	
		seq2_spf = new t_spf_array(seq2->numofbases, this->seq_man->seq2_fp, this->ppf_cli, this->seq_man->ppf_cli->fold_prior2_path, true);
		printf("%s: %lf megabytes allocated in spf2.\n", this->seq_man->seq2_fp, seq2_spf->n_bytes_alloced / MEG_BYTES);
	}
	else
	{
		seq2_spf = new t_spf_array(seq2->numofbases, this->seq_man->seq2_fp, this->ppf_cli, NULL, true);
		printf("%s: %lf megabytes allocated in spf2.\n", this->seq_man->seq2_fp, seq2_spf->n_bytes_alloced / MEG_BYTES);
	}

	t_template_pf_array* template_pf_array = new t_template_pf_array(this->seq_man, 													
														seq1_spf->folding_constraints->coinc_pointer_relocation_map,
														seq2_spf->folding_constraints->coinc_pointer_relocation_map,
														false);

	double mem_usage = template_pf_array->n_alloced_bytes;

	delete(template_pf_array);
	delete(seq1_spf);
	delete(seq2_spf);

	return(mem_usage);
}

#include "phmm_aln.h"
#include "structure/structure_object.h"
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "utils/xmath/log/xlog_math.h"
#include "utils/file/utils.h"
#include "phmm_paths.h"
#include "phmm.h"
#include "phmm_array.h"
//#include "../constrained_ppf/paths.h"

bool _DUMP_PHMM_PP_LOOPS_MESSAGES_ = false;

// Compute, verify and return foreward-backward results.
t_pp_result* t_phmm_aln::compute_posterior_probs()
{
if(_DUMP_PHMM_PP_LOOPS_MESSAGES_)
	printf("Computing posterior probabilities of alignment.\n");
	// Compute with the constraints.
	t_ML_result* ml_result = this->compute_ML_alignment();

	// Validate and load parameter files: The loop constants.	
	const char* data_dir = resolve_data_dir();
	if(data_dir == NULL)
	{
		printf("Could not resolve thermodynamics data directory.\n");
		exit(0);
	}
	char* phmm_pars_fp = (char*)malloc(sizeof(char*) * (strlen(data_dir) + strlen(PHMM_FAM_PARS_FP) + 3));
	sprintf(phmm_pars_fp, "%s/%s", data_dir, PHMM_FAM_PARS_FP);
	this->phmm = new t_phmm(phmm_pars_fp);
	free(phmm_pars_fp);

	//printf("similarity is %f\n", ml_result->ml_similarity);
	phmm->set_parameters_by_sim(ml_result->ml_similarity);
	//phmm->dump_parameters();

	//double n_fore_array_bytes_alloced = 0.0f;
	//double n_back_array_bytes_alloced = 0.0f;
	t_phmm_array* fore_array = new t_phmm_array(l1(), l2(), this->phmm_band_constraint_size, true);

if(_DUMP_PHMM_PP_LOOPS_MESSAGES_)
	printf("Allocated %lf bytes for forward array\n", fore_array->n_bytes_alloced);

	t_phmm_array* back_array = new t_phmm_array(l1(), l2(), this->phmm_band_constraint_size, true);

if(_DUMP_PHMM_PP_LOOPS_MESSAGES_)
	printf("Allocated %lf bytes for backward array\n", back_array->n_bytes_alloced);
	
	// Initialize the arrays.
	this->init_forward_array(fore_array);
	this->init_backward_array(back_array);

	// Compute the arrays.
	this->compute_forward_array(fore_array);
	this->compute_backward_array(back_array);

	// Calculate output probabilities from forward computations.
	double fore_output_prob = LOG_OF_ZERO;

	for(int cnt = 0; cnt < N_STATES; cnt++)
	{
		double last_state_prob = this->phmm->get_trans_prob(cnt, STATE_ALN); // Prob. of finishing alignment with this state.
		//printf("%d: %.25f\n", cnt, fore_array->x(l1(), l2(), cnt));
		fore_output_prob = xlog_sum(fore_output_prob, xlog_mul(fore_array->x(l1(), l2(), cnt), last_state_prob));
	}

	// Output probability from backward computations.
	double back_output_prob = back_array->x(0,0,STATE_ALN);

if(_DUMP_PHMM_PP_LOOPS_MESSAGES_)
{
	printf("Output probability from backward array computation: %.25f\n", back_output_prob);
	printf("Output probability from forward array computation: %.25f\n", fore_output_prob);
}


	//printf("Output probability of sequences is %.5f\n", output_prob);
	if(!xlog_comp(fore_output_prob, back_output_prob))
	{
		//printf("Output probabilities from forward and backward computations are not same: %.25f\n", xexp(xlog_div(xlog_sub(fore_output_prob, back_output_prob), fore_output_prob)));
//		printf("Output probabilities from forward and backward computations are not same: %.25f\n", fore_output_prob, back_output_prob);
		exit(0);
	}

if(_DUMP_PHMM_PP_LOOPS_MESSAGES_)
{
	printf("Allocating and computing 3 planes of pairwise alignment probabilities.\n");
}
	// Compute the probabilities.
	double** aln_probs = (double**)malloc(sizeof(double*) * (l1() + 3));
	for(int i = 0; i <= l1(); i++)
	{
		int min_k = fore_array->low_phmm_array_limits[i];
		int max_k = fore_array->high_phmm_array_limits[i];

		aln_probs[i] = (double*)malloc(sizeof(double) * (max_k - min_k + 1));
		aln_probs[i] -= min_k;
		for(int k = min_k; k <= max_k; k++)
		{
			aln_probs[i][k] = xlog_div(xlog_mul(fore_array->x(i,k,STATE_ALN), back_array->x(i,k,STATE_ALN)), back_output_prob);
		}
	}

	double** ins1_probs = (double**)malloc(sizeof(double*) * (l1() + 3));
	// This loop should start with i=1 because there cannot be STATE_INS1 at i=0.
	for(int i = 0; i <= l1(); i++)
	{
		int min_k = fore_array->low_phmm_array_limits[i];
		int max_k = fore_array->high_phmm_array_limits[i];

		ins1_probs[i] = (double*)malloc(sizeof(double) * (max_k - min_k + 1));
		ins1_probs[i] -= min_k;
		for(int k = min_k; k <= max_k; k++)
		{
			ins1_probs[i][k] = xlog_div(xlog_mul(fore_array->x(i,k,STATE_INS1), back_array->x(i,k,STATE_INS1)), back_output_prob);
		}
	}

	double** ins2_probs = (double**)malloc(sizeof(double*) * (l1() + 3));
	for(int i = 0; i <= l1(); i++)
	{
		int min_k = fore_array->low_phmm_array_limits[i];
		int max_k = fore_array->high_phmm_array_limits[i];

		ins2_probs[i] = (double*)malloc(sizeof(double) * (max_k - min_k + 1));
		ins2_probs[i] -= min_k;
		for(int k = min_k; k <= max_k; k++)
		{
			ins2_probs[i][k] = xlog_div(xlog_mul(fore_array->x(i,k,STATE_INS2), back_array->x(i,k,STATE_INS2)), back_output_prob);
		}
	}

if(_DUMP_PHMM_PP_LOOPS_MESSAGES_)
	printf("Verifying forward backward computation for %s.\n", this->seq1->ctlabel);

	// Verify forward-backward computations for 1st sequence.
	for(int i = 1; i <= l1(); i++)
	{
		int min_k = fore_array->low_phmm_array_limits[i];
		int max_k = fore_array->high_phmm_array_limits[i];

		// This nucleotide in 1st sequence is either inserted or it is aligned to another nucleotide.		
		double total_aln = xlog(0.0);
		for(int k = min_k; k <= max_k; k++)
		{
			total_aln = xlog_sum(total_aln, aln_probs[i][k]);
		}

		double total_ins1 = xlog(0.0);
		for(int k = min_k; k <= max_k; k++)
		{
			total_ins1 = xlog_sum(total_ins1, ins1_probs[i][k]);
		}

		// All these above should sum up to 1.
		if(!xlog_comp(xlog_sum(total_ins1, total_aln), xlog(1.0f)))
		{
//			printf("Forward-backward computation failed for nucleotide %d of 1st sequence: %.25f <-> %.25f\n", i, xlog_sum(total_ins1, total_aln), xlog(1.0f));
			//exit(0);
		}
	}

	// Verify forward-backward computations for 2nd sequence.
if(_DUMP_PHMM_PP_LOOPS_MESSAGES_)
	printf("Verifying forward backward computation for %s.\n", this->seq2->ctlabel);

	for(int k = 1; k <= l2(); k++)
	{
		// This nucleotide in 1st sequence is either inserted or it is aligned to another nucleotide.		
		double total_aln = xlog(0.0);
		for(int i = 1; i <= l1(); i++)
		{
			if(fore_array->check_phmm_boundary(i, k))
			{
				total_aln = xlog_sum(total_aln, aln_probs[i][k]);
			}
		}

		double total_ins2 = xlog(0.0);
		for(int i = 0; i <= l1(); i++)
		{
			if(fore_array->check_phmm_boundary(i, k))
			{
				total_ins2 = xlog_sum(total_ins2, ins2_probs[i][k]);
			}
		}

		// All these above should sum up to 1.
		if(!xlog_comp(xlog_sum(total_ins2, total_aln), xlog(1.0)))
		{
//			printf("Forward backward computation failed for nucleotide %d of 2nd sequence: %.5f <-> %.5f\n", k, xlog_sum(total_ins2, total_aln), back_output_prob);
			//exit(0);
		}
	}

	// Save the results in pp_results structure.
	t_pp_result* pp_result = (t_pp_result*)malloc(sizeof(t_pp_result));
	pp_result->aln_probs = aln_probs;
	pp_result->ins1_probs = ins1_probs;
	pp_result->ins2_probs = ins2_probs;
	pp_result->op_prob = back_output_prob;
	pp_result->ml_similarity = ml_result->ml_similarity;
	pp_result->fam_threshold = phmm->get_fam_threshold(pp_result->ml_similarity);
	pp_result->phmm_band_constraint_size = fore_array->phmm_band_constraint_size;
	
	// Free array memories since they are not needed any more.
	// Free anything tha does not go to the results. 
	// The object should be self-contained such that any computation that is
	// dont after this should be independent of this one.
	this->free_ML_result(ml_result);	
	delete(fore_array);
	delete(back_array);
	delete(this->phmm);

	return(pp_result);
}

void* t_phmm_aln::free_pp_result(t_pp_result* pp_result)
{
	for(int i = 0; i <= l1(); i++)
	{
		int min_k = t_phmm_array::low_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);
		pp_result->aln_probs[i] += min_k;
		free(pp_result->aln_probs[i]);
	}
	free(pp_result->aln_probs);

	//pp_result->ins1_probs = ins1_probs;
	for(int i = 0; i <= l1(); i++)
	{
		int min_k = t_phmm_array::low_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);
		pp_result->ins1_probs[i] += min_k;
		free(pp_result->ins1_probs[i]);
	}
	free(pp_result->ins1_probs);

	//pp_result->ins2_probs = ins2_probs;
	for(int i = 0; i <= l1(); i++)
	{
		int min_k = t_phmm_array::low_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);
		pp_result->ins2_probs[i] += min_k;
		free(pp_result->ins2_probs[i]);
	}
	free(pp_result->ins2_probs);

	//pp_result->op_prob = back_output_prob;
	//pp_result->ml_similarity = ml_result->ml_similarity;
	//pp_result->fam_threshold = phmm->get_fam_threshold(pp_result->ml_similarity);

	free(pp_result);

	// Return NULL on success.
	return(NULL);
}

// Parameter to following function should be a forward hmm array.
void t_phmm_aln::init_forward_array(t_phmm_array* fore_array)
{
	// Set 0,0 alignment as 1 prob.
	fore_array->x(0, 0, STATE_START) = xlog(1.0);
	fore_array->x(0, 0, STATE_INS1) = xlog(0);
	fore_array->x(0, 0, STATE_INS2) = xlog(0);
}

void t_phmm_aln::init_backward_array(t_phmm_array* back_array)
{
	int n1 = l1();
	int n2 = l2();

	// Set probability of STATE_END state at the end as 1.0.
	back_array->x(n1+1, n2+1, STATE_END) = xlog(1.0);
	back_array->x(n1+1, n2+1, STATE_INS1) = xlog(0);
	back_array->x(n1+1, n2+1, STATE_INS2) = xlog(0);
}

void t_phmm_aln::compute_forward_array(t_phmm_array* fore_array)
{
	//printf("Calculating forward variable...\n");

	// Replicate the constraints for fast processing of constraints.
	//int* seq2_aln_const = this->get_seq2_aln_const(seq1_aln_const);

	// Order of calculation is important. Do not recalculate boundaries, use them.
	int n1 = l1();
	int n2 = l2();

	for(int i = 0; i <= n1; i++)
	{
		int low_k = fore_array->low_phmm_array_limits[i];
		int high_k = fore_array->high_phmm_array_limits[i];
		for(int k = low_k; k <= high_k; k++)
		{		
			// If this is a constrained alignment position, only compute the STATE_ALN probability passing through here.
			bool forbid_STATE_INS1 = false;
			bool forbid_STATE_INS2 = false;
			bool forbid_STATE_ALN = false;

			this->get_aln_permissions(forbid_STATE_ALN, 
										forbid_STATE_INS1, 
										forbid_STATE_INS2, 
										i,
										k);

if(_DUMP_PHMM_PP_LOOPS_MESSAGES_)
			printf("fore(%d, x)\r", i);

			// This loop is for iterating over possible states in this alignment pair.
			for(int current_state = 0; current_state < N_STATES; current_state++)
			{
				if(i != 0 || k != 0)
				{
					fore_array->x(i, k, current_state) = xlog(0);
				}

				// This loop is for marginalizing over previous state.
				for(int prev_state = 0; prev_state < N_STATES; prev_state++)
				{
					if(!forbid_STATE_ALN && 
						current_state == STATE_ALN && 
						i > 0 && 
						k > 0 &&
						fore_array->check_phmm_boundary(i-1, k-1))
					{

//                                            cerr << "forward aln " << i << " " << k << "\n";
                                            //printf("prev_state = %s\n", state_names[prev_state]);
						//printf("%s -> %s\n", state_names[prev_state), state_names[current_state]);
                                            double trans_emit_prob = xlog_mul(get_trans_emit_prob(prev_state, current_state, i, k), this->get_match_prior(i,k,n1,n2));

						fore_array->x(i,k,current_state) = xlog_sum(fore_array->x(i,k,current_state), xlog_mul(fore_array->x(i-1,k-1,prev_state), trans_emit_prob));							
					}
					
					if(!forbid_STATE_INS1 && 
						current_state == STATE_INS1 && 
						i > 0 &&
						fore_array->check_phmm_boundary(i-1, k))
					{
						//printf("prev_state = %s\n", state_names[prev_state]);
						//printf("%s -> %s\n", state_names[prev_state), state_names[current_state]);
                                            double trans_emit_prob = xlog_mul(get_trans_emit_prob(prev_state, current_state, i, k), xlog(1)/*this->get_indel_prior(i,0,n1,n2)*/ );

						fore_array->x(i,k,current_state) = xlog_sum(fore_array->x(i,k,current_state), xlog_mul(fore_array->x(i-1,k,prev_state), trans_emit_prob));
					}
					
					if(!forbid_STATE_INS2 &&
						current_state == STATE_INS2 && 
						k > 0 &&
						fore_array->check_phmm_boundary(i, k-1))
					{
						//printf("prev_state = %s\n", state_names[prev_state]);
						//printf("%s -> %s\n", state_names[prev_state), state_names[current_state]);
                                            double trans_emit_prob = xlog_mul(get_trans_emit_prob(prev_state, current_state, i, k), xlog(1)/*this->get_indel_prior(0,k,n1,n2)*/ );

						fore_array->x(i,k,current_state) = xlog_sum(fore_array->x(i,k,current_state), xlog_mul(fore_array->x(i,k-1,prev_state), trans_emit_prob));
					}
					//else
					//{
					//	printf("Problem while computing forward variable, code has seen the thing that should not be @ %s(%d)\n", __FILE__, __LINE__);
					//	printf("%d, %d\n", i, k);
					//	//exit(0); // This is a fatal error.
					//}
				}
			} // State loop.
		} // k loop.
	} // i loop.
}

// Note that I am not considering emission of current state, considering emission of next state
// based on current set of observed outputs.
void t_phmm_aln::compute_backward_array(t_phmm_array* back_array)
{
	// Replicate the constraints for fast processing of constraints.
	//int* seq2_aln_const = this->get_seq2_aln_const(seq1_aln_const);
    int n1 = l1();
    int n2 = l2();


	//printf("Computing back array.\n");
	for(int i = l1(); i >= 0; i--)
	{
		int min_k = back_array->low_phmm_array_limits[i];
		int max_k = back_array->high_phmm_array_limits[i];

if(_DUMP_PHMM_PP_LOOPS_MESSAGES_)
		printf("back(%d, x)\r", i);

		for(int k = max_k; k >= min_k; k--)
		{
			//printf("back comp(%d, %d)\n", i, k);
			bool forbid_STATE_INS1 = false;
			bool forbid_STATE_INS2 = false;
			bool forbid_STATE_ALN = false;

			this->get_aln_permissions(forbid_STATE_ALN, 
										forbid_STATE_INS1, 
										forbid_STATE_INS2, 
										i,
										k);

			// This loop is for iterating over possible states in this alignment pair.
			for(int current_state = 0; current_state < N_STATES; current_state++)
			{
				//printf("------------------------------------\n");
				back_array->x(i, k, current_state) = xlog(0);

				bool forbid_current_state_comp = false;
				if((current_state == STATE_INS1 && forbid_STATE_INS1) ||
					(current_state == STATE_INS2 && forbid_STATE_INS2) ||
					(current_state == STATE_ALN && forbid_STATE_ALN))
				{
					forbid_current_state_comp = true;
				}

				// This loop is for marginalizing over next state.
				for(int next_state = 0; !forbid_current_state_comp && next_state < N_STATES; next_state++)
				{
					// Transition to next state and emission of next pair of symbols.
					if(next_state == STATE_ALN && 
						i <= l1() && 
						k <= l2() &&
						back_array->check_phmm_boundary(i+1, k+1))
					{
						// Next state is emitting starting with i + 1 and k + 1.
						//printf("next_state = %s (prob = %f)\n", state_names[next_state,  back_array->x(i+1, k+1, next_state]);
						//printf("%s -> %s\n", state_names[current_state,  state_names[next_state]);
						//printf("%s(%d)\n", __FILE__, __LINE__);
//                                            cerr << "back aln " << i << " " << k << "\n";
                                            double trans_emit_prob = xlog_mul(this->get_match_prior(i+1, k+1, n1, n2), get_trans_emit_prob(current_state, next_state, i+1, k+1));
/*
						if(i == hmm_array->n_length1)
						{
							printf("Trans-emit probability with i = l1, k = %d(%d) for transition and emission is %.5f\n", k, hmm_array->n_length2, trans_emit_prob);
							cin.get();
						}*/

						//printf("Recursing on (%d, %d) with prob. %f.\n", i+1, k+1, back_array->x(i+1, k+1, next_state]);
						back_array->x(i, k, current_state) = xlog_sum(back_array->x(i, k, current_state),  xlog_mul(back_array->x(i+1, k+1, next_state),  trans_emit_prob));
					}
					
					if(next_state == STATE_INS1 &&
						i <= l1() && 
						back_array->check_phmm_boundary(i+1, k))
					{
						// Next state is emitting starting with i + 1 and k.
						//printf("next_state = %s (prob = %f)\n", state_names[`next_state,  back_array->x(i+1, k, next_state]);
						//printf("%s -> %s\n", state_names[current_state,  state_names[next_state]);
						//printf("%s(%d)\n", __FILE__, __LINE__);
                                            double trans_emit_prob = xlog_mul(/*this->get_indel_prior(i+1, 0, n1, n2)*/ xlog(1), get_trans_emit_prob(current_state, next_state, i+1, k));

						//printf("Recursing on (%d, %d) with prob. %f.\n", i+1, k, back_array->x(i+1, k, next_state]);
						back_array->x(i, k, current_state) = xlog_sum(back_array->x(i, k, current_state),  xlog_mul(back_array->x(i+1, k, next_state),  trans_emit_prob));

						/*if(i == hmm_array->n_length1)
						{
							printf("Trans-emit probability with i = l1, k = %d for transition and emission is %.5f\n", k, trans_emit_prob);
							cin.get();
						}*/
					}
					
					if(next_state == STATE_INS2 &&
						k <= l2() && 
						back_array->check_phmm_boundary(i, k+1))
					{
						// Next state is emitting starting with i and k + 1.
						//printf("next_state = %s (prob = %f)\n\n", state_names[next_state,  back_array->x(i, k+1, next_state]);
						//printf("%s -> %s\n", state_names[current_state,  state_names[next_state]);
						//printf("%s(%d)\n", __FILE__, __LINE__);
                                            double trans_emit_prob = xlog_mul(/*this->get_indel_prior(0, k+1, n1, n2)*/ xlog(1), get_trans_emit_prob(current_state, next_state, i, k+1));

						/*if(i == 0 && k == hmm_array->n_length2)
						{
							printf("Trans-emit probability for transition and emission is %.5f\n", trans_emit_prob);
							cin.get();
						}*/

						//printf("Recursing on (%d, %d) with prob. %f.\n", i, k+1, back_array->x(i+1, k+1, next_state]);
						back_array->x(i, k, current_state) = xlog_sum(back_array->x(i, k, current_state),  xlog_mul(back_array->x(i, k+1, next_state),  trans_emit_prob));
					}	
				} // next_state loop.

				//printf("%d, %d, current state: %s, prob = %f (%f)\n", i, k, state_names[current_state,  back_array->x(i, k, current_state,  xexp(back_array->x(i, k, current_state]) );
			} // current_state loop.
		} // i index loop.
	} // k index loop.
}

// Dump everything necessary from pp_result.
//void t_phmm_aln::dump_data_files(t_pp_result* pp_result)
//{
//	char aln_prob_matrix_fn[MAX_PATH];
//	char ins1_prob_matrix_fn[MAX_PATH];
//	char ins2_prob_matrix_fn[MAX_PATH];
//
//	// Set up file names to dump alignment prioris.
//	sprintf(aln_prob_matrix_fn, "%s/%s", DATA_PATH, PHMM_STATE_ALN_PRIORI_FN);
//	sprintf(ins1_prob_matrix_fn, "%s/%s", DATA_PATH, PHMM_STATE_INS1_PRIORI_FN);
//	sprintf(ins2_prob_matrix_fn, "%s/%s", DATA_PATH, PHMM_STATE_INS2_PRIORI_FN);
//
//	// Open priori files.
//	FILE* state_aln_prior_file = open_f(aln_prob_matrix_fn, "w");
//	FILE* state_ins1_prior_file = open_f(ins1_prob_matrix_fn, "w");
//	FILE* state_ins2_prior_file = open_f(ins2_prob_matrix_fn, "w");
//
//	// Dump pairwise alignment probabilities and insertion probabilities for individual bases.
//	// Need 00 probabilities.
//	for(int i = 0; i < l1() + 1; i++)
//	{
//		for(int k = 0; k < l2() + 1; k++)
//		{
//			// Write current aln, ins1 and ins2 probabilities (scores) to files with at least 25th digit of tenths.
//			fprintf(state_aln_prior_file, "%-6.25f ", pp_result->aln_probs[i][k]);
//			fprintf(state_ins1_prior_file, "%-6.25f ", pp_result->ins1_probs[i][k]);
//			fprintf(state_ins2_prior_file, "%-6.25f ", pp_result->ins2_probs[i][k]);
//		}
//
//		fprintf(state_aln_prior_file, "\n");
//		fprintf(state_ins1_prior_file, "\n");
//		fprintf(state_ins2_prior_file, "\n");
//	}
//
//	fclose(state_aln_prior_file);
//	fclose(state_ins1_prior_file);
//	fclose(state_ins2_prior_file);
//
//	// Dump alignment envelope file. Need the family threshold for this.
//	double fam_threshold = pp_result->fam_threshold; 
//	printf("Family threshold is %.5f\n", fam_threshold);
//	bool** aln_env = (bool**)malloc(sizeof(bool*) * (l1() + 3));
//	for(int i = 1; i < l1() + 1; i++)
//	{
//		aln_env[i] = (bool*)malloc(sizeof(bool) * (l2() + 3));
//
//		for(int k = 1; k < l2() + 1; k++)
//		{
//			double coinc_prob = xlog_sum(xlog_sum(pp_result->ins1_probs[i][k], pp_result->ins2_probs[i][k]), pp_result->aln_probs[i][k]);
//			if(coinc_prob > fam_threshold)
//			{
//				aln_env[i][k] = true;
//			}
//			else
//			{
//				aln_env[i][k] = false;
//			}
//		}
//	}
//
//	// Prune alignment envelope.
//	bool** pruned_aln_env = this->prune_aln_env(aln_env);
//
//	FILE* f_aln_env = open_f("aln_env.txt", "w");
//	for(int i = 1; i < l1() + 1; i++)
//	{
//		for(int k = 1; k < l2() + 1; k++)
//		{
//			fprintf(f_aln_env, "%d ", pruned_aln_env[i][k]);
//		}
//		fprintf(f_aln_env, "\n");
//	}
//	fclose(f_aln_env);
//
//	// Compute and dump ll file.
//	// Dump loop limits.
//	char ll_fn[1000];
//	sprintf(ll_fn, "%s/%s", DATA_PATH, PHMM_LL_FN);
//	FILE* ll_file = open_f(ll_fn, "w");
//
//	// Write loop limits.
//	for(int i = 1; i <= l1(); i++)
//	{
//		fprintf(ll_file, "%d ", i); // Dump seq1 index.
//
//		for(int k = 1; k <= l2(); k++)
//		{
//			if(pruned_aln_env[i][k])
//			{
//				fprintf(ll_file, "%d ", k); // Dump low limit.
//				break;
//			}
//		}
//
//		for(int k = l2(); k >= 1; k--)
//		{
//			if(pruned_aln_env[i][k])
//			{
//				fprintf(ll_file, "%d", k); // Dump high limit.
//				break;
//			}
//		}
//
//		fprintf(ll_file, "\n");
//	}
//
//	// close files.
//	fclose(ll_file);
//}



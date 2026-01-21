#include "phmm_aln.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "structure/structure_object.h"
#include "utils/xmath/log/xlog_math.h"
#include "utils/file/utils.h"
#include "phmm.h"
#include "phmm_array.h"
#include "p_alignment.h"

#include <algorithm>
#include <iostream>

using namespace std;

bool _DUMP_PHMM_ML_LOOPS_MESSAGES_ = false;

// Parameter to following function should be a forward hmm array.
void t_phmm_aln::init_ML_array(t_phmm_array* ml_array)
{
	// Set 0,0 alignment as 1 prob.
	ml_array->x(0, 0, STATE_START) = xlog(1.0);
	ml_array->x(0, 0, STATE_INS1) = xlog(0);
	ml_array->x(0, 0, STATE_INS2) = xlog(0);
}

t_ML_result* t_phmm_aln::compute_ML_alignment()
{
	// Compute ML array and return the alignment in ML results.
	this->phmm = new t_phmm(ML_emit_probs, ML_trans_probs);

	t_ML_result* ml_res = (t_ML_result*)malloc(sizeof(t_ML_result));

	t_phmm_array* ml_array = this->compute_ML_array(ml_res);

	this->traceback_ml_array(ml_array, ml_res);

	delete(ml_array);
	delete(this->phmm); // This is not needed after the computation any more.

	return(ml_res);
}

void* t_phmm_aln::free_ML_result(t_ML_result* ml_result)
{
	ml_result->seq1_aln_line->clear();
	ml_result->seq2_aln_line->clear();

	delete(ml_result->seq1_aln_line);
	delete(ml_result->seq2_aln_line);

	delete(ml_result->ml_alignment);
	free(ml_result);

	return(NULL);
}

/*
aligned_positions is the alignmnt constraint: It is of length seq1 and for i
in seq1 that is constrained to align k in seq2, aligned_positions[i]=k is set, otw
aligned_positions[i]=0 for all other nucleotides.
*/
t_phmm_array* t_phmm_aln::compute_ML_array(t_ML_result* ml_res)
{
	// Replicate the constraints for fast processing of constraints.
	//int* seq2_aln_const = this->get_seq2_aln_const(seq1_aln_const);

	// ML array.
	t_phmm_array* ml_array = new t_phmm_array(l1(), l2(), this->phmm_band_constraint_size, true);

if(_DUMP_PHMM_ML_LOOPS_MESSAGES_)
	printf("Allocated %lf bytes for ML array\n", ml_array->n_bytes_alloced);

	// Init ML array.
	this->init_ML_array(ml_array);
	// Calculate max scores matrix and set max scoring path all along.
	for(int i = 0; i <= l1(); i++)
	{
		int min_k = ml_array->low_phmm_array_limits[i];
		int max_k = ml_array->high_phmm_array_limits[i];

if(_DUMP_PHMM_ML_LOOPS_MESSAGES_)
		printf("ML(%d, x)\r", i);

		for(int k = min_k; k <= max_k; k++)
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
			// Go over cur states.
			for(int cur_state = 0; cur_state < N_STATES; cur_state++)
			{
				// Choose max path through next states using transition and emission of next symbol.				
				double trans_emit_prob = xlog(0);

				double max_score = xlog(0);
				char max_scoring_state = 0;

				for(int prev_state = 0; prev_state < N_STATES; prev_state++)
				{
					if(!forbid_STATE_ALN && 
						cur_state == STATE_ALN && 
						i > 0 && 
						k > 0 &&
						ml_array->check_phmm_boundary(i-1, k-1))
					{
                                            trans_emit_prob = xlog_mul(this->get_match_prior(i,k,l1(),l2()), get_trans_emit_prob(prev_state, cur_state, i, k));

						//printf("trans_emit = %.3f\n", trans_emit_prob);


						// Set max, that is min for probabilities.
						if(xlog_mul(ml_array->x(i - 1, k - 1, prev_state), trans_emit_prob) > max_score)
						{
							//ml_array->x(i, k] = xlog_mul(ml_array->x(i, k), trans_emit_prob);
							max_score = xlog_mul(ml_array->x(i - 1, k - 1, prev_state), trans_emit_prob);

							max_scoring_state = prev_state;
						}
					}

					if(!forbid_STATE_INS1 && 
						cur_state == STATE_INS1 && 
						i > 0 &&
						ml_array->check_phmm_boundary(i-1, k))
					{
                                            trans_emit_prob = xlog_mul(/*this->get_indel_prior(i,0,l1(),l2())*/ xlog(1), get_trans_emit_prob(prev_state, cur_state, i, k));

						//printf("trans_emit = %.3f\n", trans_emit_prob);


						// Set max, that is min for probabilities.
						if(xlog_mul(ml_array->x(i - 1, k, prev_state), trans_emit_prob) > max_score)
						{
							//ml_array->x(i, k] = xlog_mul(ml_array->x(i, k), trans_emit_prob);
							max_score = xlog_mul(ml_array->x(i - 1, k, prev_state), trans_emit_prob);
						}
					}

					if(!forbid_STATE_INS2 &&
						cur_state == STATE_INS2 && 
						k > 0 &&
						ml_array->check_phmm_boundary(i, k-1))
					{
                                            trans_emit_prob = xlog_mul(/*this->get_indel_prior(0,k,l1(),l2())*/ xlog(1), get_trans_emit_prob(prev_state, cur_state, i, k));
//                                            cerr << "lala1\n";
						//printf("trans_emit = %.3f\n", trans_emit_prob);

						// Set max, that is min for probabilities.
						if(xlog_mul(ml_array->x(i, k - 1, prev_state), trans_emit_prob) > max_score)
						{
							//ml_array->x(i, k] = xlog_mul(ml_array->x(i, k), trans_emit_prob);
							max_score = xlog_mul(ml_array->x(i, k - 1, prev_state), trans_emit_prob);
						}
					}
				} // prev_state loop.

				// Copy max score.
				if(i != 0 || k != 0)
				{
					ml_array->x(i, k, cur_state) = max_score;
					//fprintf(f_ml, "%d %d %d %.5f\n", i,k, cur_state, ml_array->x(i, k, cur_state));
					//printf("%d %d %d %.5f\n", i,k, cur_state, ml_array->x(i, k, cur_state));
				}
			} // current_state loop
		} // k loop
	} // i loop

	return(ml_array);
}

void t_phmm_aln::traceback_ml_array(t_phmm_array* ml_array, t_ML_result* ml_res)
{
	// Calculate score for ending sequence and set ML path for ending sequence.
	double max_score = xlog(0);
	char max_scoring_state = 0;
	for(int i = 0; i < N_STATES; i++)
	{
		double log_trans_prob = this->phmm->get_trans_prob(i, STATE_ALN);
	
		if(xlog_mul(ml_array->x(l1(), l2(), i), log_trans_prob) > max_score)
		{
			 max_score = xlog_mul(ml_array->x(l1(), l2(), i), log_trans_prob);
			 max_scoring_state = i;
		}
	}

if(_DUMP_PHMM_ML_LOOPS_MESSAGES_)
	printf("Max score = %.6f\n", max_score);
	//getc(stdin);

	ml_res->ml_prob = max_score;
	ml_res->seq1_aln_line = new vector<char>();
	ml_res->seq2_aln_line = new vector<char>();

	// Trace from end.
	int i = l1();
	int k = l2();
	int cur_state = max_scoring_state;

	/* 
		Trace until all of the nucleotides 
		are assigned with a state in the 
		state sequence.
	*/ 
	while(i != 0 || k != 0)
	{
		//printf("%d, %d -> %s\n", i, k, t_phmm::state_names[cur_state]);

		// Look for the previous state.
		int prev_state = 1234; // This is an invalid value to check for errors.

		if(i > 0 && k > 0 && cur_state == STATE_ALN)
		{
			ml_res->seq1_aln_line->push_back(this->seq1->nucs[i]);
			ml_res->seq2_aln_line->push_back(this->seq2->nucs[k]);

			for(int i_prev_state = 0; i_prev_state < N_STATES; i_prev_state++)
			{
				//printf("%.100f <-> %.100f\n", xlog_mul(ml_array->x(i-1, k-1, i_prev_state), this->get_trans_emit_prob(i_prev_state, cur_state, i, k)), ml_array->x(i, k, cur_state));
				//if(ml_array->x(i, k, cur_state) == xlog_mul(ml_array->x(i-1, k-1, i_prev_state), this->get_trans_emit_prob(i_prev_state, cur_state, i, k)))
                            if(xlog_comp( ml_array->x(i, k, cur_state), xlog_mul(ml_array->x(i-1, k-1, i_prev_state), xlog_mul(this->get_match_prior(i,k,l1(),l2()), this->get_trans_emit_prob(i_prev_state, cur_state, i, k))) ))
				{
					//printf("Found! %d\n", i_prev_state);
					prev_state = i_prev_state;
				}
			} // prev_state loop
		}

		if(i > 0 && cur_state == STATE_INS1)
		{
			ml_res->seq1_aln_line->push_back(this->seq1->nucs[i]);
			ml_res->seq2_aln_line->push_back(GAP_SYM);

			for(int i_prev_state = 0; i_prev_state < N_STATES; i_prev_state++)
			{
				//printf("%.15f <-> %.15f\n", xlog_mul(ml_array->x(i-1, k, i_prev_state), this->get_trans_emit_prob(i_prev_state, cur_state, i, k)), ml_array->x(i, k, cur_state));
				//if(ml_array->x(i, k, cur_state) == xlog_mul(ml_array->x(i-1, k, i_prev_state), this->get_trans_emit_prob(i_prev_state, cur_state, i, k)))
                            if(xlog_comp( ml_array->x(i, k, cur_state), xlog_mul(ml_array->x(i-1, k, i_prev_state), xlog_mul(/*this->get_indel_prior(i,0,l1(),l2())*/ xlog(1), this->get_trans_emit_prob(i_prev_state, cur_state, i, k))) ))
				{
					prev_state = i_prev_state;
				}
			} // prev_state loop
		}

		if(k > 0 && cur_state == STATE_INS2)
		{
			ml_res->seq1_aln_line->push_back(GAP_SYM);
			ml_res->seq2_aln_line->push_back(this->seq2->nucs[k]);

			for(int i_prev_state = 0; i_prev_state < N_STATES; i_prev_state++)
			{
				//printf("%.15f <-> %.15f\n", xlog_mul(ml_array->x(i, k-1, i_prev_state), this->get_trans_emit_prob(i_prev_state, cur_state, i, k)), ml_array->x(i, k, cur_state));
				//if(ml_array->x(i, k, cur_state) == xlog_mul(ml_array->x(i, k-1, i_prev_state), this->get_trans_emit_prob(i_prev_state, cur_state, i, k)))
                            if(xlog_comp( ml_array->x(i, k, cur_state), xlog_mul(ml_array->x(i, k-1, i_prev_state), xlog_mul(/*this->get_indel_prior(0,k,l1(),l2())*/ xlog(1), this->get_trans_emit_prob(i_prev_state, cur_state, i, k))) ))
				{
					prev_state = i_prev_state;
					//printf("Found! %d\n", i_prev_state);
				}
			} // prev_state loop
		}

		// A check if previous state was traced successfully.
		if(prev_state == 1234)
		{
			printf("Could not resolve previous state (%d, %d) (%d, %d)\n", i, k, l1(), l2());
			exit(0);
		}

		// Update the indices based on the current indices. This should be the last
		// thing to do since termination check is right after here.
		if(cur_state == STATE_ALN)
		{
			i--;
			k--;
		}
		else if(cur_state == STATE_INS1)
		{
			i--;
		}
		else if(cur_state == STATE_INS2)
		{
			k--;
		}
		else
		{
			printf("An invalid state in ML traceback.\n");
			exit(0);
		}

		cur_state = prev_state;
	} // Main traceback loop of indices i and k.

	char* seq1_aln_line_str = (char*)malloc(sizeof(char) * (ml_res->seq1_aln_line->size() + 3)); 
	char* seq2_aln_line_str = (char*)malloc(sizeof(char) * (ml_res->seq2_aln_line->size() + 3));

	int l_aln = (int)ml_res->seq1_aln_line->size();

	std::reverse(ml_res->seq1_aln_line->begin(), ml_res->seq1_aln_line->end());
	std::reverse(ml_res->seq2_aln_line->begin(), ml_res->seq2_aln_line->end());

	// Reverse the strings to get alignment lines.
	for(int i = 0; i < (int)ml_res->seq1_aln_line->size(); i++)
	{
		//seq1_aln_line_str[i] = ml_res->seq1_aln_line->at(l_aln - i - 1);
		//seq2_aln_line_str[i] = ml_res->seq2_aln_line->at(l_aln - i - 1);
		seq1_aln_line_str[i] = ml_res->seq1_aln_line->at(i);
		seq2_aln_line_str[i] = ml_res->seq2_aln_line->at(i);
	}

	// Finish the strings.
	seq1_aln_line_str[l_aln] = 0;
	seq2_aln_line_str[l_aln] = 0;

	t_p_alignment* ml_aln = new t_p_alignment(seq1_aln_line_str, seq2_aln_line_str);

if(_DUMP_PHMM_ML_LOOPS_MESSAGES_)
	printf("Similarity = %lf\n", ml_aln->get_aln_similarity(GAP_SYM));

	ml_res->ml_similarity = ml_aln->get_aln_similarity(GAP_SYM);
	ml_res->ml_alignment = ml_aln;

if(_DUMP_PHMM_ML_LOOPS_MESSAGES_)
{
	FILE* f_ml_aln = open_f("ml_alignment.txt", "w");
	fprintf(f_ml_aln, "%s\n%s\n", seq1_aln_line_str, seq2_aln_line_str);
	fclose(f_ml_aln);
}

	//printf("%s\n%s\n", seq1_aln_line_str, seq2_aln_line_str);

	free(seq1_aln_line_str);
	free(seq2_aln_line_str);
}


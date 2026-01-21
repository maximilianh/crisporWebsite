#include "phmm_aln.h"
#include "structure/structure_object.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "utils/xmath/log/xlog_math.h"
#include "utils/file/utils.h"
#include "phmm.h"
#include "phmm_array.h"

bool _DUMP_ALN_ENV_UTILS_MESSAGES_ = false;

t_aln_env_result* t_phmm_aln::compute_alignment_envelope(int aln_env_type, int par)
{
	t_pp_result* pp_result = this->compute_posterior_probs();
	double log_threshold = pp_result->fam_threshold;
	return(compute_alignment_envelope(aln_env_type, pp_result, log_threshold, par));
}

// Backend function for computing alignment envelope.
t_aln_env_result* t_phmm_aln::compute_alignment_envelope(int aln_env_type, 
														 t_pp_result* _pp_result, 
														 double log_threshold, 
														 int par)
{
if(_DUMP_ALN_ENV_UTILS_MESSAGES_)
	printf("Computing alignment envelope...\n");

	// if pp_result is not supplied, recompute it.
	t_pp_result* pp_result = NULL;
	if(_pp_result == NULL)
	{
		pp_result = this->compute_posterior_probs();
	}
	else
	{
		pp_result = _pp_result;
	}

	// alignment envelope type affects how the limits are set.
	// Limit indices are 1 based.
	int* low_limits = (int*)malloc(sizeof(int) * (this->l1() + 2));
	int* high_limits = (int*)malloc(sizeof(int) * (this->l1() + 2));

	// Initialize loop limits.
	for(int i = 0; i <= this->l1(); i++)
	{
		low_limits[i] = 0;
		high_limits[i] = 0;
	}

	if(aln_env_type == PROB_ALN_ENV)
	{
		// Compute alignment envelope.
if(_DUMP_ALN_ENV_UTILS_MESSAGES_)
		printf("Allocating alignment envelope...\n");
		bool** aln_env = (bool**)malloc((this->l1() + 1) * sizeof(bool*));
		double n_aln_env_bytes = 0.0f;

		for(int i = 0; i <= this->l1(); i++)
		{
			int low_k = t_phmm_array::low_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);
			int high_k = t_phmm_array::high_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);
			aln_env[i] = (bool*)malloc((high_k - low_k + 1) * sizeof(bool));
			n_aln_env_bytes += ((high_k - low_k + 1) * sizeof(bool));
			aln_env[i] -= low_k;
		}
if(_DUMP_ALN_ENV_UTILS_MESSAGES_)
		printf("Allocated %lf bytes for alignment envelope.\n", n_aln_env_bytes);

if(_DUMP_ALN_ENV_UTILS_MESSAGES_)
		printf("Computing alignment envelope from probability planes.\n");

		for(int i = 0; i <= this->l1(); i++)
		{
			int low_k = t_phmm_array::low_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);
			int high_k = t_phmm_array::high_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);

			for(int k = low_k; k <= high_k; k++)
			{
				//printf("(%d, %d): %f, %f\n", cnt1, cnt2, xlog_div(global_aln_info.aln_probs[cnt1][cnt2], global_aln_info.op_prob), log_threshold);
				double ins1_prob = pp_result->ins1_probs[i][k];
				double ins2_prob = pp_result->ins2_probs[i][k];
				double aln_prob = pp_result->aln_probs[i][k];
				double three_plane_sum = xlog_sum(ins1_prob, xlog_sum(ins2_prob, aln_prob));

				if(three_plane_sum < log_threshold)
				{
					aln_env[i][k] = false;
				}
				else
				{
					aln_env[i][k] = true;
				}
			}
		}

		//FILE* f_aln_env = fopen("aln_env.txt", "w");
		//for(int i = 0; i <= this->l1(); i++)
		//{
		//	int low_k = t_phmm_array::low_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);
		//	int high_k = t_phmm_array::high_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);

		//	for(int k = low_k; k <= high_k; k++)
		//	{
		//		fprintf(f_aln_env, "%d ", aln_env[i][k]);
		//	} // k loop

		//	fprintf(f_aln_env, "\n");
		//} // i loop
		//fclose(f_aln_env);

if(_DUMP_ALN_ENV_UTILS_MESSAGES_)
		printf("Validating alignment envelope connectivity...\n");

		// If alignment envelope is not connected, return NULL.
		if(!this->check_connection(aln_env))
		{
			printf("Alignment envelope not connected.\n");

			// If pp_result is allocated, free it.
			if(_pp_result == NULL)
			{
				this->free_pp_result(pp_result);
			}

			// Free the limits.
			free(low_limits);
			free(high_limits);

			// Free aln. env. since it is of no use any more.
			for(int i = 0; i <= this->l1(); i++)
			{
				int low_k = t_phmm_array::low_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);
				aln_env[i] += low_k;
				free( aln_env[i] );
			}	

			free(aln_env);	

			return(NULL);
		}

if(_DUMP_ALN_ENV_UTILS_MESSAGES_)
		printf("Pruning alignment envelope...\n");
		// Calculate pruned alignment envelope and set it to global_aln_info's alignment envelope, 
		// calculate also the size of alignment envelope.
		//#define _PRUNE_ALN_
		//#ifdef _PRUNE_ALN_
		bool** pruned_aln_env = this->prune_aln_env(aln_env);
		//#else
		//	copy_aln_env(aln_env);
		//#endif

if(_DUMP_ALN_ENV_UTILS_MESSAGES_)
		printf("Releasing alignment envelope memory.\n");

		// Free aln. env. since it is of no use any more.
		for(int i = 0; i <= this->l1(); i++)
		{
			int low_k = t_phmm_array::low_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);
			aln_env[i] += low_k;
			free( aln_env[i] );
		}	

		free(aln_env);	

if(_DUMP_ALN_ENV_UTILS_MESSAGES_)
		printf("Computing loop limits.\n");
		// Compute the loop limits.
		for(int i = 1; i <= this->l1(); i++)
		{
			int low_k = t_phmm_array::low_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);
			int high_k = t_phmm_array::high_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);

			for(int k = low_k; k <= high_k; k++)
			{
				if(pruned_aln_env[i][k])
				{
					//fprintf(ll_file, "%d ", cnt2); // Dump low limit.
					low_limits[i] = k;
					break;
				}
			}

			for(int k = high_k; k >= low_k; k--)
			{
				if(pruned_aln_env[i][k])
				{
					//fprintf(ll_file, "%d", cnt2); // Dump high limit.
					high_limits[i] = k;
					break;
				}
			}
		} // loop limit computation loop.

		// Free pruned aln. env. since it is of no use any more.
if(_DUMP_ALN_ENV_UTILS_MESSAGES_)
		printf("Releasing pruned alignment envelope memory.\n");

		for(int i = 1; i <= this->l1(); i++)
		{
			int low_k = t_phmm_array::low_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);
			pruned_aln_env[i] += low_k;
			free(pruned_aln_env[i]);
		}	

		free(pruned_aln_env);
	} // PROB_ALN_ENV
	else if(aln_env_type == BANDED_ALN_ENV)
	{
		// par argument contains band size.
		const int band_size = par;
		const int n1 = this->l1();
		const int n2 = this->l2();

		// Initialize loop limits.
		for(int i = 1; i <= n1; i++)
		{
			low_limits[i]  = MIN(0,  int(i*n2/(float)n1 - band_size));
			high_limits[i] = MAX(n2, int(i*n2/(float)n1 + band_size));

if(_DUMP_ALN_ENV_UTILS_MESSAGES_)
			printf("%d -> (%d, %d)\n", i, low_limits[i], high_limits[i]);
		}

		//exit(0);

	} // BANDED_ALN_ENV
	else if(aln_env_type == FULL_ALN_ENV)
	{
		// Initialize loop limits.
		for(double i = 0.0f; i <= this->l1(); i++)
		{
			low_limits[(int)i] = 0;
			high_limits[(int)i] = this->l2();
		}
	} // FULL_ALN_ENV
	else if(aln_env_type == MANUAL_ALN_ENV)
	{
		this->load_map_limits_from_map("aln_map.txt", low_limits, high_limits);
	}
	else
	{
		printf("Invalid alignment envelope type: %d\n", aln_env_type);
		exit(0);
	} // switch according to selected alignment envelope type.

	low_limits[0] = low_limits[1];
	high_limits[0] = high_limits[1];

	// Set low limits with values 1 to 0, so that the initialized values can be recursed correctly.
	for(int i = 0; i <= this->l1(); i++)
	{
		if(low_limits[i] == 1)
		{
			low_limits[i] = 0;
		}
	}

	// Allocate and set aln_env_result.
	t_aln_env_result* aln_env_result = (t_aln_env_result*)malloc(sizeof(t_aln_env_result));
	aln_env_result->high_limits = high_limits;
	aln_env_result->low_limits = low_limits;
	//aln_env_result->pp_result = pp_result;

	// Check for alignment constraints in the alignment envelope.
	this->check_ins1_ins2(aln_env_result);

	// Dump the probability planes. (all of it)
if(_DUMP_ALN_ENV_UTILS_MESSAGES_)
{
	FILE* f_aln_probs = open_f("aln_plane_probs", "wb");
	FILE* f_ins1_probs = open_f("ins1_plane_probs", "wb");
	FILE* f_ins2_probs = open_f("ins2_plane_probs", "wb");
	for(int i1 = 1; i1 <= this->l1(); i1++)
	{
		int low_i2 = t_phmm_array::low_phmm_limit(i1, l1(), l2(), this->phmm_band_constraint_size);
		int high_i2 = t_phmm_array::high_phmm_limit(i1, l1(), l2(), this->phmm_band_constraint_size);

		for(int i2 = low_i2; i2 <= high_i2; i2++)
		{
			if(pp_result->aln_probs[i1][i2] != xlog(0.0))
			{
				double cur_aln_prob = pp_result->aln_probs[i1][i2];
				fwrite(&i1, sizeof(int), 1, f_aln_probs);
				fwrite(&i2, sizeof(int), 1, f_aln_probs);
				fwrite(&cur_aln_prob, sizeof(double), 1, f_aln_probs);
			}

			if(pp_result->ins1_probs[i1][i2] != xlog(0.0))
			{
				double cur_ins1_prob = pp_result->ins1_probs[i1][i2];
				fwrite(&i1, sizeof(int), 1, f_ins1_probs);
				fwrite(&i2, sizeof(int), 1, f_ins1_probs);
				fwrite(&cur_ins1_prob, sizeof(double), 1, f_ins1_probs);
			}

			if(pp_result->ins2_probs[i1][i2] != xlog(0.0))
			{
				double cur_ins2_prob = pp_result->ins2_probs[i1][i2];
				fwrite(&i1, sizeof(int), 1, f_ins2_probs);
				fwrite(&i2, sizeof(int), 1, f_ins2_probs);
				fwrite(&cur_ins2_prob, sizeof(double), 1, f_ins2_probs);
			}
		} // i2 loop.
	} // i1 loop.

	fclose(f_aln_probs);
	fclose(f_ins1_probs);
	fclose(f_ins2_probs);

	FILE* f_lls = open_f("loop_limits.txt", "w");

	// Dump the loop limits.
	for(int i = 0; i <= this->l1(); i++)
	{
		fprintf(f_lls, "%d %d %d\n", i, low_limits[i], high_limits[i]);
	}

	fclose(f_lls);
} // message dump check.

	//printf("Dumping alignment map.\n");
	//FILE* f_aln_map = open_f("aln_map.txt", "w");

	//for(int i = 1; i <= this->l1(); i++)
	//{
	//	for(int j = 1; j <= this->l2(); j++)
	//	{
	//		if(j < low_limits[i])
	//		{
	//			fprintf(f_aln_map, "0");
	//		}
	//		else if(j <= high_limits[i])
	//		{
	//			fprintf(f_aln_map, "1");
	//		}
	//		else
	//		{
	//			fprintf(f_aln_map, "0");
	//		}
	//	}
	//	fprintf(f_aln_map, "\n");
	//}

	//fclose(f_aln_map);

if(_DUMP_ALN_ENV_UTILS_MESSAGES_)
	printf("Computed alignment envelope.\n");

	return(aln_env_result);
}

void* t_phmm_aln::free_aln_env_result(t_aln_env_result* aln_env_result)
{	
	free(aln_env_result->high_limits);
	free(aln_env_result->low_limits);
	//this->free_pp_result(aln_env_result->pp_result);
	free(aln_env_result);

	return(NULL);
}

void t_phmm_aln::check_ins1_ins2(t_aln_env_result* aln_env_result)
{
	for(int i = 1; i < this->seq1->numofbases; i++)
	{
		if(aln_env_result->low_limits[i] > aln_env_result->high_limits[i-1])
		{
			aln_env_result->high_limits[i-1] = aln_env_result->low_limits[i];
			//printf("**ENLARGED ALIGNMENT ENVELOPE TO COVER %d, %d**\n", i-1, aln_env_result->high_limits[i-1]);
		}
	}
}

bool** t_phmm_aln::prune_aln_env(bool** aln_env)
{
if(_DUMP_ALN_ENV_UTILS_MESSAGES_)
	printf("Allocating pruned alignment envelope.\n");
	bool** pruned_aln_env = (bool**)malloc(sizeof(bool*) * (l1() + 3));
	for(int i = 1; i <= l1(); i++)
	{
		int low_k = t_phmm_array::low_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);
		int high_k = t_phmm_array::high_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);

		pruned_aln_env[i] = (bool*)malloc(sizeof(bool) * (high_k - low_k + 1));
		pruned_aln_env[i] -= low_k;

		for(int k = low_k; k <= high_k; k++)
		{
			pruned_aln_env[i][k] = aln_env[i][k];
		}
	}

if(_DUMP_ALN_ENV_UTILS_MESSAGES_)
	printf("Checking backward connections.\n");

	// Find the first 1 in the first row of aln_env or column.
	for(int i = 1; i <= l1(); i++)
	{
		int low_k = t_phmm_array::low_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);
		int high_k = t_phmm_array::high_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);

		for(int k = low_k; k <= high_k; k++)
		{
			if(check_backward_connection(i, k, pruned_aln_env))
			{
				pruned_aln_env[i][k] = 1;
			}
			else
			{
				pruned_aln_env[i][k] = 0;
			}
		}
	}

if(_DUMP_ALN_ENV_UTILS_MESSAGES_)
	printf("Checking forward connections.\n");
	for(int i = l1(); i >= 1; i--)
	{
		int low_k = t_phmm_array::low_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);
		int high_k = t_phmm_array::high_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);

		for(int k = high_k; k >= low_k; k--)
		{
			if(check_forward_connection(i, k, pruned_aln_env))
			{
				pruned_aln_env[i][k] = 1;
			}
			else
			{
				pruned_aln_env[i][k] = 0;
			}
		}
	}

	return(pruned_aln_env);
}


// Checks the connection between two points in aln_env, checks immediate neighbors.
bool t_phmm_aln::check_backward_connection(int i, int k, bool** pruned_aln_env)
{
	if(pruned_aln_env[i][k])
	{
		// If at the first row or first column, just return true.
		if(i == 1 || k == 1)
		{
			return(true);
		}
		else
		{
			if((t_phmm_array::check_phmm_boundary(i-1, k, this->l1(), this->l2(), this->phmm_band_constraint_size) && pruned_aln_env[i-1][k] == 1)
				|| (t_phmm_array::check_phmm_boundary(i, k-1, this->l1(), this->l2(), this->phmm_band_constraint_size) && pruned_aln_env[i][k-1] == 1 )
				|| (t_phmm_array::check_phmm_boundary(i-1, k-1, this->l1(), this->l2(), this->phmm_band_constraint_size) && pruned_aln_env[i-1][k-1] == 1))
			{
				return(true);
			}
			else
			{
				return(false);
			}
		}
	}
	else
	{
		return(false);
	}
}

// Checks the connection between two points in aln_env, checks immediate neighbors.
bool t_phmm_aln::check_forward_connection(int i, int k, bool** pruned_aln_env)
{
	if(pruned_aln_env[i][k])
	{
		// If at the first row or first column, just return true.
		if(i == l1() || k == l2())
		{
			return(true);
		}
		else
		{
			if((t_phmm_array::check_phmm_boundary(i+1, k, this->l1(), this->l2(), this->phmm_band_constraint_size) && pruned_aln_env[i+1][k]) || 
				(t_phmm_array::check_phmm_boundary(i, k+1, this->l1(), this->l2(), this->phmm_band_constraint_size) && pruned_aln_env[i][k+1]) || 
				(t_phmm_array::check_phmm_boundary(i+1, k+1, this->l1(), this->l2(), this->phmm_band_constraint_size) && pruned_aln_env[i+1][k+1]))
			{
				return(true);
			}
			else
			{
				return(false);
			}
		}
	}
	else
	{
		return(false);
	}
}

bool t_phmm_aln::check_connection(t_aln_env_result* aln_env_res)
{
	for(int i = 1; i < l1(); i++)
	{
		if((aln_env_res->high_limits[i] + 1) < aln_env_res->low_limits[i+1])
		{
			return(false);
		}
	}

	return(true);
}

bool t_phmm_aln::check_connection(bool** aln_env)
{
	bool** connectedness = (bool**)malloc(sizeof(bool*) * (l1() + 3));
	for(int i = 0; i <= l1(); i++)
	{
		int low_k = t_phmm_array::low_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);
		int high_k = t_phmm_array::high_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);

		connectedness[i] = (bool*)malloc(sizeof(bool) * (high_k - low_k + 1));
		connectedness[i] -= low_k;

		for(int k = low_k; k <= high_k; k++)
		{
			connectedness[i][k] = false;
		}
	}

	// Set 0,0 as connected, boundary condition.
	connectedness[0][0] = true;

	for(int i = 0; i <= l1(); i++)
	{
		int low_k = t_phmm_array::low_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);
		int high_k = t_phmm_array::high_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);

		for(int k = low_k; k <= high_k; k++)
		{
			if(connectedness[i][k])
			{
				// INS1
				if(i+1 <= l1() && 
					aln_env[i+1][k] &&
					t_phmm_array::check_phmm_boundary(i+1, k, l1(), l2(), this->phmm_band_constraint_size))
				{
					connectedness[i+1][k] = true;
				}

				// INS2
				if(k+1 <= l2() && 
					aln_env[i][k+1] &&
					t_phmm_array::check_phmm_boundary(i, k+1, l1(), l2(), this->phmm_band_constraint_size))
				{
					connectedness[i][k+1] = true;
				}

				// ALN
				if(i+1 <= l1() &&
					k+1 <= l2() &&
					 aln_env[i+1][k+1] &&
					 t_phmm_array::check_phmm_boundary(i+1, k+1, l1(), l2(), this->phmm_band_constraint_size))
				{
					connectedness[i+1][k+1] = true;
				}
			}// connectedness check for current indices.
		} // k loop
	} // i loop

	//FILE* f_aln_env = fopen("connectedness.txt", "w");
	//for(int i = 0; i <= this->l1(); i++)
	//{
	//	int low_k = t_phmm_array::low_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);
	//	int high_k = t_phmm_array::high_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);

	//	for(int k = low_k; k <= high_k; k++)
	//	{
	//		fprintf(f_aln_env, "%d ", connectedness[i][k]);
	//	} // k loop

	//	fprintf(f_aln_env, "\n");
	//} // i loop
	//fclose(f_aln_env);

	bool overall_connectedness = connectedness[l1()][l2()];

	for(int i = 0; i <= l1(); i++)
	{
		int low_k = t_phmm_array::low_phmm_limit(i, l1(), l2(), this->phmm_band_constraint_size);
		connectedness[i] += low_k;
		free(connectedness[i]);
	}
	free(connectedness);

	// If connectedness reached to end of sequences, return true.
	return(overall_connectedness);
}

void t_phmm_aln::load_map_limits_from_map(const char* aln_map_fn, int* low_limits, int* high_limits)
{
if(_DUMP_ALN_ENV_UTILS_MESSAGES_)
	printf("Setting alignment loop limits from map.\n");
	int N1 = this->l1();
	int N2 = this->l2();

	// Open alignment map file.
	FILE* aln_map_file = open_f(aln_map_fn, "r");

	if(aln_map_file == NULL)
	{
		printf("Could not find alignment map file %s @ %s(%d), exiting.\n", aln_map_fn, __FILE__, __LINE__);
		exit(0);
	}

	for(int i1 = 1; i1 <= N1; i1++)
	{
		// Reset limits.
		low_limits[i1] = -1;
		high_limits[i1] = -1;

		for(int i2 = 1; i2 <= N2; i2++)
		{
			// Current flag for current position in alignment map file.
			int cur_flag;

			// Read map file which consists of 1s and 0s for correct positions.
			fscanf(aln_map_file, "%d", &cur_flag); // Read current flag.

if(_DUMP_ALN_ENV_UTILS_MESSAGES_)
			printf("%d ", cur_flag);

			// Set low limit at the point where 1's start.
			if(low_limits[i1] == -1 && cur_flag == 1)
			{
				low_limits[i1] = i2;
			}

			// Set high limit if high limit is not already set and lomw limit is already set.
			if(high_limits[i1] == -1 && low_limits[i1] != -1 && cur_flag == 0)
			{
				high_limits[i1] = i2 - 1;
			}

			// If high limit is not set and loop hit end of alignment line, set high limit to end of 2nd sequence.
			if(high_limits[i1] == -1 && i2 == N2)
			{
				high_limits[i1] = N2;
			}
		}

		printf("\n");
	}

	fclose(aln_map_file);

	//// Have to set limits for 0th nucleotide.
	//low_limits[0] = low_limits[1];
	//high_limits[0] = high_limits[1];

	//low_limits[0] = 1;
	//high_limits[0] = N2;

//	// For i > N1, just add N2 to limits for i < N1.
//	for(int i1 = N1 + 1; i1 <= 2 * N1; i1++)
//	{
//		low_limits[i1] = low_limits[i1 - N1] + N2;
//		high_limits[i1] = high_limits[i1 - N1] + N2;
//
//if(_DUMP_LOOP_LIMIT_MESSAGES_)
//{
//		printf("low[%d] = %d, high[%d]=%d\n", i1, low_limits[i1], i1, high_limits[i1]);
//}
//	}
//
//	// Set low limits with values 1 to 0, so that the initialized values can be recursed correctly.
//	for(int i = 0; i <= N1; i++)
//	{
//		if(low_limits[i] == 1)
//		{
//			low_limits[i] = 0;
//		}

////if(_DUMP_LOOP_LIMIT_MESSAGES_)
//{
//		//printf("low[%d] = %d, high[%d]=%d\n", i, low_limits[i], i, high_limits[i]);
//}
	/*}*/
}



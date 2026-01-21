#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "parts_compilation_directives.h"

#include "ppf_w.h"
#include "ppf_w_mb.h"

#include "process_sequences.h"
#include "ppf_math.h"
#include "single_pf_array.h"
#include "phmm_parameters.h"
#include "alignment_priors.h"

#include "template_pf_array.h"
#include "ppf_progress_bar.h"
#include "ppf_cli.h"
#include "../../../src/phmm/structure/folding_constraints.h"
#include "../../../src/phmm/utils/file/utils.h"

#include <iostream>

#include <limits.h>

using namespace std;

#include <algorithm>
#define MIN(x,y) min(x,y)
#define MAX(x,y) max(x,y)

int t_template_pf_array::n_arrays = 0;

// #include "mem_tracker/memory_tracker.h" // Override standard memory functions.

bool _DUMP_TEMPLATE_PF_ARRAY_MESSAGES_ = false;

int* t_template_pf_array::high_limits = NULL;
int* t_template_pf_array::low_limits = NULL;

t_template_pf_array::t_template_pf_array(t_seq_man* seq_man, bool mallocate)
{
	// Copy sequence lengths.
	this->N1 = seq_man->get_l_seq1();
	this->N2 = seq_man->get_l_seq2();

	// Init fold envelopes, no aln envelopes with this constructor.
	this->seq1_fold_env = NULL;
	this->seq1_str_coinc_env = NULL;
	this->seq2_fold_env = NULL;
	this->seq2_str_coinc_env = NULL;

	// Now allocate and initialize pf array.
	this->pf_array = NULL;
	alloc_init_pf_array(seq_man, mallocate);

	// This is needed for returning a reference in accession functions.
	this->zero = ZERO;

	if(mallocate)
	{
		t_template_pf_array::n_arrays++;
	}
}

// This is another constructor with fold envelope,
// Note that the fold envelope is just one sequence 1 which is used
// for saving memory for array points which are certainly going to be 0 beforehand by looking
// at pairing priors. 
// Note that seq1_fold_env is used both for allocation and calculations however seq2_fold_env is
// only used for constraining calculations.
t_template_pf_array::t_template_pf_array(t_seq_man* _seq_man, 
										short** _seq1_ptr_rel_map,
										short** _seq2_ptr_rel_map,
										bool mallocate)
{
	this->seq_man = _seq_man;

	// Copy sequence lengths.
	this->N1 = seq_man->get_l_seq1();
	this->N2 = seq_man->get_l_seq2();

	// Copy fold envelope's, note that fold envelopes are coming from spf arrays and they should not be 
	// manipulated.
	//this->seq1_fold_env = _seq1_fold_env;
	//this->seq1_str_coinc_env = _seq1_str_coinc_env;
	//this->seq2_fold_env = _seq2_fold_env;
	//this->seq2_str_coinc_env = _seq2_str_coinc_env;

	this->seq1_ptr_rel_map = _seq1_ptr_rel_map;
	this->seq2_ptr_rel_map = _seq2_ptr_rel_map;

	this->pf_array = NULL;

	// Now allocate and initialize pf array.
	alloc_init_pf_array(seq_man, mallocate);

	// Init 0! Can be avoided by correct initialization of pf array.
	this->zero = ZERO;

	if(mallocate)
	{
		t_template_pf_array::n_arrays++;
	}
}

t_template_pf_array::~t_template_pf_array()
{
if(_DUMP_TEMPLATE_PF_ARRAY_MESSAGES_)
	printf("Destruct'ing a t_template_pf_array object.\n");

	if(this->pf_array != NULL)
	{
		// Allocate pf array.
		//if(mallocate)
		//{
			//this->pf_array = (double****)malloc(sizeof(double***) * (seq_man->get_l_seq1() + 3));
		//}

		//this->n_alloced_bytes += sizeof(double***) * (seq_man->get_l_seq1() + 3);

		int high_j = seq_man->get_l_seq1();

		for(int i = 1; i <= this->N1; i++)
		{
			//printf("i = %d\n", i);
			int min_j = 0x1fffffff;
			int max_j = 0;
			int n_j_allocs = 0;

			for(int j = i; j <= MIN(i+this->seq_man->ppf_cli->max_n_separation_between_nucs, this->N1); j++)
			{
				if(this->seq1_ptr_rel_map[i][j] != POS_MEM_NOT_EXIST)
				{
					if(j < min_j)
					{
						min_j = j;
					}

					if(j > max_j)
					{
						max_j = j;
					}

					n_j_allocs++;
				}
			}

			//if(mallocate && n_j_allocs > 0)
			//{
			//	this->pf_array[i] = (double***)malloc( sizeof(double**) * (n_j_allocs+1) );
			//	this->pf_array[i] -= this->seq1_ptr_rel_map[i][min_j];
			//}
			
			//if(n_j_allocs > 0)
			//{
			//	this->n_alloced_bytes += ( sizeof(double**) * (n_j_allocs+1) );
			//}

			for(int j = min_j; j <= max_j; j++)
			{
				// Is this allocatable??
				if( (this->seq1_ptr_rel_map == NULL) || 
					((this->seq1_ptr_rel_map != NULL) && seq1_ptr_rel_map[i][j] != POS_MEM_NOT_EXIST) )
				{
					//printf("j = %d\n", j);
					// Allocate kxl plane using loop limits for v.
					int max_k = t_template_pf_array::high_limits[i-1]+1;
					int min_k = t_template_pf_array::low_limits[i-1]+1;

					//if(mallocate && (max_k - min_k + 1) > 0)
					//{
					//	this->pf_array[i][this->seq1_ptr_rel_map[i][j]] = (double**)malloc( sizeof(double*) * (max_k - min_k + 1) );
					//	this->pf_array[i][this->seq1_ptr_rel_map[i][j]] -= min_k;
					//}
					
					//if((max_k - min_k + 1) > 0)
					//{
					//	this->n_alloced_bytes += ( sizeof(double**) * (max_k - min_k + 1) );
					//}

					for(int k = min_k; k <= max_k; k++)
					{
						//printf("k = %d\n", k);

						// Also use the 
						int min_l = 0x1fffffff;
						int max_l = 0;
						int n_l_allocs = 0;

						for(int l = k; l <= MIN(k+this->seq_man->ppf_cli->max_n_separation_between_nucs, this->N2); l++)
						{
							if(this->seq2_ptr_rel_map[k][l] != POS_MEM_NOT_EXIST)
							{
								if(l >= t_template_pf_array::low_limits[j] && l <= t_template_pf_array::high_limits[j])
								{
									if(min_l > l)
									{
										min_l = l;
									}

									if(max_l < l)
									{
										max_l = l;
									}

									n_l_allocs++;
								}

								//printf("l(%d): min_l: %d, max_l: %d\n", l, min_l, max_l);
							}
						}

						//if(mallocate && n_l_allocs > 0)
						//{	
						//	this->pf_array[i][this->seq1_ptr_rel_map[i][j]][k] = (double*)malloc( sizeof(double) * (n_l_allocs + 1) );
						//	this->pf_array[i][this->seq1_ptr_rel_map[i][j]][k] -= seq2_ptr_rel_map[k][min_l];
						//}
						
						//if(n_l_allocs > 0)
						//{
						//	this->n_alloced_bytes += (sizeof(double) * (n_l_allocs + 1));
						//}

						//for(int l = min_l; l <= max_l; l++)
						//{
						//	//printf("Initing %d, %d, %d, %d\n", i,j,k,l);
						//	if( (this->seq2_ptr_rel_map == NULL) || 
						//		((this->seq2_ptr_rel_map != NULL) && seq2_ptr_rel_map[k][l] != POS_MEM_NOT_EXIST) )
						//	{
						//		//printf("Initing %d, %d, %d, %d\n", i,j,k,l);
						//		if(mallocate && n_l_allocs > 0)
						//		{
						//			this->pf_array[i][this->seq1_ptr_rel_map[i][j]][k][this->seq2_ptr_rel_map[k][l]] = CONVERT_FROM_LIN(0.0f);
						//		}
						//	} // seq2_str_coinc check.
						//} // l loop

						if(n_l_allocs > 0)
						{	
							this->pf_array[i][this->seq1_ptr_rel_map[i][j]][k] += seq2_ptr_rel_map[k][min_l];
							free(this->pf_array[i][this->seq1_ptr_rel_map[i][j]][k]);
						}

					} // k loop

					if((max_k - min_k + 1) > 0)
					{
						this->pf_array[i][this->seq1_ptr_rel_map[i][j]] += min_k;
						free(this->pf_array[i][this->seq1_ptr_rel_map[i][j]]);
					}

				} // seq1_str_coinc check.
			} // j loop

			if(n_j_allocs > 0)
			{
				this->pf_array[i] += this->seq1_ptr_rel_map[i][min_j];
				free(this->pf_array[i]);
			}
		} // i loop.

		free(this->pf_array);
	} // NULL check for the array memory.
}

// Use seq1_fold_env if present.
void t_template_pf_array::alloc_init_pf_array(t_seq_man* seq_man, bool mallocate)
{
if(_DUMP_TEMPLATE_PF_ARRAY_MESSAGES_)
{
	printf("Allocating and initing a ppf array.\n");
}

	this->n_alloced_bytes = 0.0f;

	// Allocate pf array.
	if(mallocate)
	{
		this->pf_array = (double****)malloc(sizeof(double***) * (seq_man->get_l_seq1() + 3));
	}

	this->n_alloced_bytes += sizeof(double***) * (seq_man->get_l_seq1() + 3);

	int high_j = seq_man->get_l_seq1();

	t_ppf_progress_bar* tpa_alloc_progress_bar = new t_ppf_progress_bar(NULL, '=', true, high_j);

	for(int i = 1; i <= this->N1; i++)
	{
		//printf("i = %d\n", i);
		int min_j = 0x1fffffff;
		int max_j = 0;
		int n_j_allocs = 0;

		for(int j = i; j <= MIN(i+this->seq_man->ppf_cli->max_n_separation_between_nucs, this->N1); j++)
		{
			if(this->seq1_ptr_rel_map[i][j] != POS_MEM_NOT_EXIST)
			{
				if(j < min_j)
				{
					min_j = j;
				}

				if(j > max_j)
				{
					max_j = j;
				}

				n_j_allocs++;
			}
		}

		if(mallocate && n_j_allocs > 0)
		{
			this->pf_array[i] = (double***)malloc( sizeof(double**) * (n_j_allocs+1) );
			this->pf_array[i] -= this->seq1_ptr_rel_map[i][min_j];
		}
		
		if(n_j_allocs > 0)
		{
			this->n_alloced_bytes += ( sizeof(double**) * (n_j_allocs+1) );
		}

		for(int j = min_j; j <= max_j; j++)
		{
			// Is this allocatable??
			if( (this->seq1_ptr_rel_map == NULL) || 
				((this->seq1_ptr_rel_map != NULL) && seq1_ptr_rel_map[i][j] != POS_MEM_NOT_EXIST) )
			{
				//printf("j = %d\n", j);
				// Allocate kxl plane using loop limits for v.
				int max_k = t_template_pf_array::high_limits[i-1]+1;
				int min_k = t_template_pf_array::low_limits[i-1]+1;

				if(mallocate && (max_k - min_k + 1) > 0)
				{
					this->pf_array[i][this->seq1_ptr_rel_map[i][j]] = (double**)malloc( sizeof(double*) * (max_k - min_k + 1) );
					this->pf_array[i][this->seq1_ptr_rel_map[i][j]] -= min_k;
				}
				
				if((max_k - min_k + 1) > 0)
				{
					this->n_alloced_bytes += ( sizeof(double**) * (max_k - min_k + 1) );
				}

				for(int k = min_k; k <= max_k; k++)
				{
					//printf("k = %d\n", k);

					// Also use the 
					int min_l = 0x1fffffff;
					int max_l = 0;
					int n_l_allocs = 0;

					for(int l = k; l <= MIN(k+this->seq_man->ppf_cli->max_n_separation_between_nucs, this->N2); l++)
					{
						if(this->seq2_ptr_rel_map[k][l] != POS_MEM_NOT_EXIST)
						{
							if(l >= t_template_pf_array::low_limits[j] && l <= t_template_pf_array::high_limits[j])
							{
								if(min_l > l)
								{
									min_l = l;
								}

								if(max_l < l)
								{
									max_l = l;
								}

								n_l_allocs++;
							}

							//printf("l(%d): min_l: %d, max_l: %d\n", l, min_l, max_l);
						}
					}

					if(mallocate && n_l_allocs > 0)
					{	
						this->pf_array[i][this->seq1_ptr_rel_map[i][j]][k] = (double*)malloc( sizeof(double) * (n_l_allocs + 1) );
						this->pf_array[i][this->seq1_ptr_rel_map[i][j]][k] -= seq2_ptr_rel_map[k][min_l];
					}
					
					if(n_l_allocs > 0)
					{
						this->n_alloced_bytes += (sizeof(double) * (n_l_allocs + 1));
					}

					for(int l = min_l; l <= max_l; l++)
					{
						//printf("Initing %d, %d, %d, %d\n", i,j,k,l);
						if( (this->seq2_ptr_rel_map == NULL) || 
							((this->seq2_ptr_rel_map != NULL) && seq2_ptr_rel_map[k][l] != POS_MEM_NOT_EXIST) )
						{
							//printf("Initing %d, %d, %d, %d\n", i,j,k,l);
							if(mallocate && n_l_allocs > 0)
							{
								this->pf_array[i][this->seq1_ptr_rel_map[i][j]][k][this->seq2_ptr_rel_map[k][l]] = CONVERT_FROM_LIN(0.0f);
							}
						} // seq2_str_coinc check.
					} // l loop
				} // k loop
			} // seq1_str_coinc check.
		} // j loop
	} // i loop.

	delete(tpa_alloc_progress_bar);

if(_DUMP_TEMPLATE_PF_ARRAY_MESSAGES_)
	printf("template_pf_array %lf bytes\n", this->n_alloced_bytes);
}

void t_template_pf_array::alloc_init_loop_limits(t_seq_man* seq_man, int* _low_limits, int* _high_limits)
{	
if(_DUMP_TEMPLATE_PF_ARRAY_MESSAGES_)
{
	printf("Allocating and initing template pf array loop limits..\n");
}

	// Allocate loop limits for i and j from 1 to 2 x N1.
	if(t_template_pf_array::low_limits == NULL)
	{
		t_template_pf_array::low_limits = (int*)malloc( sizeof(int) * (seq_man->get_l_seq1() * 2 + 2) );
		t_template_pf_array::high_limits = (int*)malloc( sizeof(int) * (seq_man->get_l_seq1() * 2 + 2) );

		for(int i = 0; i <= seq_man->get_l_seq1(); i++)
		{
			t_template_pf_array::low_limits[i] = _low_limits[i];
			t_template_pf_array::high_limits[i] = _high_limits[i];
		}
	}
	else
	{
		return;
	}
}

void t_template_pf_array::update_loop_limits(t_seq_man* seq_man, char* lls_fp)
{
	FILE* f_lls = open_f(lls_fp, "r");
	if(f_lls == NULL)
	{
		printf("Could not read loop limits file @ %s.\n", lls_fp);
		exit(0);
	}

	// Allocate loop limits for i and j from 1 to 2 x N1.
	if(t_template_pf_array::low_limits == NULL)
	{
		t_template_pf_array::low_limits = (int*)malloc( sizeof(int) * (seq_man->get_l_seq1() * 2 + 2) );
		t_template_pf_array::high_limits = (int*)malloc( sizeof(int) * (seq_man->get_l_seq1() * 2 + 2) );

		for(int i = 0; i <= seq_man->get_l_seq1(); i++)
		{
			int cur_i = 0;
			int low_i = 0;
			int high_i = 0;
			fscanf(f_lls, "%d %d %d", &cur_i, &low_i, &high_i);

			t_template_pf_array::low_limits[cur_i] = low_i;
			t_template_pf_array::high_limits[cur_i] = high_i;
		}
	}
	else
	{
		free(t_template_pf_array::low_limits);
		free(t_template_pf_array::high_limits);

		t_template_pf_array::low_limits = (int*)malloc( sizeof(int) * (seq_man->get_l_seq1() * 2 + 2) );
		t_template_pf_array::high_limits = (int*)malloc( sizeof(int) * (seq_man->get_l_seq1() * 2 + 2) );

		for(int i = 0; i <= seq_man->get_l_seq1(); i++)
		{
			int cur_i = 0;
			int low_i = 0;
			int high_i = 0;
			fscanf(f_lls, "%d %d %d", &cur_i, &low_i, &high_i);

			t_template_pf_array::low_limits[cur_i] = low_i;
			t_template_pf_array::high_limits[cur_i] = high_i;
		}
	}

	fclose(f_lls);
}

void t_template_pf_array::update_loop_limits(t_seq_man* seq_man, int* _low_limits, int* _high_limits)
{	
if(_DUMP_TEMPLATE_PF_ARRAY_MESSAGES_)
{
	printf("Allocating and initing template pf array loop limits..\n");

	printf("Updating template pf array loop limits.\n");
}

	// Allocate loop limits for i and j from 1 to 2 x N1.
	if(t_template_pf_array::low_limits == NULL)
	{
		t_template_pf_array::low_limits = (int*)malloc( sizeof(int) * (seq_man->get_l_seq1() * 2 + 2) );
		t_template_pf_array::high_limits = (int*)malloc( sizeof(int) * (seq_man->get_l_seq1() * 2 + 2) );

		for(int i = 0; i <= seq_man->get_l_seq1(); i++)
		{
			t_template_pf_array::low_limits[i] = _low_limits[i];
			t_template_pf_array::high_limits[i] = _high_limits[i];
		}
	}
	else
	{
		free(t_template_pf_array::low_limits);
		free(t_template_pf_array::high_limits);

		t_template_pf_array::low_limits = (int*)malloc( sizeof(int) * (seq_man->get_l_seq1() * 2 + 2) );
		t_template_pf_array::high_limits = (int*)malloc( sizeof(int) * (seq_man->get_l_seq1() * 2 + 2) );

		for(int i = 0; i <= seq_man->get_l_seq1(); i++)
		{
			t_template_pf_array::low_limits[i] = _low_limits[i];
			t_template_pf_array::high_limits[i] = _high_limits[i];
		}
	}
}

bool t_template_pf_array::check_boundary(int i1, int i2)
{
	if(t_template_pf_array::low_limits[i1] <= i2 && t_template_pf_array::high_limits[i1] >= i2)
	{
		return(true);
	}
	else
	{
		//printf("Check boundary invalid: i1=%d -> i2=%d\n", i1, i2);
		return(false);
	}
}

bool t_template_pf_array::check_4D_ll(int i, int j, int k, int l)
{
	bool ll_check = (this->check_boundary(i-1, k-1)) && (this->check_boundary(j,l));
	bool str_check = (j > i && (j-i) <= this->seq_man->ppf_cli->max_n_separation_between_nucs) && 
						(l > k && (l-k) <= this->seq_man->ppf_cli->max_n_separation_between_nucs) &&
						(this->seq1_ptr_rel_map[i][j] != POS_MEM_NOT_EXIST && this->seq2_ptr_rel_map[k][l] != POS_MEM_NOT_EXIST);

	return(ll_check && str_check);
}

// Access to partition function array by reference.
double& t_template_pf_array::x(int i, int j, int k, int l)
{
	if(this->zero != ZERO)
	{
		printf("ZERO changed for %d, %d, %d, %d\n", i,j,k,l);
		int* p = NULL;
		*p = 0;
	}

	if(this->seq1_ptr_rel_map[i][j] == POS_MEM_NOT_EXIST || this->seq2_ptr_rel_map[k][l] == POS_MEM_NOT_EXIST)
	{
		printf("-1!\n");
		int* p = NULL;
		*p = 0;
		getc(stdin);
	}

	return(this->pf_array[i][this->seq1_ptr_rel_map[i][j]][k][this->seq2_ptr_rel_map[k][l]]);
}


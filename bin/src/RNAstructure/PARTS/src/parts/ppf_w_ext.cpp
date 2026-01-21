#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "parts_compilation_directives.h"
#include "ppf_w_l.h"
#include "ppf_cli.h"
#include "process_sequences.h"
#include "ppf_math.h"
#include <iostream>
#include "alignment_priors.h"
#include "phmm_parameters.h"

#include "single_pf_array.h"
#include "template_pf_array.h"

#include "ppf_v_mhe.h"
#include "ppf_w_ext.h"
#include "ppf_loops.h"
#include "ppf_scale.h"

#include <math.h>

#include "ppf_tb_stack.h"
#include "stoch_tb/stoch_sampled_str_aln_sample_set.h"
#include "stoch_tb/stoch_sampling_math.h"
#include "stoch_tb/stoch_sampled_alignment.h"
#include "stoch_tb/stoch_sampled_structures.h"

#include "ppf_tb_stack.h"
#include "map_structures.h"
#include "map_alignment.h"
#include "stoch_tb/stoch_sampling_math.h"
#include "ppf_cli.h"

#include "../../../src/phmm/structure/folding_constraints.h"

using namespace std;

#include <algorithm>
#define MIN(x,y) min(x,y)
#define MAX(x,y) max(x,y)

bool _DUMP_PPF_W_EXT_MESSAGES_ = false;

t_ppf_WEXT::t_ppf_WEXT(t_ppf_loops* _ppf_loops)
{
if(_DUMP_PPF_W_EXT_MESSAGES_)
	printf(" W_EXT ");

	this->ppf_loops = _ppf_loops;

	this->seq_man = this->ppf_loops->seq_man;

	// Copy sequence lengths.
	this->N1 = seq_man->get_l_seq1();
	this->N2 = seq_man->get_l_seq2();

	this->seq1_spf = this->ppf_loops->seq1_spf;
	this->seq2_spf = this->ppf_loops->seq2_spf;
	this->aln_priors = this->ppf_loops->aln_priors;

	// Initialize to NULL.
	this->pf_array = NULL;
	this->ext_pf_array = NULL;
	
	this->alloc_init_w_ext_array();	
}

t_ppf_WEXT::~t_ppf_WEXT()
{
	//this->pf_array = (double**)malloc(sizeof(double*) * (this->N1 + 1));

	if(this->pf_array != NULL)
	{
		for(int i = 0; i <= this->N1; i++)
		{
			int min_k = MAX(0, t_template_pf_array::low_limits[i]);
			int max_k = MIN(this->N2, t_template_pf_array::high_limits[i]);
			
			//this->pf_array[i] = (double*)malloc(sizeof(double) * (max_k - min_k + 1));
			//this->pf_array[i] -= min_k;
			//for(int k = min_k; k <= max_k; k++)
			//{
			//	this->pf_array[i][k] = CONVERT_FROM_LIN(0.0f);
			//} // k loop

			this->pf_array[i] += min_k;
			free(this->pf_array[i]);
		} // i loop

		free(this->pf_array);
	}

	if(this->ext_pf_array != NULL)
	{
		//this->ext_pf_array = (double**)malloc(sizeof(double*) * (this->N1 + 3));
		for(int i = 1; i <= this->N1+1; i++)
		{
			int min_k = MAX(1, t_template_pf_array::low_limits[i-1]+1);
			int max_k = MIN(this->N2+1, t_template_pf_array::high_limits[i-1]+1);
			
			this->ext_pf_array[i] += min_k;
			free(this->ext_pf_array[i]);		
		} // i loop

		free(this->ext_pf_array);
	}
}

void t_ppf_WEXT::alloc_init_w_ext_array()
{
	this->pf_array = (double**)malloc(sizeof(double*) * (this->N1 + 3));

	for(int i = 0; i <= this->N1; i++)
	{
		int min_k = MAX(0, t_template_pf_array::low_limits[i]);
		int max_k = MIN(this->N2, t_template_pf_array::high_limits[i]);
		
		this->pf_array[i] = (double*)malloc(sizeof(double) * (max_k - min_k + 1));
		this->pf_array[i] -= min_k;
		for(int k = min_k; k <= max_k; k++)
		{
			this->pf_array[i][k] = CONVERT_FROM_LIN(0.0f);
		} // k loop
	} // i loop

	// Initialize.
	for(int i = 0; i <= this->N1; i++)
	{
		int min_k = MAX(0, t_template_pf_array::low_limits[i]);
		int max_k = MIN(this->N2, t_template_pf_array::high_limits[i]);
		
		for(int k = min_k; k <= max_k; k++)
		{
			if(i == 0 && k != 0 && this->check_boundary(i, k-1) && k > 0)
			{
				// k is inserted.
				double k_ins_unp_score = MUL(this->aln_priors->x(i,k,STATE_INS2), seq2_spf->ux_3p(1,k));
				this->pf_array[i][k] = MUL(this->pf_array[i][k-1], k_ins_unp_score);
			}
			else if(i != 0 && k == 0 && this->check_boundary(i-1, k) && i > 0)
			{
				// i is inserted.
				double i_ins_unp_score = MUL(this->aln_priors->x(i,k,STATE_INS1), seq1_spf->ux_3p(1,i));
				this->pf_array[i][k] = MUL(this->pf_array[i-1][k], i_ins_unp_score);
			}
			else if(k == 0 && i == 0)
			{
				this->pf_array[i][k] = CONVERT_FROM_LIN(1.0f);
			}
			else // Reset all others to 0.
			{
				this->pf_array[i][k] = CONVERT_FROM_LIN(0.0f);
			}
		} // k loop
	} // i loop
} 

void t_ppf_WEXT::alloc_init_ext_w_ext_array()
{
	//_DUMP_PPF_W_EXT_MESSAGES_ = true;
if(_DUMP_PPF_W_EXT_MESSAGES_)
	printf("Initing W_ext[ext]\n");

	this->ext_pf_array = (double**)malloc(sizeof(double*) * (this->N1 + 3));

	for(int i = 1; i <= this->N1+1; i++)
	{
		int min_k = MAX(1, t_template_pf_array::low_limits[i-1]+1);
		int max_k = MIN(this->N2+1, t_template_pf_array::high_limits[i-1]+1);
		
		this->ext_pf_array[i] = (double*)malloc(sizeof(double) * (max_k - min_k + 1));
		this->ext_pf_array[i] -= min_k;
		for(int k = max_k; k >= min_k; k--)
		{
			this->ext_pf_array[i][k] = CONVERT_FROM_LIN(0.0f);
		} // k loop
	} // i loop

	// Initialize.
	for(int i = this->N1+1; i >= 1; i--)
	{
		int min_k = MAX(1, t_template_pf_array::low_limits[i-1]+1);
		int max_k = MIN(this->N2+1, t_template_pf_array::high_limits[i-1]+1);
		
		for(int k = max_k; k >= min_k; k--)
		{
			if(i == N1+1 && k != N2+1 && this->check_boundary(i-1, k) && k > 0)
			{
				// k is inserted.
				double k_ins_unp_score = MUL(this->aln_priors->x(i-1,k,STATE_INS2), seq2_spf->ux_5p(k,N2));
				this->ext_pf_array[i][k] = MUL(this->ext_pf_array[i][k+1], k_ins_unp_score);

if(_DUMP_PPF_W_EXT_MESSAGES_)
				printf("%s(%d): Wext[ext][Init](%d, %d) = %.10f x %.10f x %.10f\n", __FILE__, __LINE__, i,k, 
					this->aln_priors->x(i-1,k,STATE_INS2), this->ext_pf_array[i][k+1], seq2_spf->ux_5p(k,N2));
			}
			else if(i != N1+1 && k == N2+1 && this->check_boundary(i, k-1) && i > 0)
			{
				// i is inserted.
				double i_ins_unp_score = MUL(this->aln_priors->x(i,k-1,STATE_INS1), seq1_spf->ux_5p(i,N1));
				this->ext_pf_array[i][k] = MUL(this->ext_pf_array[i+1][k], i_ins_unp_score);

if(_DUMP_PPF_W_EXT_MESSAGES_)
				printf("%s(%d): Wext[ext][Init](%d, %d) = %.10f x %.10f x %.10f\n", __FILE__, __LINE__, i,k, 
					this->aln_priors->x(i,k-1,STATE_INS1), this->ext_pf_array[i+1][k], seq1_spf->ux_5p(i,N1));
			}
			else if(i == N1+1 && k == N2+1)
			{
				this->ext_pf_array[i][k] = CONVERT_FROM_LIN(1.0f);
			}
			else // Reset all others to 0.
			{
				this->ext_pf_array[i][k] = CONVERT_FROM_LIN(0.0f);
			}

			if(_DUMP_PPF_W_EXT_MESSAGES_ &&
				this->x_ext(i,k) != ZERO)
			{
				printf("W_ext[ext][init](%d, %d) = %.5f\n", i, k, this->x_ext(i,k));
			}
		} // k loop
	} // i loop

	//getc(stdin);
}

void t_ppf_WEXT::calculate_ext_W_ext()
{
	this->alloc_init_ext_w_ext_array();

	// Initialize.
	for(int i = this->N1; i >= 1; i--)
	{
		int min_k = MAX(1, t_template_pf_array::low_limits[i-1]+1);
		int max_k = MIN(this->N2, t_template_pf_array::high_limits[i-1]+1);
		
		for(int k = max_k; k >= min_k; k--)
		{
#ifdef _LINEAR_COMPUTATIONS_
			do 
			{
#endif

				// Initialize to 0.
				this->x_ext(i, k) = ZERO;

				// insert i
				if(k-1 >= 0 &&
					this->check_boundary(i, k-1))
				{
					double i_ins_unp_score = MUL(this->aln_priors->x(i, k-1, STATE_INS1), seq1_spf->ux_3p(1, i));
					this->x_ext(i, k) = this->ppf_loops->MAX_SUM(this->x_ext(i, k), MUL(this->x_ext(i+1, k), i_ins_unp_score));
				}

				// insert k.
				if(i-1 >= 0 &&
					this->check_boundary(i-1, k))
				{
					double k_ins_unp_score = MUL(this->aln_priors->x(i-1, k, STATE_INS2), seq2_spf->ux_3p(1, k));
					this->x_ext(i, k) = this->ppf_loops->MAX_SUM(this->x_ext(i, k), MUL(this->x_ext(i, k+1), k_ins_unp_score));
				}

				// align i and k.
				if(i >= 0 && k >= 0 &&
					this->check_boundary(i, k))
				{
					double ik_aln_unp_score = MUL3(this->aln_priors->x(i, k, STATE_ALN), seq1_spf->ux_3p(i-1, i), seq2_spf->ux_3p(k-1, k));
					this->x_ext(i, k) = this->ppf_loops->MAX_SUM(this->x_ext(i, k), MUL(this->x_ext(i+1, k+1), ik_aln_unp_score));
				}

				// add domain.
				for(int ip = i+1; ip <= MIN(N1, i+this->seq_man->ppf_cli->max_n_separation_between_nucs); ip++)
				{
					int min_kp = MAX(k+1, t_template_pf_array::low_limits[ip]);
					int max_kp = MIN(t_template_pf_array::high_limits[ip], MIN(N2, k+this->seq_man->ppf_cli->max_n_separation_between_nucs));

					for(int kp = min_kp; kp <= max_kp; kp++)
					{		
						double current_concatenation = MUL(this->ppf_loops->V_mhe->x(i, ip, k, kp), this->x_ext(ip+1, kp+1));
						this->x_ext(i,k) = this->ppf_loops->MAX_SUM(this->x_ext(i,k), current_concatenation);
						//if(_DUMP_PPF_W_EXT_MESSAGES_ &&  
						//	this->x_ext(ip, kp) != ZERO && 
						//	this->ppf_loops->V_mhe->x(ip+1, i, kp+1, k) != ZERO)
						//{
						//	printf("Concating Wext(%d, %d)-(Vmhe)(%d, %d, %d, %d): %lf x %lf\n", ip, kp, ip+1, i, kp+1, k, this->x(ip, kp), this->ppf_loops->V_mhe->x(ip+1, i, kp+1, k));
						//}
					} // kp loop
				} // ip loop

				if(_DUMP_PPF_W_EXT_MESSAGES_ &&
					this->x_ext(i,k) != ZERO)
				{
					printf("W_ext[ext][comp](%d, %d) = %.5f\n", i, k, this->x_ext(i,k));
				}
#ifdef _LINEAR_COMPUTATIONS_
			}
			while(this->ppf_loops->ppf_scaler->check_pp_ppf_external_array_rescale(1, i, 1, k));
#endif

		} // k loop
	} // i loop
}

double& t_ppf_WEXT::x(int i, int k)
{
	return(this->pf_array[i][k]);
}

double& t_ppf_WEXT::x_ext(int i, int k)
{
	return(this->ext_pf_array[i][k]);
}

bool t_ppf_WEXT::check_boundary(int i, int k)
{
	if(t_template_pf_array::low_limits[i] <= k && t_template_pf_array::high_limits[i] >= k)
	{
		return(true);
	}
	else
	{
		//printf("Check boundary invalid: i1=%d -> i2=%d\n", i1, i2);
		return(false);
	}
}

// Calculate all values for W_ext array.
void t_ppf_WEXT::calculate_W_ext()
{
	//this->alloc_init_w_ext_array();

	for(int i = 1; i <= this->N1; i++)
	{
		int min_k = MAX(1, t_template_pf_array::low_limits[i]);
		int max_k = MIN(this->N2, t_template_pf_array::high_limits[i]);
		
		for(int k = min_k; k <= max_k; k++)
		{
#ifdef _LINEAR_COMPUTATIONS_
			do 
			{
#endif
	if(_DUMP_PPF_W_EXT_MESSAGES_)
				printf("WEXT(%d, %d)\n", i, k);

				// Initialize to 0.
				this->x(i, k) = ZERO;

				// Extend to 3' side.
				// insert i.
				if(i-1 >= 0 &&
					this->check_boundary(i-1, k))
				{
					double i_ins_unp_score = MUL(this->aln_priors->x(i, k, STATE_INS1), seq1_spf->ux_3p(i-1, i));
					this->x(i, k) = this->ppf_loops->MAX_SUM(this->x(i, k), MUL(this->x(i-1, k), i_ins_unp_score));

	if(_DUMP_PPF_W_EXT_MESSAGES_)
					printf("ins i: %lf x %lf\n", this->x(i-1, k), i_ins_unp_score);
				}
				else
				{
	if(_DUMP_PPF_W_EXT_MESSAGES_)
					printf("ins i not allowed\n");
				}

				// insert k.
				if(k-1 >= 0 &&
					this->check_boundary(i, k-1))
				{
					double k_ins_unp_score = MUL(this->aln_priors->x(i, k, STATE_INS2), seq2_spf->ux_3p(k-1, k));
					this->x(i, k) = this->ppf_loops->MAX_SUM(this->x(i, k), MUL(this->x(i, k-1), k_ins_unp_score));

	if(_DUMP_PPF_W_EXT_MESSAGES_)
					printf("ins k: %lf x %lf\n", this->x(i, k-1), k_ins_unp_score);
				}
				else
				{
	if(_DUMP_PPF_W_EXT_MESSAGES_)
					printf("ins k not allowed\n");
				}

				// align i and k.
				if(i-1 >= 0 &&
					k-1 >= 0 &&
					this->check_boundary(i-1, k-1))
				{
					double ik_aln_unp_score = MUL3(this->aln_priors->x(i, k, STATE_ALN), seq1_spf->ux_3p(i-1, i), seq2_spf->ux_3p(k-1, k));
					this->x(i, k) = this->ppf_loops->MAX_SUM(this->x(i, k), MUL(this->x(i-1, k-1), ik_aln_unp_score));

	if(_DUMP_PPF_W_EXT_MESSAGES_)
					printf("align i-k: %lf x (%lf x %lf x %lf)\n", this->x(i-1, k-1), this->aln_priors->x(i, k, STATE_ALN), seq1_spf->ux_3p(i-1, i), seq2_spf->ux_3p(k-1, k));
				}
				else
				{
	if(_DUMP_PPF_W_EXT_MESSAGES_)
					printf("i-k alignment not allowed\n");
				}

				// Extend via MHR's in V and V_mhe.
				// Search via ip, kp loop.
				// It is important to set the low limit of ip and kp to 0 so that concatenations of V and V_mhe 
				// V(1, j, 1, l) are counted.
				for(int ip = i-1; ip >= MAX(0, i-this->seq_man->ppf_cli->max_n_separation_between_nucs); ip--)
				{
					int min_kp = MAX(MAX(0, t_template_pf_array::low_limits[ip]), k-this->seq_man->ppf_cli->max_n_separation_between_nucs);
					int max_kp = MIN(k-1, t_template_pf_array::high_limits[ip]);

					for(int kp = max_kp; kp >= min_kp; kp--)
					{		
						//if(V->v_pf_array->check_4D_ll(ip+1, i, kp+1, k))
						{						
							double current_concatenation = MUL(this->x(ip, kp), this->ppf_loops->V_mhe->x(ip+1, i, kp+1, k));
							this->x(i,k) = this->ppf_loops->MAX_SUM(this->x(i,k), current_concatenation);
							if(_DUMP_PPF_W_EXT_MESSAGES_ &&  
								this->x(ip, kp) != ZERO && 
								this->ppf_loops->V_mhe->x(ip+1, i, kp+1, k) != ZERO)
							{
								printf("Concating Wext(%d, %d)-(Vmhe)(%d, %d, %d, %d): %lf x %lf\n", ip, kp, ip+1, i, kp+1, k, this->x(ip, kp), this->ppf_loops->V_mhe->x(ip+1, i, kp+1, k));
							}
						} // accession check.
					} // kp loop
				} // ip loop

	if(_DUMP_PPF_W_EXT_MESSAGES_)
	{
				printf("WEXT(%d, %d): %lf\n", i,k, this->x(i,k));
				printf("------------------------------------------------\n");
	}
#ifdef _LINEAR_COMPUTATIONS_
			}
			while(this->ppf_loops->ppf_scaler->check_pp_ppf_array_rescale(1, i, 1, k));
#endif
		} // k loop
	} // i loop
} // calculate_W_ext function.

void t_ppf_WEXT::rescale_internal(bool up_scale)
{
	for(int i = 0; i <= this->N1; i++)
	{
		int min_k = MAX(0, t_template_pf_array::low_limits[i]);
		int max_k = MIN(this->N2, t_template_pf_array::high_limits[i]);
		
		for(int k = min_k; k <= max_k; k++)
		{
if(_DUMP_PPF_W_EXT_MESSAGES_)
			printf("Rescaling WEXT(%d, %d)\n", i, k);

			double factor_de_rescale = pow(this->ppf_loops->ppf_scaler->rescaling_increment_factor_per_nucleotide, (i+k));

			if(!up_scale)
			{
				factor_de_rescale = 1.0f / factor_de_rescale;
			}

			this->x(i, k) = MUL(factor_de_rescale, this->x(i, k));
		} // k loop
	} // i loop
} // calculate_W_ext function.

void t_ppf_WEXT::rescale_external(bool up_scale)
{
	// Initialize.
	for(int i = this->N1+1; i >= 1; i--)
	{
		int min_k = MAX(1, t_template_pf_array::low_limits[i-1]+1);
		int max_k = MIN(this->N2+1, t_template_pf_array::high_limits[i-1]+1);
		
		for(int k = max_k; k >= min_k; k--)
		{
if(_DUMP_PPF_W_EXT_MESSAGES_)
			printf("Rescaling WEXT[ext](%d, %d)\n", i, k);

			double factor_de_rescale = pow(this->ppf_loops->ppf_scaler->rescaling_increment_factor_per_nucleotide, (this->seq_man->get_l_seq1() + this->seq_man->get_l_seq2()) - ((i-1)+(k-1)));

			if(!up_scale)
			{
				factor_de_rescale = 1.0f / factor_de_rescale;
			}

			this->x_ext(i, k) = MUL(factor_de_rescale, this->x_ext(i, k));
		} // k loop
	} // i loop
}

void t_ppf_WEXT::backtrack(int i, int k)
{
	double random_cumulative = ZERO;
	double current_cumulative = ZERO;
	bool pushed = false;

	// Setup cumulatives.
	if(this->ppf_loops->ppf_cli->mode == PARTS_RUN_MODE_STOCH_SAMPLE)
	{
		random_cumulative = MUL(this->ppf_loops->stoch_sampler->random_double_1_0(), this->x(i,k));
	}
	else if(this->ppf_loops->ppf_cli->mode == PARTS_RUN_MODE_MAP)
	{
		random_cumulative = this->x(i,k);
	}
	else
	{
		printf("Cannot run backtracking with mode %d\n", this->ppf_loops->ppf_cli->mode);
		exit(0);
	}

if(_DUMP_PPF_W_EXT_MESSAGES_)
	printf("Stoch_tb WEXT(%d, %d)\n", i, k);

	// If there is nothing to backtrack, return.
	if(i == 0 && k == 0)
	{
		return;
	}

	// Extend to 3' side.
	// insert i.
	if(i-1 >= 0 && 
		this->check_boundary(i-1, k))
	{
		double i_ins_unp_score = MUL(this->aln_priors->x(i, k, STATE_INS1), seq1_spf->ux_3p(i-1, i));
		//this->x(i, k) = this->ppf_loops->MAX_SUM(this->x(i, k), MUL(this->x(i-1, k), i_ins_unp_score));
		current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->x(i-1, k), i_ins_unp_score));

		if(!pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
		{
			//this->ppf_loops->map_alignment->set_seq1_ins(i,k);
			this->ppf_loops->set_seq1_ins(i,k);
			this->ppf_loops->tb_stack->push_str(1, i-1, 1, k, TRACE_W_ext);
			pushed = true;
		}
if(_DUMP_PPF_W_EXT_MESSAGES_)
		printf("ins i: %lf x %lf\n", this->x(i-1, k), i_ins_unp_score);
	}
	else
	{
if(_DUMP_PPF_W_EXT_MESSAGES_)
		printf("ins i not allowed\n");
	}

	// insert k.
	if(k-1 >= 0 &&
		this->check_boundary(i, k-1))
	{
		double k_ins_unp_score = MUL(this->aln_priors->x(i, k, STATE_INS2), seq2_spf->ux_3p(k-1, k));
		//this->x(i, k) = this->ppf_loops->MAX_SUM(this->x(i, k), MUL(this->x(i, k-1), k_ins_unp_score));
		current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->x(i, k-1), k_ins_unp_score));

		if(!pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
		{
			//this->ppf_loops->map_alignment->set_seq2_ins(i,k);
			this->ppf_loops->set_seq2_ins(i,k);
			this->ppf_loops->tb_stack->push_str(1, i, 1, k-1, TRACE_W_ext);
			pushed = true;
		}

if(_DUMP_PPF_W_EXT_MESSAGES_)
		printf("ins k: %lf x %lf\n", this->x(i, k-1), k_ins_unp_score);
	}
	else
	{
if(_DUMP_PPF_W_EXT_MESSAGES_)
		printf("ins k not allowed\n");
	}

	// align i and k.
	if(i-1 >= 0 && 
		k-1 >= 0 &&
		this->check_boundary(i-1, k-1))
	{
		double ik_aln_unp_score = MUL3(this->aln_priors->x(i, k, STATE_ALN), seq1_spf->ux_3p(i-1, i), seq2_spf->ux_3p(k-1, k));
		//this->x(i, k) = this->ppf_loops->MAX_SUM(this->x(i, k), MUL(this->x(i-1, k-1), ik_aln_unp_score));
		current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->x(i-1, k-1), ik_aln_unp_score));

		if(!pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
		{
			//this->ppf_loops->map_alignment->set_aln(i,k);
			this->ppf_loops->set_aln(i,k);
			this->ppf_loops->tb_stack->push_str(1, i-1, 1, k-1, TRACE_W_ext);
			pushed = true;
		}

if(_DUMP_PPF_W_EXT_MESSAGES_)
		printf("align i-k: %lf x (%lf x %lf x %lf)\n", this->x(i-1, k-1), this->aln_priors->x(i, k, STATE_ALN), seq1_spf->ux_3p(i-1, i), seq2_spf->ux_3p(k-1, k));
	}
	else
	{
if(_DUMP_PPF_W_EXT_MESSAGES_)
		printf("i-k alignment not allowed\n");
	}	

	// Extend via MHR's in V and V_mhe.
	// Search via ip, kp loop.
	for(int ip = i-1; ip >= MAX(0, i-this->seq_man->ppf_cli->max_n_separation_between_nucs); ip--)
	{
		int min_kp = MAX(MAX(0, t_template_pf_array::low_limits[ip]), k-this->seq_man->ppf_cli->max_n_separation_between_nucs);
		int max_kp = MIN(k-1, t_template_pf_array::high_limits[ip]);

		for(int kp = max_kp; kp >= min_kp; kp--)
		{		
			//if(V->v_pf_array->check_4D_ll(ip+1, i, kp+1, k))
			{
				double current_vmhe_concatenation = MUL(this->x(ip, kp), this->ppf_loops->V_mhe->x(ip+1, i, kp+1, k));
				current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, current_vmhe_concatenation);
				if(!pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
				{
					this->ppf_loops->tb_stack->push_str(ip+1, i, kp+1, k, TRACE_V_mhe);
					this->ppf_loops->tb_stack->push_str(1, ip, 1, kp, TRACE_W_ext);
					pushed = true;
				}

			} // accession check.
		} // kp loop
	} // ip loop

	if(!pushed || !COMPARE(this->x(i,k), current_cumulative))
	{
		printf("Traceback error @ %s(%d)\n", __FILE__, __LINE__);

		if(!pushed)
		{
			printf("Not pushed: %lf, %lf, %lf\n", current_cumulative, random_cumulative, this->x(i,k));
		}

		if(!COMPARE(this->x(i,k), current_cumulative))
		{
			printf("Not matching: %lf, %lf, %lf\n", current_cumulative, random_cumulative, this->x(i,k));
		}

		exit(0);
	}

if(_DUMP_PPF_W_EXT_MESSAGES_)
{
	printf("WEXT(%d, %d): %lf\n", i,k, this->x(i,k));
	printf("------------------------------------------------\n");
}
}


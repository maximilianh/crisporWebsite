#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "parts_compilation_directives.h"
#include "ppf_w.h"
#include "ppf_w_l.h"
#include "ppf_w_mb.h"
#include "ppf_w_mbl.h"
#include "ppf_v_mhe.h"
#include "process_sequences.h"
#include "ppf_math.h"
#include "math.h"
#include "structural_parameters.h"

#include "template_pf_array.h"
#include "single_pf_array.h"
#include "ppf_loops.h"

#include "ppf_tb_stack.h"
#include "map_structures.h"
#include "map_alignment.h"
#include "stoch_tb/stoch_sampling_math.h"
#include "ppf_cli.h"

#include "../../../src/phmm/structure/folding_constraints.h"

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

bool _DUMP_PPF_WMBL_MESSAGES_ = false;

t_ppf_WMBL::t_ppf_WMBL(t_ppf_loops* _ppf_loops)
{
if(_DUMP_PPF_WMBL_MESSAGES_)
	printf(" WMBL ");

	this->ppf_loops = _ppf_loops;

	this->seq1_spf = ppf_loops->seq1_spf;
	this->seq2_spf = ppf_loops->seq2_spf;
	this->aln_priors = aln_priors;

	// Copy sequence lengths.
	this->N1 = ppf_loops->seq_man->get_l_seq1();
	this->N2 = ppf_loops->seq_man->get_l_seq2();

	this->pf_array = new t_template_pf_array(this->ppf_loops->seq_man, 													
													seq1_spf->folding_constraints->coinc_pointer_relocation_map,
													seq2_spf->folding_constraints->coinc_pointer_relocation_map,
													true);

	this->ext_pf_array = new t_template_pf_array(this->ppf_loops->seq_man, 													
													seq1_spf->folding_constraints->coinc_pointer_relocation_map,
													seq2_spf->folding_constraints->coinc_pointer_relocation_map,
													true);
}

t_ppf_WMBL::~t_ppf_WMBL()
{
if(_DUMP_PPF_WMBL_MESSAGES_)
	printf("Destruct'ing a t_ppf_WMBL object.\n");

	delete(this->pf_array);
	delete(this->ext_pf_array);
}

bool t_ppf_WMBL::check_boundary(int i1, int i2)
{
	return(this->pf_array->check_boundary(i1,i2));
}

// Access to partition function array by reference.
double t_ppf_WMBL::x(int i, int j, int k, int l)
{
	if(this->pf_array->check_4D_ll(i,j,k,l))
	{
		return(this->pf_array->x(i, j, k, l));
	}
	else
	{
		return(ZERO);	
	}
}

double& t_ppf_WMBL::x_ext(int i, int j, int k, int l)
{
	if(this->pf_array->check_4D_ll(i,j,k,l))
	{
		return(this->ext_pf_array->x(i, j, k, l));
	}
	else
	{
		printf("%s(%d)\n", __FILE__, __LINE__);
		int* p = NULL;
		*p = 0;
		exit(0);
	}
}

double& t_ppf_WMBL::x_setter(int i, int j, int k, int l)
{
	if(this->pf_array->check_4D_ll(i,j,k,l))
	{
		return(this->pf_array->x(i, j, k, l));
	}
	else
	{
		printf("%s(%d)\n", __FILE__, __LINE__);
		exit(0);
	}
}

void t_ppf_WMBL::calculate_WMBL(int i, int j, int k, int l, bool bt) // Arrays to use.
{
	if(!this->pf_array->check_4D_ll(i,j,k,l))
	{
		//printf("Returning from WMBL calculation @ %s(%d) before calculating WMBL(%d, %d, %d, %d)\n", __FILE__, __LINE__, i,j,k,l);
		return; // Do not calculate this value.
	}

if(_DUMP_PPF_WMBL_MESSAGES_)
	printf("\nWMBL(%d, %d, %d, %d):\n", i,j,k,l);

	// Do calculations from WMBL side, i.e. by extending WMBL and WMB to 3' side from WMBL.
	// Set up backtracking if requested backtracking.
	double random_cumulative = ZERO;
	double current_cumulative = ZERO;
	bool pushed = false;
	if(bt)
	{
		if(this->ppf_loops->ppf_cli->mode == PARTS_RUN_MODE_STOCH_SAMPLE)
		{
			random_cumulative = MUL(this->ppf_loops->stoch_sampler->random_double_1_0(), this->x(i,j,k,l));
		}
		else if(this->ppf_loops->ppf_cli->mode == PARTS_RUN_MODE_MAP)
		{
			random_cumulative = this->x(i,j,k,l);
		}
		else
		{
			printf("Cannot run backtracking with mode %d\n", this->ppf_loops->ppf_cli->mode);
			exit(0);
		}
	}

	for(int ip = i + 1; ip < j; ip++)
	{
		int kp_start, kp_end;

		kp_start = max(k + 1, t_template_pf_array::low_limits[ip]);
		kp_end = min(l - 1, t_template_pf_array::high_limits[ip]);

		for(int kp = kp_start; kp <= kp_end; kp++)
		{
			if(this->pf_array->check_4D_ll(i, ip, k, kp))
			{
if(_DUMP_PPF_WMBL_MESSAGES_)
{
				printf("WMBL(%d, %d, %d, %d): ip = %d, kp = %d (N1 = %d, N2 = %d)\n", i,j,k,l, ip, kp, N1, N2);
}
				// V_mhe on 5'.
				current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->ppf_loops->V_mhe->x(i, ip, k, kp), this->ppf_loops->V_mhe->x(ip+1, j, kp+1, l)));
				if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
				{
					this->ppf_loops->tb_stack->push_str(i, ip, k, kp, TRACE_V_mhe);
					this->ppf_loops->tb_stack->push_str(ip+1, j, kp+1, l, TRACE_V_mhe);
					pushed = true;
				}

				current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->ppf_loops->V_mhe->x(i, ip, k, kp), this->ppf_loops->WL->x(ip+1, j, kp+1, l)));
				if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
				{
					this->ppf_loops->tb_stack->push_str(i, ip, k, kp, TRACE_V_mhe);
					this->ppf_loops->tb_stack->push_str(ip+1, j, kp+1, l, TRACE_WL);
					pushed = true;
				}

				current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->ppf_loops->V_mhe->x(i, ip, k, kp), this->ppf_loops->WMBL->x(ip+1, j, kp+1, l)));
				if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
				{
					this->ppf_loops->tb_stack->push_str(i, ip, k, kp, TRACE_V_mhe);
					this->ppf_loops->tb_stack->push_str(ip+1, j, kp+1, l, TRACE_WMBL);
					pushed = true;
				}

				// WL on 5'.
				current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->ppf_loops->WL->x(i, ip, k, kp), this->ppf_loops->V_mhe->x(ip+1, j, kp+1, l)));
				if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
				{
					this->ppf_loops->tb_stack->push_str(i, ip, k, kp, TRACE_WL);
					this->ppf_loops->tb_stack->push_str(ip+1, j, kp+1, l, TRACE_V_mhe);
					pushed = true;
				}

				current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->ppf_loops->WL->x(i, ip, k, kp), this->ppf_loops->WL->x(ip+1, j, kp+1, l)));
				if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
				{
					this->ppf_loops->tb_stack->push_str(i, ip, k, kp, TRACE_WL);
					this->ppf_loops->tb_stack->push_str(ip+1, j, kp+1, l, TRACE_WL);
					pushed = true;
				}

				current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->ppf_loops->WL->x(i, ip, k, kp), this->ppf_loops->WMBL->x(ip+1, j, kp+1, l)));
				if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
				{
					this->ppf_loops->tb_stack->push_str(i, ip, k, kp, TRACE_WL);
					this->ppf_loops->tb_stack->push_str(ip+1, j, kp+1, l, TRACE_WMBL);
					pushed = true;
				}
			} // bound check
		} // kp loop
	} // ip loop

	// Assign or compare.
	if(bt)
	{
		if(!COMPARE(this->x(i,j,k,l), current_cumulative))
		{
			printf("Traceback error @ %s(%d)\n", __FILE__, __LINE__);
			exit(0);
		}
	}
	else
	{
		this->x_setter(i,j,k,l) = current_cumulative;
	}

}

void t_ppf_WMBL::calculate_ext_dependencies(int i, int j, int k, int l) // Arrays to use.
{
	if(!this->pf_array->check_4D_ll(i,j,k,l))
	{
		//printf("Returning from WMBL calculation @ %s(%d) before calculating WMBL(%d, %d, %d, %d)\n", __FILE__, __LINE__, i,j,k,l);
		return; // Do not calculate this value.
	}

if(_DUMP_PPF_WMBL_MESSAGES_)
	printf("\nWMBL(%d, %d, %d, %d):\n", i,j,k,l);

	for(int ip = i + 1; ip < j; ip++)
	{
		int kp_start, kp_end;

		kp_start = max(k + 1, t_template_pf_array::low_limits[ip]);
		kp_end = min(l - 1, t_template_pf_array::high_limits[ip]);

		for(int kp = kp_start; kp <= kp_end; kp++)
		{
			if(this->pf_array->check_4D_ll(i, ip, k, kp))
			{
if(_DUMP_PPF_WMBL_MESSAGES_)
{
				printf("WMBL(%d, %d, %d, %d): ip = %d, kp = %d (N1 = %d, N2 = %d)\n", i,j,k,l, ip, kp, N1, N2);
}

				// V_mhe on 5'.
				//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
				//	MUL(this->ppf_loops->V_mhe->x(i, ip, k, kp), this->ppf_loops->V_mhe->x(ip+1, j, kp+1, l)));
				this->ppf_loops->V_mhe->x_ext(i, ip, k, kp) = this->ppf_loops->MAX_SUM(this->ppf_loops->V_mhe->x_ext(i, ip, k, kp), 
					MUL(this->x_ext(i,j,k,l), this->ppf_loops->V_mhe->x(ip+1, j, kp+1, l)));


				if(this->ppf_loops->V_mhe->pf_array->check_4D_ll(ip+1, j, kp+1, l))
				{
					this->ppf_loops->V_mhe->x_ext(ip+1, j, kp+1, l) = this->ppf_loops->MAX_SUM(this->ppf_loops->V_mhe->x_ext(ip+1, j, kp+1, l), 
						MUL(this->x_ext(i,j,k,l), this->ppf_loops->V_mhe->x(i, ip, k, kp)));
				}


				//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
				//	MUL(this->ppf_loops->V_mhe->x(i, ip, k, kp), this->ppf_loops->WL->x(ip+1, j, kp+1, l)));
				this->ppf_loops->V_mhe->x_ext(i, ip, k, kp) = this->ppf_loops->MAX_SUM(this->ppf_loops->V_mhe->x_ext(i, ip, k, kp), 
					MUL(this->x_ext(i,j,k,l), this->ppf_loops->WL->x(ip+1, j, kp+1, l)));


				if(this->ppf_loops->WL->pf_array->check_4D_ll(ip+1, j, kp+1, l))
				{
					this->ppf_loops->WL->x_ext(ip+1, j, kp+1, l) = this->ppf_loops->MAX_SUM(this->ppf_loops->WL->x_ext(ip+1, j, kp+1, l), 
						MUL(this->x_ext(i,j,k,l), this->ppf_loops->V_mhe->x(i, ip, k, kp)));
				}


				//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
				//	MUL(this->ppf_loops->V_mhe->x(i, ip, k, kp), this->ppf_loops->WMBL->x(ip+1, j, kp+1, l)));
				this->ppf_loops->V_mhe->x_ext(i, ip, k, kp) = this->ppf_loops->MAX_SUM(this->ppf_loops->V_mhe->x_ext(i, ip, k, kp), 
					MUL(this->x_ext(i,j,k,l), this->ppf_loops->WMBL->x(ip+1, j, kp+1, l)));


				if(this->ppf_loops->WMBL->pf_array->check_4D_ll(ip+1, j, kp+1, l))
				{
					this->ppf_loops->WMBL->x_ext(ip+1, j, kp+1, l) = this->ppf_loops->MAX_SUM(this->ppf_loops->WMBL->x_ext(ip+1, j, kp+1, l), 
						MUL(this->x_ext(i,j,k,l), this->ppf_loops->V_mhe->x(i, ip, k, kp)));
				}

				// WL on 5'.

				//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
				//	MUL(this->ppf_loops->WL->x(i, ip, k, kp), this->ppf_loops->V_mhe->x(ip+1, j, kp+1, l)));
				this->ppf_loops->WL->x_ext(i, ip, k, kp) = this->ppf_loops->MAX_SUM(this->ppf_loops->WL->x_ext(i, ip, k, kp), 
					MUL(this->x_ext(i,j,k,l), this->ppf_loops->V_mhe->x(ip+1, j, kp+1, l)));

				if(this->ppf_loops->V_mhe->pf_array->check_4D_ll(ip+1, j, kp+1, l))
				{
					this->ppf_loops->V_mhe->x_ext(ip+1, j, kp+1, l) = this->ppf_loops->MAX_SUM(this->ppf_loops->V_mhe->x_ext(ip+1, j, kp+1, l), 
						MUL(this->x_ext(i,j,k,l), this->ppf_loops->WL->x(i, ip, k, kp)));
				}


				//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
				//	MUL(this->ppf_loops->WL->x(i, ip, k, kp), this->ppf_loops->WL->x(ip+1, j, kp+1, l)));
				this->ppf_loops->WL->x_ext(i, ip, k, kp) = this->ppf_loops->MAX_SUM(this->ppf_loops->WL->x_ext(i, ip, k, kp), 
					MUL(this->x_ext(i,j,k,l), this->ppf_loops->WL->x(ip+1, j, kp+1, l)));

				if(this->ppf_loops->WL->pf_array->check_4D_ll(ip+1, j, kp+1, l))
				{
					this->ppf_loops->WL->x_ext(ip+1, j, kp+1, l) = this->ppf_loops->MAX_SUM(this->ppf_loops->WL->x_ext(ip+1, j, kp+1, l), 
						MUL(this->x_ext(i,j,k,l), this->ppf_loops->WL->x(i, ip, k, kp)));
				}

				//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
				//	MUL(this->ppf_loops->WL->x(i, ip, k, kp), this->ppf_loops->WMBL->x(ip+1, j, kp+1, l)));
				this->ppf_loops->WL->x_ext(i, ip, k, kp) = this->ppf_loops->MAX_SUM(this->ppf_loops->WL->x_ext(i, ip, k, kp), 
					MUL(this->x_ext(i,j,k,l), this->ppf_loops->WMBL->x(ip+1, j, kp+1, l)));

				if(this->ppf_loops->WMBL->pf_array->check_4D_ll(ip+1, j, kp+1, l))
				{
					this->ppf_loops->WMBL->x_ext(ip+1, j, kp+1, l) = this->ppf_loops->MAX_SUM(this->ppf_loops->WMBL->x_ext(ip+1, j, kp+1, l), 
						MUL(this->x_ext(i,j,k,l), this->ppf_loops->WL->x(i, ip, k, kp)));
				}
			} // bound check
		} // kp loop
	} // ip loop
}


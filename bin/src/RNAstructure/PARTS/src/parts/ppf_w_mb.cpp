#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "parts_compilation_directives.h"
#include "ppf_w_mb.h"
#include "ppf_w_mbl.h"
#include "ppf_w_l.h"
#include "ppf_v_mhe.h"
#include "process_sequences.h"
#include "ppf_math.h"
#include "alignment_priors.h"
#include "phmm_parameters.h"

#include "single_pf_array.h"
#include "template_pf_array.h"
			
#include "ppf_loops.h"		

#include "ppf_tb_stack.h"
#include "map_structures.h"
#include "map_alignment.h"
#include "stoch_tb/stoch_sampling_math.h"
#include "ppf_cli.h"

#include "../../../src/phmm/structure/folding_constraints.h"

extern bool problem;
extern bool map_problem;

bool _DUMP_PPF_WMB_MESSAGES_ = false;

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

t_ppf_WMB::t_ppf_WMB(t_ppf_loops* _ppf_loops)
{
if(_DUMP_PPF_WMB_MESSAGES_)
	printf(" WMB ");

	this->ppf_loops = _ppf_loops;

	this->seq1_spf = this->ppf_loops->seq1_spf;
	this->seq2_spf = this->ppf_loops->seq2_spf;
	this->aln_priors = this->ppf_loops->aln_priors;

	// Copy sequence lengths.
	this->N1 = this->ppf_loops->seq_man->get_l_seq1();
	this->N2 = this->ppf_loops->seq_man->get_l_seq2();

	this->pf_array = new t_template_pf_array(ppf_loops->seq_man, 
											seq1_spf->folding_constraints->coinc_pointer_relocation_map,
											seq2_spf->folding_constraints->coinc_pointer_relocation_map,
											true);

	this->ext_pf_array = new t_template_pf_array(ppf_loops->seq_man, 
											seq1_spf->folding_constraints->coinc_pointer_relocation_map,
											seq2_spf->folding_constraints->coinc_pointer_relocation_map,
											true);
}

t_ppf_WMB::~t_ppf_WMB()
{
if(_DUMP_PPF_WMB_MESSAGES_)
	printf("Destruct'ing a t_ppf_WMB object.\n");

	delete(this->pf_array);
	delete(this->ext_pf_array);
}

bool t_ppf_WMB::check_boundary(int i1, int i2)
{
	return(this->pf_array->check_boundary(i1,i2));
}

// Access to partition function array by reference.
double& t_ppf_WMB::x_setter(int i, int j, int k, int l)
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

double t_ppf_WMB::x(int i, int j, int k, int l)
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

double& t_ppf_WMB::x_ext(int i, int j, int k, int l)
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

void t_ppf_WMB::calculate_WMB(int i, int j, int k, int l, bool bt) // Arrays to use.
{
	if(!this->pf_array->check_4D_ll(i,j,k,l))
	{
		//printf("Returning from WMB calculation @ %s(%d) before calculating WMB(%d, %d, %d, %d)\n", __FILE__, __LINE__, i,j,k,l);
		return; // Do not calculate this value.
	}

if(_DUMP_PPF_WMB_MESSAGES_)
	printf("\nWMB(%d, %d, %d, %d):\n", i,j,k,l);

	// Do calculations from WMBL side, i.e. by extending WMBL and WMB to 3' side from WMBL.
	// Set up backtracking if requested backtracking.
	double random_cumulative = ZERO;
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

	double current_cumulative = this->ppf_loops->WMBL->x(i,j,k,l);

	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i, j, k, l, TRACE_WMBL);
		pushed = true;
	}

	double j_l_ALN_unpair_score = MUL(seq1_spf->ux_3p(i,j), seq2_spf->ux_3p(k,l));
	double j_INS_unpair_score = seq1_spf->ux_3p(i,j);
	double l_INS_unpair_score = seq2_spf->ux_3p(k,l);

if(_DUMP_PPF_WMB_MESSAGES_)
	printf("WMB(%d, %d, %d, %d): inited: %f\n", i,j,k,l,(current_cumulative));

	// Extend multibranch loops by one nucleotide on 3' side, which can happen with 3 cases.
	// j and l conditionals make sure that there is no jumping from 3' end to 5' end.
	if(this->pf_array->check_4D_ll(i, j-1, k, l-1))
 	{
		double j_l_ALN_prior = aln_priors->x(j, l, STATE_ALN);
		current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL3(j_l_ALN_prior, j_l_ALN_unpair_score, this->x(i,j-1,k,l-1)));	

		if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
		{
			this->ppf_loops->tb_stack->push_str(i, j-1, k, l-1, TRACE_WMB);
			//this->ppf_loops->map_alignment->set_aln(j,l);
			this->ppf_loops->set_aln(j,l);
			pushed = true;
		}

if(_DUMP_PPF_WMB_MESSAGES_)
		printf("WMB_WMBL(%d, %d, %d, %d) = %f recurse on ALNed insertion: WMB(%d, %d, %d, %d) = %f\n", i,j,k,l,(current_cumulative), i,j-1,k,l-1, (this->x(i,j-1,k,l-1)));
	}

	if(this->pf_array->check_4D_ll(i, j, k, l-1))
	{
		double l_INS_prior = aln_priors->x(j, l, STATE_INS2);
		current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL3(l_INS_prior, l_INS_unpair_score, this->x(i,j,k,l-1)));

		if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
		{
			this->ppf_loops->tb_stack->push_str(i, j, k, l-1, TRACE_WMB);
			//this->ppf_loops->map_alignment->set_seq2_ins(j,l);
			this->ppf_loops->set_seq2_ins(j,l);
			pushed = true;
		}

if(_DUMP_PPF_WMB_MESSAGES_)
		printf("WMB_WMBL(%d, %d, %d, %d) = %f recurse on l insertion: WMB(%d, %d, %d, %d) = %f\n", i, j, k, l, (current_cumulative), i,j,k,l-1, (this->x(i,j,k,l-1)));
	}

	if(this->pf_array->check_4D_ll(i, j-1, k, l))
	{
		double j_INS_prior = aln_priors->x(j, l, STATE_INS1);
		current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL3(j_INS_prior, j_INS_unpair_score, this->x(i,j-1,k,l)));

		if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
		{
			this->ppf_loops->tb_stack->push_str(i, j-1, k, l, TRACE_WMB);
			//this->ppf_loops->map_alignment->set_seq1_ins(j,l);
			this->ppf_loops->set_seq1_ins(j,l);
			pushed = true;
		}

if(_DUMP_PPF_WMB_MESSAGES_)
		printf("WMB_WMBL(%d, %d, %d, %d) = %f recurse on j insertion: WMB(%d, %d, %d, %d) = %f\n", i, j, k, l, (current_cumulative), i,j-1,k,l, (this->x(i,j-1,k,l)));
	}

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

if(_DUMP_PPF_WMB_MESSAGES_)
	printf("WMB(%d, %d, %d, %d) = %f\n", i,j,k,l, (this->x(i,j,k,l)));

}

void t_ppf_WMB::calculate_ext_dependencies(int i, int j, int k, int l) // Arrays to use.
{
	if(!this->pf_array->check_4D_ll(i,j,k,l))
	{
		//printf("Returning from WMB calculation @ %s(%d) before calculating WMB(%d, %d, %d, %d)\n", __FILE__, __LINE__, i,j,k,l);
		return; // Do not calculate this value.
	}

if(_DUMP_PPF_WMB_MESSAGES_)
	printf("\nWMB(%d, %d, %d, %d):\n", i,j,k,l);


	//double current_cumulative = this->ppf_loops->WMBL->x(i,j,k,l);
	this->ppf_loops->WMBL->x_ext(i,j,k,l) = this->ppf_loops->MAX_SUM(this->ppf_loops->WMBL->x_ext(i,j,k,l),
		this->x_ext(i,j,k,l));

	double j_l_ALN_unpair_score = MUL(seq1_spf->ux_3p(i,j), seq2_spf->ux_3p(k,l));
	double j_INS_unpair_score = seq1_spf->ux_3p(i,j);
	double l_INS_unpair_score = seq2_spf->ux_3p(k,l);


	// Extend multibranch loops by one nucleotide on 3' side, which can happen with 3 cases.
	// j and l conditionals make sure that there is no jumping from 3' end to 5' end.
	if(this->pf_array->check_4D_ll(i, j-1, k, l-1))
 	{
		double j_l_ALN_prior = aln_priors->x(j, l, STATE_ALN);

		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL3(j_l_ALN_prior, j_l_ALN_unpair_score, this->x(i,j-1,k,l-1)));	
		this->ppf_loops->WMB->x_ext(i,j-1,k,l-1) = this->ppf_loops->MAX_SUM(this->ppf_loops->WMB->x_ext(i,j-1,k,l-1),
			MUL3(this->x_ext(i,j,k,l), j_l_ALN_prior, j_l_ALN_unpair_score));
	}

	if(this->pf_array->check_4D_ll(i, j, k, l-1))
	{
		double l_INS_prior = aln_priors->x(j, l, STATE_INS2);
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL3(l_INS_prior, l_INS_unpair_score, this->x(i,j,k,l-1)));
		this->ppf_loops->WMB->x_ext(i,j,k,l-1) = this->ppf_loops->MAX_SUM(this->ppf_loops->WMB->x_ext(i,j,k,l-1),
			MUL3(this->x_ext(i,j,k,l), l_INS_prior, l_INS_unpair_score));
	}

	if(this->pf_array->check_4D_ll(i, j-1, k, l))
	{
		double j_INS_prior = aln_priors->x(j, l, STATE_INS1);
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL3(j_INS_prior, j_INS_unpair_score, this->x(i,j-1,k,l)));
		this->ppf_loops->WMB->x_ext(i,j-1,k,l) = this->ppf_loops->MAX_SUM(this->ppf_loops->WMB->x_ext(i,j-1,k,l),
			MUL3(this->x_ext(i,j,k,l), j_INS_prior, j_INS_unpair_score));
	}
}


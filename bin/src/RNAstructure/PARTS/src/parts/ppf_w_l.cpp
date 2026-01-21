#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "parts_compilation_directives.h"
#include "ppf_w_l.h"
#include "process_sequences.h"
#include "ppf_math.h"
#include <iostream>
#include "alignment_priors.h"
#include "phmm_parameters.h"

#include "single_pf_array.h"
#include "template_pf_array.h"
#include "ppf_v_mhe.h"
#include "ppf_loops.h"

#include "ppf_tb_stack.h"
#include "map_structures.h"
#include "map_alignment.h"
#include "stoch_tb/stoch_sampling_math.h"
#include "ppf_cli.h"

#include "../../../src/phmm/structure/folding_constraints.h"

using namespace std;

bool _DUMP_PPF_WL_MESSAGES_ = false;

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

t_ppf_WL::t_ppf_WL(t_ppf_loops* _ppf_loops)
{
if(_DUMP_PPF_WL_MESSAGES_)
	printf(" WL ");

	this->ppf_loops = _ppf_loops;

	this->seq1_spf = this->ppf_loops->seq1_spf;
	this->seq2_spf = this->ppf_loops->seq2_spf;
	this->aln_priors = this->ppf_loops->aln_priors;

	// Copy sequence lengths.
	this->N1 = this->ppf_loops->seq_man->get_l_seq1();
	this->N2 = this->ppf_loops->seq_man->get_l_seq2();

	this->pf_array = new t_template_pf_array(this->ppf_loops->seq_man,
													seq1_spf->folding_constraints->coinc_pointer_relocation_map,
													seq2_spf->folding_constraints->coinc_pointer_relocation_map,
													true);

	this->ext_pf_array = new t_template_pf_array(this->ppf_loops->seq_man,
													seq1_spf->folding_constraints->coinc_pointer_relocation_map,
													seq2_spf->folding_constraints->coinc_pointer_relocation_map,
													true);
}

t_ppf_WL::~t_ppf_WL()
{
if(_DUMP_PPF_WL_MESSAGES_)
	printf("Destruct'ing a t_ppf_WL object.\n");

	delete(this->pf_array);
	delete(this->ext_pf_array);

}

bool t_ppf_WL::check_boundary(int i1, int i2)
{
	return(this->pf_array->check_boundary(i1, i2));
}

// Access to partition function array by reference.
double& t_ppf_WL::x_setter(int i, int j, int k, int l)
{
	if(this->pf_array->check_4D_ll(i,j,k,l))
	{
		return(this->pf_array->x(i,j,k,l));
	}
	else
	{
		printf("%s(%d)\n", __FILE__, __LINE__);
		exit(0);
	}
}

double t_ppf_WL::x(int i, int j, int k, int l)
{
	if(this->pf_array->check_4D_ll(i,j,k,l))
	{
		return(this->pf_array->x(i,j,k,l));
	}
	else
	{
		return(ZERO);
	}
}

double& t_ppf_WL::x_ext(int i, int j, int k, int l)
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

void t_ppf_WL::calculate_WL(int i, int j, int k, int l, bool bt) // Arrays to use.
{
	// Do the first check if the alignment positions are inside alignment constraint.
	// Check for i-1, k-1.
	if( !this->pf_array->check_4D_ll(i,j,k,l) )
	{
		return;
	}

	bool i_k = this->pf_array->check_boundary(i,k);
	bool i_dec_k_dec = this->pf_array->check_boundary(i-1, k-1);
	bool i_dec_k = this->pf_array->check_boundary(i-1, k);
	bool i_k_dec = this->pf_array->check_boundary(i, k-1);

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

if(_DUMP_PPF_WL_MESSAGES_)
	printf("WL(%d, %d, %d, %d) inited as %f\n", i,j,k,l, (this->x(i, j, k, l)) );

	// ik align
	double i_k_ALN_unpair_score = MUL(seq1_spf->ux_5p(i,j), seq2_spf->ux_5p(k,l));

	double i_k_ALN_prior = ZERO;
	if(i_k)
	{
		i_k_ALN_prior = aln_priors->x(i, k, STATE_ALN);
	}
	current_cumulative = this->ppf_loops->MAX_SUM( current_cumulative, MUL(this->x(i+1, j, k+1, l), MUL(i_k_ALN_prior, i_k_ALN_unpair_score)) ); // i and k inserted: 
	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i+1, j, k+1, l, TRACE_WL);
		this->ppf_loops->set_aln(i,k);
		pushed = true;
	}

	current_cumulative = this->ppf_loops->MAX_SUM( current_cumulative, MUL(this->ppf_loops->V_mhe->x(i+1, j, k+1, l), MUL(i_k_ALN_prior, i_k_ALN_unpair_score)) ); // i and k inserted: 
	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i+1, j, k+1, l, TRACE_V_mhe);
		this->ppf_loops->set_aln(i,k);
		pushed = true;
	}

if(_DUMP_PPF_WL_MESSAGES_)
{
	printf("WL(%d, %d, %d, %d): recursed on ALNed insertion with score WL(%d, %d, %d, %d) = %f\n", i,j,k,l, i+1, j,k+1,l, (this->x(i+1, j, k+1, l)) );
	printf("i_k_ALN_unpair_score = %f, i_k_ALN_prior = %f\n", i_k_ALN_unpair_score, i_k_ALN_prior);
}

	// k insert.
	double k_INS_unpair_score = seq2_spf->ux_5p(k,l);

	double k_INS_prior = ZERO;
	if(i_dec_k)
	{
		k_INS_prior = aln_priors->x(i - 1, k, STATE_INS2);
	}
	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->x(i, j, k+1, l), MUL(k_INS_prior, k_INS_unpair_score)) ); // k is inserted.
	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i, j, k+1, l, TRACE_WL);
		this->ppf_loops->set_seq2_ins(i-1,k);
		pushed = true;
	}

	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->ppf_loops->V_mhe->x(i, j, k+1, l), MUL(k_INS_prior, k_INS_unpair_score)) ); // k is inserted.
	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i, j, k+1, l, TRACE_V_mhe);
		this->ppf_loops->set_seq2_ins(i-1,k);
		pushed = true;
	}

if(_DUMP_PPF_WL_MESSAGES_)
{
	printf("WL(%d, %d, %d, %d): recursed on k insertion with score WL(%d, %d, %d, %d) = %f\n", i,j,k,l,i, j, k+1, l, (this->x(i, j, k+1, l)) );
	printf("k_INS_unpair_score = %f, k_INS_prior = %f\n", k_INS_unpair_score, k_INS_prior);
}

	// i insert
	double i_INS_unpair_score = seq1_spf->ux_5p(i,j);

	double i_INS_prior = ZERO;
	if(i_k_dec)
	{
		i_INS_prior = aln_priors->x(i, k - 1, STATE_INS1);
	}
	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->x(i+1, j, k, l), MUL(i_INS_prior, i_INS_unpair_score)) ); // i is inserted.
	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i+1, j, k, l, TRACE_WL);
		this->ppf_loops->set_seq1_ins(i,k-1);
		pushed = true;
	}

	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->ppf_loops->V_mhe->x(i+1, j, k, l), MUL(i_INS_prior, i_INS_unpair_score)) ); // i is inserted.
	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i+1, j, k, l, TRACE_V_mhe);
		this->ppf_loops->set_seq1_ins(i,k-1);
		pushed = true;
	}

if(_DUMP_PPF_WL_MESSAGES_)
{
	printf("WL(%d, %d, %d, %d): recursed on i insertion with score WL(%d, %d, %d, %d) = %f\n", i,j,k,l,i+1, j,k,l, (this->x(i+1, j, k, l)) );
	printf("i_INS_unpair_score = %f, i_INS_prior = %f\n", i_INS_unpair_score, i_INS_prior);
}

if(_DUMP_PPF_WL_MESSAGES_)
	printf("WL(%d, %d, %d, %d) = %f\n\n", i,j,k,l, (this->x(i,j,k,l)) );

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

void t_ppf_WL::calculate_ext_dependencies(int i, int j, int k, int l) // Arrays to use.
{
	// Do the first check if the alignment positions are inside alignment constraint.
	// Check for i-1, k-1.
	if( !this->pf_array->check_4D_ll(i,j,k,l) )
	{
		return;
	}

	bool i_k = this->pf_array->check_boundary(i,k);
	bool i_dec_k_dec = this->pf_array->check_boundary(i-1, k-1);
	bool i_dec_k = this->pf_array->check_boundary(i-1, k);
	bool i_k_dec = this->pf_array->check_boundary(i, k-1);

	// ik align
	double i_k_ALN_unpair_score = MUL(seq1_spf->ux_5p(i,j), seq2_spf->ux_5p(k,l));

	double i_k_ALN_prior = ZERO;
	if(i_k)
	{
		i_k_ALN_prior = aln_priors->x(i, k, STATE_ALN);
	}

	if(this->pf_array->check_4D_ll(i+1, j, k+1, l))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM( current_cumulative, MUL(this->x(i+1, j, k+1, l), MUL(i_k_ALN_prior, i_k_ALN_unpair_score)) ); // i and k inserted: 
		this->x_ext(i+1, j, k+1, l) = this->ppf_loops->MAX_SUM(this->x_ext(i+1, j, k+1, l), 
			MUL3(i_k_ALN_prior, i_k_ALN_unpair_score, this->x_ext(i,j,k,l)));
	}

	if(this->ppf_loops->V_mhe->pf_array->check_4D_ll(i+1, j, k+1, l))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM( current_cumulative, MUL(this->ppf_loops->V_mhe->x(i+1, j, k+1, l), MUL(i_k_ALN_prior, i_k_ALN_unpair_score)) ); // i and k inserted: 
		this->ppf_loops->V_mhe->x_ext(i+1, j, k+1, l) = this->ppf_loops->MAX_SUM(this->ppf_loops->V_mhe->x_ext(i+1, j, k+1, l), 
			MUL3(i_k_ALN_prior, i_k_ALN_unpair_score, this->x_ext(i,j,k,l)));
	}

	// k insert.
	double k_INS_unpair_score = seq2_spf->ux_5p(k,l);

	double k_INS_prior = ZERO;
	if(i_dec_k)
	{
		k_INS_prior = aln_priors->x(i - 1, k, STATE_INS2);
	}

	if(this->pf_array->check_4D_ll(i, j, k+1, l))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->x(i, j, k+1, l), MUL(k_INS_prior, k_INS_unpair_score)) ); // k is inserted.
		this->x_ext(i, j, k+1, l) = this->ppf_loops->MAX_SUM(this->x_ext(i, j, k+1, l), 
			MUL3(k_INS_prior, k_INS_unpair_score, this->x_ext(i,j,k,l)));
	}

	if(this->ppf_loops->V_mhe->pf_array->check_4D_ll(i, j, k+1, l))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->ppf_loops->V_mhe->x(i, j, k+1, l), MUL(k_INS_prior, k_INS_unpair_score)) ); // k is inserted.
		this->ppf_loops->V_mhe->x_ext(i, j, k+1, l) = this->ppf_loops->MAX_SUM(this->ppf_loops->V_mhe->x_ext(i, j, k+1, l), 
			MUL3(k_INS_prior, k_INS_unpair_score, this->x_ext(i,j,k,l)));
	}

	// i insert
	double i_INS_unpair_score = seq1_spf->ux_5p(i,j);

	double i_INS_prior = ZERO;
	if(i_k_dec)
	{
		i_INS_prior = aln_priors->x(i, k - 1, STATE_INS1);
	}

	if(this->pf_array->check_4D_ll(i+1, j, k, l))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->x(i+1, j, k, l), MUL(i_INS_prior, i_INS_unpair_score)) ); // i is inserted.
		this->x_ext(i+1, j, k, l) = this->ppf_loops->MAX_SUM(this->x_ext(i+1, j, k, l), 
			MUL3(i_INS_prior, i_INS_unpair_score, this->x_ext(i,j,k,l)));
	}

	if(this->ppf_loops->V_mhe->pf_array->check_4D_ll(i+1, j, k, l))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->ppf_loops->V_mhe->x(i+1, j, k, l), MUL(i_INS_prior, i_INS_unpair_score)) ); // i is inserted.
		this->ppf_loops->V_mhe->x_ext(i+1, j, k, l) = this->ppf_loops->MAX_SUM(this->ppf_loops->V_mhe->x_ext(i+1, j, k, l), 
			MUL3(i_INS_prior, i_INS_unpair_score, this->x_ext(i,j,k,l)));
	}
}


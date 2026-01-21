#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "parts_compilation_directives.h"
#include "ppf_w.h"
#include "ppf_w_l.h"
#include "ppf_v_mhe.h"
#include "process_sequences.h"
#include "ppf_math.h"
#include <iostream>
#include "alignment_priors.h"
#include "phmm_parameters.h"

#include "single_pf_array.h"
#include "template_pf_array.h"

#include "ppf_loops.h"
//#include "../../../src/phmm/nnm_energy/xlog_math.h"

#include "ppf_tb_stack.h"
#include "map_structures.h"
#include "map_alignment.h"
#include "stoch_tb/stoch_sampling_math.h"
#include "ppf_cli.h"

#include "../../../src/phmm/structure/folding_constraints.h"

using namespace std;

extern bool problem;
extern bool map_problem;

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

bool _DUMP_PPF_W_MESSAGES_ = false;

t_ppf_W::t_ppf_W(t_ppf_loops* _ppf_loops)
{
if(_DUMP_PPF_W_MESSAGES_)
	printf(" W ");

	this->ppf_loops = _ppf_loops;

	// Copy sequence lengths.
	this->N1 = this->ppf_loops->seq_man->get_l_seq1();
	this->N2 = this->ppf_loops->seq_man->get_l_seq2();

	this->seq1_spf = this->ppf_loops->seq1_spf;
	this->seq2_spf = this->ppf_loops->seq2_spf;
	this->aln_priors = this->ppf_loops->aln_priors;

	this->pf_array = new t_template_pf_array(this->ppf_loops->seq_man, 
												seq1_spf->folding_constraints->coinc_pointer_relocation_map,
												seq2_spf->folding_constraints->coinc_pointer_relocation_map,
												true);

	this->ext_pf_array = new t_template_pf_array(this->ppf_loops->seq_man, 
												seq1_spf->folding_constraints->coinc_pointer_relocation_map,
												seq2_spf->folding_constraints->coinc_pointer_relocation_map,
												true);

}

t_ppf_W::~t_ppf_W()
{
if(_DUMP_PPF_W_MESSAGES_)
	printf("Destruct'ing a t_ppf_W object.\n");

	delete(this->pf_array);
	delete(this->ext_pf_array);
}


bool t_ppf_W::check_boundary(int i1, int i2)
{
	return(this->pf_array->check_boundary(i1, i2));
}

// Access to partition function array by reference.
double& t_ppf_W::x_setter(int i, int j, int k, int l)
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

double t_ppf_W::x(int i, int j, int k, int l)
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

double& t_ppf_W::x_ext(int i, int j, int k, int l)
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

void t_ppf_W::calculate_W(int i, int j, int k, int l, bool bt) // Arrays to use.
{
	// Have to fork for j > N, do bound checking...
	if( !this->pf_array->check_4D_ll(i,j,k,l) )
	{
		//printf("Returning from W calculation @ %s(%d) before calculating W(%d, %d, %d, %d)\n", __FILE__, __LINE__, i,j,k,l);
		return;
	}

if(_DUMP_PPF_W_MESSAGES_)
	printf("\nW(%d, %d, %d, %d):\n", i,j,k,l);

	// Following are needed to check if i and k unpaired nucleotide additions are inside
	// alignment constraints.
	bool i_dec_k_dec = this->check_boundary(i - 1, k - 1);
	bool i_dec_k = this->check_boundary(i - 1, k);
	bool i_k_dec = this->check_boundary(i, k - 1);

	// Following are needed for extending W to 5' side, that is,
	// when calculating W from WR.
	bool i_inc_k_inc = this->check_boundary(i + 1, k + 1);
	bool i_inc_k = this->check_boundary(i + 1, k);
	bool i_k_inc = this->check_boundary(i, k + 1);

	// Following are needed for extending W to 3' side, that is,
	// when calculating W from WL.
	bool j_dec_l_dec = this->check_boundary(j - 1 , l - 1);
	bool j_l_dec = this->check_boundary(j, l - 1);
	bool j_dec_l = this->check_boundary(j - 1 , l);

	// Calculate W from WL.
	// Init W(i,j,k,l) as WL(i,j,k,l).
	// And extend it through 3' side.
	double j_l_ALN_unpair_score = MUL(seq1_spf->ux_3p(i,j), seq2_spf->ux_3p(k,l));
	double j_INS_unpair_score = seq1_spf->ux_3p(i,j);
	double l_INS_unpair_score = seq2_spf->ux_3p(k,l);

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

	// Initialize with WL.
	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, this->ppf_loops->WL->x(i,j,k,l));
	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i, j, k, l, TRACE_WL);
		pushed = true;
	}

if(_DUMP_PPF_W_MESSAGES_)
	printf("W_WL(%d, %d, %d, %d): inited: %f\n", i,j,k,l,(current_cumulative));

	// Note that the same can be done by initing W as WR and extending to 5' side, this is in fact a very good test case for verifying partition function.
	// align ij
	double j_l_ALN_prior = aln_priors->x(j, l, STATE_ALN);
	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->x(i, j - 1, k, l - 1), MUL(j_l_ALN_prior, j_l_ALN_unpair_score)) );
	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i, j-1, k, l-1, TRACE_W);
		this->ppf_loops->set_aln(j,l);
		pushed = true;
	}

	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->ppf_loops->V_mhe->x(i, j - 1, k, l - 1), MUL(j_l_ALN_prior, j_l_ALN_unpair_score)) );
	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i, j-1, k, l-1, TRACE_V_mhe);
		this->ppf_loops->set_aln(j,l);
		pushed = true;
	}

if(_DUMP_PPF_W_MESSAGES_)
	printf("W_WL(%d, %d, %d, %d) = %f, recursed on ALNed insertion with score W(%d, %d, %d, %d) = %f\n", i,j,k,l, (current_cumulative), i, j - 1, k, l - 1, (this->x(i, j - 1, k, l - 1)));

	// insert l
	double l_INS2_prior = aln_priors->x(j, l, STATE_INS2);
	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->x(i, j, k, l - 1), MUL(l_INS2_prior, l_INS_unpair_score)) );
	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i, j, k, l-1, TRACE_W);
		this->ppf_loops->set_seq2_ins(j,l);
		pushed = true;
	}

	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->ppf_loops->V_mhe->x(i, j, k, l - 1), MUL(l_INS2_prior, l_INS_unpair_score)) );
	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i, j, k, l-1, TRACE_V_mhe);
		this->ppf_loops->set_seq2_ins(j,l);
		pushed = true;
	}

if(_DUMP_PPF_W_MESSAGES_)
	printf("W_WL(%d, %d, %d, %d) = %f, recursed on l insertion with score W(%d, %d, %d, %d) = %f\n", i,j,k,l, (current_cumulative), i, j, k, l - 1, (this->x(i, j, k, l - 1)));

	// insert j
	double j_INS1_prior = aln_priors->x(j, l, STATE_INS1);
	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->x(i, j - 1, k, l), MUL(j_INS1_prior, j_INS_unpair_score)) );
	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i, j-1, k, l, TRACE_W);
		this->ppf_loops->set_seq1_ins(j,l);
		pushed = true;
	}

	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->ppf_loops->V_mhe->x(i, j - 1, k, l), MUL(j_INS1_prior, j_INS_unpair_score)) );
	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i, j-1, k, l, TRACE_V_mhe);
		this->ppf_loops->set_seq1_ins(j,l);
		pushed = true;
	}

if(_DUMP_PPF_W_MESSAGES_)
	printf("W_WL(%d, %d, %d, %d) = %f, recursed on j insertion with score W(%d, %d, %d, %d) = %f\n", i, j, k, l, (current_cumulative), i, j - 1, k, l, (this->x(i, j - 1, k, l)));

	if(current_cumulative < ZERO)
	{
		printf("W_WL(%d, %d, %d, %d) = %f\n", i,j,k,l,current_cumulative );
		getc(stdin);
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

if(_DUMP_PPF_W_MESSAGES_)
	printf("W(%d, %d, %d, %d) = %f\n", i, j, k, l, (this->x(i,j,k,l)));
}

void t_ppf_W::calculate_ext_dependencies(int i, int j, int k, int l) // Arrays to use.
{
	// Have to fork for j > N, do bound checking...
	if( !this->pf_array->check_4D_ll(i,j,k,l) )
	{
		//printf("Returning from W calculation @ %s(%d) before calculating W(%d, %d, %d, %d)\n", __FILE__, __LINE__, i,j,k,l);
		return;
	}

if(_DUMP_PPF_W_MESSAGES_)
	printf("\nW(%d, %d, %d, %d):\n", i,j,k,l);

	// Following are needed to check if i and k unpaired nucleotide additions are inside
	// alignment constraints.
	bool i_dec_k_dec = this->check_boundary(i - 1, k - 1);
	bool i_dec_k = this->check_boundary(i - 1, k);
	bool i_k_dec = this->check_boundary(i, k - 1);

	// Following are needed for extending W to 5' side, that is,
	// when calculating W from WR.
	bool i_inc_k_inc = this->check_boundary(i + 1, k + 1);
	bool i_inc_k = this->check_boundary(i + 1, k);
	bool i_k_inc = this->check_boundary(i, k + 1);

	// Following are needed for extending W to 3' side, that is,
	// when calculating W from WL.
	bool j_dec_l_dec = this->check_boundary(j - 1 , l - 1);
	bool j_l_dec = this->check_boundary(j, l - 1);
	bool j_dec_l = this->check_boundary(j - 1 , l);

	// Calculate W from WL.
	// Init W(i,j,k,l) as WL(i,j,k,l).
	// And extend it through 3' side.
	double W_WL = this->ppf_loops->WL->x(i,j,k,l);
	double j_l_ALN_unpair_score = MUL(seq1_spf->ux_3p(i,j), seq2_spf->ux_3p(k,l));
	double j_INS_unpair_score = seq1_spf->ux_3p(i,j);
	double l_INS_unpair_score = seq2_spf->ux_3p(k,l);

if(_DUMP_PPF_W_MESSAGES_)
	printf("W_WL(%d, %d, %d, %d): inited: %f\n", i,j,k,l,(W_WL));

	// Initialize with W_WL.	
	if(this->ppf_loops->WL->pf_array->check_4D_ll(i,j,k,l))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, this->ppf_loops->WL->x(i,j,k,l));
		this->ppf_loops->WL->x_ext(i,j,k,l) = this->ppf_loops->MAX_SUM(this->ppf_loops->WL->x_ext(i,j,k,l), this->x_ext(i,j,k,l));
	}

	// Note that the same can be done by initing W as WR and extending to 5' side, this is in fact a very good test case for verifying partition function.
	// align ij
	double j_l_ALN_prior = aln_priors->x(j, l, STATE_ALN);

	if(this->pf_array->check_4D_ll(i, j-1, k, l-1))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->x(i, j - 1, k, l - 1), MUL(j_l_ALN_prior, j_l_ALN_unpair_score)) );
		this->x_ext(i,j-1,k,l-1) = this->ppf_loops->MAX_SUM(this->x_ext(i,j-1,k,l-1), 
			MUL3(this->x_ext(i,j,k,l), j_l_ALN_prior, j_l_ALN_unpair_score));
	}

	if(this->ppf_loops->V_mhe->pf_array->check_4D_ll(i, j-1, k, l-1))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->ppf_loops->V_mhe->x(i, j - 1, k, l - 1), MUL(j_l_ALN_prior, j_l_ALN_unpair_score)) );
		this->ppf_loops->V_mhe->x_ext(i,j-1,k,l-1) = this->ppf_loops->MAX_SUM(this->ppf_loops->V_mhe->x_ext(i,j-1,k,l-1), 
			MUL3(this->x_ext(i,j,k,l), j_l_ALN_prior, j_l_ALN_unpair_score));
	}

	// insert l
	double l_INS2_prior = aln_priors->x(j, l, STATE_INS2);

	if(this->pf_array->check_4D_ll(i, j, k, l-1))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->x(i, j, k, l - 1), MUL(l_INS2_prior, l_INS_unpair_score)) );
		this->x_ext(i,j,k,l-1) = this->ppf_loops->MAX_SUM(this->x_ext(i,j,k,l-1), 
			MUL3(this->x_ext(i,j,k,l), l_INS2_prior, l_INS_unpair_score));
	}

	if(this->ppf_loops->V_mhe->pf_array->check_4D_ll(i, j, k, l-1))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->ppf_loops->V_mhe->x(i, j, k, l - 1), MUL(l_INS2_prior, l_INS_unpair_score)) );
		this->ppf_loops->V_mhe->x_ext(i,j,k,l-1) = this->ppf_loops->MAX_SUM(this->ppf_loops->V_mhe->x_ext(i,j,k,l-1), 
			MUL3(this->x_ext(i,j,k,l),  l_INS2_prior, l_INS_unpair_score));
	}

	// insert j
	double j_INS1_prior = aln_priors->x(j, l, STATE_INS1);

	if(this->pf_array->check_4D_ll(i, j-1, k, l))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->x(i, j - 1, k, l), MUL(j_INS1_prior, j_INS_unpair_score)) );
		this->x_ext(i,j-1,k,l) = this->ppf_loops->MAX_SUM(this->x_ext(i,j-1,k,l), 
			MUL3(this->x_ext(i,j,k,l), j_INS1_prior, j_INS_unpair_score));
	}

	if(this->ppf_loops->V_mhe->pf_array->check_4D_ll(i, j-1, k, l))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, MUL(this->ppf_loops->V_mhe->x(i, j - 1, k, l), MUL(j_INS1_prior, j_INS_unpair_score)) );
		this->ppf_loops->V_mhe->x_ext(i,j-1,k,l) = this->ppf_loops->MAX_SUM(this->ppf_loops->V_mhe->x_ext(i,j-1,k,l), 
			MUL3(this->x_ext(i,j,k,l),  j_INS1_prior, j_INS_unpair_score));
	}
}

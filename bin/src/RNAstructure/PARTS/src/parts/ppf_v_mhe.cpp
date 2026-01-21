#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "parts_compilation_directives.h"
#include "ppf_v_mhe.h"
#include "ppf_w_mhi.h"
#include "ppf_w.h"
#include "ppf_ss.h"
#include "ppf_w_mb.h"
#include "ppf_cli.h"
#include "stoch_tb/stoch_sampling_math.h"
#include "process_sequences.h"
#include "ppf_math.h"
#include "alignment_priors.h"
#include "phmm_parameters.h"
#include "template_pf_array.h"
#include "single_pf_array.h"
#include "ppf_loops.h"
#include "ppf_tb_stack.h"
#include "map_structures.h"
#include "map_alignment.h"
#include <iostream>

#include "../../../src/phmm/structure/folding_constraints.h"

using namespace std;

// t_ppf_V_mhe: Handles base pair insertions on V array, note that t_ppf_V_mhe inserts base pairs one by one to each sequence
// so it can handle INS1 followed by an INS2 or vice versa, so this is basically different in how insertions are handled
// by consecutive insertion handling with t_ppf_V::calculate_W_helix_ins(...). In the calculations this will go to
// initialization of W array. t_ppf_V_mhe covers structural alignments where i-j and k-l are paired however either one of those
// are inserted. And base pair insertions start on top of aligned base pairs, that is, on top of V array.

extern bool problem;

bool _DUMP_PPF_V_MHE_MESSAGES_ = false;

t_ppf_V_mhe::t_ppf_V_mhe(t_ppf_loops* _ppf_loops)
{
if(_DUMP_PPF_V_MHE_MESSAGES_)
	printf(" VMHE ");

	this->ppf_loops = _ppf_loops;

	// Copy sequence lengths.
	this->N1 = ppf_loops->seq_man->get_l_seq1();
	this->N2 = ppf_loops->seq_man->get_l_seq2();

	this->seq1_spf = ppf_loops->seq1_spf;
	this->seq2_spf = ppf_loops->seq2_spf;
	this->aln_priors = ppf_loops->aln_priors;


	// Alloc and init two template pf arrays.
	// Those are needed for keeping track of insertion places.
	this->pf_array = new t_template_pf_array(this->ppf_loops->seq_man, 													
											seq1_spf->folding_constraints->coinc_pointer_relocation_map,
											seq2_spf->folding_constraints->coinc_pointer_relocation_map,
											true);

	this->ext_pf_array = new t_template_pf_array(this->ppf_loops->seq_man, 													
											seq1_spf->folding_constraints->coinc_pointer_relocation_map,
											seq2_spf->folding_constraints->coinc_pointer_relocation_map,
											true);
}

t_ppf_V_mhe::~t_ppf_V_mhe()
{
if(_DUMP_PPF_V_MHE_MESSAGES_)
	printf("Destruct'ing a t_ppf_V_mhe object.\n");

	delete(this->pf_array);
	delete(this->ext_pf_array);
}

bool t_ppf_V_mhe::check_boundary(int i1, int i2)
{
	// It does not matter which array to use to check boundaries.
	return(this->pf_array->check_boundary(i1, i2));
}

// Access to partition function array by reference.
double& t_ppf_V_mhe::x_setter(int i, int j, int k, int l)
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

double t_ppf_V_mhe::x(int i, int j, int k, int l)
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

double& t_ppf_V_mhe::x_ext(int i, int j, int k, int l)
{
	if(this->pf_array->check_4D_ll(i,j,k,l))
	{
		return(this->ext_pf_array->x(i,j,k,l));
	}
	else
	{
		printf("%s(%d)\n", __FILE__, __LINE__);
		int* p = NULL;
		*p = 0;
		exit(0);
	}
}

void t_ppf_V_mhe::calculate_V_mhe(int i, int j, int k, int l, bool bt) // Arrays to use. 
{
	bool i_inc_k_inc = this->check_boundary(i+1, k+1);
	bool i_dec_k_dec = this->check_boundary(i-1, k-1);
	bool j_dec_l_dec = this->check_boundary(j-1, l-1);
	bool j_dec_l = this->check_boundary(j-1, l);
	bool j_l_dec = this->check_boundary(j, l-1);
	bool i_k = this->check_boundary(i, k);
	bool i_k_dec = this->check_boundary(i, k-1);
	bool i_dec_k = this->check_boundary(i-1, k);

	// alignment energies
	double i_INS_prior = ZERO;
	double j_INS_prior = ZERO;
	double k_INS_prior = ZERO;
	double l_INS_prior = ZERO;
	double i_j_k_l_ALN_prior = ZERO;

	if(i_k_dec)
	{
		i_INS_prior = aln_priors->x(i, k - 1, STATE_INS1);
	}
	if(j_dec_l)
	{
		j_INS_prior = aln_priors->x(j, l, STATE_INS1);
	}
	if(i_dec_k)
	{
		k_INS_prior = aln_priors->x(i - 1, k, STATE_INS2);
	}

	if(j_l_dec)
	{
		l_INS_prior = aln_priors->x(j, l, STATE_INS2);
	}

	if(i_k && j_dec_l_dec)
		i_j_k_l_ALN_prior = MUL(aln_priors->x(i,k, STATE_ALN), aln_priors->x(j,l, STATE_ALN)); // Alignment score for this addition.

	// pairing energies
	double ij_pairing_score = seq1_spf->px(i, j);
	double kl_pairing_score = seq2_spf->px(k, l);
	double i_unpairing_score = seq1_spf->ux_5p(i,j);
	double j_unpairing_score = seq1_spf->ux_3p(i,j);
	double k_unpairing_score = seq2_spf->ux_5p(k,l);
	double l_unpairing_score = seq2_spf->ux_3p(k,l);

	double ij_bpi_score = MUL3(i_INS_prior, j_INS_prior, ij_pairing_score);
	double kl_bpi_score = MUL3(k_INS_prior, l_INS_prior, kl_pairing_score);
	double ij_bau_score = MUL4(i_j_k_l_ALN_prior, ij_pairing_score, k_unpairing_score, l_unpairing_score);
	double kl_bau_score = MUL4(i_j_k_l_ALN_prior, kl_pairing_score, i_unpairing_score, j_unpairing_score);
	double abp_score = MUL3(i_j_k_l_ALN_prior, ij_pairing_score, kl_pairing_score);

	double current_cumulative = ZERO;

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

	// bpi/bau's on top of Wmhi's
	// ij bpi
	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
													MUL(this->ppf_loops->W_mhi->x_kl(i+1, j-1, k, l), ij_bpi_score));

	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i+1, j-1, k, l, TRACE_W_kl_mhi);
		this->ppf_loops->set_seq1_ins(i,k-1);
		this->ppf_loops->set_seq1_ins(j,l);

		this->ppf_loops->add_ct1_bp(i, j);
		pushed = true;
	}

	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
													MUL(this->ppf_loops->V_mhe->x(i+1, j-1, k, l), ij_bpi_score));

	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i+1, j-1, k, l, TRACE_V_mhe);
		this->ppf_loops->set_seq1_ins(i,k-1);
		this->ppf_loops->set_seq1_ins(j,l);

		this->ppf_loops->add_ct1_bp(i, j);
		pushed = true;
	}

	// kl bpi
	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
													MUL(this->ppf_loops->W_mhi->x_ij(i, j, k+1, l-1), kl_bpi_score));

	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i, j, k+1, l-1, TRACE_W_ij_mhi);
		this->ppf_loops->set_seq2_ins(i-1,k);
		this->ppf_loops->set_seq2_ins(j,l);

		this->ppf_loops->add_ct2_bp(k, l);
		pushed = true;
	}

	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
													MUL(this->ppf_loops->V_mhe->x(i, j, k+1, l-1), kl_bpi_score));

	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i, j, k+1, l-1, TRACE_V_mhe);
		this->ppf_loops->set_seq2_ins(i-1,k);
		this->ppf_loops->set_seq2_ins(j,l);

		this->ppf_loops->add_ct2_bp(k, l);
		pushed = true;
	}

	// ij bau
	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
													MUL(this->ppf_loops->W_mhi->x_kl(i+1, j-1, k+1, l-1), ij_bau_score));


	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i+1, j-1, k+1, l-1, TRACE_W_kl_mhi);
		this->ppf_loops->set_aln(i,k);
		this->ppf_loops->set_aln(j,l);

		this->ppf_loops->add_ct1_bp(i, j);
		pushed = true;
	}

	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
													MUL(this->ppf_loops->V_mhe->x(i+1, j-1, k+1, l-1), ij_bau_score));

	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i+1, j-1, k+1, l-1, TRACE_V_mhe);
		this->ppf_loops->set_aln(i,k);
		this->ppf_loops->set_aln(j,l);

		this->ppf_loops->add_ct1_bp(i, j);
		pushed = true;
	}

	// kl bau
	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
													MUL(this->ppf_loops->W_mhi->x_ij(i+1, j-1, k+1, l-1), kl_bau_score));

	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i+1, j-1, k+1, l-1, TRACE_W_ij_mhi);
		this->ppf_loops->set_aln(i,k);
		this->ppf_loops->set_aln(j,l);

		this->ppf_loops->add_ct2_bp(k, l);
		pushed = true;
	}

	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
													MUL(this->ppf_loops->V_mhe->x(i+1, j-1, k+1, l-1), kl_bau_score));

	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i+1, j-1, k+1, l-1, TRACE_V_mhe);
		this->ppf_loops->set_aln(i,k);
		this->ppf_loops->set_aln(j,l);

		this->ppf_loops->add_ct2_bp(k, l);
		pushed = true;
	}

	// abp
	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
													MUL(this->ppf_loops->W_mhi->x_ij(i+1, j-1, k+1, l-1), abp_score));

	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i+1, j-1, k+1, l-1, TRACE_W_ij_mhi);
		this->ppf_loops->set_aln(i,k);
		this->ppf_loops->set_aln(j,l);

		this->ppf_loops->add_ct1_bp(i, j);
		this->ppf_loops->add_ct2_bp(k, l);
		pushed = true;
	}

	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
													MUL(this->ppf_loops->W_mhi->x_kl(i+1, j-1, k+1, l-1), abp_score));

	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i+1, j-1, k+1, l-1, TRACE_W_kl_mhi);
		this->ppf_loops->set_aln(i,k);
		this->ppf_loops->set_aln(j,l);

		this->ppf_loops->add_ct1_bp(i, j);
		this->ppf_loops->add_ct2_bp(k, l);
		pushed = true;
	}

	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
													MUL(this->ppf_loops->V_mhe->x(i+1, j-1, k+1, l-1), abp_score));

	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i+1, j-1, k+1, l-1, TRACE_V_mhe);
		this->ppf_loops->set_aln(i,k);
		this->ppf_loops->set_aln(j,l);

		this->ppf_loops->add_ct1_bp(i, j);
		this->ppf_loops->add_ct2_bp(k, l);
		pushed = true;
	}

	// abp also recurses on non Vmhe.
	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
													MUL(abp_score, this->ppf_loops->SS->x(i+1,j-1,k+1,l-1)));

	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i+1, j-1, k+1, l-1, TRACE_SS);
		this->ppf_loops->set_aln(i,k);
		this->ppf_loops->set_aln(j,l);

		this->ppf_loops->add_ct1_bp(i, j);
		this->ppf_loops->add_ct2_bp(k, l);
		pushed = true;
	}


	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
													MUL(abp_score, this->ppf_loops->W->x(i+1,j-1,k+1,l-1)));

	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i+1, j-1, k+1, l-1, TRACE_W);
		this->ppf_loops->set_aln(i,k);
		this->ppf_loops->set_aln(j,l);

		this->ppf_loops->add_ct1_bp(i, j);
		this->ppf_loops->add_ct2_bp(k, l);
		pushed = true;
	}


	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
													MUL(abp_score, this->ppf_loops->WMB->x(i+1,j-1,k+1,l-1)));

	if(bt && !pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i+1, j-1, k+1, l-1, TRACE_WMB);
		this->ppf_loops->set_aln(i,k);
		this->ppf_loops->set_aln(j,l);

		this->ppf_loops->add_ct1_bp(i, j);
		this->ppf_loops->add_ct2_bp(k, l);
		pushed = true;
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
}

void t_ppf_V_mhe::calculate_ext_dependencies(int i, int j, int k, int l) // Arrays to use. 
{
	bool i_inc_k_inc = this->check_boundary(i+1, k+1);
	bool i_dec_k_dec = this->check_boundary(i-1, k-1);
	bool j_dec_l_dec = this->check_boundary(j-1, l-1);
	bool j_dec_l = this->check_boundary(j-1, l);
	bool j_l_dec = this->check_boundary(j, l-1);
	bool i_k = this->check_boundary(i, k);
	bool i_k_dec = this->check_boundary(i, k-1);
	bool i_dec_k = this->check_boundary(i-1, k);

	// alignment energies
	double i_INS_prior = ZERO;
	double j_INS_prior = ZERO;
	double k_INS_prior = ZERO;
	double l_INS_prior = ZERO;
	double i_j_k_l_ALN_prior = ZERO;

	if(i_k_dec)
	{
		i_INS_prior = aln_priors->x(i, k - 1, STATE_INS1);
	}
	if(j_dec_l)
	{
		j_INS_prior = aln_priors->x(j, l, STATE_INS1);
	}
	if(i_dec_k)
	{
		k_INS_prior = aln_priors->x(i - 1, k, STATE_INS2);
	}

	if(j_l_dec)
	{
		l_INS_prior = aln_priors->x(j, l, STATE_INS2);
	}

	if(i_k && j_dec_l_dec)
		i_j_k_l_ALN_prior = MUL(aln_priors->x(i,k, STATE_ALN), aln_priors->x(j,l, STATE_ALN)); // Alignment score for this addition.

	// pairing energies
	double ij_pairing_score = seq1_spf->px(i, j);
	double kl_pairing_score = seq2_spf->px(k, l);
	double i_unpairing_score = seq1_spf->ux_5p(i,j);
	double j_unpairing_score = seq1_spf->ux_3p(i,j);
	double k_unpairing_score = seq2_spf->ux_5p(k,l);
	double l_unpairing_score = seq2_spf->ux_3p(k,l);

	double ij_bpi_score = MUL3(i_INS_prior, j_INS_prior, ij_pairing_score);
	double kl_bpi_score = MUL3(k_INS_prior, l_INS_prior, kl_pairing_score);
	double ij_bau_score = MUL4(i_j_k_l_ALN_prior, ij_pairing_score, k_unpairing_score, l_unpairing_score);
	double kl_bau_score = MUL4(i_j_k_l_ALN_prior, kl_pairing_score, i_unpairing_score, j_unpairing_score);
	double abp_score = MUL3(i_j_k_l_ALN_prior, ij_pairing_score, kl_pairing_score);

	double current_cumulative = ZERO;

	// bpi/bau's on top of Wmhi's
	// ij bpi
	if(this->ppf_loops->W_mhi->kl_pf_array->check_4D_ll(i+1, j-1, k, l))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
		//												MUL(this->ppf_loops->W_mhi->x_kl(i+1, j-1, k, l), ij_bpi_score));
		this->ppf_loops->W_mhi->x_ext_kl(i+1, j-1, k, l) = this->ppf_loops->MAX_SUM(this->ppf_loops->W_mhi->x_ext_kl(i+1, j-1, k, l), 
			MUL(this->x_ext(i,j,k,l), ij_bpi_score));
	}

	if(this->ppf_loops->V_mhe->pf_array->check_4D_ll(i+1, j-1, k, l))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
		//												MUL(this->ppf_loops->V_mhe->x(i+1, j-1, k, l), ij_bpi_score));
		this->ppf_loops->V_mhe->x_ext(i+1, j-1, k, l) = this->ppf_loops->MAX_SUM(this->ppf_loops->V_mhe->x_ext(i+1, j-1, k, l), 
			MUL(this->x_ext(i,j,k,l), ij_bpi_score));
	}

	// kl bpi
	if(this->ppf_loops->W_mhi->ij_pf_array->check_4D_ll(i, j, k+1, l-1))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
		//												MUL(this->ppf_loops->W_mhi->x_ij(i, j, k+1, l-1), kl_bpi_score));
		this->ppf_loops->W_mhi->x_ext_ij(i, j, k+1, l-1) = this->ppf_loops->MAX_SUM(this->ppf_loops->W_mhi->x_ext_ij(i, j, k+1, l-1), 
			MUL(this->x_ext(i,j,k,l), kl_bpi_score));
	}

	if(this->ppf_loops->V_mhe->pf_array->check_4D_ll(i, j, k+1, l-1))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
		//												MUL(this->ppf_loops->V_mhe->x(i, j, k+1, l-1), kl_bpi_score));
		this->ppf_loops->V_mhe->x_ext(i, j, k+1, l-1) = this->ppf_loops->MAX_SUM(this->ppf_loops->V_mhe->x_ext(i, j, k+1, l-1), 
			MUL(this->x_ext(i,j,k,l), kl_bpi_score));
	}

	// ij bau
	if(this->ppf_loops->W_mhi->kl_pf_array->check_4D_ll(i+1, j-1, k+1, l-1))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
		//												MUL(this->ppf_loops->W_mhi->x_kl(i+1, j-1, k+1, l-1), ij_bau_score));
		this->ppf_loops->W_mhi->x_ext_kl(i+1, j-1, k+1, l-1) = this->ppf_loops->MAX_SUM(this->ppf_loops->W_mhi->x_ext_kl(i+1, j-1, k+1, l-1), 
			MUL(this->x_ext(i,j,k,l), ij_bau_score));
	}

	if(this->ppf_loops->V_mhe->pf_array->check_4D_ll(i+1, j-1, k+1, l-1))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
		//												MUL(this->ppf_loops->V_mhe->x(i+1, j-1, k+1, l-1), ij_bau_score));
		this->ppf_loops->V_mhe->x_ext(i+1, j-1, k+1, l-1) = this->ppf_loops->MAX_SUM(this->ppf_loops->V_mhe->x_ext(i+1, j-1, k+1, l-1), 
			MUL(this->x_ext(i,j,k,l), ij_bau_score));
	}

	// kl bau
	if(this->ppf_loops->W_mhi->ij_pf_array->check_4D_ll(i+1, j-1, k+1, l-1))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
		//												MUL(this->ppf_loops->W_mhi->x_ij(i+1, j-1, k+1, l-1), kl_bau_score));
		this->ppf_loops->W_mhi->x_ext_ij(i+1, j-1, k+1, l-1) = this->ppf_loops->MAX_SUM(this->ppf_loops->W_mhi->x_ext_ij(i+1, j-1, k+1, l-1), 
			MUL(this->x_ext(i,j,k,l), kl_bau_score));
	}

	if(this->ppf_loops->V_mhe->pf_array->check_4D_ll(i+1, j-1, k+1, l-1))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
		//												MUL(this->ppf_loops->V_mhe->x(i+1, j-1, k+1, l-1), kl_bau_score));
		this->ppf_loops->V_mhe->x_ext(i+1, j-1, k+1, l-1) = this->ppf_loops->MAX_SUM(this->ppf_loops->V_mhe->x_ext(i+1, j-1, k+1, l-1), 
			MUL(this->x_ext(i,j,k,l), kl_bau_score));
	}

	// abp
	if(this->ppf_loops->W_mhi->ij_pf_array->check_4D_ll(i+1, j-1, k+1, l-1))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
		//												MUL(this->ppf_loops->W_mhi->x_ij(i+1, j-1, k+1, l-1), abp_score));
		this->ppf_loops->W_mhi->x_ext_ij(i+1, j-1, k+1, l-1) = this->ppf_loops->MAX_SUM(this->ppf_loops->W_mhi->x_ext_ij(i+1, j-1, k+1, l-1), 
			MUL(this->x_ext(i,j,k,l), abp_score));
	}

	if(this->ppf_loops->W_mhi->kl_pf_array->check_4D_ll(i+1, j-1, k+1, l-1))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
		//												MUL(this->ppf_loops->W_mhi->x_kl(i+1, j-1, k+1, l-1), abp_score));
		this->ppf_loops->W_mhi->x_ext_kl(i+1, j-1, k+1, l-1) = this->ppf_loops->MAX_SUM(this->ppf_loops->W_mhi->x_ext_kl(i+1, j-1, k+1, l-1), 
			MUL(this->x_ext(i,j,k,l), abp_score));
	}

	if(this->ppf_loops->V_mhe->pf_array->check_4D_ll(i+1, j-1, k+1, l-1))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
		//												MUL(this->ppf_loops->V_mhe->x(i+1, j-1, k+1, l-1), abp_score));
		this->ppf_loops->V_mhe->x_ext(i+1, j-1, k+1, l-1) = this->ppf_loops->MAX_SUM(this->ppf_loops->V_mhe->x_ext(i+1, j-1, k+1, l-1), 
			MUL(this->x_ext(i,j,k,l), abp_score));
	}

	// abp also recurses on non Vmhe.
	if(this->ppf_loops->SS->check_str_coinc_ll(i+1, j-1, k+1, l-1))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
		//												MUL(abp_score, this->ppf_loops->SS->x(i+1,j-1,k+1,l-1)));
		this->ppf_loops->SS->x_ext(i+1, j-1, k+1, l-1) = this->ppf_loops->MAX_SUM(this->ppf_loops->SS->x_ext(i+1, j-1, k+1, l-1), 
			MUL(this->x_ext(i,j,k,l), abp_score));
	}

	if(this->ppf_loops->W->pf_array->check_4D_ll(i+1, j-1, k+1, l-1))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
		//												MUL(abp_score, this->ppf_loops->W->x(i+1,j-1,k+1,l-1)));
		this->ppf_loops->W->x_ext(i+1, j-1, k+1, l-1) = this->ppf_loops->MAX_SUM(this->ppf_loops->W->x_ext(i+1, j-1, k+1, l-1), 
			MUL(this->x_ext(i,j,k,l), abp_score));
	}

	if(this->ppf_loops->WMB->pf_array->check_4D_ll(i+1, j-1, k+1, l-1))
	{
		//current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, 
		//												MUL(abp_score, this->ppf_loops->WMB->x(i+1,j-1,k+1,l-1)));
		this->ppf_loops->WMB->x_ext(i+1, j-1, k+1, l-1) = this->ppf_loops->MAX_SUM(this->ppf_loops->WMB->x_ext(i+1, j-1, k+1, l-1), 
			MUL(this->x_ext(i,j,k,l), abp_score));
	}
}


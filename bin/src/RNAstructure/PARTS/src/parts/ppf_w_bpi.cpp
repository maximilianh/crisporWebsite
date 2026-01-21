#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "parts_compilation_directives.h"
#include "ppf_w_bpi.h"
#include "ppf_v.h"
#include "ppf_w.h"
#include "ppf_w_mb.h"
#include "ppf_v_bpi.h"
#include "ss_str.h"
#include "process_sequences.h"
#include "ppf_math.h"
#include "alignment_priors.h"
#include "phmm_parameters.h"
#include "template_pf_array.h"
#include "single_pf_array.h"
#include "structural_parameters.h"
#include "ppf_w_bp_aln_up.h"
#include "ppf_v_bp_aln_up.h"
#include "ppf_v_mhe.h"
#include "ppf_w_mhi.h"
#include "ppf_loops.h"
#include <iostream>
#include "loop_limits.h"

using namespace std;
 
#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

// t_ppf_W_bpi: Handles base pair insertions on V array, note that t_ppf_W_bpi inserts base pairs one by one to each sequence
// so it can handle INS1 followed by an INS2 or vice versa, so this is basically different in how insertions are handled
// by consecutive insertion handling with t_ppf_V::calculate_W_helix_ins(...). In the calculations this will go to
// initialization of W array. t_ppf_W_bpi covers structural alignments where i-j and k-l are paired however either one of those
// are inserted. And base pair insertions start on top of aligned base pairs, that is, on top of V array.

extern bool problem;

bool _DUMP_PPF_OS_BPI_MESSAGES_ = false;

t_ppf_W_bpi::t_ppf_W_bpi(t_ppf_loops* _ppf_loops)
{
	this->ppf_loops = _ppf_loops;

	// Copy sequence lengths.
	this->N1 = this->ppf_loops->seq_man->get_l_seq1();
	this->N2 = this->ppf_loops->seq_man->get_l_seq2();

	// Alloc and init two template pf arrays.
	// Those are needed for keeping track of insertion places.
/*
	this->os_ij_bpi_pf_array = new t_template_pf_array(seq_man);
	this->os_kl_bpi_pf_array = new t_template_pf_array(seq_man);
*/
	
	this->seq1_spf = this->ppf_loops->seq1_spf;
	this->seq2_spf = this->ppf_loops->seq2_spf;
	this->aln_priors = this->ppf_loops->aln_priors;
}

t_ppf_W_bpi::~t_ppf_W_bpi()
{
	printf("Destruct'ing a t_ppf_W_bpi object.\n");
}

bool t_ppf_W_bpi::check_boundary(int i1, int i2)
{
	if(t_template_pf_array::low_limits[i1] <= i2 && t_template_pf_array::high_limits[i1] >= i2)
	{
		return(true);
	}
	else
	{
		return(false);
	}
}

// Access to partition function array by reference.
double& t_ppf_W_bpi::x_ij_bpi(int i, int j, int k, int l)
{
	if(this->os_ij_bpi_pf_array->check_4D_ll(i,j,k,l))
	{
		return(this->os_ij_bpi_pf_array->x(i,j,k,l));
	}
	else
	{
		return(this->os_ij_bpi_pf_array->zero);
	}
}

// Access to partition function array by reference.
double& t_ppf_W_bpi::x_kl_bpi(int i, int j, int k, int l)
{
	return(this->os_kl_bpi_pf_array->x(i,j,k,l));
}

// Return SUM of possible base pair insertions.
double t_ppf_W_bpi::x_bpi(int i, int j, int k, int l)
{
//if(_DUMP_PPF_OS_BPI_MESSAGES_)
	//printf("Returning SUM!!!\n");

	return(SUM(this->os_kl_bpi_pf_array->x(i,j,k,l), this->os_ij_bpi_pf_array->x(i,j,k,l)));
}

// Return max of possible base pair insertions.
double t_ppf_W_bpi::x_map_bpi(int i, int j, int k, int l)
{
	return(max(this->os_kl_bpi_pf_array->x(i,j,k,l), this->os_ij_bpi_pf_array->x(i,j,k,l)));
}

double t_ppf_W_bpi::calculate_os_ij_bpi(int i, int j, int k, int l) // Arrays to use.
{
	// Have to fork for j > N, do bound checking...
	if( !this->ppf_loops->W_mhi->w_ij_mhi_pf_array->check_4D_ll(i, j, k, l) )
	{
		printf("Returning from w_bpi calculation @ %s(%d) before calculating W_bi(%d, %d, %d, %d)\n", __FILE__, __LINE__, i,j,k,l);
		return(ZERO);
	}

if(_DUMP_PPF_OS_BPI_MESSAGES_)
	printf("\nw_bpi(%d, %d, %d, %d):\n", i,j,k,l);

	// Following checks are for making sure that insertions are valid.
	bool i_inc_k_inc = this->check_boundary(i + 1, k + 1);
	bool i_inc_k = this->check_boundary(i + 1, k);
	bool i_k_inc = this->check_boundary(i, k + 1);
	bool i_dec_k_dec = this->check_boundary(i - 1, k - 1);
	bool i_dec_k = this->check_boundary(i - 1, k);
	bool i_k_dec = this->check_boundary(i, k - 1);

	bool j_dec_l_dec = this->check_boundary(j - 1 , l - 1);
	bool j_l_dec = this->check_boundary(j, l - 1);
	bool j_dec_l = this->check_boundary(j - 1 , l);

	// Insert i-j base pair on top of (W+WMB+os_ij_bpi+ss_str)(i+1, j-1, k, l)
	// Note that we are taking k-l as an open structure, which is for ss_str, an open ss. 
	double i_j_bp_ins_score = CONVERT_FROM_LIN(0.0);
	double i_INS_prior = CONVERT_FROM_LIN(0.0);
	double j_INS_prior = CONVERT_FROM_LIN(0.0);
	double i_j_pairing_score = CONVERT_FROM_LIN(0.0);

	if((j <= N1 && (j-i) > MIN_LOOP)) // This check makes sure that MIN_LOOP constraint is applied correctly.
		i_j_pairing_score = seq1_spf->px(i, j);

	// An i insertion implies a transition from alignment position i-1, k-1 to alignment position i, k-1.
	//if(i_dec_k_dec && i_k_dec)
	if(this->ppf_loops->W->w_pf_array->check_4D_ll(i+1, j-1, k, l))
	{
		i_INS_prior = aln_priors->x(i, k - 1, STATE_INS1);
		j_INS_prior = aln_priors->x(j, l, STATE_INS1);
	}

	i_j_bp_ins_score = MUL(i_j_pairing_score, MUL(i_INS_prior, j_INS_prior));

	double int_ppf_W_WMB_os_score = CONVERT_FROM_LIN(0.0);

	if(((j-i) > MIN_LOOP && j <= N1))
	{
		//double closed_W_score = SUM(V->x(i+1, j-1, k, l), SUM(V_bpi->x_bpi(i+1, j-1, k, l), v_bp_aln_up->x(i+1, j-1, k, l)));
		//double closed_W_score = SUM(V->x(i+1, j-1, k, l), V_bpi->x_bpi(i+1, j-1, k, l));
		double closed_W_score = this->ppf_loops->MAX_SUM(this->ppf_loops->V->x(i+1, j-1, k, l), this->ppf_loops->V_mhe->x_mhe(i+1, j-1, k, l));
		double open_W_score = SUB(this->ppf_loops->W->x(i+1, j-1, k, l), closed_W_score);
		int_ppf_W_WMB_os_score = this->ppf_loops->MAX_SUM(open_W_score, this->ppf_loops->WMB->x(i+1, j-1, k, l));

		//int_ppf_W_WMB_os_score = SUM(int_ppf_W_WMB_os_score, this->x_ij_bpi(i+1, j-1, k, l));
		//int_ppf_W_WMB_os_score = SUM(int_ppf_W_WMB_os_score, w_bp_aln_up->x_ij_aln_up(i+1, j-1, k,l));
		int_ppf_W_WMB_os_score = SUM(int_ppf_W_WMB_os_score, this->ppf_loops->W_mhi->x_ij_mhi(i+1, j-1, k, l));
	}

	// Must take ij exclusive of scores since ij is emitted by base pair insertion.
	double int_ss_score = this->ppf_loops->SS->x(i+1, j-1, k, l); 

	//this->x_ij_bpi(i, j, k, l) = MUL(i_j_bp_ins_score, SUM(int_ppf_W_WMB_os_score, int_ss_score));
	return( MUL(i_j_bp_ins_score, SUM(int_ppf_W_WMB_os_score, int_ss_score)) );
}

double t_ppf_W_bpi::calculate_os_kl_bpi(int i, int j, int k, int l) // Arrays to use.
{
	// Have to fork for j > N, do bound checking...
	if( !this->ppf_loops->W_mhi->w_kl_mhi_pf_array->check_4D_ll(i, j, k, l) )
	{
		//printf("Returning from w_bpi calculation @ %s(%d) before calculating W_bi(%d, %d, %d, %d)\n", __FILE__, __LINE__, i,j,k,l);
		return(ZERO);
	}

if(_DUMP_PPF_OS_BPI_MESSAGES_)
	printf("\nw_bpi(%d, %d, %d, %d):\n", i,j,k,l);

	// Following checks are for making sure that insertions are valid.
	bool i_inc_k_inc = this->check_boundary(i + 1, k + 1);
	bool i_inc_k = this->check_boundary(i + 1, k);
	bool i_k_inc = this->check_boundary(i, k + 1);
	bool i_dec_k_dec = this->check_boundary(i - 1, k - 1);
	bool i_dec_k = this->check_boundary(i - 1, k);
	bool i_k_dec = this->check_boundary(i, k - 1);

	bool j_dec_l_dec = this->check_boundary(j - 1 , l - 1);
	bool j_l_dec = this->check_boundary(j, l - 1);
	bool j_dec_l = this->check_boundary(j - 1 , l);
	// Insert k-l base pair on top of V(i, j, k+1, l-1) and w_bpi(i, j, k+1, l-1).
	double k_l_bp_ins_score = CONVERT_FROM_LIN(0.0);
	double k_INS_prior = CONVERT_FROM_LIN(0.0);
	double l_INS_prior = CONVERT_FROM_LIN(0.0);
	double k_l_pairing_score = CONVERT_FROM_LIN(0.0);

	if((l-k) > MIN_LOOP && l <= N2) // This check makes sure that MIN_LOOP constraint is applied correctly.
		k_l_pairing_score = seq2_spf->px(k, l);

	// k insertion implies a transition from (i-1, k-1) into (i-1, k).
	if(this->ppf_loops->W->w_pf_array->check_4D_ll(i, j, k+1, l-1))
	{
		k_INS_prior = aln_priors->x(i-1, k, STATE_INS2);
		l_INS_prior = aln_priors->x(j, l, STATE_INS2);
	}

	k_l_bp_ins_score = MUL(k_l_pairing_score, MUL(k_INS_prior, l_INS_prior));

	double int_ppf_W_WMB_os_score = CONVERT_FROM_LIN(0.0);

	// Recurse on i+1, j-1, k, l of w_bpi, W and WMB for inserting i-j pair.
	if(((k > (l-N2) && l > N2) || ((l-k) > MIN_LOOP && l <= N2)) 
		&& k != N2 && l != N2+1)
	{
		double closed_W_score = this->ppf_loops->MAX_SUM(this->ppf_loops->V->x(i, j, k+1, l-1), this->ppf_loops->V_mhe->x_mhe(i, j, k+1, l-1));

		double open_W_score = SUB(this->ppf_loops->W->x(i, j, k+1, l-1), closed_W_score);
		int_ppf_W_WMB_os_score = SUM(open_W_score, this->ppf_loops->WMB->x(i, j, k+1, l-1));

		int_ppf_W_WMB_os_score = SUM(int_ppf_W_WMB_os_score, this->ppf_loops->W_mhi->x_kl_mhi(i, j, k+1, l-1));

if(_DUMP_PPF_OS_BPI_MESSAGES_)
{
		printf("open_W_score = %f - %f = %f\n", this->ppf_loops->W->x(i, j, k+1, l-1), closed_W_score, open_W_score);
		printf("int_ppf_W_WMB_os_score = %f + %f + %f = %f\n", open_W_score, this->ppf_loops->WMB->x(i, j, k+1, l-1), this->ppf_loops->W_mhi->x_kl_mhi(i, j, k+1, l-1), int_ppf_W_WMB_os_score);
}
	}

	// Determine ss containing structure scores.
	double int_ss_score = this->ppf_loops->SS->x(i, j, k+1, l-1); // Must take inclusive ss scores.

	return( MUL(k_l_bp_ins_score, SUM(int_ppf_W_WMB_os_score, int_ss_score)) );
}

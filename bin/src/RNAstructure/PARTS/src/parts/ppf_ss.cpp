#include <string.h>
#include <limits.h>
#include "parts_compilation_directives.h"
#include "ppf_ss.h"
#include "ppf_cli.h"
#include "ppf_w.h"

#include "ppf_w_mb.h"
#include <stdio.h>
#include <stdlib.h>
#include "../../../src/phmm/structure/folding_constraints.h"
#include "template_pf_array.h"
#include "process_sequences.h"
#include "ppf_math.h"
#include "ppf_scale.h"
#include "single_pf_array.h"
#include "alignment_priors.h"
#include "phmm_parameters.h"
#include "ppf_loops.h"
#include "ppf_progress_bar.h"


#include "ppf_tb_stack.h"
#include "map_structures.h"
#include "map_alignment.h"
#include "stoch_tb/stoch_sampling_math.h"
#include "ppf_cli.h"
// #include "../../../src/phmm/phmm.h"
#include "../../../src/phmm/structure/folding_constraints.h"

#include <limits.h>

#include <algorithm>
#define MIN(x,y) std::min(x,y)
#define MAX(x,y) std::max(x,y)

bool _DUMP_PPF_SS_MESSAGES_ = false;

// C'tor.
//t_ppf_SS::t_ppf_SS(t_seq_man* seq_man, t_ppf_W* W, t_ppf_WMB* WMB, t_frag_aln_enum_array* frag_aln_enumer)
t_ppf_SS::t_ppf_SS(t_ppf_loops* _ppf_loops) // Constructor.
{
	this->ppf_loops = _ppf_loops;

	// Copy sequence lengths.
	this->N1 = ppf_loops->seq_man->get_l_seq1();
	this->N2 = ppf_loops->seq_man->get_l_seq2();

	// Copy alignment and folding priors.
	this->seq1_spf = this->ppf_loops->seq1_spf;
	this->seq2_spf = this->ppf_loops->seq2_spf;
	this->aln_priors = this->ppf_loops->aln_priors;

	this->seq1_ptr_reloc_map = seq1_spf->folding_constraints->coinc_pointer_relocation_map;
	this->seq2_ptr_reloc_map = seq2_spf->folding_constraints->coinc_pointer_relocation_map;

	this->seq_man = this->ppf_loops->seq_man;

	this->alloc_init_pf_array();
	this->alloc_init_ext_pf_array();
}

void t_ppf_SS::alloc_init_pf_array()
{
if(_DUMP_PPF_SS_MESSAGES_)
	printf("Initing SS array.\n");

	this->n_alloced_bytes = 0.0f;

	// Allocate pf array.
	this->pf_array = (double****)malloc(sizeof(double***) * (this->ppf_loops->seq_man->get_l_seq1() + 3));

	this->n_alloced_bytes += sizeof(double***) * (this->ppf_loops->seq_man->get_l_seq1() + 3);

	int high_j = this->ppf_loops->seq_man->get_l_seq1();

	t_ppf_progress_bar* tpa_alloc_progress_bar = new t_ppf_progress_bar(NULL, '=', true, high_j);

	for(int i = 1; i <= this->N1; i++)
	{
		//printf("i = %d\n", i);
		int min_j = i-1;
		int max_j = this->N1;

if(_DUMP_PPF_SS_MESSAGES_)
		printf("Computing relocation map for i=%d:\n", i);

		this->pf_array[i] = (double***)malloc( sizeof(double**) * (max_j - min_j + 1) );
		this->pf_array[i] -= min_j;

		for(int j = min_j; j <= max_j; j++)
		{
			// Check if the remaining indices are allocateable.
			if(i == j+1 || 
				this->seq1_ptr_reloc_map[i][j] != POS_MEM_NOT_EXIST)
			{
				// Check if the remaining indices are allocateable.
				// Allocate kxl plane using loop limits for v.
				int max_k = t_template_pf_array::high_limits[i-1]+1;
				int min_k = MAX(1, t_template_pf_array::low_limits[i-1]+1);

				if((max_k - min_k + 1) > 0)
				{
					this->pf_array[i][j] = (double**)malloc( sizeof(double*) * (max_k - min_k + 1) );
					this->pf_array[i][j] -= min_k;
				}
				
				if((max_k - min_k + 1) > 0)
				{
					this->n_alloced_bytes += ( sizeof(double**) * (max_k - min_k + 1) );
				}

				for(int k = min_k; k <= max_k; k++)
				{
					//printf("k = %d\n", k);

					// Also use the 
					int min_l = MAX(k-1, t_template_pf_array::low_limits[j]);
					int max_l = MIN(N2, t_template_pf_array::high_limits[j]);

					this->pf_array[i][j][k] = (double*)malloc( sizeof(double) * (max_l - min_l + 1) );
					this->pf_array[i][j][k] -= min_l;
					
					for(int l = min_l; l <= max_l; l++)
					{
						this->pf_array[i][j][k][l] = CONVERT_FROM_LIN(0.0f);
						//printf("Initing %d, %d, %d, %d\n", i,j,k,l);
					} // l loop
				} // k loop
			} // str coinc check for i-j.
		} // j loop
	} // i loop.

	delete(tpa_alloc_progress_bar);

if(_DUMP_PPF_SS_MESSAGES_)
	printf("SS array allocated %lf bytes\n", this->n_alloced_bytes);

	// Initialize pf_array[i][i-1][k][l] and pf_array[i][j][k][k-1].
	for(int i = 1; i <= this->N1; i++)
	{
		//printf("i = %d\n", i);
		int min_j = i-1;
		int max_j = this->N1;

		for(int j = min_j; j <= max_j; j++)
		{
			int max_k = t_template_pf_array::high_limits[i-1]+1;
			int min_k = MAX(1, t_template_pf_array::low_limits[i-1]+1);

			for(int k = min_k; k <= max_k; k++)
			{
				int min_l = MAX(k-1, t_template_pf_array::low_limits[j]);
				int max_l = MIN(N2, t_template_pf_array::high_limits[j]);

				for(int l = min_l; l <= max_l; l++)
				{
					// Check the validity of the indices first.
					if(this->check_str_coinc_ll(i,j,k,l))
					{
						// Check for indices to be initialized.
						if(j == i-1 && l == k-1)
						{
							this->x_setter(i,j,k,l) = CONVERT_FROM_LIN(1.0f);
						}
						else if(j == i-1)
						{
							this->x_setter(i,j,k,l) = MUL3(this->aln_priors->x(j, l, STATE_INS2), this->x(i,j,k,l-1), this->seq2_spf->ux_3p(k, l));
						}
						else if(l == k-1)
						{
							this->x_setter(i,j,k,l) = MUL3(this->aln_priors->x(j, l, STATE_INS1), this->x(i,j-1,k,l), this->seq1_spf->ux_3p(i, j));
						}
						else
						{
							this->x_setter(i,j,k,l) = CONVERT_FROM_LIN(0.0f);
						}					

						if(_DUMP_PPF_SS_MESSAGES_ &&
							this->x(i,j,k,l) != ZERO)
						{
							printf("SS[Init](%d, %d, %d, %d) = %.10f\n", i,j,k,l, this->x(i,j,k,l));
						}
					}
				} // l loop
			} // k loop
		} // j loop
	} // i loop.
}

void t_ppf_SS::alloc_init_ext_pf_array()
{
if(_DUMP_PPF_SS_MESSAGES_)
	printf("Initing SS array.\n");

	this->n_alloced_bytes = 0.0f;

	// Allocate pf array.
	this->ext_pf_array = (double****)malloc(sizeof(double***) * (this->ppf_loops->seq_man->get_l_seq1() + 3));

	this->n_alloced_bytes += sizeof(double***) * (this->ppf_loops->seq_man->get_l_seq1() + 3);

	int high_j = this->ppf_loops->seq_man->get_l_seq1();

	t_ppf_progress_bar* tpa_alloc_progress_bar = new t_ppf_progress_bar(NULL, '=', true, high_j);

	for(int i = 1; i <= this->N1; i++)
	{
		//printf("i = %d\n", i);
		int min_j = i-1;
		int max_j = this->N1;

if(_DUMP_PPF_SS_MESSAGES_)
		printf("Computing relocation map for i=%d:\n", i);

		this->ext_pf_array[i] = (double***)malloc( sizeof(double**) * (max_j - min_j + 1) );
		this->ext_pf_array[i] -= min_j;

		for(int j = min_j; j <= max_j; j++)
		{
			// Check if the remaining indices are allocateable.
			if(i == j+1 || 
				this->seq1_ptr_reloc_map[i][j] != POS_MEM_NOT_EXIST)
			{
				// Allocate kxl plane using loop limits for v.
				int max_k = t_template_pf_array::high_limits[i-1]+1;
				int min_k = MAX(1, t_template_pf_array::low_limits[i-1]+1);

				if((max_k - min_k + 1) > 0)
				{
					this->ext_pf_array[i][j] = (double**)malloc( sizeof(double*) * (max_k - min_k + 1) );
					this->ext_pf_array[i][j] -= min_k;
				}
				
				if((max_k - min_k + 1) > 0)
				{
					this->n_alloced_bytes += ( sizeof(double**) * (max_k - min_k + 1) );
				}

				for(int k = min_k; k <= max_k; k++)
				{
					//printf("k = %d\n", k);

					// Also use the 
					int min_l = MAX(k-1, t_template_pf_array::low_limits[j]);
					int max_l = MIN(N2, t_template_pf_array::high_limits[j]);

					this->ext_pf_array[i][j][k] = (double*)malloc( sizeof(double) * (max_l - min_l + 1) );
					this->ext_pf_array[i][j][k] -= min_l;
					
					for(int l = min_l; l <= max_l; l++)
					{
						this->ext_pf_array[i][j][k][l] = CONVERT_FROM_LIN(0.0f);
						//printf("Initing %d, %d, %d, %d\n", i,j,k,l);
					} // l loop
				} // k loop
			} // seq1_str_coinc check.
		} // j loop
	} // i loop.

	delete(tpa_alloc_progress_bar);

if(_DUMP_PPF_SS_MESSAGES_)
	printf("SS[ext] array allocated %lf bytes\n", this->n_alloced_bytes);
}

bool t_ppf_SS::check_str_coinc_ll(int i, int j, int k, int l)
{
	bool ll_check = (this->check_boundary(i-1, k-1)) && (this->check_boundary(j,l));
	bool i_limits_check = (i >= 1 && k >= 1 && (j >= i-1) && (l >= k-1));
	bool str_coinc_check = (i == j+1 || this->seq1_ptr_reloc_map[i][j] != POS_MEM_NOT_EXIST) && (k == l+1 || this->seq2_ptr_reloc_map[k][l] != POS_MEM_NOT_EXIST);

	//if(i == 1 &&
	//	j == 10 &&
	//	k == 1 &&
	//	l == 9)
	//{
	//	printf("(%d, %d, %d, %d): ll_check = %d\ni_limits_check=%d\nstr_coinc_check=%d\n", i,j,k,l,ll_check, i_limits_check, str_coinc_check);
	//}

	return(ll_check && i_limits_check && str_coinc_check);
}

bool t_ppf_SS::check_boundary(int i1, int i2)
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

double& t_ppf_SS::x_setter(int i, int j, int k, int l)
{
	if(this->check_str_coinc_ll(i,j,k,l))
	{
		return(this->pf_array[i][j][k][l]);
	}
	else
	{
		printf("Referred to inaccessible index @ %s(%d)\n", __FILE__, __LINE__);
		exit(0);
	}
}

double t_ppf_SS::x(int i, int j, int k, int l)
{
	if(this->check_str_coinc_ll(i,j,k,l))
	{
		return(this->pf_array[i][j][k][l]);
	}
	else
	{
		return(ZERO);
	}
}

double& t_ppf_SS::x_ext(int i, int j, int k, int l)
{
	if(this->check_str_coinc_ll(i,j,k,l))
	{
		return(this->ext_pf_array[i][j][k][l]);
	}
	else
	{
		printf("%s(%d)\n", __FILE__, __LINE__);
		exit(0);
	}
}

void t_ppf_SS::compute() // Constructor.
{
	// Initialize pf_array[i][i-1][k][l] and pf_array[i][j][k][k-1].
	for(int i = 1; i <= this->N1; i++)
	{
		//printf("i = %d\n", i);
		int min_j = i;
		int max_j = this->N1;

if(_DUMP_PPF_SS_MESSAGES_)
		printf("Computing relocation map for i=%d:\n", i);

		for(int j = min_j; j <= max_j; j++)
		{
			int max_k = t_template_pf_array::high_limits[i-1]+1;
			int min_k = MAX(1, t_template_pf_array::low_limits[i-1]+1);

			for(int k = min_k; k <= max_k; k++)
			{
				int min_l = MAX(k, t_template_pf_array::low_limits[j]);
				int max_l = MIN(N2, t_template_pf_array::high_limits[j]);

				for(int l = min_l; l <= max_l; l++)
				{		
					//getc(stdin);
					//if(!this->check_str_coinc_ll(1,10,1,9))
					//{
					//	printf("Cannot access 1,5,1,4\n");
					//	getc(stdin);
					//}
					//else
					//{
					//	printf("Can access!\n");
					//	getc(stdin);
					//}

					// Check the validity of the indices first.
					if(this->check_str_coinc_ll(i,j,k,l))
					{
						//this->x_setter(i,j,k,l) = ZERO;

						//if(i == 1 &&
						//	j == 5 &&
						//	k == 1 &&
						//	l == 4)
						//{
						//	printf("Computing 1,5,1,4\n");
						//	_DUMP_PPF_SS_MESSAGES_ = true;
						//}

						// Check for indices to be initialized.
if(_DUMP_PPF_SS_MESSAGES_)
						printf("---------------------------------------------\nSS[Compute](%d, %d, %d, %d)\n", i,j,k,l, this->x(i,j,k,l));

						// Align j-l.
						double j_unpair_score = this->seq1_spf->ux_3p(i, j); // Inserted nuc is s2+1.
						double l_unpair_score = this->seq2_spf->ux_3p(k, l); // Inserted nuc is s2+1.
						double ALN_rec_score = MUL4(this->aln_priors->x(j, l, STATE_ALN), 
													j_unpair_score, 
													l_unpair_score, 
													this->x(i, j-1, k, l-1));

						this->x_setter(i,j,k,l) = this->ppf_loops->MAX_SUM(this->x(i,j,k,l), ALN_rec_score);

						if(_DUMP_PPF_SS_MESSAGES_ &&
							this->x(i,j,k,l) != ZERO)
						{
							printf("SS[Compute](%d, %d, %d, %d): ALN j-l: %.5f x %.5f x %.5f x %.5f = %.5f (%.10f)\n", i,j,k,l, 
																														this->aln_priors->x(j, l, STATE_ALN), 
																														j_unpair_score, 
																														l_unpair_score, 
																														this->x(i, j-1, k, l-1),
																														ALN_rec_score,
																														this->x(i,j,k,l));
						}

						// ins j.
						double INS1_rec_score = MUL3(this->aln_priors->x(j, l, STATE_INS1), j_unpair_score, this->x(i,j-1, k,l));
						this->x_setter(i,j,k,l) = this->ppf_loops->MAX_SUM(this->x(i,j,k,l), INS1_rec_score);

						if(_DUMP_PPF_SS_MESSAGES_ &&
							this->x(i,j,k,l) != ZERO)
						{
							printf("SS[Compute](%d, %d, %d, %d): INS j: %.5f x %.5f x %.5f = %.5f (%.10f)\n", i,j,k,l, 
																														this->aln_priors->x(j, l, STATE_INS1), 
																														j_unpair_score, 
																														this->x(i,j-1, k,l),
																														INS1_rec_score,
																														this->x(i,j,k,l));
						}

						// ins l.
						double INS2_rec_score = MUL3(this->aln_priors->x(j, l, STATE_INS2), l_unpair_score, this->x(i,j,k,l-1));
						this->x_setter(i,j,k,l) = this->ppf_loops->MAX_SUM(this->x(i,j,k,l), INS2_rec_score);

						if(_DUMP_PPF_SS_MESSAGES_ &&
							this->x(i,j,k,l) != ZERO)
						{
							printf("SS[Compute](%d, %d, %d, %d): INS l: %.5f x %.5f x %.5f = %.5f (%.10f)\n", i,j,k,l, 
																												this->aln_priors->x(j, l, STATE_INS2), 
																												l_unpair_score, 
																												this->x(i,j,k,l-1),
																												INS2_rec_score,
																												this->x(i,j,k,l));
						}

						if(_DUMP_PPF_SS_MESSAGES_ &&
							this->x(i,j,k,l) != ZERO)
						{
							printf("SS[Compute](%d, %d, %d, %d) = %.10f\n", i,j,k,l, this->x(i,j,k,l));
						}

						if(i == 1 &&
							j == 5 &&
							k == 1 &&
							l == 4)
						{
							_DUMP_PPF_SS_MESSAGES_ = false;
						}
					}
				} // l loop
			} // k loop
		} // j loop
	} // i loop.
}

void t_ppf_SS::rescale(bool up_scale) // Constructor.
{
	// Initialize pf_array[i][i-1][k][l] and pf_array[i][j][k][k-1].
	for(int i = 1; i <= this->N1; i++)
	{
		//printf("i = %d\n", i);
		int min_j = i-1;
		int max_j = this->N1;		

		for(int j = min_j; j <= max_j; j++)
		{
			int max_k = t_template_pf_array::high_limits[i-1]+1;
			int min_k = MAX(1, t_template_pf_array::low_limits[i-1]+1);

			for(int k = min_k; k <= max_k; k++)
			{
				int min_l = MAX(k-1, t_template_pf_array::low_limits[j]);
				int max_l = MIN(N2, t_template_pf_array::high_limits[j]);

				for(int l = min_l; l <= max_l; l++)
				{
					// Check the validity of the indices first.
					if(this->check_str_coinc_ll(i,j,k,l))
					{
						//printf("SS rescale(%d, %d, %d, %d)\n", i,j,k,l);

						// Check for indices to be initialized.
						double factor_de_rescale = this->ppf_loops->ppf_scaler->cumulative_rescale_factor(i,j,k,l);

						if(up_scale)
						{
							this->x_setter(i,j,k,l) *= factor_de_rescale;
						}
						else
						{
							this->x_setter(i,j,k,l) /= factor_de_rescale;
						}
					}
				} // l loop
			} // k loop
		} // j loop
	} // i loop.
}

void t_ppf_SS::backtrack(int i, int j, int k, int l)
{
	// There is nothing left to backtrack.
	if(i > j && k > l)
	{
		return;
	}

	double random_cumulative = ZERO;
	double current_cumulative = ZERO;
	bool pushed = false;

	// Setup cumulatives.
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

	// Check for indices to be initialized.
if(_DUMP_PPF_SS_MESSAGES_)
	printf("---------------------------------------------\nSS[Compute](%d, %d, %d, %d)\n", i,j,k,l, this->x(i,j,k,l));

	// Align j-l.
	double j_unpair_score = this->seq1_spf->ux_3p(i, j); // Inserted nuc is s2+1.
	double l_unpair_score = this->seq2_spf->ux_3p(k, l); // Inserted nuc is s2+1.
	double ALN_rec_score = MUL4(this->aln_priors->x(j, l, STATE_ALN), 
								j_unpair_score, 
								l_unpair_score, 
								this->x(i, j-1, k, l-1));

	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, ALN_rec_score);
	if(!pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i, j-1, k, l-1, TRACE_SS);
		this->ppf_loops->set_aln(j,l);
		pushed = true;
	}

	// ins j.
	double INS1_rec_score = MUL3(this->aln_priors->x(j, l, STATE_INS1), j_unpair_score, this->x(i,j-1, k,l));
	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, INS1_rec_score);
	if(!pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i, j-1, k, l, TRACE_SS);
		this->ppf_loops->set_seq1_ins(j,l);
		pushed = true;
	}

	// ins l.
	double INS2_rec_score = MUL3(this->aln_priors->x(j, l, STATE_INS2), l_unpair_score, this->x(i,j,k,l-1));
	current_cumulative = this->ppf_loops->MAX_SUM(current_cumulative, INS2_rec_score);
	if(!pushed && this->ppf_loops->BT_OP(current_cumulative, random_cumulative))
	{
		this->ppf_loops->tb_stack->push_str(i, j, k, l-1, TRACE_SS);
		this->ppf_loops->set_seq2_ins(j,l);
		pushed = true;
	}

	if(!pushed ||
		!COMPARE(this->x(i,j,k,l), current_cumulative))
	{
		printf("Traceback error @ %s(%d)\n", __FILE__, __LINE__);
		exit(0);
	}
	this->ppf_loops->ppf_scaler->get_unscaled_log_value_by_indices(i, j, k, l, this->x(i,j,k,l));
}

t_ppf_SS::~t_ppf_SS()
{
if(_DUMP_PPF_SS_MESSAGES_)
	printf("Initing SS array.\n");

	//this->n_alloced_bytes = 0.0f;

	// Allocate pf array.
	//this->pf_array = (double****)malloc(sizeof(double***) * (this->ppf_loops->seq_man->get_l_seq1() + 3));
	//this->n_alloced_bytes += sizeof(double***) * (this->ppf_loops->seq_man->get_l_seq1() + 3);

	int high_j = this->ppf_loops->seq_man->get_l_seq1();

	for(int i = 1; i <= this->N1; i++)
	{
		//printf("i = %d\n", i);
		int min_j = i-1;
		int max_j = this->N1;

if(_DUMP_PPF_SS_MESSAGES_)
		printf("Computing relocation map for i=%d:\n", i);

		//this->pf_array[i] = (double***)malloc( sizeof(double**) * (max_j - min_j + 1) );
		//this->pf_array[i] -= min_j;

		for(int j = min_j; j <= max_j; j++)
		{
			// Check if the remaining indices are allocateable.
			if(i == j + 1 || 
				this->seq1_ptr_reloc_map[i][j] != POS_MEM_NOT_EXIST)
			{
				// Allocate kxl plane using loop limits for v.
				int max_k = t_template_pf_array::high_limits[i-1]+1;
				int min_k = MAX(1, t_template_pf_array::low_limits[i-1]+1);

				//if((max_k - min_k + 1) > 0)
				//{
				//	this->pf_array[i][j] = (double**)malloc( sizeof(double*) * (max_k - min_k + 1) );
				//	this->pf_array[i][j] -= min_k;
				//}
				
				//if((max_k - min_k + 1) > 0)
				//{
				//	this->n_alloced_bytes += ( sizeof(double**) * (max_k - min_k + 1) );
				//}

				for(int k = min_k; k <= max_k; k++)
				{
					//printf("k = %d\n", k);

					// Also use the 
					int min_l = MAX(k-1, t_template_pf_array::low_limits[j]);
					int max_l = MIN(N2, t_template_pf_array::high_limits[j]);

					this->pf_array[i][j][k] += min_l;
					free(this->pf_array[i][j][k]);

					this->ext_pf_array[i][j][k] += min_l;
					free(this->ext_pf_array[i][j][k]);
					
					//for(int l = min_l; l <= max_l; l++)
					//{
					//	this->pf_array[i][j][k][l] = CONVERT_FROM_LIN(0.0f);
					//	//printf("Initing %d, %d, %d, %d\n", i,j,k,l);
					//} // l loop
				} // k loop

				if((max_k - min_k + 1) > 0)
				{
					this->pf_array[i][j] += min_k;
					free(this->pf_array[i][j]);

					this->ext_pf_array[i][j] += min_k;
					free(this->ext_pf_array[i][j]);
				}
			} // seq1_str_coinc check.
		} // j loop

		this->pf_array[i] += min_j;
		free(this->pf_array[i]);

		this->ext_pf_array[i] += min_j;
		free(this->ext_pf_array[i]);
	} // i loop.

	free(this->pf_array);
	free(this->ext_pf_array);
}



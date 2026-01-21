#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ppf_scale.h"
#include "process_sequences.h"
#include "single_pf_array.h"
#include "ppf_math.h"
#include "ppf_cli.h"

#include "math.h"
#include "template_pf_array.h"

// Include array headers.
#include "ppf_w_l.h"
#include "ppf_w.h"
#include "ppf_w_mbl.h"
#include "ppf_w_mb.h"
#include "ppf_w_mhi.h"
#include "ppf_v_mhe.h"
#include "ppf_loops.h"
#include "ppf_w_ext.h"

#include "ppf_math.h"

#include <iostream>

using namespace std;

/*
class t_ppf_scale
{
public:
	t_seq_man* seq_man;
	double rescale_factor; // Initial rescale factor.
	double* seq1_scales; // Scaling factors of each nucleotide in 1st sequence.
	double* seq2_scales; // Scaling factors of each nucleotide in 2nd sequence.

	t_ppf_scale(t_seq_man* seq_man);
	t_ppf_scale(t_seq_man* seq_man, double init_rescale_factor);

	// Rescale check function, calls each array to rescale if going out of precision for
	// one of the arrays.
	void check_rescale(int i, int j, int k, int l, double value);

	// Read scaler values by nucleotide index and 
	// i1 is nucleotide index of sequence 1 and i2 is ...
	double get_scale_value(int i1, int i2); 

	~t_ppf_scale();
}

#endif // _PPF_SCALE_
*/

/*
ALL CALCULATIONS ARE LINEAR!!!
*/

bool _DUMP_PPF_SCALE_MESSAGES_ = false;

//t_ppf_scale::t_ppf_scale(t_seq_man* seq_man, t_ppf_cli* ppf_cli, t_spf_array* this->ppf_loops->seq1_spf, t_spf_array* this->ppf_loops->seq2_spf)
t_ppf_scale::t_ppf_scale(t_ppf_loops* _ppf_loops)
{
	this->ppf_loops = _ppf_loops;
	this->seq_man = ppf_loops->seq_man; 
	this->ppf_cli = ppf_loops->ppf_cli; // Need this just for panicking, bad design.

	// Allocate scales.
	this->seq1_scales = (double*)malloc((seq_man->get_l_seq1() + 2) * sizeof(double));
	this->seq2_scales = (double*)malloc((seq_man->get_l_seq2() + 2) * sizeof(double));

	double seq1_cum_scale = 1.0;
	double seq2_cum_scale = 1.0;

	// Calculate scale of each nucleotide.
	// Use (score_pair(i) * score_unpair(i))^(-.5)
	for(int cnt1 = 1; cnt1 <= seq_man->get_l_seq1(); cnt1++)
	{
		this->seq1_scales[cnt1] = sqrt(this->ppf_loops->seq1_spf->ind_pairing_array[cnt1] * this->ppf_loops->seq1_spf->ind_unpairing_array[cnt1]);
        //if( COMPARE(this->seq1_scales[cnt1], 0.0))
		if( COMPARE(this->ppf_loops->seq1_spf->ind_pairing_array[cnt1], 0.0) || COMPARE(this->ppf_loops->seq1_spf->ind_pairing_array[cnt1], 1.0))
        {
			this->seq1_scales[cnt1] = 1.0;
        }
        else
        {
			this->seq1_scales[cnt1] = 1 / this->seq1_scales[cnt1];
		}

		seq1_cum_scale *= this->seq1_scales[cnt1];
if(_DUMP_PPF_SCALE_MESSAGES_)
		printf("seq1 %d. scale: %.20f (%.20f, %.20f)\n", cnt1, this->seq1_scales[cnt1], this->ppf_loops->seq1_spf->ind_pairing_array[cnt1], this->ppf_loops->seq1_spf->ind_unpairing_array[cnt1]);
	}

	for(int cnt2 = 1; cnt2 <= seq_man->get_l_seq2(); cnt2++)
	{
		this->seq2_scales[cnt2] = sqrt(this->ppf_loops->seq2_spf->ind_pairing_array[cnt2] * this->ppf_loops->seq2_spf->ind_unpairing_array[cnt2]);
		//if( COMPARE(this->seq2_scales[cnt2], ZERO))
		if( COMPARE(this->ppf_loops->seq2_spf->ind_pairing_array[cnt2], 0.0) || COMPARE(this->ppf_loops->seq2_spf->ind_pairing_array[cnt2], 1.0))
		{
			 this->seq2_scales[cnt2] = 1.0;
		}
		else
		{
			this->seq2_scales[cnt2] = 1 / this->seq2_scales[cnt2];
		}
if(_DUMP_PPF_SCALE_MESSAGES_)
                printf("seq2 %d. scale: %.20f (%.20f, %.20f)\n", cnt2, this->seq2_scales[cnt2], this->ppf_loops->seq2_spf->ind_pairing_array[cnt2], this->ppf_loops->seq2_spf->ind_unpairing_array[cnt2]);

		seq2_cum_scale *= this->seq2_scales[cnt2];
	}

if(_DUMP_PPF_SCALE_MESSAGES_)
{
	printf("Cumulative scale for seq1 is %G, for seq2 is %G\n", seq1_cum_scale, seq2_cum_scale);
	cin.get();
}

	// Set the rescaling factor per nucleotide, set it to default value if it is not overriden from command line.
	if(ppf_cli->cmd_rescaling_increment_factor_per_nucleotide != 0.0)
	{
		this->rescaling_increment_factor_per_nucleotide = ppf_cli->cmd_rescaling_increment_factor_per_nucleotide; // Override.
	}
	else
	{
		this->rescaling_increment_factor_per_nucleotide = DEFAULT_PER_NUCLEOTIDE_RESCALE_INCREMENT; // Set to default.
	}
}

// Destructor.
t_ppf_scale::~t_ppf_scale()
{
	free(this->seq1_scales);
	free(this->seq2_scales);
}

// Return MULtiplication of scaling factors at two indices in arguments.
double t_ppf_scale::get_scale_value(int i1, int i2)
{
	return(this->seq1_scales[i1] * this->seq2_scales[i2]);
}

// Return partial scale value for x_1(s1, e1) and x2_(s2, e2).
double t_ppf_scale::get_partial_scale_factor(int s1, int e1, int s2, int e2)
{
	double partial_scale = 1.0;

	for(int cnt1 = s1; cnt1 <= e1; cnt1++)
	{
		if(cnt1 > this->seq_man->get_l_seq1())
		{
			partial_scale *= this->seq1_scales[cnt1 - this->seq_man->get_l_seq1()];
		}
		else
		{
			partial_scale *= this->seq1_scales[cnt1];
		}
	}

	for(int cnt2 = s2; cnt2 <= e2; cnt2++)
	{
		if(cnt2 > this->seq_man->get_l_seq2())
		{
			partial_scale *= this->seq2_scales[cnt2 - this->seq_man->get_l_seq2()];
		}
		else
		{
			partial_scale *= this->seq2_scales[cnt2];
		}
	}

	return(partial_scale);
}

// Return partial scale value for x_1(s1, e1) and x2_(s2, e2).
double t_ppf_scale::get_log_partial_scale_factor(int s1, int e1, int s2, int e2)
{
	double log_partial_scale = 0.0;

	for(int cnt1 = s1; cnt1 <= e1; cnt1++)
	{
		if(cnt1 > this->seq_man->get_l_seq1())
		{
			log_partial_scale += log(this->seq1_scales[cnt1 - this->seq_man->get_l_seq1()]);
		}
		else
		{
			log_partial_scale += log(this->seq1_scales[cnt1]);
		}
	}

	for(int cnt2 = s2; cnt2 <= e2; cnt2++)
	{
		if(cnt2 > this->seq_man->get_l_seq2())
		{
			log_partial_scale += log(this->seq2_scales[cnt2 - this->seq_man->get_l_seq2()]);
		}
		else
		{
			log_partial_scale += log(this->seq2_scales[cnt2]);
		}
	}

	return(log_partial_scale);
}


double t_ppf_scale::get_unscaled_log_value_by_indices(int i, int j, int k, int l, double scaled_value)
{
	// Can do linear multiplication comfortably.
	//printf("%.20f\n", log(scaled_value));
	//printf("%.20f\n", this->get_partial_scale_factor(i,j,k,l));
	//return( log(scaled_value) - log(this->get_partial_scale_factor(i,j,k,l)) ); 
	if(scaled_value == 0.0)
	{
		return(ZERO);
	}
	else
	{
		return( log(scaled_value) - this->get_log_partial_scale_factor(i,j,k,l) );
	}
}

void t_ppf_scale::copy_scales_to_spf_by_id(int id, double* spf_scales)
{
	// Copy scales from ppf_scaler for fast accession to arrays.
	if(id == 0) // first sequence?
	{
		for(int cnt = 1; cnt <= this->seq_man->get_l_seq1(); cnt++)
		{
			spf_scales[cnt] = this->seq1_scales[cnt];
		}
	}
	else if(id == 1) // second sequence?
	{
		for(int cnt = 1; cnt <= this->seq_man->get_l_seq2(); cnt++)
		{
			spf_scales[cnt] = this->seq2_scales[cnt];
		}
	}
	else
	{
		printf("WTF @ %s(%d)?\n", __FILE__, __LINE__);
		exit(0);
	}
}

// Check if a rescale is needed for the arrays in arguments.
// This function interrupts main loop calculations in pp mode.
bool t_ppf_scale::check_pp_ppf_array_rescale(int i, int j, int k, int l)
{
	// Check all array values against MAX_ARRAY_VALUE_LIMIT,
	// if there is no possible overflow yet, return false and return.
	// otherwise rescale all arrays over again with rescale_factor.
	if(this->ppf_loops->WL->x(i,j,k,l) > MAX_ARRAY_VALUE_LIMIT ||
		this->ppf_loops->W->x(i,j,k,l) > MAX_ARRAY_VALUE_LIMIT ||
		this->ppf_loops->WMBL->x(i,j,k,l) > MAX_ARRAY_VALUE_LIMIT ||
        this->ppf_loops->WMB->x(i,j,k,l) > MAX_ARRAY_VALUE_LIMIT ||
		this->ppf_loops->W_mhi->x_ij(i,j,k,l) > MAX_ARRAY_VALUE_LIMIT ||
		this->ppf_loops->W_mhi->x_kl(i,j,k,l) > MAX_ARRAY_VALUE_LIMIT ||
		this->ppf_loops->V_mhe->x(i,j,k,l) > MAX_ARRAY_VALUE_LIMIT ||
		this->ppf_loops->W_ext->x(j,l) > MAX_ARRAY_VALUE_LIMIT)
	{
		char panic_msg[500];
		sprintf(panic_msg, "One of the PPF arrays of indices (%d, %d, %d, %d) need down scaling @ %s(%d)", i,j,k,l, __FILE__, __LINE__);
		//ppf_panic(panic_msg, this->ppf_cli);

		// Update all scaling factors.
		// This is done for once and the rescaling 
		this->increment_scales_by_rescale(false);

		// Need a rescaling update and rescaling here.
		this->ppf_loops->rescale_all_ppf_arrays_by_indices(i, j, k, l, false);

		return(true);
	}
	else
	{
		return false;
	}
}

bool t_ppf_scale::check_pp_ppf_external_array_rescale(int i, int j, int k, int l)
{
	bool need_rescaling = false;
	// Check all array values against MAX_ARRAY_VALUE_LIMIT,
	// if there is no possible overflow yet, return false and return.
	// otherwise rescale all arrays over again with rescale_factor.
	if(this->ppf_loops->WL->pf_array->check_4D_ll(i,j,k,l))
	{
		if(this->ppf_loops->WL->x_ext(i,j,k,l) > MAX_ARRAY_VALUE_LIMIT)
		{
			need_rescaling = true;
		}
	}

	if(this->ppf_loops->W->pf_array->check_4D_ll(i,j,k,l))
	{
		if(this->ppf_loops->W->x_ext(i,j,k,l) > MAX_ARRAY_VALUE_LIMIT)
		{
			need_rescaling = true;
		}
	}

	if(this->ppf_loops->WMBL->pf_array->check_4D_ll(i,j,k,l))
	{
		if(this->ppf_loops->WMBL->x_ext(i,j,k,l) > MAX_ARRAY_VALUE_LIMIT)
		{
			need_rescaling = true;
		}
	}

	if(this->ppf_loops->WMB->pf_array->check_4D_ll(i,j,k,l))
	{
		if(this->ppf_loops->WMB->x_ext(i,j,k,l) > MAX_ARRAY_VALUE_LIMIT)
		{
			need_rescaling = true;
		}
	}

	if(this->ppf_loops->WMB->pf_array->check_4D_ll(i,j,k,l))
	{
		if(this->ppf_loops->WMB->x_ext(i,j,k,l) > MAX_ARRAY_VALUE_LIMIT)
		{
			need_rescaling = true;
		}
	}

	if(this->ppf_loops->W_mhi->kl_pf_array->check_4D_ll(i,j,k,l))
	{
		if(this->ppf_loops->W_mhi->x_ext_kl(i,j,k,l) > MAX_ARRAY_VALUE_LIMIT)
		{
			need_rescaling = true;
		}
	}

	if(this->ppf_loops->W_mhi->ij_pf_array->check_4D_ll(i,j,k,l))
	{
		if(this->ppf_loops->W_mhi->x_ext_ij(i,j,k,l) > MAX_ARRAY_VALUE_LIMIT)
		{
			need_rescaling = true;
		}
	}

	if(this->ppf_loops->W_ext->check_boundary(j-1,l-1))
	{
		if(this->ppf_loops->W_ext->x_ext(j,l) > MAX_ARRAY_VALUE_LIMIT)
		{
			need_rescaling = true;
		}
	}

	if(need_rescaling)
	{
		char panic_msg[4096];
		sprintf(panic_msg, "One of the PPF arrays of indices (%d, %d, %d, %d) need down scaling @ %s(%d)", i,j,k,l, __FILE__, __LINE__);

		this->increment_scales_by_rescale(false);

		this->ppf_loops->rescale_all_ppf_arrays_by_indices(i, j, k, l, false);
		this->ppf_loops->rescale_all_external_ppf_arrays_by_indices(i, j, k, l, false);

		return(true);
	}
	else
	{
		return false;
	}
}

// Update all scales.
// Note that if up_scale is true, this means that we need to increase 
// scaling values, that is, an up_scale is needed because of an underflow.
void t_ppf_scale::increment_scales_by_rescale(bool up_scale)
{


	for(int cnt1 = 1; cnt1 <= this->seq_man->get_l_seq1(); cnt1++)
	{
		if(up_scale)
		{
			this->seq1_scales[cnt1] *= this->rescaling_increment_factor_per_nucleotide;
			this->ppf_loops->seq1_spf->ppf_scales[cnt1] *= this->rescaling_increment_factor_per_nucleotide;
		}
		else
		{
			this->seq1_scales[cnt1] /= this->rescaling_increment_factor_per_nucleotide;
			this->ppf_loops->seq1_spf->ppf_scales[cnt1] /= this->rescaling_increment_factor_per_nucleotide;
		}
		
if(_DUMP_PPF_SCALE_MESSAGES_)
{
		printf("Rescaled seq1 scale %d: %.10f\n", cnt1, this->seq1_scales[cnt1]);
}
	}

	for(int cnt2 = 1; cnt2 <= this->seq_man->get_l_seq2(); cnt2++)
	{
		if(up_scale)
		{
			this->seq2_scales[cnt2] *= this->rescaling_increment_factor_per_nucleotide;
			this->ppf_loops->seq2_spf->ppf_scales[cnt2] *= this->rescaling_increment_factor_per_nucleotide;
		}
		else
		{
			this->seq2_scales[cnt2] /= this->rescaling_increment_factor_per_nucleotide;
			this->ppf_loops->seq2_spf->ppf_scales[cnt2] /= this->rescaling_increment_factor_per_nucleotide;
		}

if(_DUMP_PPF_SCALE_MESSAGES_)
{
		printf("Rescaled seq2 scale %d: %.10f\n", cnt2, this->seq2_scales[cnt2]);
}
	}
}

// Return one cumulative rescaling factor for one round of rescaling to be used for rescaling all arrays.
// Call this function for rescaling of an array at indices (i,j,k,l). Since the rescaling factor for each nucleotide
// is constant, it is calculated by this->rescaling_increment_factor_per_nucleotide^(j-i+l-k+2).
double t_ppf_scale::cumulative_rescale_factor(int i, int j, int k, int l)
{
	return(pow(this->rescaling_increment_factor_per_nucleotide, (j - i + l - k + 2)));
}

double t_ppf_scale::cumulative_rescale_external_factor(int i, int j, int k, int l)
{
	return(pow(this->rescaling_increment_factor_per_nucleotide, (this->seq_man->get_l_seq1() + this->seq_man->get_l_seq2()) - (j - i + l - k + 2)));
}

// Needed for recaling factor of frag aln enum array, since it does not emit i and k.
double t_ppf_scale::cumulative_rescale_factor_exc_ik(int i, int j, int k, int l)
{
	return(pow(this->rescaling_increment_factor_per_nucleotide, (j - i + l - k)));
}


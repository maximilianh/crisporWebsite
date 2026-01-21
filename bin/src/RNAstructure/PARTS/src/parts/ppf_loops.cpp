#include <string.h>
#include <limits.h>
#include "process_sequences.h"
#include "single_pf_array.h"
#include <stdio.h>
#include <stdlib.h>
#include "ppf_loops.h"
#include "parts_compilation_directives.h"
#include "alignment_priors.h"
#include "parts_paths.h"
#include "ppf_cli.h"
#include "template_pf_array.h"
#include "array_mem_manager.h"

#include "../../../src/phmm/structure/structure_object.h"
#include "../../../src/phmm/phmm_aln.h"
#include "../../../src/phmm/phmm.h"
#include "ppf_scale.h"

#include "ppf_operators.h"
#include "ppf_loops.h"
#include <stdio.h>
#include <stdlib.h>

#include "parts_compilation_directives.h"
#include "ppf_loops.h"
#include "process_sequences.h"
#include "single_pf_array.h"
#include "ppf_w.h"
#include "ppf_w_l.h"
#include "ppf_w_mb.h"
#include "ppf_w_mbl.h"
#include "ppf_w_ext.h"
#include "process_sequences.h"
#include "structural_parameters.h"
#include "ppf_math.h"
#include "pf_alignment.h"
#include "alignment_priors.h"
#include "posterior_pairing_prob_loops.h"
#include "parts_paths.h"
#include "ppf_ss.h"
#include "ppf_scale.h"
#include "ppf_cli.h"
#include "ppf_timer.h"
#include "ppf_v_mhe.h"
#include "ppf_w_mhi.h"
#include "ppf_scale.h"

#include "ppf_progress_bar.h"
#include "ppf_tb_stack.h"
#include "../../../src/phmm/structure/folding_constraints.h"
#include "../../../src/phmm/utils/file/utils.h"

#include "array_file_manager.h"

// Stochastic traceback code.
#include "stoch_tb/stoch_sampled_str_aln_sample_set.h"
#include "stoch_tb/stoch_sampling_math.h"
#include "stoch_tb/stoch_sampled_alignment.h"
#include "stoch_tb/stoch_sampled_structures.h"

#include "map_structures.h"
#include "map_alignment.h"

#include "pp_results.h"

#include "template_pf_array.h"

#include <iostream>

using namespace std;

bool _DUMP_PPF_LOOPS_MESSAGES_ = false;

t_ppf_loops::t_ppf_loops(int argc, char* argv[])
{
	this->last_error_code = NO_PPF_ERROR; 
	t_ppf_cli* ppf_cli = NULL;
	if(argc == 2)
	{
		this->ppf_cli = new t_ppf_cli(argv[1]);

		if(this->ppf_cli->GetErrorCode() != NO_CLI_ERROR)
		{
			this->last_error_code = ERR_CLI_ERROR;
			return;
		}
		else
		{
			this->init_loops();
		}
	}
	else
	{
		this->ppf_cli = new t_ppf_cli(argc, argv);

		if(this->ppf_cli->GetErrorCode() != NO_CLI_ERROR)
		{
			this->last_error_code = ERR_CLI_ERROR;
			return;
		}
		else
		{
			this->init_loops();
		}
	}
}

void t_ppf_loops::set_output_prefixes(const char* seq1_prefix, const char* seq2_prefix){
  ppf_cli->set_output_prefixes(seq1_prefix, seq2_prefix);
  return;
}

t_ppf_loops::t_ppf_loops(char* conf_file)
{
	this->last_error_code = NO_PPF_ERROR; 

	this->ppf_cli = new t_ppf_cli(conf_file);

	if(this->ppf_cli->GetErrorCode() != NO_CLI_ERROR)
	{
		this->last_error_code = ERR_CLI_ERROR;
	}

	this->init_loops();
}

t_ppf_loops::t_ppf_loops(char* seq1_fp,
						char* seq2_fp, 
						int mode)
{
	this->last_error_code = NO_PPF_ERROR; 

	this->ppf_cli = new t_ppf_cli(seq1_fp, seq2_fp, mode);

	if(this->ppf_cli->GetErrorCode() != NO_CLI_ERROR)
	{
		this->last_error_code = ERR_CLI_ERROR;
	}
	else
	{
		this->init_loops();
	}
}

t_ppf_loops::t_ppf_loops(char* seq1_fp,
						char* seq2_fp)
{
	this->last_error_code = NO_PPF_ERROR; 

	this->ppf_cli = new t_ppf_cli(seq1_fp, seq2_fp, PARTS_RUN_MODE_INIT);

	if(this->ppf_cli->GetErrorCode() != NO_CLI_ERROR)
	{
		this->last_error_code = ERR_CLI_ERROR;
	}
	else
	{
		this->init_loops();
	}
}

t_ppf_loops::~t_ppf_loops()
{
	// Delete arrays.
	delete(this->seq1_spf);
	delete(this->seq2_spf);
	delete(this->aln_priors);

	if(this->seq1_pp != NULL)
	{	
		for(int i = 1; i <= this->seq_man->get_l_seq1(); i++)
		{
			this->seq1_pp[i] += (i);
			free(this->seq1_pp[i]);

		} // i loop
		free(this->seq1_pp);
	}

	// Allocate and compute seq2 pairing probabilities.
	if(this->seq2_pp != NULL)
	{	
		for(int k = 1; k <= this->seq_man->get_l_seq2(); k++)
		{
			this->seq2_pp[k] += (k);
			free(this->seq2_pp[k]);
		} // k allocation loop
		free(this->seq2_pp);
	}

#ifdef _LINEAR_COMPUTATIONS_
	delete(this->ppf_scaler);
#endif

	if(this->map_alignment != NULL)
	{
		delete(this->map_alignment);
	}

	if(this->map_strs != NULL)
	{
		delete(this->map_strs);
	}

	delete(this->ppf_cli);

	delete(this->seq_man);	
}

int t_ppf_loops::GetErrorCode()
{
	return(this->last_error_code);
}

char* t_ppf_loops::GetErrorMessage(const int error_code)
{
	if(error_code == ERR_CLI_ERROR)
	{
		return(ppf_cli_error_msgs[this->ppf_cli->GetErrorCode()]);
	}
	else
	{
		return(ppf_loops_error_msgs[error_code]);
	}
}

int t_ppf_loops::mallocate_arrays()
{
	//this->SS = NULL;
	//this->V_mhe = NULL;
	//this->WL = NULL;
	//this->W = NULL;
	//this->WMB = NULL;
	//this->WMBL = NULL;
	//this->W_mhi = NULL;
	//this->W_ext = NULL;

	// Add error checking in the following.

if(_DUMP_PPF_LOOPS_MESSAGES_)
	printf("[SS]\n");

	this->SS = new t_ppf_SS(this);

if(_DUMP_PPF_LOOPS_MESSAGES_)
	printf("[Vmhe]\n");

	this->V_mhe = new t_ppf_V_mhe(this);

if(_DUMP_PPF_LOOPS_MESSAGES_)
	printf("[WL]\n");

	this->WL = new t_ppf_WL(this);

if(_DUMP_PPF_LOOPS_MESSAGES_)
	printf("[W]\n");

	this->W = new t_ppf_W(this);

if(_DUMP_PPF_LOOPS_MESSAGES_)
	printf("[WMB]\n");

	this->WMB = new t_ppf_WMB(this);

if(_DUMP_PPF_LOOPS_MESSAGES_)
	printf("[WMBL]\n");

	this->WMBL = new t_ppf_WMBL(this);

if(_DUMP_PPF_LOOPS_MESSAGES_)
	printf("[Wmhi]\n");

	this->W_mhi = new t_ppf_W_mhi(this);

if(_DUMP_PPF_LOOPS_MESSAGES_)
	printf("[Wext]\n");

	this->W_ext = new t_ppf_WEXT(this);

	return(0);
}

void t_ppf_loops::delete_arrays()
{
	delete(this->SS);
	delete(this->V_mhe);
	delete(this->WL);
	delete(this->W);
	delete(this->WMB);
	delete(this->WMBL);
	delete(this->W_mhi);
	delete(this->W_ext);
}


void t_ppf_loops::init_loops()
{
	printf("Initing loops.\n");

	this->seq_man = new t_seq_man(ppf_cli);	

	// Compute spf's and aln_priors.	
	// Compute/Read the loop limits and allocate alignment priors. 
	this->aln_priors = NULL;
	if(ppf_cli->loop_limits_fn != NULL)
	{
		// Load loop limits.
		t_template_pf_array::update_loop_limits(seq_man, ppf_cli->loop_limits_fn);

		// Copy the posterior probabilities.
		aln_priors = new t_aln_priors(seq_man, ppf_cli, true);
	}
	else
	{	
		// Run phmm on sequences.
		//t_phmm_aln* phmm_aln = new t_phmm_aln(ppf_cli->seq1_path, ppf_cli->seq2_path, ppf_cli->phmm_band_constraint_size);
		t_phmm_aln* phmm_aln = new t_phmm_aln(ppf_cli->seq1_path, ppf_cli->seq2_path);

		t_aln_env_result* aln_env_result = NULL;
		t_pp_result* pp_result = phmm_aln->compute_posterior_probs();
		if(ppf_cli->array_mem_limit_in_megs != 0.0f)
		{
			t_array_mem_manager* array_mem_manager = new t_array_mem_manager(seq_man, ppf_cli);
			aln_env_result = array_mem_manager->get_env_threshold_per_mem_limit(ppf_cli->array_mem_limit_in_megs * MEG_BYTES);

			if(aln_env_result == NULL)
			{
				printf("Alignment envelope computation failed.\n");
				exit(0);
			}

			t_template_pf_array::update_loop_limits(seq_man, aln_env_result->low_limits, aln_env_result->high_limits);

			//exit(0);
		}
		else // Utilize the threshold from phmm data.
		{
			// No alignment constraints.
			aln_env_result = phmm_aln->compute_alignment_envelope(PROB_ALN_ENV, 0);

			if(aln_env_result == NULL)
			{
				printf("Alignment envelope computation failed.\n");
				exit(0);
			}

			t_template_pf_array::update_loop_limits(seq_man, aln_env_result->low_limits, aln_env_result->high_limits);
		}

		//printf("before t_aln_priors\n");
		//getc(stdin);

		// Copy the posterior probabilities.
		aln_priors = new t_aln_priors(seq_man, 
										ppf_cli,
										pp_result->aln_probs,
										pp_result->ins1_probs,
										pp_result->ins2_probs, 
										true);

		//printf("after t_aln_priors\n");
		//getc(stdin);

		phmm_aln->free_aln_env_result(aln_env_result);
		delete(phmm_aln);
	}

	//printf("before t_spf_array\n");
	//getc(stdin);

	// At this point the loop limits and alignment priors should be set with all extra memory deleted.

	// Following ifdef's are particulalryl important since they determine how program treats
	// the sequences.
	// Determine fold prior 1 from spf.
	if(ppf_cli->fold_prior1_path != NULL)
	{
		this->seq1_spf = new t_spf_array(seq_man->get_l_seq1(), ppf_cli->seq1_path, ppf_cli, ppf_cli->fold_prior1_path, true);
	}
	else
	{
		this->seq1_spf = new t_spf_array(seq_man->get_l_seq1(), ppf_cli->seq1_path, ppf_cli, NULL, true);
	}

	// Determine fold prior 2 from spf.
	if(ppf_cli->fold_prior2_path != NULL)
	{
		this->seq2_spf = new t_spf_array(seq_man->get_l_seq2(), ppf_cli->seq2_path, ppf_cli, ppf_cli->fold_prior2_path, true);
	}
	else
	{
		this->seq2_spf = new t_spf_array(seq_man->get_l_seq2(), ppf_cli->seq2_path, ppf_cli, NULL, true);
	}

#ifdef _LINEAR_COMPUTATIONS_
	// Inistantiate a ppf scaler if doing the calaultions in linear space.
	this->ppf_scaler = new t_ppf_scale(this);

	// Set scalers for each sequence.
	seq1_spf->set_scaler(ppf_scaler->seq1_scales);
	seq2_spf->set_scaler(ppf_scaler->seq2_scales);
#else
	t_ppf_scale* ppf_scaler = NULL;
#endif // _LINEAR_COMPUTATIONS_

	this->seq1_pp = NULL;
	this->seq2_pp = NULL;

	this->stoch_sampled_set = NULL;
	this->stoch_sampler = NULL;
	this->tb_stack = NULL;

	this->map_alignment = NULL;
	this->map_strs = NULL;
	this->seq1_pp = NULL;
	this->seq2_pp = NULL;
}

// Following adds the samples base pairs/aligned positions to the appropriate structure/alignment.
void t_ppf_loops::set_aln(int i, int k)
{
	if(this->ppf_cli->mode == PARTS_RUN_MODE_STOCH_SAMPLE)
	{
		//random_cumulative = MUL(this->ppf_loops->stoch_sampler->random_double_1_0(), this->x(i,j,k,l));
		this->stoch_sampled_set->current_sampled_alignment()->set_aln(i, k);
	}
	else if(this->ppf_cli->mode == PARTS_RUN_MODE_MAP)
	{
		//random_cumulative = this->x(i,j,k,l);
		this->map_alignment->set_aln(i, k);
	}
}

void t_ppf_loops::set_seq1_ins(int i, int k)
{
	if(this->ppf_cli->mode == PARTS_RUN_MODE_STOCH_SAMPLE)
	{
		this->stoch_sampled_set->current_sampled_alignment()->set_seq1_ins(i, k);
	}
	else if(this->ppf_cli->mode == PARTS_RUN_MODE_MAP)
	{
		this->map_alignment->set_seq1_ins(i, k);
	}
}

void t_ppf_loops::set_seq2_ins(int i, int k)
{
	if(this->ppf_cli->mode == PARTS_RUN_MODE_STOCH_SAMPLE)
	{
		this->stoch_sampled_set->current_sampled_alignment()->set_seq2_ins(i, k);
	}
	else if(this->ppf_cli->mode == PARTS_RUN_MODE_MAP)
	{
		this->map_alignment->set_seq2_ins(i, k);
	}
}

signed int t_ppf_loops::Get_MAP_Alignment_Length()
{
	if(this->map_alignment == NULL)
	{
		this->last_error_code = ERR_MAP_COMPUTATIONS_NOT_PERFORMED;
		return(-1);
	}
	else
	{
		return(this->map_alignment->get_l_aln());
	}
}

//signed int t_ppf_loops::Get_Seq1_Alignment(const int i)
//{
//	if(this->map_alignment == NULL)
//	{
//		this->last_error_code = ERR_MAP_COMPUTATIONS_NOT_PERFORMED;
//		return(-1);
//	}
//	else
//	{
//		return(this->map_alignment->seq1_alns[i][0]);
//	}
//}
//
//signed int t_ppf_loops::Get_Seq2_Alignment(const int k)
//{
//	if(this->map_alignment == NULL)
//	{
//		this->last_error_code = ERR_MAP_COMPUTATIONS_NOT_PERFORMED;
//		return(-1);
//	}
//	else
//	{
//		return(this->map_alignment->seq2_alns[k][0]);
//	}
//}

signed int t_ppf_loops::GetAlignedNucleotideSeq1(const int nucleotide1)
{
        if(this->map_alignment == NULL)
        {
                this->last_error_code = ERR_MAP_COMPUTATIONS_NOT_PERFORMED;
                return(-1);
        }
        else if(nucleotide1 > this->GetLengthSequence1() || nucleotide1 < 1)
        {
                this->last_error_code = ERR_INVALID_ARGUMENT;
                return(-1);
        }
        else
        {
                int i_coinc_nuc = this->map_alignment->seq1_alns[nucleotide1][0];
                if(nucleotide1 == this->map_alignment->seq2_alns[i_coinc_nuc][0])
                {
                        return(i_coinc_nuc);
                }
                else
                {
                        return(0);
                }
        }
}


signed int t_ppf_loops::GetAlignedNucleotideSeq2(const int nucleotide2)
{
    if(this->map_alignment == NULL)
    {
        this->last_error_code = ERR_MAP_COMPUTATIONS_NOT_PERFORMED;
        return(-1);
    }
    else if(nucleotide2 > this->GetLengthSequence2() || nucleotide2 < 1)
    {
        this->last_error_code = ERR_INVALID_ARGUMENT;
        return(-1);
    }
    else
    {
        int k_coinc_nuc = this->map_alignment->seq2_alns[nucleotide2][0];
        if(nucleotide2 == this->map_alignment->seq1_alns[k_coinc_nuc][0])
        {
            return(k_coinc_nuc);
        }
        else
        {
            return(0);
        }
    }
}

signed int t_ppf_loops::GetAlignmentSeq1(int alignmentcolumn)
{
	if(this->map_alignment == NULL)
	{
		this->last_error_code = ERR_MAP_COMPUTATIONS_NOT_PERFORMED;
		return(-1);
	}

	int l_aln = this->map_alignment->get_l_aln();

	if(alignmentcolumn > l_aln || alignmentcolumn < 1)
	{
		this->last_error_code = ERR_INVALID_ARGUMENT;
		return(-1);
	}
	else
	{
		return(this->map_alignment->aln_index_line1[alignmentcolumn]);
	}
}

signed int t_ppf_loops::GetAlignmentSeq2(int alignmentcolumn)
{
	if(this->map_alignment == NULL)
	{
		this->last_error_code = ERR_MAP_COMPUTATIONS_NOT_PERFORMED;
		return(-1);
	}

	int l_aln = this->map_alignment->get_l_aln();

	if(alignmentcolumn > l_aln || alignmentcolumn < 1)
	{
		this->last_error_code = ERR_INVALID_ARGUMENT;
		return(-1);
	}
	else
	{
		return(this->map_alignment->aln_index_line2[alignmentcolumn]);
	}
}

int t_ppf_loops::GetLengthSequence1()
{
	return(this->seq_man->get_l_seq1());
}

int t_ppf_loops::GetLengthSequence2()
{
	return(this->seq_man->get_l_seq2());
}


double t_ppf_loops::GetPairProbabilitySeq1(const int i, const int j)
{
	if(this->seq1_pp == NULL)
	{
		this->last_error_code = ERR_PPF_COMPUTATIONS_NOT_PERFORMED;
		return(-1);
	}
	else
	{
		return(this->seq1_pp[i][j]);
	}
}

double t_ppf_loops::GetPairProbabilitySeq2(const int k, const int l)
{
	if(this->seq2_pp == NULL)
	{
		this->last_error_code = ERR_PPF_COMPUTATIONS_NOT_PERFORMED;
		return(-1);
	}
	else
	{
		return(this->seq2_pp[k][l]);
	}
}

void t_ppf_loops::add_ct1_bp(int i, int j)
{
	if(this->ppf_cli->mode == PARTS_RUN_MODE_STOCH_SAMPLE)
	{
		this->stoch_sampled_set->current_sampled_structures()->add_bp_sampled_ct1(i, j);
	}
	else if(this->ppf_cli->mode == PARTS_RUN_MODE_MAP)
	{
		this->map_strs->add_bp_to_MAP_ct1(i, j);
	}
}

void t_ppf_loops::add_ct2_bp(int k, int l)
{
	if(this->ppf_cli->mode == PARTS_RUN_MODE_STOCH_SAMPLE)
	{
		this->stoch_sampled_set->current_sampled_structures()->add_bp_sampled_ct2(k, l);
	}
	else if(this->ppf_cli->mode == PARTS_RUN_MODE_MAP)
	{
		this->map_strs->add_bp_to_MAP_ct2(k, l);
	}
}

int t_ppf_loops::compute_pseudo_free_energy_loops()
{
	if((ppf_cli->mode & PARTS_RUN_MODE_MAP)== PARTS_RUN_MODE_MAP)
		printf("PARTS (MAP)\n%s\n%s\n", ppf_cli->seq1_path, ppf_cli->seq2_path); // This is changed to MAP.
	else if((ppf_cli->mode & PARTS_RUN_MODE_PP)== PARTS_RUN_MODE_PP)
		printf("PARTS (PP)\n%s\n%s\n", ppf_cli->seq1_path, ppf_cli->seq2_path);
	else if((ppf_cli->mode & PARTS_RUN_MODE_STOCH_SAMPLE) == PARTS_RUN_MODE_STOCH_SAMPLE)
		printf("PARTS Sampling\n%s\n%s\n", ppf_cli->seq1_path, ppf_cli->seq2_path);
	
	// Set max_sum function.
	if(this->ppf_cli->mode == PARTS_RUN_MODE_PP || 
		this->ppf_cli->mode == PARTS_RUN_MODE_STOCH_SAMPLE)
	{
		this->MAX_SUM = ppf_sum_callback;
		this->BT_OP = ppf_geq_callback;
	}
	else if(this->ppf_cli->mode == PARTS_RUN_MODE_MAP)
	{
		this->MAX_SUM = ppf_max_callback;
		this->BT_OP = ppf_compare_callback;
	}
	else
	{
		printf("Cannot set pf operators for PARTS mode %d\n", this->ppf_cli->mode);
		exit(0);
	}

if(_DUMP_PPF_LOOPS_MESSAGES_)
	printf("PPF_LOOPS[MAX_SUM](5.0, 3.0) = %.5f\n", this->MAX_SUM(5.0, 3.0));

	// Allocate arrays.
	this->mallocate_arrays();

	// Compute SS first.
	this->SS->compute();

	// Compute internal loops.
	this->compute_internal_pseudo_free_energy_loops();	

	// Compute W_ext.
	W_ext->calculate_W_ext();

	// Tracebacks or pp computations based on mode of running.
	if(this->ppf_cli->mode == PARTS_RUN_MODE_MAP ||
		this->ppf_cli->mode == PARTS_RUN_MODE_STOCH_SAMPLE)
	{
		this->backtrack_pseudo_free_energy_loops();
	}
	else if(this->ppf_cli->mode == PARTS_RUN_MODE_PP)
	{
		W_ext->calculate_ext_W_ext();

if(_DUMP_PPF_LOOPS_MESSAGES_)
{
		printf("\nPPF via Wext(%d, %d) = %.10f\nPPF via Wext[ext](%d, %d) = %.10f\nIndirect: %.10f + %.10f + %.10f = %.10f\n", 
			this->seq_man->get_l_seq1(), 
			this->seq_man->get_l_seq2(), 
			W_ext->x(this->seq_man->get_l_seq1(), this->seq_man->get_l_seq2()), 
			1, 
			1, 
			W_ext->x_ext(1, 1),
			this->V_mhe->x(1, this->GetLengthSequence1(), 1, this->GetLengthSequence2()), 
			this->W->x(1, this->GetLengthSequence1(), 1, this->GetLengthSequence2()), 
			this->WMB->x(1, this->GetLengthSequence1(), 1, this->GetLengthSequence2()),
			SUM3(this->V_mhe->x(1, this->GetLengthSequence1(), 1, this->GetLengthSequence2()), 
			this->W->x(1, this->GetLengthSequence1(), 1, this->GetLengthSequence2()), 
			this->WMB->x(1, this->GetLengthSequence1(), 1, this->GetLengthSequence2())));
		getc(stdin);
}		

		this->compute_external_pseudo_free_energy_loops();

		// Compute posterior probabilities of pairing.
		this->compute_posterior_pp();
	}
	
	// Free array memory.
	this->delete_arrays();
	return(0);
}

int t_ppf_loops::compute_MAP()
{
	ppf_cli->mode = PARTS_RUN_MODE_MAP;
	printf("PARTS (MAP)\n%s\n%s\n", ppf_cli->seq1_path, ppf_cli->seq2_path); // This is changed to MAP.
	
	// Set operator.
	this->MAX_SUM = ppf_max_callback;
	this->BT_OP = ppf_compare_callback;

if(_DUMP_PPF_LOOPS_MESSAGES_)
	printf("PPF_LOOPS[MAX_SUM](5.0, 3.0) = %.5f\n", this->MAX_SUM(5.0, 3.0));

	// Allocate arrays.
	this->mallocate_arrays();

	// Compute SS first.
	this->SS->compute();

	// Compute internal loops.
	this->compute_internal_pseudo_free_energy_loops();

	// Compute W_ext.
	W_ext->calculate_W_ext();
	
	// Backtrack MAP structural alignment
	this->backtrack_pseudo_free_energy_loops();	

	// Free array memory.
	this->delete_arrays();

	return(this->last_error_code);
}

int t_ppf_loops::compute_PPF()
{
	ppf_cli->mode = PARTS_RUN_MODE_PP;

	printf("PARTS (PP)\n%s\n%s\n", ppf_cli->seq1_path, ppf_cli->seq2_path); // This is changed to MAP.
	
	// Set operators.
	this->MAX_SUM = ppf_sum_callback;
	this->BT_OP = ppf_geq_callback;

	printf("PPF_LOOPS[MAX_SUM](5.0, 3.0) = %.5f\n", this->MAX_SUM(5.0, 3.0));

	// Allocate arrays.
	this->mallocate_arrays();

	// Compute SS first.
	this->SS->compute();

	// Compute internal loops.
	this->compute_internal_pseudo_free_energy_loops();

	// Compute internal arrays.
        W_ext->calculate_W_ext();

	// Compute external arrays.
	W_ext->calculate_ext_W_ext();

if(_DUMP_PPF_LOOPS_MESSAGES_)
{
	printf("\nPPF via Wext(%d, %d) = %.10f\nPPF via Wext[ext](%d, %d) = %.10f\n", this->seq_man->get_l_seq1(), 
		this->seq_man->get_l_seq2(), W_ext->x(this->seq_man->get_l_seq1(), this->seq_man->get_l_seq2()), 
		1, 1, W_ext->x_ext(1, 1));	
}

	this->compute_external_pseudo_free_energy_loops();

	// Compute posterior probabilities of pairing.
	this->compute_posterior_pp();

if(_DUMP_PPF_LOOPS_MESSAGES_)
{
	for(int i = 1; i <= this->seq_man->get_l_seq1(); i++)
	{
		for(int j = i+1; j <= this->seq_man->get_l_seq1(); j++)
		{
			//printf("seq1(%d, %d) -> %lf\n", i, j, map_loops->GetPairProbabilitySeq1(i, j));
			printf("seq1(%d, %d) -> %lf\n", i, j, this->seq1_pp[i][j]);
		}
	}

	for(int k = 1; k <= this->seq_man->get_l_seq2(); k++)
	{
		for(int l = k+1; l <= this->seq_man->get_l_seq2(); l++)
		{
			printf("seq2(%d, %d) -> %lf\n", k, l, this->GetPairProbabilitySeq2(k, l));
		}
	}
}

	// Free array memory.
	this->delete_arrays();

	return(this->last_error_code);
}

void t_ppf_loops::compute_internal_pseudo_free_energy_loops()
{
	int r_i = 0;
	int r_j = 0;
	int r_k = 0;
	int r_l = 0;

	int N1 = seq_man->get_l_seq1();
	int N2 = seq_man->get_l_seq2();

	int low_k, high_k, low_l, high_l;

	// progress bar to print a nice progress bar.
	t_ppf_progress_bar* ppf_progress_bar = new t_ppf_progress_bar(ppf_cli, '=', true, this->seq_man->get_l_seq1());

	// Main loop calculation:
	for(int j = 1; j <= N1; j++)
	{
		ppf_progress_bar->update_bar(j);

		for( int i = j-1; i >= MAX(j - ppf_cli->max_n_separation_between_nucs, 1); i--)
		{
			bool ij_coinc = true;
			if(seq1_spf->folding_constraints->coinc_pointer_relocation_map != NULL && seq1_spf->folding_constraints->coinc_pointer_relocation_map[i][j] == POS_MEM_NOT_EXIST)
			{
				ij_coinc = false;
			}

			low_l = MAX( t_template_pf_array::low_limits[j], 1 );
			high_l = MIN( t_template_pf_array::high_limits[j], N2 );

			for ( int l = low_l; ij_coinc && l <= high_l; l++ )
			{
				low_k = MAX(MAX(l - ppf_cli->max_n_separation_between_nucs, 1), t_template_pf_array::low_limits[i-1] + 1);
				high_k = MIN(l - 1, t_template_pf_array::high_limits[i - 1]+1);

				for ( int k = high_k; k >= low_k; k-- )
				{
					bool kl_coinc = true;
					if(seq2_spf->folding_constraints->coinc_pointer_relocation_map != NULL && seq2_spf->folding_constraints->coinc_pointer_relocation_map[k][l] == POS_MEM_NOT_EXIST)
					{
						kl_coinc = false;
					}

					// Make sure that there is at least MIN_LOOP of distance between i-j and k-l.
					if(kl_coinc)
					{

if(_DUMP_PPF_LOOPS_MESSAGES_)
{
						printf("----------------------------------------------\n");
						printf("Calculating (%d, %d, %d, %d):\n", i, j, k, l);
}

#ifdef _LINEAR_COMPUTATIONS_
						do 
						{
#endif		
							if( ( (j > N1 || (j-i) > MIN_LOOP) &&  (l>N2 || (l-k) > MIN_LOOP) ) 
								&& ( (l <= N2 || k - (l - N2) > MIN_LOOP) && (j <= N1 || i - (j - N1) > MIN_LOOP) ) // To impose constraint on the length of interior loop when exterior loop is calculated.
								)
							{

								V_mhe->calculate_V_mhe(i, j, k, l); // Arrays to use.
							}

							W_mhi->calculate_ij_W_mhi(i, j, k, l); // Arrays to use. 
							W_mhi->calculate_kl_W_mhi(i, j, k, l); // Arrays to use. 
							WL->calculate_WL(i, j, k, l); // Arrays to use.
							W->calculate_W(i, j, k, l); // Arrays to use.
							WMBL->calculate_WMBL(i, j, k, l); // Arrays to use.
							WMB->calculate_WMB(i, j, k, l); // Arrays to use.

#ifdef _LINEAR_COMPUTATIONS_
						}
						while(ppf_scaler->check_pp_ppf_array_rescale(i, j, k, l));
#endif

extern bool _DUMP_PPF_V_MHE_MESSAGES_;

						//if(i == 1 && j == 73 && k == 1 && l == 72)
						//{
						//	printf("Vmhe(2, 72, 2, 71) = %lf\n", this->ppf_scaler->get_unscaled_log_value_by_indices(2,72,2,71,this->V_mhe->x(2,72,2,71)));
						//}

						//if(i == 1 && k == 1 && j >= 73 && l >= 72)
						//if(i == 1 && j == 73 && k == 1 && l == 72)						
						//{
						//	_DUMP_PPF_V_MHE_MESSAGES_ = true;

						//	//this->ppf_loops->V_mhe->calculate_V_mhe(i,j,k,l, false);

						//	//printf("V_mhe(%d, %d, %d, %d) = %.15f\n", i,j,k,l,this->ppf_loops->V_mhe->x(i,j,k,l));							

						//	//printf("Doing spurious Vmhe backtrack.\n");
						//	//printf("Vmhe(2, 72, 2, 71) = %lf\n", this->ppf_scaler->get_unscaled_log_value_by_indices(2,72,2,71,this->V_mhe->x(2,72,2,71)));

						//	// Allocate map structure objects for traceback.
						//	this->map_strs = new t_MAP_structures(seq_man, ppf_cli);

						//	// Alignment information.
						//	this->map_alignment = new t_MAP_alignment(seq_man, ppf_cli);

						//	this->tb_stack = new t_ppf_tb_stack();

						//	//this->ppf_loops->V_mhe->calculate_V_mhe(1, 73, 1, 72, true);
						//	this->V_mhe->calculate_V_mhe(2,72,2,71, true);
						//	printf("V_mhe(%d, %d, %d, %d) = %.15f\n", 2,72,2,71,this->V_mhe->x(2,72,2,71));	

						//	delete(this->map_strs);
						//	delete(this->map_alignment);
						//	delete(this->tb_stack);

						//	this->map_strs = NULL;
						//	this->map_alignment = NULL;
						//	this->tb_stack = NULL;
						//	exit(0);

						//	//_DUMP_PPF_V_MHE_MESSAGES_ = false;
						//}

					} // Loop limit if for internal and external checks on i-j and k-l.				
				} // k loop
			} // l loop
		} // i loop
	} // j loop

	delete(ppf_progress_bar);

	//this->dump_internal_arrays();

	//exit(0);
}

void t_ppf_loops::load_internal_arrays(char* arrays_dir)
{
	FILE* f_vmhe = fopen("vmhe.array", "rb");
	FILE* f_wl = fopen("wl.array", "rb");
	FILE* f_w = fopen("w.array", "rb");
	FILE* f_wmbl = fopen("wmbl.array", "rb");
	FILE* f_wmb = fopen("wmb.array", "rb");
	FILE* f_wmhi_ij = fopen("wmhi_ij.array", "rb");
	FILE* f_wmhi_kl = fopen("wmhi_kl.array", "rb");

	// Load vmhe.
	int i;
	int j;
	int k; 
	int l;
	double val;
	while(read_array_indices_value(f_vmhe, i, j, k, l, val))
	{		
		if(this->V_mhe->pf_array->check_4D_ll(i,j,k,l))
		{
			this->V_mhe->x_setter(i,j,k,l) = val;
		}
	}

	while(read_array_indices_value(f_wl, i, j, k, l, val))
	{
		if(this->WL->pf_array->check_4D_ll(i,j,k,l))
		{
			this->WL->x_setter(i,j,k,l) = val;
		}
	}

	while(read_array_indices_value(f_w, i, j, k, l, val))
	{
		if(this->W->pf_array->check_4D_ll(i,j,k,l))
		{
			this->W->x_setter(i,j,k,l) = val;
		}
	}

	while(read_array_indices_value(f_wmbl, i, j, k, l, val))
	{
		if(this->WMBL->pf_array->check_4D_ll(i,j,k,l))
		{
			this->WMBL->x_setter(i,j,k,l) = val;
		}
	}

	while(read_array_indices_value(f_wmb, i, j, k, l, val))
	{
		if(this->WMB->pf_array->check_4D_ll(i,j,k,l))
		{
			this->WMB->x_setter(i,j,k,l) = val;
		}
	}

	while(read_array_indices_value(f_wmhi_ij, i, j, k, l, val))
	{
		if(this->W_mhi->ij_pf_array->check_4D_ll(i,j,k,l))
		{
			this->W_mhi->x_ij_setter(i,j,k,l) = val;
		}
	}

	while(read_array_indices_value(f_wmhi_kl, i, j, k, l, val))
	{
		if(this->W_mhi->kl_pf_array->check_4D_ll(i,j,k,l))
		{
			this->W_mhi->x_kl_setter(i,j,k,l) = val;
		}
	}

	fclose(f_vmhe);
	fclose(f_wl);
	fclose(f_w);
	fclose(f_wmbl);
	fclose(f_wmb);
	fclose(f_wmhi_ij);
	fclose(f_wmhi_kl);
}

bool t_ppf_loops::read_array_indices_value(FILE* f_array, int& i, int& j, int& k, int& l, double& val)
{
	int cur_i = 0;
	int cur_j = 0;
	int cur_k = 0;
	int cur_l = 0;
	double cur_val = 0.0f;

	if(fread(&cur_i, sizeof(int), 1, f_array) != 1)
	{
		// EOF, return false, nothing is read.
		return(false);
	}

	if(fread(&cur_j, sizeof(int), 1, f_array) != 1)
	{
		printf("Could not read j!\n");
		exit(0);
		//break;
	}

	if(fread(&cur_k, sizeof(int), 1, f_array) != 1)
	{
		printf("Could not read k!\n");
		exit(0);
	}

	if(fread(&cur_l, sizeof(int), 1, f_array) != 1)
	{
		printf("Could not read l!\n");
		exit(0);
	}


	if(fread(&cur_val, sizeof(double), 1, f_array) != 1)
	{
		printf("Could not read cur_val!\n");
		exit(0);
	}
	else
	{
	}

	// Finished reading one set of values.
	printf("Read %d, %d, %d, %d, %lf\n", cur_i,cur_j,cur_k,cur_l,cur_val);
	i = cur_i;
	j = cur_j;
	k = cur_k;
	l = cur_l;
	val = cur_val;		

	return(true);
}

void t_ppf_loops::dump_internal_arrays()
{
	int r_i = 0;
	int r_j = 0;
	int r_k = 0;
	int r_l = 0;

	int N1 = seq_man->get_l_seq1();
	int N2 = seq_man->get_l_seq2();

	int low_k, high_k, low_l, high_l;

	FILE* f_vmhe = fopen("vmhe.array", "wb");
	FILE* f_wl = fopen("wl.array", "wb");
	FILE* f_w = fopen("w.array", "wb");
	FILE* f_wmbl = fopen("wmbl.array", "wb");
	FILE* f_wmb = fopen("wmb.array", "wb");
	FILE* f_wmhi_ij = fopen("wmhi_ij.array", "wb");
	FILE* f_wmhi_kl = fopen("wmhi_kl.array", "wb");

	// Main loop calculation:
	for(int j = 1; j <= N1; j++)
	{
		for( int i = j-1; i >= MAX(j - ppf_cli->max_n_separation_between_nucs, 1); i--)
		{
			bool ij_coinc = true;
			if(seq1_spf->folding_constraints->coinc_pointer_relocation_map != NULL && seq1_spf->folding_constraints->coinc_pointer_relocation_map[i][j] == POS_MEM_NOT_EXIST)
			{
				ij_coinc = false;
			}

			low_l = MAX( t_template_pf_array::low_limits[j], 1 );
			high_l = MIN( t_template_pf_array::high_limits[j], N2 );

			for ( int l = low_l; ij_coinc && l <= high_l; l++ )
			{
				low_k = MAX(MAX(l - ppf_cli->max_n_separation_between_nucs, 1), t_template_pf_array::low_limits[i-1] + 1);
				high_k = MIN(l - 1, t_template_pf_array::high_limits[i - 1]+1);

				for ( int k = high_k; k >= low_k; k-- )
				{
					bool kl_coinc = true;
					if(seq2_spf->folding_constraints->coinc_pointer_relocation_map != NULL && seq2_spf->folding_constraints->coinc_pointer_relocation_map[k][l] == POS_MEM_NOT_EXIST)
					{
						kl_coinc = false;
					}

					// Make sure that there is at least MIN_LOOP of distance between i-j and k-l.
					if(kl_coinc)
					{

if(_DUMP_PPF_LOOPS_MESSAGES_)
{
						printf("----------------------------------------------\n");
						printf("Dumping (%d, %d, %d, %d):\n", i, j, k, l);
}

#ifdef _LINEAR_COMPUTATIONS_
						//do 
						{
#endif		
							double cur_val = ZERO;
							if( ( (j > N1 || (j-i) > MIN_LOOP) &&  (l>N2 || (l-k) > MIN_LOOP) ) 
								&& ( (l <= N2 || k - (l - N2) > MIN_LOOP) && (j <= N1 || i - (j - N1) > MIN_LOOP) ) // To impose constraint on the length of interior loop when exterior loop is calculated.
								)
							{

								//V_mhe->calculate_V_mhe(i, j, k, l); // Arrays to use.
								cur_val = V_mhe->x(i,j,k,l);
								if(cur_val != ZERO)
								{
									//printf("Dumping %d, %d, %d, %d, %lf\n", i,j,k,l, cur_val);
									fwrite((void*)&i, sizeof(int), 1, f_vmhe);
									fwrite((void*)&j, sizeof(int), 1, f_vmhe);
									fwrite((void*)&k, sizeof(int), 1, f_vmhe);
									fwrite((void*)&l, sizeof(int), 1, f_vmhe);
									fwrite((void*)&cur_val, sizeof(double), 1, f_vmhe);
								}
							}

							//W_mhi->calculate_ij_W_mhi(i, j, k, l); // Arrays to use. 

							cur_val = W_mhi->x_ij(i,j,k,l);
							if(cur_val != ZERO)
							{
							fwrite(&i, sizeof(int), 1, f_wmhi_ij);
							fwrite(&j, sizeof(int), 1, f_wmhi_ij);
							fwrite(&k, sizeof(int), 1, f_wmhi_ij);
							fwrite(&l, sizeof(int), 1, f_wmhi_ij);
							fwrite(&cur_val, sizeof(double), 1, f_wmhi_ij);
							}

							//W_mhi->calculate_kl_W_mhi(i, j, k, l); // Arrays to use. 
							cur_val = W_mhi->x_kl(i,j,k,l);
							if(cur_val != ZERO)
							{
							fwrite(&i, sizeof(int), 1, f_wmhi_kl);
							fwrite(&j, sizeof(int), 1, f_wmhi_kl);
							fwrite(&k, sizeof(int), 1, f_wmhi_kl);
							fwrite(&l, sizeof(int), 1, f_wmhi_kl);
							fwrite(&cur_val, sizeof(double), 1, f_wmhi_kl);
							}

							//WL->calculate_WL(i, j, k, l); // Arrays to use.
							cur_val = WL->x(i,j,k,l);
							if(cur_val != ZERO)
							{
							fwrite(&i, sizeof(int), 1, f_wl);
							fwrite(&j, sizeof(int), 1, f_wl);
							fwrite(&k, sizeof(int), 1, f_wl);
							fwrite(&l, sizeof(int), 1, f_wl);
							fwrite(&cur_val, sizeof(double), 1, f_wl);
							}

							//W->calculate_W(i, j, k, l); // Arrays to use.
							cur_val = W->x(i,j,k,l);
							if(cur_val != ZERO)
							{
							fwrite(&i, sizeof(int), 1, f_w);
							fwrite(&j, sizeof(int), 1, f_w);
							fwrite(&k, sizeof(int), 1, f_w);
							fwrite(&l, sizeof(int), 1, f_w);
							fwrite(&cur_val, sizeof(double), 1, f_w);
							}

							//WMBL->calculate_WMBL(i, j, k, l); // Arrays to use.
							cur_val = WMBL->x(i,j,k,l);
							if(cur_val != ZERO)
							{
							fwrite(&i, sizeof(int), 1, f_wmbl);
							fwrite(&j, sizeof(int), 1, f_wmbl);
							fwrite(&k, sizeof(int), 1, f_wmbl);
							fwrite(&l, sizeof(int), 1, f_wmbl);
							fwrite(&cur_val, sizeof(double), 1, f_wmbl);
							}

							//WMB->calculate_WMB(i, j, k, l); // Arrays to use.
							cur_val = WMB->x(i,j,k,l);
							if(cur_val != ZERO)
							{
							fwrite(&i, sizeof(int), 1, f_wmb);
							fwrite(&j, sizeof(int), 1, f_wmb);
							fwrite(&k, sizeof(int), 1, f_wmb);
							fwrite(&l, sizeof(int), 1, f_wmb);
							fwrite(&cur_val, sizeof(double), 1, f_wmb);
							}

#ifdef _LINEAR_COMPUTATIONS_
						}
						//while(ppf_scaler->check_pp_ppf_array_rescale(i, j, k, l));
#endif

					} // Loop limit if for internal and external checks on i-j and k-l.				
				} // k loop
			} // l loop
		} // i loop
	} // j loop

	fclose(f_vmhe);
	fclose(f_wl);
	fclose(f_w);
	fclose(f_wmbl);
	fclose(f_wmb);
	fclose(f_wmhi_ij);
	fclose(f_wmhi_kl);
}

void t_ppf_loops::compute_posterior_pp()
{
	this->seq1_pp = (double**)malloc(sizeof(double*) * (this->seq_man->get_l_seq1() + 1));
	for(int i = 1; i <= this->seq_man->get_l_seq1(); i++)
	{
		this->seq1_pp[i] = (double*)malloc(sizeof(double) * (this->seq_man->get_l_seq1() - i + 4));
		this->seq1_pp[i] -= (i);

		for(int j = i; j <= this->seq_man->get_l_seq1(); j++)
		{
			this->seq1_pp[i][j] = ZERO;
		}
	}

	// Allocate and compute seq2 pairing probabilities.
	this->seq2_pp = (double**)malloc(sizeof(double*) * (this->seq_man->get_l_seq2() + 1));
	for(int k = 1; k <= this->seq_man->get_l_seq2(); k++)
	{
		this->seq2_pp[k] = (double*)malloc(sizeof(double) * (this->seq_man->get_l_seq2() - k + 4));
		this->seq2_pp[k] -= (k);

		for(int l = k; l <= this->seq_man->get_l_seq2(); l++)
		{
			this->seq2_pp[k][l] = ZERO;
		} // l allocation loop
	} // k allocation loop

	// 4D loop for marginalization.
	for(int i = 1; i <= this->seq_man->get_l_seq1(); i++)
	{
		for(int j = i+1; j <= this->seq_man->get_l_seq1(); j++)
		{
			int low_k = MAX(1, t_template_pf_array::low_limits[i-1] + 1);
			int high_k = MIN(this->seq_man->get_l_seq2(), t_template_pf_array::high_limits[i - 1]+1);

			for(int k = low_k; k <= high_k; k++)
			{
				int low_l = MAX( t_template_pf_array::low_limits[j], k+1 );
				int high_l = MIN( t_template_pf_array::high_limits[j], this->seq_man->get_l_seq2() );

				for(int l = low_l; l <= high_l; l++)
				{
					bool i_inc_k_inc = this->V_mhe->check_boundary(i+1, k+1);
					bool i_dec_k_dec = this->V_mhe->check_boundary(i-1, k-1);
					bool j_dec_l_dec = this->V_mhe->check_boundary(j-1, l-1);
					bool j_dec_l = this->V_mhe->check_boundary(j-1, l);
					bool j_l_dec = this->V_mhe->check_boundary(j, l-1);
					bool i_k = this->V_mhe->check_boundary(i, k);
					bool i_k_dec = this->V_mhe->check_boundary(i, k-1);
					bool i_dec_k = this->V_mhe->check_boundary(i-1, k);

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

					double internal_energies = ZERO;

					// Update seq1_pp
					// .. -> semiM1
					if(this->W_mhi->ij_pf_array->check_4D_ll(i,j,k,l))
					{
						// semiM1 -> semiM1
						double internal_energies = this->W_mhi->x_ij(i, j, k, l);
						this->seq1_pp[i][j] = SUM(this->seq1_pp[i][j], MUL(internal_energies, this->W_mhi->x_ext_ij(i,j,k,l)));
					}
					
					// .. -> M
					if(this->V_mhe->pf_array->check_4D_ll(i,j,k,l))
					{
						// M -> M
						// abp
						internal_energies = MUL(abp_score, this->V_mhe->x(i+1, j-1, k+1, l-1));
						this->seq1_pp[i][j] = SUM(this->seq1_pp[i][j], MUL(internal_energies, this->V_mhe->x_ext(i,j,k,l)));

						// bau1
						internal_energies = MUL(ij_bau_score, this->V_mhe->x(i+1, j-1, k+1, l-1));
						this->seq1_pp[i][j] = SUM(this->seq1_pp[i][j], MUL(internal_energies, this->V_mhe->x_ext(i,j,k,l)));

						// bpi1
						internal_energies = MUL(ij_bpi_score, this->V_mhe->x(i+1, j-1, k, l));
						this->seq1_pp[i][j] = SUM(this->seq1_pp[i][j], MUL(internal_energies, this->V_mhe->x_ext(i,j,k,l)));

						// semiM1 -> M
						// abp
						internal_energies = MUL(abp_score, this->W_mhi->x_ij(i+1, j-1, k+1, l-1));
						this->seq1_pp[i][j] = SUM(this->seq1_pp[i][j], MUL(internal_energies, this->V_mhe->x_ext(i,j,k,l)));

						// semiM2 -> M
						// abp
						internal_energies = MUL(abp_score, this->W_mhi->x_kl(i+1, j-1, k+1, l-1));
						this->seq1_pp[i][j] = SUM(this->seq1_pp[i][j], MUL(internal_energies, this->V_mhe->x_ext(i,j,k,l)));

						// bau1
						internal_energies = MUL(ij_bau_score, this->W_mhi->x_kl(i+1, j-1, k+1, l-1));
						this->seq1_pp[i][j] = SUM(this->seq1_pp[i][j], MUL(internal_energies, this->V_mhe->x_ext(i,j,k,l)));

						// bpi1
						internal_energies = MUL(ij_bpi_score, this->W_mhi->x_kl(i+1, j-1, k, l));
						this->seq1_pp[i][j] = SUM(this->seq1_pp[i][j], MUL(internal_energies, this->V_mhe->x_ext(i,j,k,l)));

						// nonMHE -> M
						// abp
						internal_energies = MUL(abp_score, this->SS->x(i+1, j-1, k+1, l-1));
						this->seq1_pp[i][j] = SUM(this->seq1_pp[i][j], MUL(internal_energies, this->V_mhe->x_ext(i,j,k,l)));

						internal_energies = MUL(abp_score, this->W->x(i+1, j-1, k+1, l-1));
						this->seq1_pp[i][j] = SUM(this->seq1_pp[i][j], MUL(internal_energies, this->V_mhe->x_ext(i,j,k,l)));

						internal_energies = MUL(abp_score, this->WMB->x(i+1, j-1, k+1, l-1));
						this->seq1_pp[i][j] = SUM(this->seq1_pp[i][j], MUL(internal_energies, this->V_mhe->x_ext(i,j,k,l)));
					}

					// Update seq2_pp
					// .. -> semiM2
					if(this->W_mhi->kl_pf_array->check_4D_ll(i,j,k,l))
					{
						// semiM2 -> semiM2
						internal_energies = this->W_mhi->x_kl(i, j, k, l);
						this->seq2_pp[k][l] = SUM(this->seq2_pp[k][l], MUL(internal_energies, this->W_mhi->x_ext_kl(i, j, k, l)));
					}

					// .. -> M
					if(this->V_mhe->pf_array->check_4D_ll(i,j,k,l))
					{
						// M -> M
						// abp
						internal_energies = MUL(abp_score, this->V_mhe->x(i+1, j-1, k+1, l-1));
						this->seq2_pp[k][l] = SUM(this->seq2_pp[k][l], MUL(internal_energies, this->V_mhe->x_ext(i, j, k, l)));

						//bau2
						internal_energies = MUL(kl_bau_score, this->V_mhe->x(i+1, j-1, k+1, l-1));
						this->seq2_pp[k][l] = SUM(this->seq2_pp[k][l], MUL(internal_energies, this->V_mhe->x_ext(i, j, k, l)));

						//bpi2
						internal_energies = MUL(kl_bpi_score, this->V_mhe->x(i, j, k+1, l-1));
						this->seq2_pp[k][l] = SUM(this->seq2_pp[k][l], MUL(internal_energies, this->V_mhe->x_ext(i, j, k, l)));

						// semiM2 -> M
						// abp
						internal_energies = MUL(abp_score, this->W_mhi->x_kl(i+1, j-1, k+1, l-1));
						this->seq2_pp[k][l] = SUM(this->seq2_pp[k][l], MUL(internal_energies, this->V_mhe->x_ext(i, j, k, l)));

						// semiM1 -> M
						// abp
						internal_energies = MUL(abp_score, this->W_mhi->x_ij(i+1, j-1, k+1, l-1));
						this->seq2_pp[k][l] = SUM(this->seq2_pp[k][l], MUL(internal_energies, this->V_mhe->x_ext(i, j, k, l)));

						// bau2
						internal_energies = MUL(kl_bau_score, this->W_mhi->x_ij(i+1, j-1, k+1, l-1));
						this->seq2_pp[k][l] = SUM(this->seq2_pp[k][l], MUL(internal_energies, this->V_mhe->x_ext(i, j, k, l)));

						// bpi2
						internal_energies = MUL(kl_bpi_score, this->W_mhi->x_ij(i, j, k+1, l-1));
						this->seq2_pp[k][l] = SUM(this->seq2_pp[k][l], MUL(internal_energies, this->V_mhe->x_ext(i, j, k, l)));

						// nonM -> M
						// abp
						internal_energies = MUL(abp_score, this->SS->x(i+1, j-1, k+1, l-1));
						this->seq2_pp[k][l] = SUM(this->seq2_pp[k][l], MUL(internal_energies, this->V_mhe->x_ext(i, j, k, l)));

						internal_energies = MUL(abp_score, this->W->x(i+1, j-1, k+1, l-1));
						this->seq2_pp[k][l] = SUM(this->seq2_pp[k][l], MUL(internal_energies, this->V_mhe->x_ext(i, j, k, l)));

						internal_energies = MUL(abp_score, this->WMB->x(i+1, j-1, k+1, l-1));
						this->seq2_pp[k][l] = SUM(this->seq2_pp[k][l], MUL(internal_energies, this->V_mhe->x_ext(i, j, k, l)));
					} // .. -> M
				} // l loop
			} // k loop
		} // j loop
	} // i loop

	char seq1_pp_dump_fp[4096];
	if(this->ppf_cli->seq1_pp_op == NULL)
	{
		sprintf(seq1_pp_dump_fp, "%s_pp.txt", this->ppf_cli->seq1_op_file_prefix);
	}
	else
	{
		strcpy(seq1_pp_dump_fp, this->ppf_cli->seq1_pp_op);
	}
	FILE* f_pp = open_f(seq1_pp_dump_fp, "w");
	for(int i = 1; i <= this->seq_man->get_l_seq1(); i++)
	{
		for(int j = 1; j <= this->seq_man->get_l_seq1(); j++)
		{
			if(j > i)
			{
				this->seq1_pp[i][j] = DIV(this->seq1_pp[i][j], this->W_ext->x_ext(1,1));
				if(GT(this->seq1_pp[i][j], CONVERT_FROM_LIN(1.0f)))
				{
					printf("Probability greater than 1 for (%d, %d), %.15f @ %s(%d)\n", i, j, this->seq1_pp[i][j], __FILE__, __LINE__);
					getc(stdin);
				}

				fprintf(f_pp, "%.5f ", this->seq1_pp[i][j]);
			}
			else if(j < i)
			{
				fprintf(f_pp, "%.5f ", this->seq1_pp[j][i]);
			}
			else
			{
				fprintf(f_pp, "%.5f ", ZERO);
			}
		}
		fprintf(f_pp, "\n");
	}
	fclose(f_pp);

	char seq2_pp_dump_fp[4096];
	if(this->ppf_cli->seq2_pp_op == NULL)
	{
		sprintf(seq2_pp_dump_fp, "%s_pp.txt", this->ppf_cli->seq2_op_file_prefix);
	}
	else
	{
		strcpy(seq2_pp_dump_fp, this->ppf_cli->seq2_pp_op);
	}

	f_pp = open_f(seq2_pp_dump_fp, "w");
	for(int k = 1; k <= this->seq_man->get_l_seq2(); k++)
	{
		for(int l = 1; l <= this->seq_man->get_l_seq2(); l++)
		{
			if(l > k)
			{
				this->seq2_pp[k][l] = DIV(this->seq2_pp[k][l], this->W_ext->x_ext(1,1));
				if(GT(this->seq2_pp[k][l], CONVERT_FROM_LIN(1.0f)))
				{
					printf("Probability greater than 1 for (%d, %d), %.15f @ %s(%d)\n", k, l, this->seq2_pp[k][l], __FILE__, __LINE__);
					getc(stdin);
				}

				fprintf(f_pp, "%.5f ", this->seq2_pp[k][l]);
			}
			else if(l < k)
			{
				fprintf(f_pp, "%.5f ", this->seq2_pp[l][k]);
			}
			else
			{
				fprintf(f_pp, "%.5f ", ZERO);
			}
		}
		fprintf(f_pp, "\n");
	}
	fclose(f_pp);
}



void t_ppf_loops::compute_external_pseudo_free_energy_loops()
{
	int r_i = 0;
	int r_j = 0;
	int r_k = 0;
	int r_l = 0;

	int N1 = seq_man->get_l_seq1();
	int N2 = seq_man->get_l_seq2();

	int low_k, high_k, low_l, high_l;

	// progress bar to print a nice progress bar.
	t_ppf_progress_bar* ppf_progress_bar = new t_ppf_progress_bar(ppf_cli, '=', true, this->seq_man->get_l_seq1());

	// Boundary conditions.
	//this->V_mhe->x_ext(1, N1, 1, N2) = CONVERT_FROM_LIN(1.0f);
	//this->W->x_ext(1, N1, 1, N2) = CONVERT_FROM_LIN(1.0f);
	//this->WMB->x_ext(1, N1, 1, N2) = CONVERT_FROM_LIN(1.0f);

	// Main loop calculation:
	for( int j = N1; j >= 1; j--)
	{
		ppf_progress_bar->update_bar(j);

		//for( int i = j-1; i >= MAX(j - ppf_cli->max_n_separation_between_nucs, 1); i--)
		for( int i = MAX(j - ppf_cli->max_n_separation_between_nucs, 1); i <= j-1 ; i++)
		{
			bool ij_coinc = true;
			if(seq1_spf->folding_constraints->coinc_pointer_relocation_map != NULL && seq1_spf->folding_constraints->coinc_pointer_relocation_map[i][j] == POS_MEM_NOT_EXIST)
			{
				ij_coinc = false;
			}

			low_l = MAX( t_template_pf_array::low_limits[j], 1 );
			high_l = MIN( t_template_pf_array::high_limits[j], N2 );

			//for ( int l = low_l; ij_coinc && l <= high_l; l++ )
			for(int l = high_l; ij_coinc && l >= low_l; l--)
			{
				low_k = MAX(MAX(l - ppf_cli->max_n_separation_between_nucs, 1), t_template_pf_array::low_limits[i-1] + 1);
				high_k = MIN(l - 1, t_template_pf_array::high_limits[i - 1]+1);

				//for ( int k = high_k; k >= low_k; k-- )
				for ( int k = low_k; k <= high_k; k++ )
				{
					bool kl_coinc = true;
					if(seq2_spf->folding_constraints->coinc_pointer_relocation_map != NULL && seq2_spf->folding_constraints->coinc_pointer_relocation_map[k][l] == POS_MEM_NOT_EXIST)
					{
						kl_coinc = false;
					}

					// Make sure that there is at least MIN_LOOP of distance between i-j and k-l.
					if(kl_coinc)
					{

if(_DUMP_PPF_LOOPS_MESSAGES_)
{
						printf("----------------------------------------------\n");
						printf("Calculating[ext] (%d, %d, %d, %d):\n", i, j, k, l);
}

#ifdef _LINEAR_COMPUTATIONS_
						while(ppf_scaler->check_pp_ppf_external_array_rescale(i, j, k, l)){}
						{
#endif		
							WMB->calculate_ext_dependencies(i, j, k, l); // Arrays to use.
							WMBL->calculate_ext_dependencies(i, j, k, l); // Arrays to use.
							W->calculate_ext_dependencies(i, j, k, l); // Arrays to use.
							WL->calculate_ext_dependencies(i, j, k, l); // Arrays to use.
							W_mhi->calculate_kl_ext_dependencies(i, j, k, l); // Arrays to use. 
							W_mhi->calculate_ij_ext_dependencies(i, j, k, l); // Arrays to use. 

							if( ( (j > N1 || (j-i) > MIN_LOOP) &&  (l>N2 || (l-k) > MIN_LOOP) ) 
								&& ( (l <= N2 || k - (l - N2) > MIN_LOOP) && (j <= N1 || i - (j - N1) > MIN_LOOP) ) // To impose constraint on the length of interior loop when exterior loop is calculated.
								)
							{
								// Initialize with Wext, this sets the boundary condition.
								V_mhe->x_ext(i,j,k,l) = MAX_SUM(V_mhe->x_ext(i,j,k,l), 
									MUL(this->W_ext->x(i-1, k-1), this->W_ext->x_ext(j+1, l+1)));

								//if(this->W_ext->x(i-1, k-1) < this->SS->x(1, i-1, 1, k-1))
								//{
								//	printf("%s(%d), %.15f, %.15f, %d, %d (%d, %d)\n", __FILE__, __LINE__, this->W_ext->x(i-1, k-1), this->SS->x(1, i-1, 1, k-1),i,k, N1, N2);
								//	getc(stdin);
								//}

								//if(this->W_ext->x_ext(j+1, l+1) < this->SS->x(j+1, N1, l+1, N2))
								//{
								//	printf("%s(%d), %.15f\n%.15f\n%d, %d (%d, %d)\n", __FILE__, __LINE__, this->W_ext->x_ext(j+1, l+1), this->SS->x(j+1, N1, l+1, N2),j,l, N1, N2);
								//	getc(stdin);
								//}

								//V_mhe->x_ext(i,j,k,l) = MAX_SUM(V_mhe->x_ext(i,j,k,l), 
								//	MUL(this->SS->x(1, i-1, 1, k-1), this->SS->x(j+1, N1, l+1, N2)));

								V_mhe->calculate_ext_dependencies(i, j, k, l); // Arrays to use.
							}

#ifdef _LINEAR_COMPUTATIONS_
						}
#endif
					} // Loop limit if for internal and external checks on i-j and k-l.				
				} // k loop
			} // l loop
		} // i loop
	} // j loop

	delete(ppf_progress_bar);

	printf("\n");
}

void t_ppf_loops::rescale_all_external_ppf_arrays_by_indices(int i_lim, int j_lim, int k_lim, int l_lim, bool up_scale)
{
	int r_i = 0;
	int r_j = 0;
	int r_k = 0;
	int r_l = 0;

	int N1 = seq_man->get_l_seq1();
	int N2 = seq_man->get_l_seq2();

	int low_k, high_k, low_l, high_l;

	// Boundary conditions.
	//this->V_mhe->x_ext(1, N1, 1, N2) = CONVERT_FROM_LIN(1.0f);
	//this->W->x_ext(1, N1, 1, N2) = CONVERT_FROM_LIN(1.0f);
	//this->WMB->x_ext(1, N1, 1, N2) = CONVERT_FROM_LIN(1.0f);

	// Main loop calculation:
	for( int j = N1; j >= 1; j--)
	{
		//for( int i = j-1; i >= MAX(j - ppf_cli->max_n_separation_between_nucs, 1); i--)
		for( int i = MAX(j - ppf_cli->max_n_separation_between_nucs, 1); i <= j-1 ; i++)
		{
			bool ij_coinc = true;
			if(seq1_spf->folding_constraints->coinc_pointer_relocation_map != NULL && seq1_spf->folding_constraints->coinc_pointer_relocation_map[i][j] == POS_MEM_NOT_EXIST)
			{
				ij_coinc = false;
			}

			low_l = MAX( t_template_pf_array::low_limits[j], 1 );
			high_l = MIN( t_template_pf_array::high_limits[j], N2 );

			//for ( int l = low_l; ij_coinc && l <= high_l; l++ )
			for(int l = high_l; ij_coinc && l >= low_l; l--)
			{
				low_k = MAX(MAX(l - ppf_cli->max_n_separation_between_nucs, 1), t_template_pf_array::low_limits[i-1] + 1);
				high_k = MIN(l - 1, t_template_pf_array::high_limits[i - 1]+1);

				//for ( int k = high_k; k >= low_k; k-- )
				for ( int k = low_k; k <= high_k; k++ )
				{
					bool kl_coinc = true;
					if(seq2_spf->folding_constraints->coinc_pointer_relocation_map != NULL && seq2_spf->folding_constraints->coinc_pointer_relocation_map[k][l] == POS_MEM_NOT_EXIST)
					{
						kl_coinc = false;
					}

					// Make sure that there is at least MIN_LOOP of distance between i-j and k-l.
					if(kl_coinc)
					{

if(_DUMP_PPF_LOOPS_MESSAGES_)
{
						printf("----------------------------------------------\n");
						printf("Rescaling[ext] (%d, %d, %d, %d):\n", i, j, k, l);
}

						// Make sure that scaling factor 
						double factor_de_rescale = this->ppf_scaler->cumulative_rescale_external_factor(i,j,k,l);

						//WMB->calculate_ext_dependencies(i, j, k, l); // Arrays to use.
						if(WMB->pf_array->check_4D_ll(i,j,k,l))
						{
							scale_pp_value(WMB->x_ext(i,j,k,l), factor_de_rescale, up_scale);
						}

						//WMBL->calculate_ext_dependencies(i, j, k, l); // Arrays to use.
						if(WMBL->pf_array->check_4D_ll(i,j,k,l))
						{
							scale_pp_value(WMBL->x_ext(i,j,k,l), factor_de_rescale, up_scale);
						}

						//W->calculate_ext_dependencies(i, j, k, l); // Arrays to use.
						if(W->pf_array->check_4D_ll(i,j,k,l))
						{
							scale_pp_value(W->x_ext(i,j,k,l), factor_de_rescale, up_scale);
						}

						//WL->calculate_ext_dependencies(i, j, k, l); // Arrays to use.
						if(WL->pf_array->check_4D_ll(i,j,k,l))
						{
							scale_pp_value(WL->x_ext(i,j,k,l), factor_de_rescale, up_scale);
						}

						//W_mhi->calculate_kl_ext_dependencies(i, j, k, l); // Arrays to use. 
						if(W_mhi->ij_pf_array->check_4D_ll(i,j,k,l))
						{
							scale_pp_value(W_mhi->x_ext_ij(i,j,k,l), factor_de_rescale, up_scale);
						}

						//W_mhi->calculate_ij_ext_dependencies(i, j, k, l); // Arrays to use. 
						if(W_mhi->kl_pf_array->check_4D_ll(i,j,k,l))
						{
							scale_pp_value(W_mhi->x_ext_kl(i,j,k,l), factor_de_rescale, up_scale);
						}

						if( ( (j > N1 || (j-i) > MIN_LOOP) &&  (l>N2 || (l-k) > MIN_LOOP) ) 
							&& ( (l <= N2 || k - (l - N2) > MIN_LOOP) && (j <= N1 || i - (j - N1) > MIN_LOOP) ) // To impose constraint on the length of interior loop when exterior loop is calculated.
							)
						{
							// Initialize with Wext, this sets the boundary condition.
							//V_mhe->x_ext(i,j,k,l) = MAX_SUM(V_mhe->x_ext(i,j,k,l), 
							//	MUL(this->W_ext->x(i-1, k-1), this->W_ext->x_ext(j+1, l+1)));

							//V_mhe->calculate_ext_dependencies(i, j, k, l); // Arrays to use.

							if(V_mhe->pf_array->check_4D_ll(i,j,k,l))
							{
								scale_pp_value(V_mhe->x_ext(i,j,k,l), factor_de_rescale, up_scale);
							}
						}
					} // Loop limit if for internal and external checks on i-j and k-l.				
				} // k loop
			} // l loop
		} // i loop
	} // j loop

	this->W_ext->rescale_external(up_scale);

	printf("\n");
}

// Following loops are determined by ppf loops, there is also an map version of this rescaling loop.
// This loop rescales all arrays with ppf loops.
void t_ppf_loops::rescale_all_ppf_arrays_by_indices(int i_lim, int j_lim, int k_lim, int l_lim, bool up_scale)
{
	// Update rescaling factor.
	int N1 = this->seq_man->get_l_seq1();
	int N2 = this->seq_man->get_l_seq2();

	printf("Rescaling: %d, %d, %d, %d\n", i_lim, j_lim, k_lim, l_lim);
	
	int low_k, high_k, low_l, high_l;

if(_DUMP_PPF_LOOPS_MESSAGES_)
{
	printf("\n--------------------------------\nPPF RESCALING LOOP @ (%d, %d, %d, %d):\n", i_lim, j_lim, k_lim, l_lim);
}

	// Rescale SS array.
	this->SS->rescale(up_scale);

	// Main loop calculation:
	// Rescaling starts with whole limits and 
	for(int j = 1; j <= N1; j++)
	{
if(_DUMP_PPF_LOOPS_MESSAGES_)
{
		printf("j = %d(N1 = %d, 2xN1 = %d)\r", j, N1, 2*N1);
		fflush(stdout);
}

		for( int i = j-1; i >= MAX(j - ppf_cli->max_n_separation_between_nucs, 1); i--)
		{
			bool ij_coinc = true;
			if(seq1_spf->folding_constraints->coinc_pointer_relocation_map != NULL && seq1_spf->folding_constraints->coinc_pointer_relocation_map[i][j] == POS_MEM_NOT_EXIST)
			{
				ij_coinc = false;
			}

			low_l = MAX( t_template_pf_array::low_limits[j], 1 );
			high_l = MIN( t_template_pf_array::high_limits[j], N2 );

			for ( int l = low_l; ij_coinc && l <= high_l; l++ )
			{
				low_k = MAX(MAX(l - ppf_cli->max_n_separation_between_nucs, 1), t_template_pf_array::low_limits[i-1] + 1);
				high_k = MIN(l - 1, t_template_pf_array::high_limits[i - 1]+1);

				for ( int k = high_k; k >= low_k; k-- )
				{
					// If the loop limits are at limits, 
					// return from scaling.
					if(i == i_lim && j == j_lim && k == k_lim && l == l_lim)
					{
if(_DUMP_PPF_LOOPS_MESSAGES_)
{
						printf("Rescaling ok.\n--------------------------------\n");
}
						//return;
					}

					bool kl_coinc = true;
					if(seq2_spf->folding_constraints->coinc_pointer_relocation_map != NULL && seq2_spf->folding_constraints->coinc_pointer_relocation_map[k][l] == POS_MEM_NOT_EXIST)
					{
						kl_coinc = false;
					}

					// Make sure that there is at least MIN_LOOP of distance between i-j and k-l.
					if(kl_coinc)
					{
						double factor_de_rescale = this->ppf_scaler->cumulative_rescale_factor(i,j,k,l);

						if( ( (j > N1 || (j-i) > MIN_LOOP) &&  (l>N2 || (l-k) > MIN_LOOP) ) 
							&& ( (l <= N2 || k - (l - N2) > MIN_LOOP) && (j <= N1 || i - (j - N1) > MIN_LOOP) ) // To impose constraint on the length of interior loop when exterior loop is calculated.
							)
						{
							if(V_mhe->pf_array->check_4D_ll(i,j,k,l))
							{
								//printf("Before scaling Vmhe: %.15f\n", this->ppf_scaler->get_unscaled_log_value_by_indices(i,j,k,l, V_mhe->x(i,j,k,l)));
								scale_pp_value(V_mhe->x_setter(i,j,k,l), factor_de_rescale, up_scale);
								//printf("After scaling Vmhe: %.15f\n", this->ppf_scaler->get_unscaled_log_value_by_indices(i,j,k,l, V_mhe->x(i,j,k,l)));

								if(i == 2 && j == 72 && k == 2 && l == 71)
								{
									printf("After Rescaling 2, 72, 2, 71: %lf\n", this->ppf_scaler->get_unscaled_log_value_by_indices(i,j,k,l, V_mhe->x(i,j,k,l)));
								}
							}
						}

						if(W_mhi->ij_pf_array->check_4D_ll(i,j,k,l))
						{
							scale_pp_value(W_mhi->x_ij_setter(i,j,k,l), factor_de_rescale, up_scale);
						}

						if(W_mhi->kl_pf_array->check_4D_ll(i,j,k,l))
						{
							scale_pp_value(W_mhi->x_kl_setter(i,j,k,l), factor_de_rescale, up_scale);
						}

						if(WL->pf_array->check_4D_ll(i,j,k,l))
						{
							scale_pp_value(WL->x_setter(i,j,k,l), factor_de_rescale, up_scale); 
						}

						if(W->pf_array->check_4D_ll(i,j,k,l))
						{
							scale_pp_value(W->x_setter(i,j,k,l), factor_de_rescale, up_scale);
						}

						if(WMBL->pf_array->check_4D_ll(i,j,k,l))
						{
							scale_pp_value(WMBL->x_setter(i,j,k,l), factor_de_rescale, up_scale);
						}

						if(WMB->pf_array->check_4D_ll(i,j,k,l))
						{
							scale_pp_value(WMB->x_setter(i,j,k,l), factor_de_rescale, up_scale);
						}
					} // Loop limit if					
				} // l loop
			} // k loop
		} // i loop
	} // j loop

	// Rescale W_ext.
	this->W_ext->rescale_internal(up_scale);
}

// Following function scales a value correctly based on scaling flag.
// Divide or multiple with cumulative rescaling factor based on the flag.
void t_ppf_loops::scale_pp_value(double& value, double factor_de_scale, bool up_scale)
{
    if(up_scale)
    {
            value *= factor_de_scale;
    }
    else
    {
            value /= factor_de_scale;
    }
}

void* t_ppf_loops::backtrack_pseudo_free_energy_loops()
{
	//int i = 77;
	//int k = 76;
	//printf("Before sampling loops: %lf x %lf\n", this->seq1_spf->ux_3p(i-1, i), this->seq2_spf->ux_3p(k-1, k));
	//getc(stdin);

if(_DUMP_PPF_LOOPS_MESSAGES_)
	printf("Backtracking structures...\n");

	if(this->ppf_cli->mode == PARTS_RUN_MODE_MAP)
	{
		if(this->map_strs != NULL)
		{			
			delete(this->map_strs);
			this->map_strs = NULL;
		}

		if(this->map_alignment != NULL)
		{
			delete(this->map_alignment);
			this->map_alignment = NULL;
		}

		// Allocate map structure objects for traceback.
		this->map_strs = new t_MAP_structures(seq_man, ppf_cli);

		// Alignment information.
		this->map_alignment = new t_MAP_alignment(seq_man, ppf_cli);
	}
	else
	{
		if(this->stoch_sampler != NULL)
		{
			delete(this->stoch_sampler);
			this->stoch_sampler = NULL;
		}

		this->stoch_sampler = new t_stoch_sampling_math(this);

		// Create sampling set.
		if(this->stoch_sampled_set != NULL)
		{
			this->stoch_sampled_set = NULL;
		}
		this->stoch_sampled_set = new t_stoch_sampled_str_aln_sample_set(this);
	}

	// Do backtracking once for MAP computations.
	int n_backtracks = 1;
	if(this->ppf_cli->mode == PARTS_RUN_MODE_STOCH_SAMPLE)
	{
		n_backtracks = stoch_sampled_set->n_samples;
	}

	for(int cnt = 0; cnt < n_backtracks; cnt++)
	{
		// Allocate stack for traceback.
		this->tb_stack = new t_ppf_tb_stack();

		int N1 = this->seq_man->get_l_seq1();
		int N2 = this->seq_man->get_l_seq2();

		//_DUMP_PPF_LOOPS_MESSAGES_ = true;

		tb_stack->push_str(1, N1, 1, N2, TRACE_W_ext);

		// Indices to be pulled in traceback.
		int i, j, k, l;
		char array_id;

		tb_stack->pop_str(i, j, k, l, array_id);

		// Start traceback.
		while(array_id != TRACE_stack_empty)
		{
	if(_DUMP_PPF_LOOPS_MESSAGES_)
			printf("Pulled %s(%d, %d, %d, %d)\n", trace_names[array_id], i, j, k, l); 

			switch(array_id)
			{
			case TRACE_SS:
				SS->backtrack(i,j,k,l);
				break;

			case TRACE_W:
				W->calculate_W(i,j,k,l, true);
				break;

			case TRACE_WL:
				WL->calculate_WL(i, j, k, l, true);
				break;

			case TRACE_WMB:
				WMB->calculate_WMB(i,j,k,l, true);
				break;

			case TRACE_WMBL:
				WMBL->calculate_WMBL(i,j,k,l, true);
				break;

			case TRACE_V_mhe:
				V_mhe->calculate_V_mhe(i, j, k, l, true);
				break;

			case TRACE_W_ij_mhi:
				W_mhi->calculate_ij_W_mhi(i, j, k, l, true);
				break;

			case TRACE_W_kl_mhi:
				W_mhi->calculate_kl_W_mhi(i, j, k, l, true);
				break;

			case TRACE_W_ext:
				W_ext->backtrack(j, l);
				break;

			default:
				printf("Defaulted with %d\n", array_id);
				exit(0);
			};

			// Pop a new partial structure to trace.
			tb_stack->pop_str(i, j, k, l, array_id);
		};

		// Dump structures and alignments.
		if(this->ppf_cli->mode == PARTS_RUN_MODE_STOCH_SAMPLE)
		{
			stoch_sampled_set->dump_strs_aln(cnt);
			stoch_sampled_set->reset_strs_aln();
		}

		delete(this->tb_stack);
		this->tb_stack = NULL;
	} // str sample size counter.

	if(this->ppf_cli->mode == PARTS_RUN_MODE_MAP)
	{
	}
	else
	{
		delete(stoch_sampled_set);
		this->stoch_sampled_set = NULL;
		delete(this->stoch_sampler);
		this->stoch_sampler = NULL;
	}

	// Dump structures.
	if(this->ppf_cli->mode == PARTS_RUN_MODE_MAP)
	{
		// Dump ct's
		map_strs->dump_map_cts();

		// Dump alignment.
		map_alignment->dump_map_alignment();
	}

	//// Determine mhr info.
	//t_MAP_mhr_info* mhr_info = new t_MAP_mhr_info(ppf_cli, map_strs, map_alignment);

	//t_map_results* map_results = new t_map_results(seq_man, ppf_cli, map_strs, map_alignment, mhr_info);
	//return(map_results);

	printf("\n");

	return(NULL);
}

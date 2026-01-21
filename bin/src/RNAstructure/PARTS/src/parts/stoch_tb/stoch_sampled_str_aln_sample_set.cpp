#include <string.h>
#include <limits.h>
#include "../parts_compilation_directives.h"
#include "../parts_paths.h"
#include "stoch_sampled_str_aln_sample_set.h"
#include "stoch_sampled_structures.h"
#include "stoch_sampled_alignment.h"
#include <stdio.h>
#include <stdlib.h>
#include "../ppf_cli.h"
#include "../ppf_math.h"
#include "../process_sequences.h"
#include "../phmm_parameters.h" 
#include "../alignment_priors.h" 
#include "../single_pf_array.h" 
#include "../../../../src/phmm/utils/xmath/log/xlog_math.h" 
#include "../ppf_loops.h" 
#include <math.h>

t_stoch_sampled_str_aln_sample_set::t_stoch_sampled_str_aln_sample_set(t_ppf_loops* _ppf_loops)
{
	this->ppf_loops = _ppf_loops;

	this->seq_man = this->ppf_loops->seq_man;
	this->ppf_cli = this->ppf_loops->ppf_cli;
	this->n_samples = this->ppf_loops->ppf_cli->n_samples; 

	this->seq1_spf = this->ppf_loops->seq1_spf;
	this->seq2_spf = this->ppf_loops->seq2_spf;
	this->aln_priors = this->ppf_loops->aln_priors;

	// Trace structural alignments.
	//this->sampled_structures = new t_stoch_sampled_structures(this->seq_man, this->ppf_cli);
	this->sampled_structures = new t_stoch_sampled_structures(seq_man, ppf_cli);
	this->sampled_alignment = new t_stoch_sampled_alignment(seq_man, ppf_cli);
}

// Free this object.
t_stoch_sampled_str_aln_sample_set::~t_stoch_sampled_str_aln_sample_set()
{
	delete(this->sampled_structures);
	delete(this->sampled_alignment);
}

void t_stoch_sampled_str_aln_sample_set::dump_strs_aln(int _sample_id)
{	
	// Create individual structure file names and paths.
	this->sampled_structures->dump_sampled_cts(_sample_id);

	// Generate file path for alignment sample and dump it.
	this->sampled_alignment->dump_sampled_alignment(_sample_id);
}

t_stoch_sampled_structures* t_stoch_sampled_str_aln_sample_set::current_sampled_structures()
{
	return(this->sampled_structures);
}

t_stoch_sampled_alignment* t_stoch_sampled_str_aln_sample_set::current_sampled_alignment()
{
	return(this->sampled_alignment);
}

void t_stoch_sampled_str_aln_sample_set::reset_strs_aln()
{
	this->sampled_structures->reset_bps();
	this->sampled_alignment->reset_aln();
}




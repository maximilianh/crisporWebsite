#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "parts_compilation_directives.h"
#include "array_file_manager.h"

#include "ppf_w.h"
#include "ppf_w_mb.h"


#include "process_sequences.h"
#include "ppf_math.h"
#include "single_pf_array.h"
#include "phmm_parameters.h"
#include "alignment_priors.h"

#include "template_pf_array.h"
#include "ppf_ss.h"

#include "ppf_w_mhi.h"
#include "ppf_cli.h"

#include "../../../src/phmm/utils/file/utils.h"

t_array_file_manager::t_array_file_manager(t_seq_man* _seq_man)
{
	this->seq_man = _seq_man;

	// Open the input files.
	this->in_array_files = new vector<FILE*>();
	for(int i_in_file = V_FILE; this->seq_man->ppf_cli->use_array_files && i_in_file < N_ARRAYS_LOGGED; i_in_file++)
	{
		char cur_array_fp[1000];
		sprintf(cur_array_fp, "%s/%s", ARRAY_IP_DIR, array_file_names[i_in_file]);
		FILE* cur_file = open_f(cur_array_fp, "r");
		if(cur_file == NULL)
		{
			printf("Could not open %s for reading.\n", cur_array_fp);
			exit(0);
		}
		this->in_array_files->push_back(cur_file);
	}

	// Initialize the files to write
	this->out_array_files = new vector<FILE*>();
	for(int i_out_file = V_FILE; this->seq_man->ppf_cli->save_array_files && i_out_file < N_ARRAYS_LOGGED; i_out_file++)
	{
		char cur_array_fp[1000];
		sprintf(cur_array_fp, "%s/%s", ARRAY_OP_DIR, array_file_names[i_out_file]);
		FILE* cur_file = open_f(cur_array_fp, "w");
		if(cur_file == NULL)
		{
			printf("Could not open %s for writing.\n", cur_array_fp);
		}

		this->out_array_files->push_back(cur_file);
	}
}

t_array_file_manager::~t_array_file_manager()
{
}

bool t_array_file_manager::read_array_value(int in_array, double& val)
{
	if(in_array < this->in_array_files->size() && 
		this->in_array_files->at(in_array) == NULL)
	{
		return(false);
	}

	double read_val = 0.0f;
	if(!fscanf(this->in_array_files->at(in_array), "%lf", &val) != 1)
	{
		return(false);
	}

	val = read_val;

	return(true);
}

bool t_array_file_manager::write_array_value(int out_array, double val)
{
	if(out_array < this->out_array_files->size() && 
		this->out_array_files->at(out_array) == NULL)
	{
		return(false);
	}

	if(fprintf(this->out_array_files->at(out_array), "%.10f ", val) > 0)
	{
		return(true);
	}
	else
	{
		return(false);
	}
}

#include <string.h>
#include <limits.h>
#include "process_sequences.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../../src/phmm/structure/structure_object.h"
#include "../../../src/phmm/utils/file/utils.h"
#include "ppf_cli.h"

bool _DUMP_PROCESS_SEQUENCES_MESSAGES_ = false;

t_seq_man::t_seq_man(t_ppf_cli* _ppf_cli)
{
	this->ppf_cli = _ppf_cli;

	// Check 1st sequence file.
	FILE* test_seq_file = open_f(this->ppf_cli->seq1_path, "r");
	if(test_seq_file == NULL)
	{
		printf("Could not open sequence file %s @ %s(%d), make sure it exists.\n", this->ppf_cli->seq1_path, __FILE__, __LINE__);
		exit(0);
	}
	else
	{
		fclose(test_seq_file);
	}

	// Check 2nd sequence file.
	test_seq_file = open_f(this->ppf_cli->seq2_path, "r");
	if(test_seq_file == NULL)
	{
		printf("Could not open sequence file %s @ %s(%d), make sure it exists.\n", this->ppf_cli->seq2_path, __FILE__, __LINE__);
		exit(0);
	}
	else
	{
		fclose(test_seq_file);
	}


	// Copy file names to file paths.
	strcpy(this->seq1_fp, this->ppf_cli->seq1_path);
	strcpy(this->seq2_fp, this->ppf_cli->seq2_path);

	// Allocate and instantiate structure objects.
	this->seq1 = new t_structure(this->seq1_fp);
	this->seq2 = new t_structure(this->seq2_fp);
/*
	// Load sequences.
	if(openseq(this->seq1, seq1_fn) == 0)
	{
		printf("Could not open sequence file %s @ %s(%d), make sure it exists.\n", seq1_fn, __FILE__, __LINE__);
		exit(0);
	}

	if(openseq(this->seq2, seq2_fn) == 0)
	{
		printf("Could not open sequence file %s @ %s(%d), make sure it exists.\n", seq2_fn, __FILE__, __LINE__);
		exit(0);
	}
	*/
}

int t_seq_man::get_l_seq1()
{
	return(this->seq1->numofbases);
}

int t_seq_man::get_l_seq2()
{
	return(this->seq2->numofbases);
}

char t_seq_man::get_nuc_seq1(int index)
{
	return(this->seq1->numseq[index]);
}

char t_seq_man::get_nuc_seq2(int index)
{
	return(this->seq2->numseq[index]);
}

t_seq_man::~t_seq_man()
{
	delete(this->seq1);
	delete(this->seq2);
}

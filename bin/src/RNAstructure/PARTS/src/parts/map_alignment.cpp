#include <string.h>
#include <limits.h>
#include "map_alignment.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "process_sequences.h"
#include "ppf_cli.h"
#include <math.h>
// For state names and enumerations, those should only be made in phmm_parameters file.
#include "phmm_parameters.h" 

#include "parts_paths.h"

#include "../../../src/phmm/utils/file/utils.h"

#define max(a,b) (((a)>(b))?(a):(b))

bool _DUMP_MAP_ALIGNMENT_MESSAGES_ = false;

t_MAP_alignment::t_MAP_alignment(t_seq_man* _seq_man, t_ppf_cli* _ppf_cli)
{
	this->seq_man = _seq_man;
	this->ppf_cli = _ppf_cli;

	// Allocate alignment arrays.
	this->seq1_alns = (int**)malloc(sizeof(int*) * (seq_man->get_l_seq1() + 2));
	this->seq2_alns = (int**)malloc(sizeof(int*) * (seq_man->get_l_seq2() + 2));

	for(int cnt = 0; cnt <= seq_man->get_l_seq1(); cnt++)
	{
		this->seq1_alns[cnt] = (int*)malloc(sizeof(int) * 2);
		this->seq1_alns[cnt][0] = 0;
		this->seq1_alns[cnt][1] = STATE_ALN;
	}

	for(int cnt = 0; cnt <= seq_man->get_l_seq2(); cnt++)
	{
		this->seq2_alns[cnt] = (int*)malloc(sizeof(int) * 2);
		this->seq2_alns[cnt][0] = 0;
		this->seq2_alns[cnt][1] = STATE_ALN;
	}

	this->aln_str1 = NULL;
	this->aln_str2 = NULL;

	this->aln_index_line1 = NULL;
	this->aln_index_line2 = NULL;

	this->l_aln = 0;
}

t_MAP_alignment::~t_MAP_alignment()
{
	for(int cnt = 0; cnt <= seq_man->get_l_seq1(); cnt++)
	{
		free(this->seq1_alns[cnt]);
	}

	for(int cnt = 0; cnt <= seq_man->get_l_seq2(); cnt++)
	{
		free(this->seq2_alns[cnt]);
	}
	free(this->seq1_alns);
	free(this->seq2_alns);

	if(this->aln_str1 != NULL)
	{
		free(this->aln_str1);
	}

	if(this->aln_str1 != NULL)
	{
		free(this->aln_str2);
	}

	if(this->aln_index_line1 != NULL)
	{
		free(this->aln_index_line1);
	}

	if(this->aln_index_line2 != NULL)
	{
		free(this->aln_index_line2);
	}
}

// Following function are for setting alignment values.
void t_MAP_alignment::set_seq1_ins(int seq1_index, int seq2_index)
{
	this->seq1_alns[seq1_index][0] = seq2_index;
	this->seq1_alns[seq1_index][1] = STATE_INS1;
}

void t_MAP_alignment::set_seq2_ins(int seq1_index, int seq2_index)
{
	this->seq2_alns[seq2_index][0] = seq1_index;
	this->seq2_alns[seq2_index][1] = STATE_INS2;
}

void t_MAP_alignment::set_aln(int seq1_index, int seq2_index)
{
	this->seq1_alns[seq1_index][0] = seq2_index;
	this->seq1_alns[seq1_index][1] = STATE_ALN;

	this->seq2_alns[seq2_index][0] = seq1_index;
	this->seq2_alns[seq2_index][1] = STATE_ALN;
}

int t_MAP_alignment::get_l_aln()
{
	if(this->l_aln != 0)
	{
		return(this->l_aln);
	}
	else
	{
		int aln_str_index = 0;

		// Following points to last alignment position in coincidence map.
		int last_i1 = 0;
		int last_i2 = 0;

		while(last_i1 != seq_man->get_l_seq1() || 
			last_i2 != seq_man->get_l_seq2())
		{
			//printf("%d(%d), %d(%d)\n", last_i1, seq_man->get_l_seq1(), last_i2, seq_man->get_l_seq2());

			// Check for alignment case.
			if(last_i1 != seq_man->get_l_seq1() && 
				last_i2 != seq_man->get_l_seq2() &&
				this->seq1_alns[last_i1 + 1][1] == STATE_ALN && 
				this->seq1_alns[last_i1 + 1][0] == last_i2 + 1)
			{
				// If next nuc. in sequence 1 is aligned, is it aligned to
				// next nuc. in sequence 2?
				last_i1++;
				last_i2++;

				aln_str_index++;

	if(_DUMP_MAP_ALIGNMENT_MESSAGES_)
				printf("Align %d, %d\n", last_i1, last_i2);
			}
			// Check for alignment case.
			else if(last_i1 != seq_man->get_l_seq1() && 
					this->seq1_alns[last_i1 + 1][1] == STATE_INS1 && 
					this->seq1_alns[last_i1 + 1][0] == last_i2)
			{
				// If next nuc. in sequence 1 is inserted, is it inserted on top of current nuc. in sequence 2?
				last_i1++;

				aln_str_index++;

	if(_DUMP_MAP_ALIGNMENT_MESSAGES_)
				printf("Insert1 %d, %d\n", last_i1, last_i2);

			}
			// Check for alignment case.
			else if(last_i2 != seq_man->get_l_seq2() && 
					this->seq2_alns[last_i2 + 1][1] == STATE_INS2 && 
					this->seq2_alns[last_i2 + 1][0] == last_i1)
			{
				// If next nuc. in sequence 2 is inserted, is it inserted on top of current nuc. in sequence 1?
				last_i2++;

				aln_str_index++;

	if(_DUMP_MAP_ALIGNMENT_MESSAGES_)
				printf("Insert2 %d, %d\n", last_i1, last_i2);
			}
		} // map alignment string formation loop.

		this->l_aln = aln_str_index;
		return(aln_str_index);
	}
}

// Dump map alignment.
void t_MAP_alignment::dump_map_alignment()
{
if(_DUMP_MAP_ALIGNMENT_MESSAGES_)
{
	FILE* map_aln_file = open_f("ppf_map_alignment.txt", "w");
	for(int cnt1 = 1; cnt1 <= this->seq_man->get_l_seq1(); cnt1++)
	{
		fprintf(map_aln_file, "%d %d %s\n", cnt1, this->seq1_alns[cnt1][0], state_names[this->seq1_alns[cnt1][1]]);
	}

	fprintf(map_aln_file, "\n\n");

	for(int cnt2 = 1; cnt2 <= this->seq_man->get_l_seq2(); cnt2++)
	{
		fprintf(map_aln_file, "%d %d %s\n", cnt2, this->seq2_alns[cnt2][0], state_names[this->seq2_alns[cnt2][1]]);
	}

	fclose(map_aln_file);
}

	// Both alignment arrays correspond to same coincidence path.
	// In order to represent those alignment arrays, have to trace them correctly
	// into alignment strings.
	//int aln_str_length = seq_man->get_l_seq1() + seq_man->get_l_seq2();
	int l_aln = this->get_l_aln();
	this->aln_str1 = (char*)malloc(sizeof(char) * (l_aln + 2));
	this->aln_str2 = (char*)malloc(sizeof(char) * (l_aln + 2));
	
	this->aln_index_line1 = (int*)malloc(sizeof(int) * (l_aln + 3));
	this->aln_index_line2 = (int*)malloc(sizeof(int) * (l_aln + 3));

	// Following points to last alignment position in coincidence map.
	int last_i1 = 0;
	int last_i2 = 0;

	char nucs[] = "NACGUI";

	// Problem is determining if next state is an event in first seq (ins1) or an event in second sequence (ins2)
	// or if it is an event in both sequences (aln). So check seq1_alns[last_i1 + 1] and seq1_alns[last_i2 + 1]
	// indices and states; see if the state and indices are corectly adding up on last_i1 and last_i2.
	// e.g. if seq1_alns[1][0] = 0 and seq1_alns[1][1] = STATE_INS1, then this means that there is an insertion
	// in first sequence which will be over 0, 0. 
	int aln_str_index = 0;
	while(last_i1 != seq_man->get_l_seq1() || 
		last_i2 != seq_man->get_l_seq2())
	{
if(_DUMP_MAP_ALIGNMENT_MESSAGES_)
		printf("%d(%d), %d(%d)\n", last_i1, seq_man->get_l_seq1(), last_i2, seq_man->get_l_seq2());

		// Check for alignment case.
		if((last_i1+1) <= seq_man->get_l_seq1() && 
			(last_i2+1) <= seq_man->get_l_seq2() &&
			this->seq1_alns[last_i1 + 1][1] == STATE_ALN && 
			this->seq1_alns[last_i1 + 1][0] == last_i2 + 1)
		{
			// If next nuc. in sequence 1 is aligned, is it aligned to
			// next nuc. in sequence 2?
			last_i1++;
			last_i2++;

			aln_str1[aln_str_index] = nucs[this->seq_man->get_nuc_seq1(last_i1)];
			aln_str2[aln_str_index] = nucs[this->seq_man->get_nuc_seq2(last_i2)];

			this->aln_index_line1[aln_str_index+1] = last_i1;
			this->aln_index_line2[aln_str_index+1] = last_i2;

			aln_str_index++;

if(_DUMP_MAP_ALIGNMENT_MESSAGES_)
			printf("Align %d, %d\n", last_i1, last_i2);
		}
		// Check for alignment case.
		else if((last_i1+1) <= seq_man->get_l_seq1() && 
				this->seq1_alns[last_i1 + 1][1] == STATE_INS1 && 
				this->seq1_alns[last_i1 + 1][0] == last_i2)
		{
			// If next nuc. in sequence 1 is inserted, is it inserted on top of current nuc. in sequence 2?
			last_i1++;

			aln_str1[aln_str_index] = nucs[this->seq_man->get_nuc_seq1(last_i1)];
			aln_str2[aln_str_index] = '.';

			this->aln_index_line1[aln_str_index+1] = last_i1;
			this->aln_index_line2[aln_str_index+1] = 0;

			aln_str_index++;

if(_DUMP_MAP_ALIGNMENT_MESSAGES_)
			printf("Insert1 %d, %d\n", last_i1, last_i2);

		}
		// Check for alignment case.
		else if((last_i2+1) <= seq_man->get_l_seq2() && 
				this->seq2_alns[last_i2 + 1][1] == STATE_INS2 && 
				this->seq2_alns[last_i2 + 1][0] == last_i1)
		{
			// If next nuc. in sequence 2 is inserted, is it inserted on top of current nuc. in sequence 1?
			last_i2++;

			aln_str1[aln_str_index] = '.';
			aln_str2[aln_str_index] = nucs[this->seq_man->get_nuc_seq2(last_i2)];

			this->aln_index_line1[aln_str_index+1] = 0;
			this->aln_index_line2[aln_str_index+1] = last_i2;

			aln_str_index++;

if(_DUMP_MAP_ALIGNMENT_MESSAGES_)
			printf("Insert2 %d, %d\n", last_i1, last_i2);
		}
		else
		{
			printf("Could not decode next coincidence position in alignment at (%d, %d) @ %s(%d).\n", last_i1, last_i2, __FILE__, __LINE__);
			exit(0);
		}
	} // map alignment string formation loop.

	// Finish alignment strings.
	aln_str1[aln_str_index] = 0;
	aln_str2[aln_str_index] = 0;

if(_DUMP_MAP_ALIGNMENT_MESSAGES_)
	printf("MAP Alignment:\n%s\n%s\n", aln_str1, aln_str2);

	char aln_fp[4096];
	if(this->ppf_cli->map_aln_op == NULL)
	{
		sprintf(aln_fp, "%s_%s_map_aln.aln", this->ppf_cli->seq1_op_file_prefix, this->ppf_cli->seq2_op_file_prefix);
	}
	else
	{
		strcpy(aln_fp, this->ppf_cli->map_aln_op);
	}
	FILE* aln_file = open_f(aln_fp, "w");
	fprintf(aln_file, "%s-%s MAP Alignment:\n%s\n%s \n", this->ppf_cli->seq1_op_file_prefix, this->ppf_cli->seq2_op_file_prefix, aln_str1, aln_str2);
	fclose(aln_file);
}


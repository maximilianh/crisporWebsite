#include <string.h>
#include <limits.h>
#include "stoch_sampled_alignment.h"
#include <stdio.h>
#include <stdlib.h>
#include "../process_sequences.h"
#include "../ppf_cli.h"
#include <math.h>
#include <string.h>
// For state names and enumerations, those should only be made in phmm_parameters file.
#include "../phmm_parameters.h" 

#include "../parts_paths.h"

#include <vector>
#include <algorithm>

#include "../../../../src/phmm/utils/file/utils.h"

using namespace std;

#define max(a,b) (((a)>(b))?(a):(b))

bool _DUMP_STOCH_SAMPLED_ALIGNMENT_MESSAGES_ = false;

t_stoch_sampled_alignment::t_stoch_sampled_alignment(t_seq_man* seq_man, t_ppf_cli* ppf_cli)
{
	this->seq_man = seq_man;
	this->ppf_cli = ppf_cli;

	// Allocate alignment arrays.
	//this->seq1_alns = (int**)malloc(sizeof(int*) * (seq_man->get_l_seq1() + 2));
	//this->seq2_alns = (int**)malloc(sizeof(int*) * (seq_man->get_l_seq2() + 2));

	this->alignment_path = new vector<t_aln_site*>();

	//for(int cnt = 0; cnt <= seq_man->get_l_seq1(); cnt++)
	//{
	//	this->seq1_alns[cnt] = (int*)malloc(sizeof(int) * 2);

	//	if(cnt <= seq_man->get_l_seq1())
	//	{
	//		this->seq1_alns[cnt][0] = 0;
	//		this->seq1_alns[cnt][1] = STATE_ALN;
	//	}
	//}

	//for(int cnt = 0; cnt <= seq_man->get_l_seq2() + 1; cnt++)
	//{
	//	this->seq2_alns[cnt] = (int*)malloc(sizeof(int) * 2);
	//	this->seq2_alns[cnt][0] = 0;
	//	this->seq2_alns[cnt][1] = STATE_ALN;
	//}
}

t_stoch_sampled_alignment::~t_stoch_sampled_alignment()
{
	for(int i = 0; i < this->alignment_path->size(); i++)
	{
		free(this->alignment_path->at(i));
	}
	this->alignment_path->clear();
	delete(this->alignment_path);
}

void t_stoch_sampled_alignment::reset_aln()
{
	for(int i_reset = 0; i_reset < this->alignment_path->size(); i_reset++)
	{
		free(this->alignment_path->at(i_reset));
	}

	this->alignment_path->clear();
}

// Following function are for setting alignment values.
void t_stoch_sampled_alignment::set_seq1_ins(int seq1_index, int seq2_index)
{
	t_aln_site* aln_site = (t_aln_site*)malloc(sizeof(t_aln_site));
	aln_site->i1 = seq1_index;
	aln_site->i2 = seq2_index;
	this->alignment_path->push_back(aln_site);

if(_DUMP_STOCH_SAMPLED_ALIGNMENT_MESSAGES_)
	printf("Adding STATE_INS1(%d, %d)\n", seq1_index, seq2_index);
}

void t_stoch_sampled_alignment::set_seq2_ins(int seq1_index, int seq2_index)
{
	t_aln_site* aln_site = (t_aln_site*)malloc(sizeof(t_aln_site));
	aln_site->i1 = seq1_index;
	aln_site->i2 = seq2_index;
	this->alignment_path->push_back(aln_site);


if(_DUMP_STOCH_SAMPLED_ALIGNMENT_MESSAGES_)
	printf("Adding STATE_INS2(%d, %d)\n", seq1_index, seq2_index);
}

void t_stoch_sampled_alignment::set_aln(int seq1_index, int seq2_index)
{
	t_aln_site* aln_site = (t_aln_site*)malloc(sizeof(t_aln_site));
	aln_site->i1 = seq1_index;
	aln_site->i2 = seq2_index;
	this->alignment_path->push_back(aln_site);


if(_DUMP_STOCH_SAMPLED_ALIGNMENT_MESSAGES_)
	printf("Adding STATE_ALN(%d, %d)\n", seq1_index, seq2_index);
}

// Dump map alignment.
void t_stoch_sampled_alignment::dump_sampled_alignment(int _sample_id)
{
	// Both alignment arrays correspond to same coincidence path.
	// In order to represent those alignment arrays, have to trace them correctly
	// into alignment strings.
	int aln_str_length = this->alignment_path->size() + 3;
	this->aln_str1 = (char*)malloc(sizeof(char) * aln_str_length);
	this->aln_str2 = (char*)malloc(sizeof(char) * aln_str_length);

	char nucs[] = "NACGUI";

	// Sort the alignment paths.
	vector<t_aln_site*> aln_path = *this->alignment_path;
	sort(aln_path.begin(), aln_path.end(), cmp_aln_site);
	//for(vector<t_aln_site*>::iterator i = aln_path.begin(); i < aln_path.end(); i++)
	//{
	//	t_aln_site* cur_aln_site = (t_aln_site*)*i;
	//	printf("%d %d\n", cur_aln_site->i1, cur_aln_site->i2);
	//}
	//getc(stdin);

	if(aln_path.back()->i1 == this->seq_man->get_l_seq1() && aln_path.back()->i2 == this->seq_man->get_l_seq2())
	{
	}
	else
	{
		printf("The alignment path returned from stochastic sampling is not complete @ %s(%d)\n", __FILE__, __LINE__);
		exit(0);
	}

	int last_i1 = 0;
	int last_i2 = 0;

	// Trace all the path.
	for(int aln_str_index = 0; aln_str_index < aln_path.size(); aln_str_index++)
	{
		int new_i1 = aln_path[aln_str_index]->i1;
		int new_i2 = aln_path[aln_str_index]->i2;

		// Check sequence 1 alignment state.
		if(new_i1 == last_i1 + 1)
		{
			aln_str1[aln_str_index] = nucs[this->seq_man->get_nuc_seq1(new_i1)];
		}
		else if(new_i1 == last_i1)
		{
			aln_str1[aln_str_index] = '.';

			if(new_i2 != last_i2 + 1)
			{
				printf("Alignment path not connected @ %s(%d)\n", __FILE__, __LINE__);
				exit(0);
			}
		}
		else
		{
			printf("Alignment path not connected @ %s(%d)\n", __FILE__, __LINE__);
			exit(0);
		}

		// Check sequence 2 alignment state.
		if(new_i2 == last_i2 + 1)
		{
			aln_str2[aln_str_index] = nucs[this->seq_man->get_nuc_seq2(new_i2)];
		}
		else if(new_i2 == last_i2)
		{
			aln_str2[aln_str_index] = '.';

			if(new_i1 != last_i1 + 1)
			{
				printf("Alignment path not connected @ %s(%d)\n", __FILE__, __LINE__);
				exit(0);
			}
		}
		else
		{
			printf("Alignment path not connected @ %s(%d)\n", __FILE__, __LINE__);
			exit(0);
		}

		// Update last indices.
		last_i1 = new_i1;
		last_i2 = new_i2;
	}

	// Finish alignment strings.
	aln_str1[aln_path.size()] = 0;
	aln_str2[aln_path.size()] = 0;

if(_DUMP_STOCH_SAMPLED_ALIGNMENT_MESSAGES_)
	printf("Sampled Alignment:\n%s\n%s\n", aln_str1, aln_str2);

	char aln_fp[4096];
	if(this->ppf_cli->sample_aln_op == NULL)
	{
		sprintf(aln_fp, "%s_%s_sampled_alns.aln", this->ppf_cli->seq1_op_file_prefix, this->ppf_cli->seq2_op_file_prefix);
	}
	else
	{
		strcpy(aln_fp, this->ppf_cli->sample_aln_op);
	}

	FILE* aln_file = open_f(aln_fp, "a");
	fprintf(aln_file, "%s-%s Sampled Alignment %d:\n%s\n%s \n", this->ppf_cli->seq1_op_file_prefix, this->ppf_cli->seq2_op_file_prefix, _sample_id, aln_str1, aln_str2);
	fclose(aln_file);

	free(this->aln_str1);
	free(this->aln_str2);
}

//void t_stoch_sampled_alignment::dump_sampled_alignment(double probability, char* extra_id)
//{
//	// Both alignment arrays correspond to same coincidence path.
//	// In order to represent those alignment arrays, have to trace them correctly
//	// into alignment strings.
//	int aln_str_length = this->alignment_path->size() + 3;
//	this->aln_str1 = (char*)malloc(sizeof(char) * aln_str_length);
//	this->aln_str2 = (char*)malloc(sizeof(char) * aln_str_length);
//
//	char nucs[] = "NACGUI";
//
//	// Sort the alignment paths.
//	vector<t_aln_site*> aln_path = *this->alignment_path;
//	sort(aln_path.begin(), aln_path.end(), cmp_aln_site);
//	for(int i = 0; i < aln_path.size(); i++)
//	{
//		printf("%d %d\n", aln_path[i]->i1, aln_path[i]->i2);
//	}
//
//	if(aln_path.back()->i1 == this->seq_man->get_l_seq1() && aln_path.back()->i2 == this->seq_man->get_l_seq2())
//	{
//	}
//	else
//	{
//		printf("The alignment path returned from stochastic sampling is not complete @ %s(%d)\n", __FILE__, __LINE__);
//		exit(0);
//	}
//
//	int last_i1 = 0;
//	int last_i2 = 0;
//
//	// Trace all the path.
//	for(int aln_str_index = 0; aln_str_index < aln_path.size(); aln_str_index++)
//	{
//		int new_i1 = aln_path[aln_str_index]->i1;
//		int new_i2 = aln_path[aln_str_index]->i2;
//
//		// Check sequence 1 alignment state.
//		if(new_i1 == last_i1 + 1)
//		{
//			aln_str1[aln_str_index] = nucs[this->seq_man->get_nuc_seq1(new_i1)];
//		}
//		else if(new_i1 == last_i1)
//		{
//			aln_str1[aln_str_index] = '.';
//
//			if(new_i2 != last_i2 + 1)
//			{
//				printf("Alignment path not connected @ %s(%d)\n", __FILE__, __LINE__);
//				exit(0);
//			}
//		}
//		else
//		{
//			printf("Alignment path not connected @ %s(%d)\n", __FILE__, __LINE__);
//			exit(0);
//		}
//
//		// Check sequence 2 alignment state.
//		if(new_i2 == last_i2 + 1)
//		{
//			aln_str2[aln_str_index] = nucs[this->seq_man->get_nuc_seq2(new_i2)];
//		}
//		else if(new_i2 == last_i2)
//		{
//			aln_str2[aln_str_index] = '.';
//
//			if(new_i1 != last_i1 + 1)
//			{
//				printf("Alignment path not connected @ %s(%d)\n", __FILE__, __LINE__);
//				exit(0);
//			}
//		}
//		else
//		{
//			printf("Alignment path not connected @ %s(%d)\n", __FILE__, __LINE__);
//			exit(0);
//		}
//
//		// Update last indices.
//		last_i1 = new_i1;
//		last_i2 = new_i2;
//	}
//
//	// Finish alignment strings.
//        aln_str1[aln_path.size()] = 0;
//        aln_str2[aln_path.size()] = 0;
//
//if(_DUMP_STOCH_SAMPLED_ALIGNMENT_MESSAGES_)
//	printf("Sampled Alignment:\n%s\n%s\n", aln_str1, aln_str2);
//
//	char sampled_alns_dir[500];
//	char sampled_aln_fn[100];
//	char sampled_aln_fp[500];
//
//	if(strcmp(extra_id, "max_prob") != 0)
//	{
//		concat_path_fn(sampled_alns_dir, SAMPLED_ALNS_DIR, ppf_cli->aln_id);
//		sprintf(sampled_aln_fn, "sample_%d_%s_aln.txt", sample_id, extra_id);
//		concat_path_fn(sampled_aln_fp, sampled_alns_dir, sampled_aln_fn);
//	}	
//	else
//	{
//		// Write max prob ct's to map_cts directory.
//		sprintf(sampled_aln_fn, "%s_max_prob_aln.txt", this->ppf_cli->aln_id);
//		concat_path_fn(sampled_aln_fp, SAMPLED_ALNS_DIR, sampled_aln_fn);
//	}
//
//	FILE* aln_file = open_f(sampled_aln_fp, "w");
//	fprintf(aln_file, "%s-%s Sampled Alignment with Probability=%.5f:\n%s\n%s \n", ppf_cli->seq1_path, ppf_cli->seq2_path, probability, aln_str1, aln_str2);
//	fclose(aln_file);
//}

// For sorting alignment paths.
bool cmp_aln_site(t_aln_site* aln_site1, t_aln_site* aln_site2)
{
	if(aln_site1->i1 < aln_site2->i1)
	{
		return(true);
	}
	else if(aln_site1->i1 > aln_site2->i1)
	{
		return(false);
	}
	else if(aln_site1->i1 == aln_site2->i1 && aln_site1->i2 < aln_site2->i2)
	{
		return(true);
	}
	else if(aln_site1->i1 == aln_site2->i1 && aln_site1->i2 > aln_site2->i2)
	{
		return(false);
	}
}


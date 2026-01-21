#include <string.h>
#include <limits.h>
#include "../map_structures.h"
#include "../process_sequences.h"
#include <stdio.h>
#include <stdlib.h>
#include "../../../../src/phmm/structure/structure_object.h"
#include "../../../../src/phmm/utils/file/utils.h"
#include "../ppf_cli.h"
#include "../parts_paths.h"
#include "stoch_sampled_structures.h"
#include "../ppf_math.h"
#include <string.h>

/*
class t_stoch_sampled_structures
{
public:
	t_stoch_sampled_structures(t_seq_man* seq_man, t_ppf_cli* ppf_cli);
	~t_stoch_sampled_structures();

	int* seq1_sampled_ct_bps;
	int* seq2_sampled_ct_bps;
	int N1;
	int N2;
	t_seq_man* seq_man;
	t_ppf_cli* ppf_cli;

	void add_bp_sampled_ct1(int i, int j);
	void add_bp_sampled_ct2(int k, int l);

	void dump_sampled_cts();

};
*/

bool _DUMP_STOCH_SAMPLED_STRUCTURES_MESSAGES_ = false;

t_stoch_sampled_structures::t_stoch_sampled_structures(t_seq_man* seq_man, t_ppf_cli* ppf_cli)
{
	this->seq_man = seq_man;
	this->ppf_cli = ppf_cli;

	// Copy sequence lengths.
	this->N1 = seq_man->get_l_seq1();
	this->N2 = seq_man->get_l_seq2();

	// Allocate base pair arrays.
	this->seq1_sampled_ct_bps = (int*)malloc(sizeof(int) * (N1 + 2));
	this->seq2_sampled_ct_bps = (int*)malloc(sizeof(int) * (N2 + 2));

	this->reset_bps();
/*
	// Init base pairing info.
	for(int cnt = 0; cnt <= N1; cnt++)
	{
		this->seq1_sampled_ct_bps[cnt] = 0;
	}

	// Init base pairing info.
	for(int cnt = 0; cnt <= N2; cnt++)
	{
		this->seq2_sampled_ct_bps[cnt] = 0;
	}
	*/
}

t_stoch_sampled_structures::~t_stoch_sampled_structures()
{
	free(this->seq1_sampled_ct_bps);
	free(this->seq2_sampled_ct_bps);
}

void t_stoch_sampled_structures::reset_bps()
{
	// Init base pairing info.
	for(int cnt = 0; cnt <= N1; cnt++)
	{
		this->seq1_sampled_ct_bps[cnt] = 0;
	}

	// Init base pairing info.
	for(int cnt = 0; cnt <= N2; cnt++)
	{
		this->seq2_sampled_ct_bps[cnt] = 0;
	}
}

void t_stoch_sampled_structures::add_bp_sampled_ct1(int i, int j)
{
	// Set i-j base pair, symmetrically.
	if(this->seq1_sampled_ct_bps[i] != 0 || this->seq1_sampled_ct_bps[j] != 0)
	{
		printf("%d paired %d and %d paired %d\n", i, this->seq1_sampled_ct_bps[i], j, this->seq1_sampled_ct_bps[j]);
		exit(0);
	}
	else
	{
	}

	this->seq1_sampled_ct_bps[i] = j;
	this->seq1_sampled_ct_bps[j] = i;
}

void t_stoch_sampled_structures::add_bp_sampled_ct2(int k, int l)
{
	// Set k-l base pair, symmetrically.
	if(this->seq2_sampled_ct_bps[k] != 0 || this->seq2_sampled_ct_bps[l] != 0)
	{
		printf("%d paired %d and %d paired %d\n", k, this->seq2_sampled_ct_bps[k], l, this->seq2_sampled_ct_bps[l]);
		exit(0);
	}
	else
	{
if(_DUMP_STOCH_SAMPLED_STRUCTURES_MESSAGES_)
		printf("Adding %d pair %d in ct2.\n", k,l);
	}

	this->seq2_sampled_ct_bps[k] = l;
	this->seq2_sampled_ct_bps[l] = k;
}

// The sampled ct's are usually a set of cts.
// Must enumerate the cts.
//void t_stoch_sampled_structures::dump_sampled_cts(double probability, char* extra_file_id)
//{
//if(_DUMP_STOCH_SAMPLED_STRUCTURES_MESSAGES_)
//	printf("Dumping first structure.\n");
//
//	// Create the sampling output directory.
//	char sampled_cts_dir[500];
//	concat_path_fn(sampled_cts_dir, SAMPLED_CTS_DIR, ppf_cli->aln_id);
//
//	char sampled_ct1_fn[100];
//	char sampled_ct2_fn[100];
//	char sampled_ct1_fp[500];
//	char sampled_ct2_fp[500];
//
//	if(strcmp(extra_file_id, "max_prob") != 0)
//	{
//		sprintf(sampled_ct1_fn, "sample_%d_ct1_%s.ct", sample_id, extra_file_id);
//		sprintf(sampled_ct2_fn, "sample_%d_ct2_%s.ct", sample_id, extra_file_id);
//		concat_path_fn(sampled_ct1_fp, sampled_cts_dir, sampled_ct1_fn);
//		concat_path_fn(sampled_ct2_fp, sampled_cts_dir, sampled_ct2_fn);
//	}
//	else
//	{
//		// Write max prob ct's to map_cts directory.
//		sprintf(sampled_ct1_fn, "%s_max_prob_sampled_ct1.ct", this->ppf_cli->aln_id);
//		sprintf(sampled_ct2_fn, "%s_max_prob_sampled_ct2.ct", this->ppf_cli->aln_id);
//		concat_path_fn(sampled_ct1_fp, MAP_CTS_DIR, sampled_ct1_fn);
//		concat_path_fn(sampled_ct2_fp, MAP_CTS_DIR, sampled_ct2_fn);
//	}
//
//	FILE* ct1_file = open_f(sampled_ct1_fp, "w");
//	fprintf(ct1_file, "%d\tct1_sampled Probability=%.5f\n", this->N1, probability);
//
//	// Dump all base pairing info for ct1.
//	for(int cnt = 1; cnt <= N1; cnt++)
//	{
//		// Sth. like following:
//		//     1 G       0    2   73    1
//		fprintf(ct1_file, "%d %c\t%d\t%d\t%d\t%d\n", cnt, this->seq_man->seq1->nucs[cnt], cnt-1, cnt+1, this->seq1_sampled_ct_bps[cnt], cnt);
//	}
//
//	fclose(ct1_file);
//
//if(_DUMP_STOCH_SAMPLED_STRUCTURES_MESSAGES_)
//	printf("Dumping second sampled structure.\n");
//
//	FILE* ct2_file = open_f(sampled_ct2_fp, "w");
//	fprintf(ct2_file, "%d\tct2_sampled Probability=%.5f\n", this->N2, probability);
//	// Dump all base pairing info for ct1.
//	for(int cnt = 1; cnt <= N2; cnt++)
//	{
//		// Sth. like following:
//		//     1 G       0    2   73    1
//		fprintf(ct2_file, "%d %c\t%d\t%d\t%d\t%d\n", cnt, this->seq_man->seq2->nucs[cnt], cnt-1, cnt+1, this->seq2_sampled_ct_bps[cnt], cnt);
//	}
//
//	fclose(ct2_file);
//
//}

// The sampled ct's are usually a set of cts.
// Must enumerate the cts.
void t_stoch_sampled_structures::dump_sampled_cts(int _sample_id)
{
if(_DUMP_STOCH_SAMPLED_STRUCTURES_MESSAGES_)
	printf("Dumping first structure.\n");

	// Create the sampling output directory.
	char ct1_fp[4096];
	if(this->ppf_cli->seq1_sample_ct_op == NULL)
	{
		sprintf(ct1_fp, "%s_sampled_cts.ct", this->ppf_cli->seq1_op_file_prefix);
	}
	else
	{
		strcpy(ct1_fp, this->ppf_cli->seq1_sample_ct_op);
	}
	FILE* ct1_file = NULL;
	ct1_file = open_f(ct1_fp, "a");

	fprintf(ct1_file, "%d\t %s_%d\t Energy=0.0\n", this->N1, this->ppf_cli->seq1_op_file_prefix, _sample_id);

	// Dump all base pairing info for ct1.
	for(int cnt = 1; cnt <= N1; cnt++)
	{
		// Sth. like following:
		//     1 G       0    2   73    1
		fprintf(ct1_file, "%d %c\t%d\t%d\t%d\t%d\n", cnt, this->seq_man->seq1->nucs[cnt], cnt-1, cnt+1, this->seq1_sampled_ct_bps[cnt], cnt);
	}

	fclose(ct1_file);

if(_DUMP_STOCH_SAMPLED_STRUCTURES_MESSAGES_)
	printf("Dumping second sampled structure.\n");

	char ct2_fp[4096];
	if(this->ppf_cli->seq2_sample_ct_op == NULL)
	{
		sprintf(ct2_fp, "%s_sampled_cts.ct", this->ppf_cli->seq2_op_file_prefix);
	}
	else
	{
		strcpy(ct2_fp, this->ppf_cli->seq2_sample_ct_op);
	}

	FILE* ct2_file = open_f(ct2_fp, "a");
	fprintf(ct2_file, "%d\t %s_%d\t Energy=0.0\n", this->N2, this->ppf_cli->seq2_op_file_prefix, _sample_id);

	// Dump all base pairing info for ct1.
	for(int cnt = 1; cnt <= N2; cnt++)
	{
		// Sth. like following:
		//     1 G       0    2   73    1
		fprintf(ct2_file, "%d %c\t%d\t%d\t%d\t%d\n", cnt, this->seq_man->seq2->nucs[cnt], cnt-1, cnt+1, this->seq2_sampled_ct_bps[cnt], cnt);
	}

	fclose(ct2_file);

}

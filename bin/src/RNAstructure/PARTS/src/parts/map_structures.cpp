#include <string.h>
#include <limits.h>
#include "map_structures.h"
#include "process_sequences.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../../src/phmm/structure/structure_object.h"
#include "../../../src/phmm/utils/file/utils.h"
#include "ppf_cli.h"
#include "parts_paths.h"

/*
class t_MAP_structures
{
public:
	int N1;
	int N2;
	int* seq1_map_ct_bps;
	int* seq2_map_ct_bps;	

	// Constructors and destructors.
	t_MAP_structures(int N1, int N2);
	~t_MAP_structures();

	// Add base pairs to structures while doing MAP decoding.
	void add_bp_to_MAP_ct1(int i, int j);
	void add_bp_to_MAP_ct2(int k, int l);

	// Dump MAP decoded ct's in ct format.
	void dump_map_cts();
};

#endif
*/

bool _DUMP_MAP_STRUCTURES_MESSAGES_ = false;

t_MAP_structures::t_MAP_structures(t_seq_man* seq_man, t_ppf_cli* ppf_cli)
{
	this->seq_man = seq_man;
	this->ppf_cli = ppf_cli;

	// Copy sequence lengths.
	this->N1 = seq_man->get_l_seq1();
	this->N2 = seq_man->get_l_seq2();

	// Allocate base pair arrays.
	this->seq1_map_ct_bps = (int*)malloc(sizeof(int) * (N1 + 2));
	this->seq2_map_ct_bps = (int*)malloc(sizeof(int) * (N2 + 2));

	// Init base pairing info.
	for(int cnt = 0; cnt <= N1; cnt++)
	{
		this->seq1_map_ct_bps[cnt] = 0;
	}

	// Init base pairing info.
	for(int cnt = 0; cnt <= N2; cnt++)
	{
		this->seq2_map_ct_bps[cnt] = 0;
	}
}

t_MAP_structures::~t_MAP_structures()
{
	free(this->seq1_map_ct_bps);
	free(this->seq2_map_ct_bps);
}

void t_MAP_structures::add_bp_to_MAP_ct1(int i, int j)
{
	// Set i-j base pair, symmetrically.
	this->seq1_map_ct_bps[i] = j;
	this->seq1_map_ct_bps[j] = i;
}

void t_MAP_structures::add_bp_to_MAP_ct2(int k, int l)
{
	// Set k-l base pair, symmetrically.
	this->seq2_map_ct_bps[k] = l;
	this->seq2_map_ct_bps[l] = k;
}

void t_MAP_structures::dump_map_cts()
{
if(_DUMP_MAP_STRUCTURES_MESSAGES_)
	printf("Dumping first structure.\n");

	char ct1_fp[4096];
	if(this->ppf_cli->seq1_map_ct_op == NULL)
	{
		sprintf(ct1_fp, "%s_map_ct.ct", this->ppf_cli->seq1_op_file_prefix);
	}
	else
	{
		strcpy(ct1_fp, this->ppf_cli->seq1_map_ct_op);
	}

	FILE* ct1_file = open_f(ct1_fp, "w");

	fprintf(ct1_file, "%d\t%s_predicted\n", this->N1, this->ppf_cli->seq1_op_file_prefix);

	// Dump all base pairing info for ct1.
	for(int cnt = 1; cnt <= N1; cnt++)
	{
		// Sth. like following:
		//     1 G       0    2   73    1
		fprintf(ct1_file, "%d %c\t%d\t%d\t%d\t%d\n", cnt, this->seq_man->seq1->nucs[cnt], cnt-1, cnt+1, this->seq1_map_ct_bps[cnt], cnt);
	}

	fclose(ct1_file);

if(_DUMP_MAP_STRUCTURES_MESSAGES_)
	printf("Dumping second structure.\n");

	char ct2_fp[4096];
	if(this->ppf_cli->seq2_map_ct_op == NULL)
	{
		sprintf(ct2_fp, "%s_map_ct.ct", this->ppf_cli->seq2_op_file_prefix);
	}
	else
	{
		strcpy(ct2_fp, this->ppf_cli->seq2_map_ct_op);
	}

	FILE* ct2_file = open_f(ct2_fp, "w");
	fprintf(ct2_file, "%d\t%s_predicted\n", this->N2, this->ppf_cli->seq2_op_file_prefix);

	// Dump all base pairing info for ct1.
	for(int cnt = 1; cnt <= N2; cnt++)
	{
		// Sth. like following:
		//     1 G       0    2   73    1
		fprintf(ct2_file, "%d %c\t%d\t%d\t%d\t%d\n", cnt, this->seq_man->seq2->nucs[cnt], cnt-1, cnt+1, this->seq2_map_ct_bps[cnt], cnt);
	}

	fclose(ct2_file);

}

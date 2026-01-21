#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "../../../src/phmm/structure/structure_object.h"
#include "map_structures.h"
#include "map_alignment.h"
#include "map_mhr_info.h"
#include "map_results.h"
#include "ppf_cli.h"
#include "process_sequences.h"
#include <string.h>

/*
	// Construct the results from map structures and map alignment.
	t_map_results(t_MAP_structures* _map_strs, t_MAP_alignment* _map_alignment, t_MAP_mhr_info* _map_mhr_info);
	t_MAP_structures* map_strs;
	t_MAP_alignment* map_alignment;
	t_MAP_mhr_info* map_mhr_info;

	~t_map_results();

	// Mutate map_str's into two structure classes.
	t_structure* str1;
	t_structure* str2;
	void mutate_map_strs();*/

t_map_results::t_map_results(t_seq_man* _seq_man, 
							 t_ppf_cli* _ppf_cli, 
							 t_MAP_structures* _map_strs, 
							 t_MAP_alignment* _map_alignment, 
							 t_MAP_mhr_info* _map_mhr_info)
{
	this->map_strs = _map_strs;
	this->map_alignment = _map_alignment;
	this->map_mhr_info = _map_mhr_info;

	this->ppf_cli = _ppf_cli;
	this->seq_man = _seq_man;

	// Determine individual t_structure's.
	mutate_map_strs();
}

// Create t_structure's.
void t_map_results::mutate_map_strs()
{
	this->str1 = new t_structure();
	this->str2 = new t_structure();

	// Set str1.
	this->str1->numofbases = this->map_strs->N1;
	this->str1->ctlabel = (char*)malloc((strlen(this->ppf_cli->seq1_path)) + 3);
	sprintf(this->str1->ctlabel, this->ppf_cli->seq1_path);

	// Set str2.
	this->str2->numofbases = this->map_strs->N2;
	this->str2->ctlabel = (char*)malloc((strlen(this->ppf_cli->seq2_path)) + 3);
	sprintf(this->str2->ctlabel, this->ppf_cli->seq2_path);

	// Set file pointer.
	//this->str1->fp = NULL;
	//this->str2->fp = NULL;

	// Set nucs and numseqs.
	this->str1->nucs = (char*)malloc( this->map_strs->N1 + 3);
	this->str1->numseq = (int*)malloc( sizeof(int) * (this->map_strs->N1 + 3) );
	this->str1->basepr = (int*)malloc( sizeof(int) * (this->map_strs->N1 + 3) );
	for(int cnt = 1; cnt <= this->map_strs->N1; cnt++)
	{
		this->str1->nucs[cnt] = this->seq_man->seq1->nucs[cnt];
		this->str1->numseq[cnt] = this->seq_man->seq1->numseq[cnt];
		this->str1->basepr[cnt] = this->map_strs->seq1_map_ct_bps[cnt];
	}

	this->str2->nucs = (char*)malloc( this->map_strs->N2 + 3);
	this->str2->numseq = (int*)malloc( sizeof(int) * (this->map_strs->N2 + 3) );
	this->str2->basepr = (int*)malloc( sizeof(int) * (this->map_strs->N2 + 3) );
	for(int cnt = 1; cnt <= this->map_strs->N2; cnt++)
	{		
		this->str2->nucs[cnt] = this->seq_man->seq2->nucs[cnt];
		this->str2->numseq[cnt] = this->seq_man->seq2->numseq[cnt];
		this->str2->basepr[cnt] = this->map_strs->seq2_map_ct_bps[cnt];
	}
	
	// Everything is set.
}


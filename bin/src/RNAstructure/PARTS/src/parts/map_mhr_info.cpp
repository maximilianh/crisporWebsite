#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "map_structures.h"
#include "map_alignment.h"
#include "ppf_cli.h"
#include "map_mhr_info.h"
#include <string.h>
#include "phmm_parameters.h" 

bool _DUMP_MAP_MFHR_INFO_MSGS_ = false;

t_MAP_mhr_info::t_MAP_mhr_info(t_ppf_cli* _ppf_cli, t_MAP_structures* _map_structures, t_MAP_alignment* _map_alignment)
{
	// Copy pointers.
	this->map_alignment = _map_alignment;
	this->map_structures = _map_structures;
	this->ppf_cli = _ppf_cli;

	// Do the traceback.
	this->mhrs1 = (int*)malloc(sizeof(int) * (this->map_structures->N1 + 2));
	this->mhrs2 = (int*)malloc(sizeof(int) * (this->map_structures->N2 + 2));
	this->mhr1_sides = (int*)malloc(sizeof(int) * (this->map_structures->N1 + 2));
	this->mhr2_sides = (int*)malloc(sizeof(int) * (this->map_structures->N2 + 2));

	// Initialize mhr arrays.
	for(int cnt = 0; cnt < this->map_structures->N1 + 2; cnt++)
	{
		this->mhrs1[cnt] = 0;
		this->mhr1_sides[cnt] = -1;
	}

	for(int cnt = 0; cnt < this->map_structures->N2 + 2; cnt++)
	{
		this->mhrs2[cnt] = 0;
		this->mhr2_sides[cnt] = -1;
	}

	// Trace he mhr's and fill mhr arrays.
	this->trace_MAP_mhrs();
}

t_MAP_mhr_info::~t_MAP_mhr_info()
{
}

// Implementation of mhr tracing using whole alignment strings.
// 1. Create the base pairing mask, which is as long as the alignment strings.
// 2. Go over the pairing mask and to determine the mhr's over the matched helical region.
// 3. Using the individual alignment strings of sequences, assign the mhr's for both sequences.
void t_MAP_mhr_info::trace_MAP_mhrs()
{
	// Copy alignment strings.
	char* aln_str1 = this->map_alignment->aln_str1;
	char* aln_str2 = this->map_alignment->aln_str2;

	// The index maps.
	// Note that alignment strings are indexed starting from 0,
	// sequences are indexed starting from 1.
	int* seq1_to_aln_str1 = (int*)malloc(sizeof(int) * (map_structures->N1 + 3));
	int* seq2_to_aln_str2 = (int*)malloc(sizeof(int) * (map_structures->N2 + 3));

	int* aln_str1_to_seq1 = (int*)malloc(sizeof(int) * (strlen(aln_str1) + 3));
	int* aln_str2_to_seq2 = (int*)malloc(sizeof(int) * (strlen(aln_str1) + 3));

	// Determine the mappings:
	// seqi_to_aln_stri corresponds to 
	int seq1_cnt = 1;
	for(int aln_i = 0; aln_i < strlen(aln_str1); aln_i++)
	{
		// Is there a nucleotide in the index in alignment string?
		if(aln_str1[aln_i] != GAP_CHAR)
		{
			aln_str1_to_seq1[aln_i] = seq1_cnt;
			seq1_to_aln_str1[seq1_cnt] = aln_i;
			seq1_cnt++;
		}
		else
		{
			// This point is a gap, set it 0.
			aln_str1_to_seq1[aln_i] = 0;
		}
	}

	int seq2_cnt = 1;
	for(int aln_i = 0; aln_i < strlen(aln_str2); aln_i++)
	{
		// Is there a nucleotide in the index in alignment string?
		if(aln_str2[aln_i] != GAP_CHAR)
		{
			aln_str2_to_seq2[aln_i] = seq2_cnt;
			seq2_to_aln_str2[seq2_cnt] = aln_i;
			seq2_cnt++;
		}
		else
		{
			// This point is a gap, set it 0.
			aln_str2_to_seq2[aln_i] = 0;
		}
	}

/*
	printf("\n%s\n", aln_str2);
	for(int aln_i = 1; aln_i <= map_structures->N2; aln_i++)
	{
		printf("%d ", seq2_to_aln_str2[aln_i]);
	}
	
	for(int aln_i = 0; aln_i < strlen(aln_str2); aln_i++)
	{
		printf("%d ", aln_str2_to_seq2[aln_i]);
	}
*/
	// Allocate base pairing mask.
	// Base pairing mask represents the base pairing in alignment strings.
	int* pairing_mask = (int*)malloc(sizeof(int) * (strlen(aln_str1) + 3) );
	char* pairing_mask_str = (char*)malloc(sizeof(char) * (strlen(aln_str1) + 3) );
	int* pairing_mask_mhrs = (int*)malloc(sizeof(int) * (strlen(aln_str1) + 3) );

	// Initialize pairing mask.
	for(int aln_i = 0; aln_i <= strlen(aln_str1); aln_i++)
	{
		// Initialize pairing mask.
		// No pairing.
		pairing_mask[aln_i] = 0;
		pairing_mask_str[aln_i] = 0;
		pairing_mask_mhrs[aln_i] = 0;
	}

	// Set pairing mask.
	for(int aln_i = 0; aln_i < strlen(aln_str1); aln_i++)
	{
		// Check both sequences if ther is one base pair in one of the sequences, set the base pair.
		int seq1_index = aln_str1_to_seq1[aln_i];
		int seq2_index = aln_str2_to_seq2[aln_i];

		int seq1_pair_index = this->map_structures->seq1_map_ct_bps[seq1_index];
		int seq2_pair_index = this->map_structures->seq2_map_ct_bps[seq2_index];

if(_DUMP_MAP_MFHR_INFO_MSGS_)
		printf("aln_i = %d, seq1_index = %d, seq2_index = %d\n", aln_i, seq1_index, seq2_index);

		if(seq1_pair_index != 0)
		{
if(_DUMP_MAP_MFHR_INFO_MSGS_)
			printf("SEQ1 PAIR\n");

			pairing_mask[aln_i] = seq1_to_aln_str1[seq1_pair_index]; // Set the pair of the nucleotide in pairing mask.

			if(seq1_pair_index > seq1_index)
				pairing_mask_str[aln_i] = '<';
			else
				pairing_mask_str[aln_i] = '>';

if(_DUMP_MAP_MFHR_INFO_MSGS_)
			printf("%d <-> %d\n", aln_i, seq1_to_aln_str1[seq1_pair_index]);
		}
		else if(seq2_pair_index != 0)
		{
if(_DUMP_MAP_MFHR_INFO_MSGS_)
			printf("SEQ2 PAIR\n");

			pairing_mask[aln_i] = seq2_to_aln_str2[seq2_pair_index]; // Set the pair of the nucleotide in pairing mask.

			if(seq2_pair_index > seq2_index)
				pairing_mask_str[aln_i] = '<';
			else
				pairing_mask_str[aln_i] = '>';

if(_DUMP_MAP_MFHR_INFO_MSGS_)
			printf("%d <-> %d\n", aln_i, seq2_to_aln_str2[seq2_pair_index]);
		}
		else
		{
			// No paired base here or inserts in alignment string.
if(_DUMP_MAP_MFHR_INFO_MSGS_)
			printf("EMPTY\n");

			pairing_mask_str[aln_i] = '.';
		}
	}

if(_DUMP_MAP_MFHR_INFO_MSGS_)
{
	printf("\n%s\n", aln_str1);
	printf("%s\n", aln_str2);
	printf("%s\n", pairing_mask_str);
}

	// Have the pairing mask. 
	// Start determining the mhr's,
	// inside the pairing mask.
	int mhr_cnt = 1;
	for(int aln_i = 0; aln_i < strlen(aln_str1); aln_i++)
	{
		// Put this base pair to an mhr.
		int aln_posn_index = aln_i;
		int aln_posn_pair_index = pairing_mask[aln_i];
		int exterior_mhr = 0;

		if(pairing_mask[aln_i] != 0 && aln_posn_index < aln_posn_pair_index)
		{
			// Determine the exterior mhr.
			if(aln_i != 0)
			{
				if(aln_posn_index < aln_posn_pair_index && pairing_mask_mhrs[aln_posn_index - 1] == pairing_mask_mhrs[aln_posn_pair_index + 1])
				{
					exterior_mhr = pairing_mask_mhrs[aln_posn_pair_index + 1];
				}
			}
			else
			{
				exterior_mhr = 0;
			}


			// Determine the mhr index of this base paired position.
			if(exterior_mhr != 0)
			{
				pairing_mask_mhrs[aln_posn_index] = exterior_mhr;
				pairing_mask_mhrs[aln_posn_pair_index] = exterior_mhr;
			}
			else
			{
				pairing_mask_mhrs[aln_posn_index] = mhr_cnt;
				pairing_mask_mhrs[aln_posn_pair_index] = mhr_cnt;
				mhr_cnt++;
			} // exterior mhr check.
		} // If this aln position is a base pair in pairing mask.	

if(_DUMP_MAP_MFHR_INFO_MSGS_)
		printf("%d", pairing_mask_mhrs[aln_i]);
	} // Go over all the alignment positions in pairing mask.

	// Copy the mhr information based on the template mhr array to sequences themselves.
	for(int aln_i = 0; aln_i < strlen(aln_str1); aln_i++)
	{
		int seq1_index = aln_str1_to_seq1[aln_i];
		int seq2_index = aln_str2_to_seq2[aln_i];
		if(seq1_index != 0)
		{
			this->mhrs1[seq1_index] = pairing_mask_mhrs[aln_i];
		}

		if(seq2_index != 0)
		{
			this->mhrs2[seq2_index] = pairing_mask_mhrs[aln_i];
		}
	}

	fflush(stdout);
}



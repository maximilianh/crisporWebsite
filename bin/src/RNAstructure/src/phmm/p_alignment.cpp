#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "p_alignment.h"

t_p_alignment::t_p_alignment(char* _seq1_aln_line, char* _seq2_aln_line)
{
	if(strlen(_seq1_aln_line) != strlen(_seq2_aln_line))
	{
		printf("alignment lines are not of same length, exiting at %s(%d)\n", __FILE__, __LINE__);
		this->seq1_aln_line = NULL;
		this->seq2_aln_line = NULL;
	}
	else
	{
		this->seq1_aln_line = (char*)malloc(sizeof(char) * (strlen(_seq1_aln_line) + 3));
		this->seq2_aln_line = (char*)malloc(sizeof(char) * (strlen(_seq2_aln_line) + 3));

		strcpy(this->seq1_aln_line, _seq1_aln_line);
		strcpy(this->seq2_aln_line, _seq2_aln_line);
	}
}

t_p_alignment::~t_p_alignment()
{
	free(this->seq1_aln_line);
	free(this->seq2_aln_line);
}

double t_p_alignment::get_aln_similarity(char gap_symbol)
{
	if(seq1_aln_line == NULL || seq2_aln_line == NULL)
	{
		return(-1);
	}

	int n_match_pos = 0;

	for(int cnt = 0; cnt < (int)strlen(seq1_aln_line); cnt++)
	{
		if(seq1_aln_line[cnt] != gap_symbol && seq1_aln_line[cnt] == seq2_aln_line[cnt])
		{
			n_match_pos++;
		}
	}

	int n_total_pos = 0;

	for(int cnt = 0; cnt < (int)strlen(seq1_aln_line); cnt++)
	{
		if(seq1_aln_line[cnt] == gap_symbol && seq2_aln_line[cnt] == gap_symbol)
		{
		}
		else
		{
			n_total_pos++;
		}
	}

	return( (double)n_match_pos / n_total_pos );
}


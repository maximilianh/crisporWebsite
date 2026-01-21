#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "parts_compilation_directives.h"
#include "pf_alignment.h"
#include "ppf_w.h"
#include "ppf_w_mb.h"
#include "process_sequences.h"
#include "ppf_math.h"
#include <iostream>
using namespace std;

// Constructor: Allocate and init state prob array.
t_pf_alignment::t_pf_alignment(t_seq_man* seq_man)
{
	// Copy sequence lengths.
	this->N1 = seq_man->get_l_seq1();
	this->N2 = seq_man->get_l_seq2();

	printf("Alloc'ing and init'ing alignment states probabilities.\n");

	// Allocate and initialize alignment states probabilities.
	// Accession to proabilities array is as following:
	// 1st index is seq1 index, 2nd index is seq2 index, 3rd index is state index.
	this->state_probs_array = (double***)malloc(sizeof(double**) * (this->N1 + 3) * (this->N2 + 3) * 3); // 3 for 3 states.

	for(int i1 = 1; i1 <= this->N1; i1++)
	{
		this->state_probs_array[i1] = (double**)malloc(sizeof(double*) * (this->N2 + 3));

		for(int i2 = 1; i2 <= this->N2; i2++)
		{
			this->state_probs_array[i1][i2] = (double*)malloc(sizeof(double) * 3); // 3 for 3 states.

			// Initialize all 3 state probabilities.
			this->state_probs_array[i1][i2][STATE_ALN] = CONVERT_FROM_LIN(0);
			this->state_probs_array[i1][i2][STATE_INS1] = CONVERT_FROM_LIN(0); 
			this->state_probs_array[i1][i2][STATE_INS2] = CONVERT_FROM_LIN(0);
		}
	}
}

t_pf_alignment::~t_pf_alignment()
{
	printf("Destruct'ing a t_pf_alignment object.\n");
}

// Main function to fill state probabilities.
void t_pf_alignment::infer_seq_aln_probs(t_ppf_W* W, t_ppf_WMB* WMB)
{
	// Call calculation of all 3 state probabilities over all alignment positions.
	infer_seq_ALN_probs(W, WMB);
	infer_seq_INS1_probs(W, WMB);
	infer_seq_INS2_probs(W, WMB);
}

// Calculate ALN state proabilities over all alignment positions, based on leaving 
// current positions, alone and marginalizing over all possible open and closed structures.
// Calculation depends on leaving i and k asa aligned and marginalizing over all possible
// open and closed interior and exterior substructures possible.
void t_pf_alignment::infer_seq_ALN_probs(t_ppf_W* W, t_ppf_WMB* WMB)
{
	printf("Calculating ALN state probabilities.\n");
	
	/* 
	ATTACK PLAN:
	Current alignment positions where ALN state probability is calculated is indexed using i and k.
	j and l are used for marginalizing over all possible interior and exterior structures. So i and k i constant
	in outer loop and j and l are changing in inner loop.
	For any value;
	j starts from 1 to N1 and skips i, 
	l starts from 1 to N2 and skips k,
	if j < i
		Take interior open andclosed substructures of j to i-1, take exterior open and closed substructures of j-1 to i+1
	else
		Take interior ... substructures of i+1 to j and exterior ... substructures of i-1 to j.

	Similarly choose k and l substructures similarly by skipping k at correct places.

	Do a proper bound check in 4 loops:
	The loop limits must be satisfied, note that this will make loops not calculate probabilities of states of some of the alignment positions,

	*/

	// Outer loops is i, k loop.
	for(int i = 1; i <= this->N1; i++)
	{
		for(int k = 1; k <= this->N2; k++)
		{
			// Inner loops are j, l loops.
			for(int j = 1; j <= this->N1; j++)
			{
				for(int l = 1; l <= this->N2; l++)
				{
					/*
					if(W->check_boundary(i,k) && W->check_boundary(j,l))
					{
						// Following are needed in order to satisfy "knot" condition for alignment:
						// i-k are ALN AND j-l are coincident, so j and l MUST BE on the same side with respect
						// to i and k.
						if(j < i && l < k)
						{
							
						}
						// LEFT HERE: N1 and N2 checks for i and k...
						else if(j > i && l > k)
						{
							double cur_int_score = SUM(W->x(i+1, j, k+1, l), WMB->x(i+1, j, k+1, l));
							double cur_ext_score = SUM(W->x(j+1, i - 1 + N1, l+1, k - 1 + N2), WMB->x(W->x(j+1, i - 1 + N1, l+1, k - 1 + N2));

							double current_ALN_score = MUL(current_interior_score, current_exterior_score);
						}
					}
					*/
				}
			}
		}
	}
}

void t_pf_alignment::infer_seq_INS1_probs(t_ppf_W* W, t_ppf_WMB* WMB)
{
	printf("Calculating INS1 state probabilities.\n");
}

void t_pf_alignment::infer_seq_INS2_probs(t_ppf_W* W, t_ppf_WMB* WMB)
{
	printf("Calculating INS2 state probabilities.\n");
}

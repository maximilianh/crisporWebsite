#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include "folding_constraints.h"
#include "structure_object.h"
#include "../utils/file/utils.h"

bool _DUMP_FOLDING_CONSTRAINTS_MESSAGES_ = false;

int t_folding_constraints::pairable[5][5] = {	{ 0, 0, 0, 0, 0 },  // -
												{ 0, 0, 0, 0, 1 },  // A
												{ 0, 0, 0, 1, 0 },  // C
												{ 0, 0, 1, 0, 1 },  // G
												{ 0, 1, 0, 1, 0 }}; // U

t_folding_constraints::t_folding_constraints(t_structure* _rna_seq)
{
	//this->energy_loops = _energy_loops;
	this->rna_seq = new t_structure(_rna_seq);

	this->alloc_init_maps();

	// This is the only way to compute pointer relocation maps.
	this->mallocate_ptr_reloc_maps();

	// Compute the relocation maps.
	this->compute_ptr_reloc_maps(1.0f, NULL);
}

t_folding_constraints::t_folding_constraints(char* _rna_seq_fp)
{
	//this->energy_loops = _energy_loops;
	this->rna_seq = new t_structure(_rna_seq_fp);

	this->alloc_init_maps();

	// This is the only way to compute pointer relocation maps.
	this->mallocate_ptr_reloc_maps();

	// Compute the relocation maps.
	this->compute_ptr_reloc_maps(1.0f, NULL);
}

// Generate all the maps from pairing array and threshold in the parameters.
t_folding_constraints::t_folding_constraints(t_structure* _rna_seq, double** linear_pp, double linear_pp_threshold)
{
	this->rna_seq = new t_structure(_rna_seq);

	// Allocate and init maps.
	this->alloc_init_maps();

	// This is the only way to compute pointer relocation maps.
	this->mallocate_ptr_reloc_maps();

	// Compute the relocation maps.
	this->compute_ptr_reloc_maps(linear_pp_threshold, linear_pp);
}

// Generate all the maps from pairing array and threshold in the parameters.
t_folding_constraints::t_folding_constraints(char* _seq_fp, double** linear_pp, double linear_pp_threshold)
{
	this->rna_seq = new t_structure(_seq_fp);

	// Allocate and init maps.
	this->alloc_init_maps();

	// This is the only way to compute pointer relocation maps.
	this->mallocate_ptr_reloc_maps();

	// Compute the relocation maps.
	this->compute_ptr_reloc_maps(linear_pp_threshold, linear_pp);
}

// Copy constructor.
t_folding_constraints::t_folding_constraints(t_folding_constraints* _folding_constraints)
{
	// Copy the maps.
	this->rna_seq = new t_structure(_folding_constraints->rna_seq);

	// Allocate and initialize maps.
	this->same_loop_map = (bool**)malloc(sizeof(bool*) * (this->rna_seq->numofbases + 2));
	this->str_coinc_map = (bool**)malloc(sizeof(bool*) * (this->rna_seq->numofbases + 2));
	this->pairing_map = (bool**)malloc(sizeof(bool*) * (this->rna_seq->numofbases + 2));
	this->forbid_non_v = (bool*)malloc(sizeof(bool) * (this->rna_seq->numofbases + 2));

	for(int i = 1; i <= this->rna_seq->numofbases; i++)
	{
		this->same_loop_map[i] = (bool*)malloc(sizeof(bool) * ((this->rna_seq->numofbases + 4 - (i))));
		this->str_coinc_map[i] = (bool*)malloc(sizeof(bool) * ((this->rna_seq->numofbases + 4 - (i))));
		this->pairing_map[i] = (bool*)malloc(sizeof(bool) * ((this->rna_seq->numofbases + 4 - (i))));
		this->forbid_non_v[i] = false;

		this->same_loop_map[i] -= (i);
		this->str_coinc_map[i] -= (i);
		this->pairing_map[i] -= (i);

		// Set all to true, everything is allowed.
		for(int j = i; j <= this->rna_seq->numofbases; j++)
		{
			this->same_loop_map[i][j] = _folding_constraints->same_loop_map[i][j];
			this->str_coinc_map[i][j] = _folding_constraints->str_coinc_map[i][j];
			this->pairing_map[i][j] = _folding_constraints->pairing_map[i][j];
		} // j loop.
	} // i loop.	

	// Copy ptr relocation maps, too, if necessary.
	if(_folding_constraints->coinc_pointer_relocation_map != NULL)
	{
		this->coinc_pointer_relocation_map = (short**)malloc(sizeof(short*) * (rna_seq->numofbases + 3));
		this->paired_pointer_relocation_map = (short**)malloc(sizeof(short*) * (rna_seq->numofbases + 3));

		for(int i = 1; i <= rna_seq->numofbases; i++)
		{
			coinc_pointer_relocation_map[i] = (short*)malloc(sizeof(short) * (rna_seq->numofbases - i + 3));
			coinc_pointer_relocation_map[i] -= i; // Do pointer shift for fold envelope.

			paired_pointer_relocation_map[i] = (short*)malloc(sizeof(short) * (rna_seq->numofbases - i + 3));
			paired_pointer_relocation_map[i] -= i; // Do pointer shift for fold envelope.

			for(int j = i; j <= rna_seq->numofbases; j++)
			{
				coinc_pointer_relocation_map[i][j] = _folding_constraints->coinc_pointer_relocation_map[i][j];
				paired_pointer_relocation_map[i][j] = _folding_constraints->paired_pointer_relocation_map[i][j];
			} // Initialization loop
		}
	}
	else
	{
		this->coinc_pointer_relocation_map = NULL;
		this->paired_pointer_relocation_map = NULL;
	}
}

void t_folding_constraints::alloc_init_maps()
{
	this->coinc_pointer_relocation_map = NULL;
	this->paired_pointer_relocation_map = NULL;
	this->saturated_str_bps = NULL;

	// Allocate and initialize maps.
	this->same_loop_map = (bool**)malloc(sizeof(bool*) * (this->rna_seq->numofbases + 2));
	this->str_coinc_map = (bool**)malloc(sizeof(bool*) * (this->rna_seq->numofbases + 2));
	this->pairing_map = (bool**)malloc(sizeof(bool*) * (this->rna_seq->numofbases + 2));
	this->forbid_non_v = (bool*)malloc(sizeof(bool) * (this->rna_seq->numofbases + 2));

	for(int i = 1; i <= this->rna_seq->numofbases; i++)
	{
		this->same_loop_map[i] = (bool*)malloc(sizeof(bool) * ((this->rna_seq->numofbases + 4 - (i))));
		this->str_coinc_map[i] = (bool*)malloc(sizeof(bool) * ((this->rna_seq->numofbases + 4 - (i))));
		this->pairing_map[i] = (bool*)malloc(sizeof(bool) * ((this->rna_seq->numofbases + 4 - (i))));
		this->forbid_non_v[i] = false;

		this->same_loop_map[i] -= (i);
		this->str_coinc_map[i] -= (i);
		this->pairing_map[i] -= (i);

		// Set all to true, everything is allowed.
		for(int j = i; j <= this->rna_seq->numofbases; j++)
		{
			this->same_loop_map[i][j] = true;
			this->str_coinc_map[i][j] = true;

			// Check pairability of i and j for setting pairing map.
			if(t_folding_constraints::pairable[this->rna_seq->numseq[i]][this->rna_seq->numseq[j]])
			{
				this->pairing_map[i][j] = true;
			}
			else
			{
				this->pairing_map[i][j] = false;
			}
		} // j loop.
	} // i loop.
}

#define SATURATED_STRUCTURE_BP_DELTA_PROB (0.05)
void t_folding_constraints::compute_saturated_structure(double** pairing_probs)
{
	if(this->saturated_str_bps != NULL)
	{
		free(this->saturated_str_bps);
		this->saturated_str_bps = (int*)malloc(sizeof(int) * (this->rna_seq->numofbases + 4));
	}
	else
	{
		this->saturated_str_bps = (int*)malloc(sizeof(int) * (this->rna_seq->numofbases + 4));
	}

	// Starting from 0.50, find the lowest threshold for which the thresholding generates a valid structure.
	bool is_structure_valid = true;
	double cur_bp_threshold = 0.50f;
	while(cur_bp_threshold >= SATURATED_STRUCTURE_BP_DELTA_PROB && 
		is_structure_valid)
	{
		// (Re)Initialize base pairs.
		for(int i = 1; i <= this->rna_seq->numofbases; i++)
		{
			saturated_str_bps[i] = 0;
		} // i loop.

		// Get the base pairs for current threshold.
		for(int i = 1; is_structure_valid && (i <= this->rna_seq->numofbases); i++)
		{
			for(int j = 1; is_structure_valid && (j <= this->rna_seq->numofbases); j++)
			{
				// Probability thresholding.
				if(pairing_probs[i][j] > cur_bp_threshold)
				{
					// If this pair does not already have a base pair, add the base pair.
					if(saturated_str_bps[i] == 0 || saturated_str_bps[i] == j)
					{
						saturated_str_bps[i] = j;
						saturated_str_bps[j] = i;
					}
					else
					{
						printf("%lf: (%d, %d) @ %lf and (%d, %d) @ %lf conflicting.\n", cur_bp_threshold,
							i, j, pairing_probs[i][j], i, saturated_str_bps[i], pairing_probs[i][saturated_str_bps[i]]);

						is_structure_valid = false;
					}
				}
			} // j loop
		} // i loop.

		// Check validity of the structure.
		for(int i = 1; is_structure_valid && (i <= this->rna_seq->numofbases); i++)
		{
			// Do the check if this is paired.
			if(saturated_str_bps[i] > i)
			{
				int j = saturated_str_bps[i];

				for(int ip = i+1; ip < j; ip++)
				{
					if(saturated_str_bps[ip] > j)
					{
						printf("%lf: (%d, %d) @ %lf and (%d, %d) @ %lf are pseudo-knotted.\n", cur_bp_threshold,
							i, j, pairing_probs[i][j], ip, saturated_str_bps[ip], pairing_probs[ip][saturated_str_bps[ip]]);

						is_structure_valid = false;
					}
				}
			} // base pairing check for i.
		} // i loop.

		if(is_structure_valid)
		{
			cur_bp_threshold -= SATURATED_STRUCTURE_BP_DELTA_PROB;
		}
	} // bp_threshold loop.

	// Take the smallest base pairing probability threshold for valid structure generation.
	cur_bp_threshold += SATURATED_STRUCTURE_BP_DELTA_PROB;

	// (Re)Initialize base pairs.
	for(int i = 1; i <= this->rna_seq->numofbases; i++)
	{
		saturated_str_bps[i] = 0;
	} // i loop.

	printf("Smallest threshold for valid structure is %lf.\n", cur_bp_threshold);

	// Get the base pairs for the smallest threshold.
	is_structure_valid = true;
	for(int i = 1; is_structure_valid && (i <= this->rna_seq->numofbases); i++)
	{
		for(int j = 1; is_structure_valid && (j <= this->rna_seq->numofbases); j++)
		{
			// Probability thresholding.
			if(pairing_probs[i][j] > cur_bp_threshold)
			{
				// If this pair does not already have a base pair, add the base pair.
				if(saturated_str_bps[i] == 0 || saturated_str_bps[i] == j)
				{
					saturated_str_bps[i] = j;
					saturated_str_bps[j] = i;
				}
				else
				{
					is_structure_valid = false;
				}
			}
		} // j loop
	} // i loop.

	if(!is_structure_valid)
	{
		printf("Structure validation failed for lowest probability threshold of %lf @ %s(%d)\n", cur_bp_threshold, __FILE__, __LINE__);
		exit(0);
	}
}

void t_folding_constraints::free_maps()
{
	// Allocate and initialize maps.
	//this->same_loop_map = (bool**)malloc(sizeof(bool*) * (this->rna_seq->numofbases + 2));
	//this->str_coinc_map = (bool**)malloc(sizeof(bool*) * (this->rna_seq->numofbases + 2));
	//this->pairing_map = (bool**)malloc(sizeof(bool*) * (this->rna_seq->numofbases + 2));
	//this->forbid_non_v = (bool*)malloc(sizeof(bool) * (this->rna_seq->numofbases + 2));

	for(int i = 1; i <= this->rna_seq->numofbases; i++)
	{
		this->same_loop_map[i] += (i);
		this->str_coinc_map[i] += (i);
		this->pairing_map[i] += (i);
		free(this->same_loop_map[i]);
		free(this->str_coinc_map[i]);
		free(this->pairing_map[i]);
	}

	free(this->same_loop_map);
	free(this->str_coinc_map);
	free(this->pairing_map);
	free(this->forbid_non_v);
}

bool t_folding_constraints::forbid_non_v_emission(int i)
{
	return(this->forbid_non_v[i]);
}

bool t_folding_constraints::forbid_non_v_emission(int i, int j)
{
	return(this->forbid_non_v[i] || this->forbid_non_v[j]);
}

void t_folding_constraints::dump_constraints()
{
	char str_coinc_map_fp[1000];
	char same_loop_map_fp[1000];
	char pairing_fp[1000];

	sprintf(str_coinc_map_fp, "%s_str_coinc_map.txt", this->rna_seq->ctlabel);
	sprintf(same_loop_map_fp, "%s_same_loop_map.txt", this->rna_seq->ctlabel);
	sprintf(pairing_fp, "%s_pairing_map.txt", this->rna_seq->ctlabel);

	FILE* f_str_coinc_map = open_f(str_coinc_map_fp, "w");
	FILE* f_same_loop_map = open_f(same_loop_map_fp, "w");
	FILE* f_pairing_map = open_f(pairing_fp, "w");

	for(int i = 1; i <= this->rna_seq->numofbases; i++)
	{
		// Set all to true, everything is allowed.
		for(int j = 1; j <= this->rna_seq->numofbases; j++)
		{
			if(j > i)
			{
				fprintf(f_str_coinc_map, "%d", this->str_coinc_map[i][j]);
				fprintf(f_same_loop_map, "%d", this->same_loop_map[i][j]);
				fprintf(f_pairing_map, "%d", this->pairing_map[i][j]);
			}
			else
			{
				fprintf(f_str_coinc_map, "%d", this->str_coinc_map[j][i]);
				fprintf(f_same_loop_map, "%d", this->same_loop_map[j][i]);
				fprintf(f_pairing_map, "%d", this->pairing_map[j][i]);
			}
		} // j loop

		fprintf(f_str_coinc_map, "\n");
		fprintf(f_same_loop_map, "\n");
		fprintf(f_pairing_map, "\n");
	} // i loop

	fclose(f_str_coinc_map);
	fclose(f_same_loop_map);
	fclose(f_pairing_map);
}

bool t_folding_constraints::check_hairpin_loop(int i, int j)
{
	//printf("Checking hairpin loop %d, %d\n", i,j);
	if(i > j)
	{
		int temp = i;
		j = i;
		i = temp;
	}

	for(int _i = i; _i <= j; _i++)
	{
		if(_i != i && _i != j &&
			this->forbid_non_v_emission(_i))
		{
			//printf("forbidding hairpin loop %d, %d\n", i,j);
			//getc(stdin);
			return(false);
		}

		for(int _j = _i+1; _j <= j; _j++)
		{
			if(!this->same_loop_map[_i][_j])
			{
				//printf("Invalid hairpin loop %d, %d\n", i,j);
				//getc(stdin);
				return(false);
			} 
		} // _j loop
	} // _i loop

	//printf("Checked hairpin loop %d, %d\n", i,j);

	return(true);
}

bool t_folding_constraints::check_internal_loop(int i, int j, int inner_i, int inner_j)
{
	if(i <= inner_i && inner_i < inner_j && inner_j <= j)
	{
	}
	else
	{
		printf("Order is not right!\n");
		exit(0);
	}


	//if(this->forbid_non_v_emission(i, j) ||
	//	this->forbid_non_v_emission(inner_i, inner_j))
	//{
	//	printf("forbidding internal loop %d, %d, %d, %d\n", i,j, inner_i, inner_j);
	//	return(false);
	//}


	//printf("Checking internal loop %d, %d, %d, %d\n", i,j, inner_i, inner_j);
	// Go over all the nucleotides in the loop.
	for(int _i = i; _i <= j; _i++)
	{
		// If _i or _j hit the 
		if(_i > inner_i && _i < inner_j)
		{
			_i = inner_j;
		}

		if(_i != i && _i != j && _i != inner_i && _i != inner_j &&
			this->forbid_non_v_emission(_i))
		{
			//printf("forbidding internal loop %d, %d, %d, %d\n", i,j, inner_i, inner_j);
			//getc(stdin);
			return(false);
		}

		for(int _j = _i+1; _j <= j; _j++)
		{
			if(_j > inner_i && _j < inner_j)
			{
				_j = inner_j;
			}

			if(!this->same_loop_map[_i][_j])
			{
				//printf("Invalid internal loop %d, %d, %d, %d\n", i,j, inner_i, inner_j);
				////getc(stdin);
				return(false);
			} 
		} // _j loop
	} // _i loop

	//printf("Checked internal loop %d, %d, %d, %d\n", i,j, inner_i, inner_j);

	return(true);
}

t_folding_constraints::~t_folding_constraints()
{
	// Delete memory allocated for maps.
	this->free_maps();

	this->free_ptr_reloc_maps();

	if(this->saturated_str_bps != NULL)
	{
		free(this->saturated_str_bps);
	}

	delete(this->rna_seq);
}

bool t_folding_constraints::force_pairing(int i, int j)
{
	// Swap values if they are not in expected order.
	if(i > j)
	{
		int temp = j;
		j = i;
		i = temp;
	}

	if(i == j ||
		!t_folding_constraints::pairable[this->rna_seq->numseq[i]][this->rna_seq->numseq[j]])
	{
		printf("Cannot force pairing of non-canonical base pair between %c%d and %c%d\n", 
			this->rna_seq->nucs[i], this->rna_seq->numseq[i],
			this->rna_seq->nucs[j], this->rna_seq->numseq[j]);
		//getc(stdin);
		return(false);
	}

	// Following code is borrowed from pfunction.
	int before = 0;
	if (i > 1 && j < this->rna_seq->numofbases)
	{
		before = t_folding_constraints::pairable[this->rna_seq->numseq[i-1]][this->rna_seq->numseq[j+1]];
	}


	int after = 0;
	if((((j-i) > MIN_LOOP + 2) &&
		(j<=this->rna_seq->numofbases))&&(i < this->rna_seq->numofbases))
	{
		after = t_folding_constraints::pairable[this->rna_seq->numseq[i+1]][this->rna_seq->numseq[j-1]];
	}
	else after = 0;

	//if there are no stackable pairs to i.j then don't allow a pair i,j
	if ((before==0) &&
		(after==0)) 
	{
		printf("Cannot enforce pairing of an isolated base pair @ (%d, %d)\n", i, j);
		//getc(stdin);
		return(false);
	}

	// Forbif emission of these indices by arrays other than V.
	this->forbid_non_v[i] = true;
	this->forbid_non_v[j] = true;

	// The pairing of nucleotides seem reasonable.
	for(int _i = 1; _i <= this->rna_seq->numofbases; _i++)
	{
		for(int _j = _i+1; _j <= this->rna_seq->numofbases; _j++)
		{
			if((_i == i && _j > j) ||
				(_i < i && _j == j) ||
				(_i == i && _j == j) ||
				(_i < i && _j > j) ||
				(_i > i && _j < j) ||
				(1 <= _i && _j < i) ||
				(j < _i && _i < _j))
			{
			}
			else
			{
				this->str_coinc_map[_i][_j] = false;
			}

			if((_i >= i && _j <= j) || 
				(_i <= i && _j >= j) ||
				(1 <= _i && _j <= i) ||
				(j <= _i && _i < _j))
			{
			}
			else
			{
				this->same_loop_map[_i][_j] = false;
			}

			//if((_i == i && _j == j)||
			//	(_i < i && _j > j) ||
			//	(_i > i && _j < j))
			if((_i == i && _j == j) ||
				(_i < i && _j > j) ||
				(_i > i && _j < j) ||
				(1 <= _i && _j < i) ||
				(j < _i && _i < _j))
			{
			}
			else
			{
				this->pairing_map[_i][_j] = false;
			}
		}// _j loop
	} // _i loop	

	return(true);
}

void t_folding_constraints::mallocate_ptr_reloc_maps()
{
	if(this->coinc_pointer_relocation_map != NULL)
	{
		this->free_ptr_reloc_maps();
	}

	this->coinc_pointer_relocation_map = (short**)malloc(sizeof(short*) * (rna_seq->numofbases + 3));
	this->paired_pointer_relocation_map = (short**)malloc(sizeof(short*) * (rna_seq->numofbases + 3));

	for(int i = 1; i <= rna_seq->numofbases; i++)
	{
		coinc_pointer_relocation_map[i] = (short*)malloc(sizeof(short) * (rna_seq->numofbases - i + 3));
		coinc_pointer_relocation_map[i] -= i; // Do pointer shift for fold envelope.

		paired_pointer_relocation_map[i] = (short*)malloc(sizeof(short) * (rna_seq->numofbases - i + 3));
		paired_pointer_relocation_map[i] -= i; // Do pointer shift for fold envelope.

		for(int j = i; j <= rna_seq->numofbases; j++)
		{
			coinc_pointer_relocation_map[i][j] = POS_MEM_NOT_EXIST;
			paired_pointer_relocation_map[i][j] = POS_MEM_NOT_EXIST;
		} // Initialization loop
	}
}

void t_folding_constraints::free_ptr_reloc_maps()
{
	if(this->coinc_pointer_relocation_map == NULL)
	{
		return;
	}

	//this->coinc_pointer_relocation_map = (short**)malloc(sizeof(short*) * (rna_seq->numofbases + 3));
	//this->paired_pointer_relocation_map = (short**)malloc(sizeof(short*) * (rna_seq->numofbases + 3));

	for(int i = 1; i <= rna_seq->numofbases; i++)
	{
		this->coinc_pointer_relocation_map[i] += i; // Do pointer shift for fold envelope.
		free(this->coinc_pointer_relocation_map[i]);

		this->paired_pointer_relocation_map[i] += i; // Do pointer shift for fold envelope.
		free(this->paired_pointer_relocation_map[i]);
	}

	free(this->coinc_pointer_relocation_map);
	free(this->paired_pointer_relocation_map);

	this->coinc_pointer_relocation_map = NULL;
	this->paired_pointer_relocation_map = NULL;
}

void t_folding_constraints::compute_ptr_reloc_maps(double pp_threshold, double** pairing_probs)
{
if(_DUMP_FOLDING_CONSTRAINTS_MESSAGES_)
	printf("Generating ptr reloc maps with threshold %lf\n", pp_threshold);

	// Reallocate pointer relocation maps in case they are already allocated.
	if(this->coinc_pointer_relocation_map != NULL ||
		this->paired_pointer_relocation_map != NULL)
	{
		this->free_ptr_reloc_maps();
		this->mallocate_ptr_reloc_maps();
	}

	// Generate pairing map.
	int* sparse_str_bps = (int*)malloc(sizeof(int) * (rna_seq->numofbases + 3));
	for(int i = 1; i <= rna_seq->numofbases; i++)
	{
		sparse_str_bps[i] = 0;
	}

	// Set the ptr relocations for paired nucleotides.	
	//t_folding_constraints* folding_constraints = new t_folding_constraints((t_structure*)rna_seq);
	for(int i = 1; i <= rna_seq->numofbases; i++)
	{
		for(int j = i+1; j <= rna_seq->numofbases; j++)
		{
			double current_pp = 0.0f;
			if(pairing_probs != NULL)
			{
				current_pp = pairing_probs[i][j];
			}
			
			//if(GEQ(current_pp, CONVERT_FROM_LIN(this->cli->str_coinc_pp_threshold)))
			if(current_pp >= pp_threshold)
			{
				// Set sparse structure base pairs.
				if(pp_threshold > 0.50f)
				{
					sparse_str_bps[i] = j;
					sparse_str_bps[j] = i;
					this->force_pairing(i,j);
				}
				else
				{
					printf("Threshold is set to < 0.5, cannot compute a sparse structure with this threshold.\n");
					exit(0);
				}
			}
		} // j loop
	} // i loop

	//// Dump sparse str.
	//char str_fp[MAX_PATH];
	//sprintf(str_fp, "%s_high_pp.ct", this->sequence_id());

	//FILE* f_str = open_f(str_fp, "w");
	//fprintf(f_str, "%d\tENERGY 0\n", this->rna_seq->numofbases);
	//for(int cnt = 1; cnt <= this->rna_seq->numofbases; cnt++)
	//{
	//	// Sth. like following:
	//	//     1 G       0    2   73    1
	//	fprintf(f_str, "%d %c\t%d\t%d\t%d\t%d\n", cnt, this->rna_seq->nucs[cnt], 0, 0, sparse_str_bps[cnt], cnt);
	//}
	//fclose(f_str);

	// Free sparse structure base pairs.
	free(sparse_str_bps);

	// Set the ptr relocations for same loop nucleotides.
	//char same_loop_map_via_MAP_bps_fp[1000];
	//sprintf(same_loop_map_via_MAP_bps_fp, "%s_same_loop_map_via_MAP_bps.txt", this->rna_seq->ctlabel);
	//FILE* f_same_loop_map = open_f(same_loop_map_via_MAP_bps_fp, "w");
	for(int i = 1; i <= rna_seq->numofbases; i++)
	{
if(_DUMP_FOLDING_CONSTRAINTS_MESSAGES_)
		printf("Same loop ptr relocations for i = %d:\n", i);

		int cur_mem_j = 0;

		for(int j = 1; j <= rna_seq->numofbases; j++)
		{
			if(j < i)
			{
				//fprintf(f_same_loop_map, "0");
			}
			else
			{
				if(j == i)
				{
if(_DUMP_FOLDING_CONSTRAINTS_MESSAGES_)
					printf("%d -> %d:\n", j, cur_mem_j);

					//fprintf(f_same_loop_map, "1");
					coinc_pointer_relocation_map[i][j] = cur_mem_j;
					cur_mem_j++;
				}
				else if(this->same_loop_map[i][j])
				{

if(_DUMP_FOLDING_CONSTRAINTS_MESSAGES_)
					printf("%d -> %d:\n", j, cur_mem_j);

					//fprintf(f_same_loop_map, "1");
					coinc_pointer_relocation_map[i][j] = cur_mem_j;
					cur_mem_j++;
				}
				else
				{
					//fprintf(f_same_loop_map, "0");
					coinc_pointer_relocation_map[i][j] = POS_MEM_NOT_EXIST;
				}
			}
		} // j loop

		//fprintf(f_same_loop_map, "\n");
	} // i loop

	// Compute paired ptr relocation maps.
	for(int i = 1; i <= rna_seq->numofbases; i++)
	{
if(_DUMP_FOLDING_CONSTRAINTS_MESSAGES_)
		printf("Paired ptr relocations for i = %d:\n", i);

		int cur_mem_j = 0;

		for(int j = 1; j <= rna_seq->numofbases; j++)
		{
			if(j < i)
			{
				//fprintf(f_pairing_map, "0");
			}
			else
			{
				if(j == i)
				{
if(_DUMP_FOLDING_CONSTRAINTS_MESSAGES_)
					printf("%d -> %d:\n", j, cur_mem_j);

					//fprintf(f_pairing_map, "0");
					paired_pointer_relocation_map[i][j] = cur_mem_j;
					cur_mem_j++;
				}
				else if(this->pairing_map[i][j])
				{

if(_DUMP_FOLDING_CONSTRAINTS_MESSAGES_)
					printf("%d -> %d:\n", j, cur_mem_j);

					//fprintf(f_pairing_map, "1");
					paired_pointer_relocation_map[i][j] = cur_mem_j;
					cur_mem_j++;
				}
				else
				{
					//fprintf(f_pairing_map, "0");
					paired_pointer_relocation_map[i][j] = POS_MEM_NOT_EXIST;
				}
			}
		} // j loop
		//fprintf(f_pairing_map, "\n");
	} // i loop
}

//short** t_folding_constraints::compute_coinc_ptr_reloc_map(char* seq_fp, double pp_threshold, double** pairing_probs)
//{
//	t_structure* rna_seq = new t_structure(seq_fp);
//	short** coinc_ptr_reloc_map = t_folding_constraints::compute_coinc_ptr_reloc_map(rna_seq, pp_threshold, pairing_probs);
//	delete(rna_seq);
//	return(coinc_ptr_reloc_map);
//}
//
//short** t_folding_constraints::compute_paired_ptr_reloc_map(char* seq_fp, double pp_threshold, double** pairing_probs)
//{
//	t_structure* rna_seq = new t_structure(seq_fp);
//	short** paired_ptr_reloc_map = t_folding_constraints::compute_paired_ptr_reloc_map(rna_seq, pp_threshold, pairing_probs);
//	delete(rna_seq);
//	return(paired_ptr_reloc_map);
//}


//// Assumes that the pairing_probs are in linear space.
//short** t_folding_constraints::compute_paired_ptr_reloc_map(t_structure* rna_seq, double pp_threshold, double** pairing_probs)
//{
//	printf("Generating ptr reloc maps with threshold %lf\n", pp_threshold);
//
//	// Generate pointer relocation maps.
//	short** paired_pointer_relocation_map = (short**)malloc(sizeof(short*) * (rna_seq->numofbases + 3));
//
//	for(int i = 1; i <= rna_seq->numofbases; i++)
//	{
//		paired_pointer_relocation_map[i] = (short*)malloc(sizeof(short) * (rna_seq->numofbases - i + 3));
//		paired_pointer_relocation_map[i] -= i; // Do pointer shift for fold envelope.
//
//		for(int j = i; j <= rna_seq->numofbases; j++)
//		{
//			paired_pointer_relocation_map[i][j] = POS_MEM_NOT_EXIST;
//		} // Initialization loop
//	}
//
//	// Generate pairing map.
//	int* sparse_str_bps = (int*)malloc(sizeof(int) * (rna_seq->numofbases + 3));
//	for(int i = 1; i <= rna_seq->numofbases; i++)
//	{
//		sparse_str_bps[i] = 0;
//	}
//
//	// Set the ptr relocations for paired nucleotides.
//	t_folding_constraints* folding_constraints = new t_folding_constraints((t_structure*)rna_seq);
//	for(int i = 1; i <= rna_seq->numofbases; i++)
//	{
//		for(int j = i+1; j <= rna_seq->numofbases; j++)
//		{
//			double current_pp = pairing_probs[i][j];
//			
//			//if(GEQ(current_pp, CONVERT_FROM_LIN(this->cli->str_coinc_pp_threshold)))
//			if(current_pp >= pp_threshold)
//			{
//				// Set sparse structure base pairs.
//				if(pp_threshold > 0.50f)
//				{
//					sparse_str_bps[i] = j;
//					sparse_str_bps[j] = i;
//					folding_constraints->force_pairing(i,j);
//				}
//				else
//				{
//					printf("Threshold is set to < 0.5, cannot compute a sparse structure with this threshold.\n");
//					exit(0);
//				}
//			}
//		} // j loop
//	} // i loop
//
//	//// Dump sparse str.
//	//char str_fp[MAX_PATH];
//	//sprintf(str_fp, "%s_high_pp.ct", this->sequence_id());
//
//	//FILE* f_str = open_f(str_fp, "w");
//	//fprintf(f_str, "%d\tENERGY 0\n", this->rna_seq->numofbases);
//	//for(int cnt = 1; cnt <= this->rna_seq->numofbases; cnt++)
//	//{
//	//	// Sth. like following:
//	//	//     1 G       0    2   73    1
//	//	fprintf(f_str, "%d %c\t%d\t%d\t%d\t%d\n", cnt, this->rna_seq->nucs[cnt], 0, 0, sparse_str_bps[cnt], cnt);
//	//}
//	//fclose(f_str);
//
//	// Free sparse structure base pairs.
//	free(sparse_str_bps);
//
//	// Set the ptr relocations for same loop nucleotides.
//	//char pairing_map_via_MAP_bps_fp[1000];
//	//sprintf(pairing_map_via_MAP_bps_fp, "%s_pairing_map_via_MAP_bps.txt", this->rna_seq->ctlabel);
//	//FILE* f_pairing_map = open_f(pairing_map_via_MAP_bps_fp, "w");
//	for(int i = 1; i <= rna_seq->numofbases; i++)
//	{
//if(_DUMP_FOLDING_CONSTRAINTS_MESSAGES_)
//		printf("Paired ptr relocations for i = %d:\n", i);
//
//		int cur_mem_j = 0;
//
//		for(int j = 1; j <= rna_seq->numofbases; j++)
//		{
//			if(j < i)
//			{
//				//fprintf(f_pairing_map, "0");
//			}
//			else
//			{
//				if(j == i)
//				{
//if(_DUMP_FOLDING_CONSTRAINTS_MESSAGES_)
//					printf("%d -> %d:\n", j, cur_mem_j);
//
//					//fprintf(f_pairing_map, "0");
//					paired_pointer_relocation_map[i][j] = cur_mem_j;
//					cur_mem_j++;
//				}
//				else if(folding_constraints->pairing_map[i][j])
//				{
//
//if(_DUMP_FOLDING_CONSTRAINTS_MESSAGES_)
//					printf("%d -> %d:\n", j, cur_mem_j);
//
//					//fprintf(f_pairing_map, "1");
//					paired_pointer_relocation_map[i][j] = cur_mem_j;
//					cur_mem_j++;
//				}
//				else
//				{
//					//fprintf(f_pairing_map, "0");
//					paired_pointer_relocation_map[i][j] = POS_MEM_NOT_EXIST;
//				}
//			}
//		} // j loop
//		//fprintf(f_pairing_map, "\n");
//	} // i loop
//
//	delete(folding_constraints);
//
//	return(paired_pointer_relocation_map);
//}
//
//short** t_folding_constraints::compute_coinc_ptr_reloc_map(t_structure* rna_seq, double pp_threshold, double** pairing_probs)
//{
//	printf("Generating ptr reloc maps with threshold %lf\n", pp_threshold);
//
//	// Generate pointer relocation maps.
//	short** coinc_pointer_relocation_map = (short**)malloc(sizeof(short*) * (rna_seq->numofbases + 3));
//
//	for(int i = 1; i <= rna_seq->numofbases; i++)
//	{
//		coinc_pointer_relocation_map[i] = (short*)malloc(sizeof(short) * (rna_seq->numofbases - i + 3));
//		coinc_pointer_relocation_map[i] -= i; // Do pointer shift for fold envelope.
//
//		for(int j = i; j <= rna_seq->numofbases; j++)
//		{
//			coinc_pointer_relocation_map[i][j] = POS_MEM_NOT_EXIST;
//		} // Initialization loop
//	}
//
//	// Generate pairing map.
//	int* sparse_str_bps = (int*)malloc(sizeof(int) * (rna_seq->numofbases + 3));
//	for(int i = 1; i <= rna_seq->numofbases; i++)
//	{
//		sparse_str_bps[i] = 0;
//	}
//
//	// Set the ptr relocations for paired nucleotides.	
//	t_folding_constraints* folding_constraints = new t_folding_constraints((t_structure*)rna_seq);
//	for(int i = 1; i <= rna_seq->numofbases; i++)
//	{
//		for(int j = i+1; j <= rna_seq->numofbases; j++)
//		{
//			double current_pp = pairing_probs[i][j];
//			
//			//if(GEQ(current_pp, CONVERT_FROM_LIN(this->cli->str_coinc_pp_threshold)))
//			if(current_pp >= pp_threshold)
//			{
//				// Set sparse structure base pairs.
//				if(pp_threshold > 0.50f)
//				{
//					sparse_str_bps[i] = j;
//					sparse_str_bps[j] = i;
//					folding_constraints->force_pairing(i,j);
//				}
//				else
//				{
//					printf("Threshold is set to < 0.5, cannot compute a sparse structure with this threshold.\n");
//					exit(0);
//				}
//			}
//		} // j loop
//	} // i loop
//
//	//// Dump sparse str.
//	//char str_fp[MAX_PATH];
//	//sprintf(str_fp, "%s_high_pp.ct", this->sequence_id());
//
//	//FILE* f_str = open_f(str_fp, "w");
//	//fprintf(f_str, "%d\tENERGY 0\n", this->rna_seq->numofbases);
//	//for(int cnt = 1; cnt <= this->rna_seq->numofbases; cnt++)
//	//{
//	//	// Sth. like following:
//	//	//     1 G       0    2   73    1
//	//	fprintf(f_str, "%d %c\t%d\t%d\t%d\t%d\n", cnt, this->rna_seq->nucs[cnt], 0, 0, sparse_str_bps[cnt], cnt);
//	//}
//	//fclose(f_str);
//
//	// Free sparse structure base pairs.
//	free(sparse_str_bps);
//
//	// Set the ptr relocations for same loop nucleotides.
//	//char same_loop_map_via_MAP_bps_fp[1000];
//	//sprintf(same_loop_map_via_MAP_bps_fp, "%s_same_loop_map_via_MAP_bps.txt", this->rna_seq->ctlabel);
//	//FILE* f_same_loop_map = open_f(same_loop_map_via_MAP_bps_fp, "w");
//	for(int i = 1; i <= rna_seq->numofbases; i++)
//	{
//if(_DUMP_FOLDING_CONSTRAINTS_MESSAGES_)
//		printf("Same loop ptr relocations for i = %d:\n", i);
//
//		int cur_mem_j = 0;
//
//		for(int j = 1; j <= rna_seq->numofbases; j++)
//		{
//			if(j < i)
//			{
//				//fprintf(f_same_loop_map, "0");
//			}
//			else
//			{
//				if(j == i)
//				{
//if(_DUMP_FOLDING_CONSTRAINTS_MESSAGES_)
//					printf("%d -> %d:\n", j, cur_mem_j);
//
//					//fprintf(f_same_loop_map, "1");
//					coinc_pointer_relocation_map[i][j] = cur_mem_j;
//					cur_mem_j++;
//				}
//				else if(folding_constraints->same_loop_map[i][j])
//				{
//
//if(_DUMP_FOLDING_CONSTRAINTS_MESSAGES_)
//					printf("%d -> %d:\n", j, cur_mem_j);
//
//					//fprintf(f_same_loop_map, "1");
//					coinc_pointer_relocation_map[i][j] = cur_mem_j;
//					cur_mem_j++;
//				}
//				else
//				{
//					//fprintf(f_same_loop_map, "0");
//					coinc_pointer_relocation_map[i][j] = POS_MEM_NOT_EXIST;
//				}
//			}
//		} // j loop
//
//		//fprintf(f_same_loop_map, "\n");
//	} // i loop
//
//	//fclose(f_same_loop_map);
//	delete(folding_constraints);
//
//	return(coinc_pointer_relocation_map);
//}






#include <stdio.h>
#include <stdlib.h>
#include "structure_object.h"
#include <string.h>
#include "../utils/ansi_string/ansi_string.h"
#include "../utils/file/utils.h"
#include <ctype.h>
#include <iostream>

char* _fgets(char* buffer, int n_max, FILE* f);

// Default constructor.
t_structure::t_structure()
{
    this->numofbases = 0;
	this->numseq = NULL;
	this->nucs = NULL;
	this->basepr = NULL;
	this->ctlabel = NULL;
	this->unpaired_forced = NULL;

	this->danglings_on_branch = NULL;
	this->danglings_on_mb_closure = NULL;
	this->stackings_on_branch = NULL;
	this->stackings_on_mb_closure = NULL;
}

void t_structure::check_set_label()
{
	//printf("Validating %s\n", this->ctlabel);
	char invalid_label_chars[] = "\"\\/ '?|<>%%^&@#!*+\n\t\r,";
	for(int i_valid = 0; i_valid < (int)strlen(this->ctlabel); i_valid++)
	{
		for(int i_check = 0; i_check < (int)strlen(invalid_label_chars); i_check++)
		{
			// Check if this character needs replacement.
			if(this->ctlabel[i_valid] == invalid_label_chars[i_check])
			{
				this->ctlabel[i_valid] = '_';
			}
		} // checking loop
	} // validation loop
}

// Constructor with known ct_label and nucleotides.
t_structure::t_structure(const char* ct_label, vector<char>* nuc_vector, bool fix_seq_chars)
{
    this->numofbases = 0;
	this->numseq = NULL;
	this->nucs = NULL;
	this->basepr = NULL;
	this->ctlabel = NULL;
	this->unpaired_forced = NULL;

	this->danglings_on_branch = NULL;
	this->danglings_on_mb_closure = NULL;
	this->stackings_on_branch = NULL;
	this->stackings_on_mb_closure = NULL;

	//printf("Length of sequence is %d\n", this->numofbases);
	this->numofbases = (int)nuc_vector->size();
	this->numseq = (int*)malloc(sizeof(int) * (this->numofbases + 1));
	this->nucs = (char*)malloc(sizeof(char) * (this->numofbases + 2));
	this->basepr = (int*)malloc(sizeof(int) * (this->numofbases + 1));
	this->ctlabel = (char*)malloc(sizeof(char) * (strlen(ct_label) + 2));
	this->unpaired_forced = (bool*)malloc(sizeof(bool) * (this->numofbases + 2));

	this->danglings_on_branch = (int*)malloc(sizeof(int) * (this->numofbases + 3));
	this->danglings_on_mb_closure = (int*)malloc(sizeof(int) * (this->numofbases + 3));
	this->stackings_on_branch = (int*)malloc(sizeof(int) * (this->numofbases + 3));
	this->stackings_on_mb_closure = (int*)malloc(sizeof(int) * (this->numofbases + 3));
	
	// Copy sequence and base pairing information.
	for(int i = 0; i <= this->numofbases; i++)
	{
		this->basepr[i] = 0;
		this->danglings_on_branch[i] = 0;
		this->danglings_on_mb_closure[i] = 0;
		this->stackings_on_branch[i] = 0;
		this->stackings_on_mb_closure[i] = 0;
	}


	strcpy(this->ctlabel, ct_label);
	if (fix_seq_chars) this->check_set_label();

	this->nucs[0] = '#';
	this->numseq[0] = 0;

	int i_nuc = 1;

	for(int i_vect = 0; i_vect < (int)nuc_vector->size(); i_vect++)
	{	
		// Process this nuc.
		if(nuc_vector->at(i_vect) != '\n' && nuc_vector->at(i_vect) != ' ' && nuc_vector->at(i_vect) != '1')
		{
			map_nuc_IUPAC_code(nuc_vector->at(i_vect), this->nucs[i_nuc], this->numseq[i_nuc], this->unpaired_forced[i_nuc]);

			this->basepr[i_nuc] = 0; // No base pairing information.

			i_nuc++;
		}
		else
		{
			// This is an invalid character, do not process.
		}
	}

	// Make these a string.
	//this->nucs[nuc_vector->size()+1] = 0;
	this->nucs[i_nuc] = 0;
}

/*
Copy constructor.
*/
t_structure::t_structure(t_structure* str_2_copy)
{
    this->numofbases = 0;
	this->numseq = NULL;
	this->nucs = NULL;
	this->basepr = NULL;
	this->ctlabel = NULL;
	this->unpaired_forced = NULL;

	this->ctlabel = (char*)malloc( (strlen(str_2_copy->ctlabel) + 3) * sizeof(char) );
	strcpy(this->ctlabel, str_2_copy->ctlabel);

	this->numofbases = str_2_copy->numofbases;
	this->nucs = (char*)malloc(sizeof(char) * (this->numofbases + 3));
	this->numseq = (int*)malloc(sizeof(int) * (this->numofbases + 3));
	this->basepr = (int*)malloc(sizeof(int) * (this->numofbases + 3));

	this->danglings_on_branch = (int*)malloc(sizeof(int) * (this->numofbases + 3));
	this->danglings_on_mb_closure = (int*)malloc(sizeof(int) * (this->numofbases + 3));
	this->stackings_on_branch = (int*)malloc(sizeof(int) * (this->numofbases + 3));
	this->stackings_on_mb_closure = (int*)malloc(sizeof(int) * (this->numofbases + 3));
	this->unpaired_forced = (bool*)malloc(sizeof(bool) * (this->numofbases + 3));
	
	// Copy sequence and base pairing information.
	for(int i = 0; i <= this->numofbases; i++)
	{
		// Following are always allocated in a t_structure.
		this->nucs[i] = str_2_copy->nucs[i];
		this->numseq[i] = str_2_copy->numseq[i];
		this->basepr[i] = str_2_copy->basepr[i];
		this->unpaired_forced[i] = str_2_copy->unpaired_forced[i];

		// If one of these exists, all of them exists.
		if(str_2_copy->danglings_on_branch != NULL)
		{
			this->danglings_on_branch[i] = str_2_copy->danglings_on_branch[i];
			this->danglings_on_mb_closure[i] = str_2_copy->danglings_on_mb_closure[i];
			this->stackings_on_branch[i] = str_2_copy->stackings_on_branch[i];
			this->stackings_on_mb_closure[i] = str_2_copy->stackings_on_mb_closure[i];
		}
		else
		{
			this->danglings_on_branch[i] = 0;
			this->danglings_on_mb_closure[i] = 0;
			this->stackings_on_branch[i] = 0;
			this->stackings_on_mb_closure[i] = 0;
		}
	}

	// End nucs as a string.
	this->nucs[this->numofbases + 1] = 0;
}

// Read a fasta file that contains sequence information for multiple fasta files and return all of them in a vector.
vector<t_structure*>* t_structure::read_multi_seq(const char* const multi_seq_fp, bool fix_seq_chars)
{
	vector<t_structure*>* seqs = new vector<t_structure*>();

	FILE* f_multi_seq = open_f(multi_seq_fp, "r");
	if(f_multi_seq == NULL)
	{
		printf("Could not find the input file @ %s.\n", multi_seq_fp);
		exit(0);
	}

	// Read the file and load information. 
	vector<char>* cur_nucs = new vector<char>();
	char cur_label[MAX_HEADER_LENGTH];
	char cur_line[MAX_HEADER_LENGTH];
	while(1)
	{
		// Read current line.
		if(fgets(cur_line, MAX_HEADER_LENGTH, f_multi_seq) == NULL)
		{
			// Save the last sequence in the sequence list and initiate a new sequence.
			if(cur_nucs->size() > 0)
			{
				t_structure* new_seq = new t_structure(cur_label, cur_nucs, fix_seq_chars);
				seqs->push_back(new_seq);
			}

			delete(cur_nucs);
			break;
		}

		// Get rid of the new line, if there is one.
		if(strlen(cur_line) > 0 && cur_line[strlen(cur_line) - 1] == '\n')
		{
			cur_line[strlen(cur_line) - 1] = 0;
		}

		if(strlen(cur_line) > 0)
		{
			// if starts with a '>', then a new sequence is initiated.
			if(cur_line[0] == '>')
			{
				// Save the last sequence in the sequence list and initiate a new sequence.
				if(cur_nucs->size() > 0)
				{
					t_structure* new_seq = new t_structure(cur_label, cur_nucs, fix_seq_chars);
					seqs->push_back(new_seq);
				}

				// Read the label from the remaining of the line.
				strcpy(cur_label, &cur_line[1]);

				// Empty current nucleotides for loading next sequence, if there is any.
				cur_nucs->clear();
			}
			else if(cur_line[0] == ';')
			{
				// Save the last sequence in the sequence list and initiate a new sequence.
				if(cur_nucs->size() > 0)
				{
					//printf("instantiating with new label: %s\n", cur_label);
					t_structure* new_str = new t_structure(cur_label, cur_nucs, fix_seq_chars);
					seqs->push_back(new_str);
				}

				// Read the label from the next line.
				fgets(cur_label, MAX_HEADER_LENGTH, f_multi_seq);
				if(cur_label[strlen(cur_label)-1] == '\n')
				{
					cur_label[strlen(cur_label)-1] = 0;
				}

				//printf("Read new label: %s\n", cur_label);

				// Empty current nucleotides for loading next sequence, if there is any.
				cur_nucs->clear();
			}
			else
			{
				// This is sequence data, copy the sequence data and continue, no input validation here.
				for(int i_cpy = 0; i_cpy < (int)strlen(cur_line); i_cpy++)
				{
					// This is a necessity coming from .seq file specifications. All .seq files end with a '1' character.
					if(cur_line[i_cpy] != '1' &&
						cur_line[i_cpy] != ' ' &&
						cur_line[i_cpy] != '\n' &&
						cur_line[i_cpy] != '\t')
					{
						cur_nucs->push_back(cur_line[i_cpy]);
					}
				} // Copy the nucleotides.
			} // label/nuc data check.
		} // Length check for current line.
	}
	fclose(f_multi_seq);
	return(seqs);
}

bool t_structure::cmp_seq(t_structure *str1, t_structure *str2)
{
	if(str1->numofbases != str2->numofbases)
	{
		return(false);
	}

	for(int i_nuc = 1; i_nuc <= str1->numofbases; i_nuc++)
	{
		if(str1->nucs[i_nuc] != str2->nucs[i_nuc])
		{
			return(false);
		}
	}

	return(true);
}

t_structure::t_structure(const char* _fp)
{
	//char* temp_fp = (char*)malloc( (strlen(_fp) + 3) * sizeof(char));
	//strcpy(temp_fp, _fp);
	t_string* file_path = new t_string(_fp);
	t_string_tokens* file_path_tokens = file_path->tokenize_by_chars(".");

	t_string* last_token = file_path_tokens->back();

	//printf("last token is: %s\n", last_token);

	char seq_str[] = "seq";
	char ct_str[] = "ct";
	char fasta_str[] = "fasta";
	if(last_token->length() == 3)
	{
		bool is_seq = true;
		if(!last_token->compare_ci(seq_str))
		{
			is_seq = false;
		}

		if(is_seq)
		{
			this->openseq(_fp);
		}
	}
	else if(last_token->length() == 2)
	{
		bool is_ct = true;
		if(!last_token->compare_ci(ct_str))
		{
			is_ct = false;
		}

		if(is_ct)
		{
			this->openct(_fp);
		}
	}
	else if(last_token->length() == strlen(fasta_str))
	{
		//printf("A fasta file!\n");
		bool is_seq = true;

		if(!last_token->compare_ci(fasta_str))
		{
			is_seq = false;
		}

		if(is_seq)
		{
			this->openfasta(_fp);
		}

	}
	else
	{
		printf("Could not determine file type of input for %s @ %s(%d).\n", _fp, __FILE__, __LINE__);
		exit(0);
	}

	// Initialize dangling/stacking information.
	if(this->danglings_on_branch == NULL)
	{
		this->danglings_on_branch = (int*)malloc(sizeof(int) * (this->numofbases + 3));
		this->danglings_on_mb_closure = (int*)malloc(sizeof(int) * (this->numofbases + 3));
		this->stackings_on_branch = (int*)malloc(sizeof(int) * (this->numofbases + 3));
		this->stackings_on_mb_closure = (int*)malloc(sizeof(int) * (this->numofbases + 3));
		
		// Copy sequence and base pairing information.
		for(int i = 0; i <= this->numofbases; i++)
		{
			this->danglings_on_branch[i] = 0;
			this->danglings_on_mb_closure[i] = 0;
			this->stackings_on_branch[i] = 0;
			this->stackings_on_mb_closure[i] = 0;
		}
	}

	//free(temp_fp);
	file_path->clean_tokens(file_path_tokens);
	delete(file_path);
}

t_structure::~t_structure()
{
	if(this->numseq != NULL)
	{
		free(this->numseq);
	}

	if(this->nucs != NULL)
	{
		free(this->nucs);
	}

	if(this->basepr != NULL)
	{
		free(this->basepr);
	}

	if(this->ctlabel != NULL)
	{
		free(this->ctlabel);
	}

	if(this->danglings_on_branch != NULL)
	{
		free(this->danglings_on_branch);
	}

	if(this->danglings_on_mb_closure != NULL)
	{
		free(this->danglings_on_mb_closure);
	}

	if(this->stackings_on_branch != NULL)
	{
		free(this->stackings_on_branch);
	}

	if(this->stackings_on_mb_closure != NULL)
	{
		free(this->stackings_on_mb_closure);
	}

	if(this->unpaired_forced != NULL)
	{
		free(this->unpaired_forced);
	}
}

/*

A 
A adenine 

C 
C cytosine 

G 
G guanine 

T 
T thymine 

U 
U uracil 
   

R 
A or G purine 

Y 
C or T (U) pyrimidine 
   

M 
A or C amino 

K 
G or T (U) keto 

S 
C or G strong (3 H bonds) 

W 
A or T (U) weak (2 H bonds) 
   

B 
C or G or T (U) not A 

D 
A or G or T (U) not C 

H 
A or C or T (U) not G 

V 
A or C or G not T (U) 
   

N 
A or C or G or T (U) any nucleotide 

*/
void t_structure::map_nuc_IUPAC_code(char raw_nuc, 
									 char &trans_nuc, 
									 int &num, 
									 bool& force_unpaired)
{
	if(raw_nuc == 'a' || raw_nuc == 'c' || raw_nuc == 'g' || raw_nuc == 'u' || raw_nuc == 't')
	{
		force_unpaired = true;
	}
	else
	{
		force_unpaired = false;
	}

	if (toupper(raw_nuc) == 'A') 
	{
		trans_nuc = raw_nuc;
		num=1;
	}
	else if(toupper(raw_nuc) == 'B')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'C')
	{
		trans_nuc = raw_nuc;
		num = 2;
	}
	else if(toupper(raw_nuc) == 'D')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'G')
	{
		trans_nuc = raw_nuc;
		num = 3;
	}
	else if(toupper(raw_nuc) == 'H')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'I')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'K')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'M')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'N')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'R')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'S')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'T')
	{
		trans_nuc = raw_nuc;
		num = 4;
	}
	else if(toupper(raw_nuc) == 'U')
	{
		trans_nuc = raw_nuc;
		num = 4;
	}
	else if(toupper(raw_nuc) == 'V')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'W')
	{
		trans_nuc = 'N';
		num = 0;
	}
        else if(toupper(raw_nuc) == 'X')
        {
                trans_nuc = 'N';
                num = 0;
        }
	else if(toupper(raw_nuc) == 'Y')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else
	{
		trans_nuc = 'N';
		num = 0;
	}

	if(num == 0)
	{
		printf("Found %c\n", raw_nuc);
	}
}

void t_structure::openct(const char* ct_fp)
{
	FILE* ct_file = open_f(ct_fp, "r");
	if(ct_file == NULL)
	{
		printf("ct file %s does not exist @ %s(%d).\n", ct_fp, __FILE__, __LINE__);
		exit(1);
	}

	// Allocate header buffer.
	this->ctlabel = (char*)malloc(sizeof(char) * MAX_HEADER_LENGTH);

	// Read first line
	fscanf(ct_file, "%d", &this->numofbases);

	// Read remaining of the line, contains new line character at the end of label.
	fgets(this->ctlabel, MAX_HEADER_LENGTH, ct_file); // Read remaining of the line.
	if(this->ctlabel[strlen(this->ctlabel) - 1] == '\n')
	{
		this->ctlabel[strlen(this->ctlabel) - 1] = 0;
	}
	this->check_set_label();

	//printf("ct label: %s\n", this->ctlabel);

	this->numseq = (int*)malloc(sizeof(int) * (this->numofbases + 3));
	this->nucs = (char*)malloc(sizeof(char) * (this->numofbases + 3));
	this->basepr = (int*)malloc(sizeof(int) * (this->numofbases + 3));
	this->danglings_on_branch = (int*)malloc(sizeof(int) * (this->numofbases + 3));
	this->danglings_on_mb_closure = (int*)malloc(sizeof(int) * (this->numofbases + 3));
	this->stackings_on_branch = (int*)malloc(sizeof(int) * (this->numofbases + 3));
	this->stackings_on_mb_closure = (int*)malloc(sizeof(int) * (this->numofbases + 3));
	this->unpaired_forced = (bool*)malloc(sizeof(bool) * (this->numofbases + 2));

	for(int i = 0; i <= this->numofbases; i++)
	{
		this->basepr[i] = 0;
		this->danglings_on_branch[i] = 0;
		this->danglings_on_mb_closure[i] = 0;
		this->stackings_on_branch[i] = 0;
		this->stackings_on_mb_closure[i] = 0;
	}

	int* dangles = (int*)malloc(sizeof(int) * (this->numofbases + 3));
	int* stacks = (int*)malloc(sizeof(int) * (this->numofbases + 3));

	// Read sequence data.
	// Must read base pairing before dangles/stacks can be resolved from file.
	for(int i = 1; i <= this->numofbases; i++)
	{
		int index;
		int some_val1;
		char raw_nuc;

		//                1  G 0  2  120 1
		fscanf(ct_file, "%d %c %d %d %d %d", &index, &raw_nuc, &dangles[i], &stacks[i], &this->basepr[i], &some_val1);

		//if(this->nucs[i] == 'a' ||
		//	this->nucs[i] == 'c' ||
		//	this->nucs[i] == 'g' ||
		//	this->nucs[i] == 'u' ||
		//	this->nucs[i] == 't')
		//{
		//	this->unpaired_forced[i] = true;
		//}
		//else
		//{
		//	this->unpaired_forced[i] = false;
		//}

		//printf("%c", this->nucs[i]);

		/*
		The danglings on external loop branches are buffered as 
		danglings on branch. Note that there cannot be a stacking on external loop closure
		because by definition external loop is not closed.
		*/

		// Convert nucleotide symbols into indices: XACGUI -> 012345 
		// refer to IUPAC nucleotide symbols for more information:
		// http://www.mun.ca/biochem/courses/3107/symbols.html
		//if (toupper(this->nucs[i]) == 'A' || toupper(this->nucs[i]) == 'B') 
		//	this->numseq[i]=1;
		//else if (toupper(this->nucs[i]) == 'C' || toupper(this->nucs[i]) == 'Z') 
		//	this->numseq[i]=2;
		//else if (toupper(this->nucs[i]) == 'G' || toupper(this->nucs[i]) == 'H') 
		//	this->numseq[i]=3;
		//else if (toupper(this->nucs[i]) == 'U' || toupper(this->nucs[i]) == 'T' || toupper(this->nucs[i]) == 'V' || toupper(this->nucs[i]) == 'W' ) 
		//	this->numseq[i]=4;
		//else if (toupper(this->nucs[i]) == 'I') 
		//	this->numseq[i]=5;
		//else 
		//	this->numseq[i]=0;

		map_nuc_IUPAC_code(raw_nuc, this->nucs[i], this->numseq[i], this->unpaired_forced[i]);

		//printf("%d\n", this->basepr[i]);
	}

#undef _USE_STACKING_INFO_
#ifdef _USE_STACKING_INFO_
	// Resolve stacks and dangles.
	for(int i = 1; i <= this->numofbases; i++)
	{
		// Dangling?
		if(dangles[i] != 0)
		{
			// Dangle on branch?
			if(dangles[i] == i+1)
			{
				if(this->basepr[i+1] == 0)
				{
					printf("Dangling of %d on unpaired nucleotide %d.\n", i, i+1);
					exit(0);
				}

				if(this->basepr[i+1] > i+1)
				{
					this->danglings_on_branch[i] = i+1;
				}
				else
				{
					this->danglings_on_mb_closure[i] = i+1;
				}
			}
			// Dangle on mbl closure?
			if(dangles[i] == i-1)
			{
				if(this->basepr[i-1] == 0)
				{
					printf("Dangling of %d on unpaired nucleotide %d.\n", i, i-1);
					exit(0);
				}

				if(this->basepr[i-1] > i-1)
				{
					this->danglings_on_mb_closure[i] = i-1;
				}
				else
				{
					this->danglings_on_branch[i] = i-1;
				}
			}
		}

		// Stacking?
		if(stacks[i] != 0)
		{
			// stack on branch?
			if(stacks[i] == i+1)
			{
				if(this->basepr[i+1] == 0)
				{
					printf("Stacking of %d on unpaired nucleotide %d.\n", i, i+1);
					exit(0);
				}

				if(this->basepr[i+1] > i+1)
				{
					this->stackings_on_branch[i] = i+1;
				}
				else
				{
					this->stackings_on_mb_closure[i] = i+1;
				}
			}
			// stack on mbl closure?
			if(stacks[i] == i-1)
			{
				if(this->basepr[i-1] == 0)
				{
					printf("Stacking of %d on unpaired nucleotide %d.\n", i, i-1);
					exit(0);
				}

				if(this->basepr[i-1] > i-1)
				{
					this->stackings_on_mb_closure[i] = i-1;
				}
				else
				{
					this->stackings_on_branch[i] = i-1;
				}
			}
		}
	}

	// Do a sanity check on dangles and stacks.
	for(int i = 1; i < this->numofbases; i++)
	{
		if(this->stackings_on_branch[i] == i+1)
		{
			int current_j = this->basepr[i+1];
			if(current_j == 0 || this->stackings_on_branch[current_j+1] != current_j)
			{
				printf("Stacking check failed for stacking of %d on %d\n", i, i+1);
			}
		}

		if(this->stackings_on_mb_closure[i] == i+1)
		{
			int current_j = this->basepr[i+1];
			if(current_j == 0 || this->stackings_on_mb_closure[current_j+1] != current_j)
			{
				printf("Stacking check failed for stacking of %d on %d\n", i, i+1);
			}
		}
	}
#endif // _USE_STACKING_INFO_

	free(dangles);
	free(stacks);

	fclose(ct_file);
}

// It is very important to make sure that a seq file is in following format:
// ; ...
// [ THIS LINE SHOULD NOT CONTAIN SEQUENCE DATA, ITS SHOULD BE A LABEL OR EMPTY LINE]
// [ EMPTY LINE OR SEQUENCE DATA]
void t_structure::openseq(const char* seq_fp)
{
	// Very strict measure: Exit is sequence file is not verifiable.
	if(!this->verify_seq(seq_fp))
	{
		printf("Could not verify sequence file %s @ %s(%d)\n", seq_fp, __FILE__, __LINE__);
		exit(1);
	}

	FILE* seq_file = open_f(seq_fp, "r");
	if(seq_file == NULL)
	{
		printf("seq file %s does not exist @ %s(%d).\n", seq_fp, __FILE__, __LINE__);
		exit(1);
	}

	this->numseq = NULL;
	this->nucs = NULL;
	this->basepr = NULL;
	this->danglings_on_branch = NULL;
	this->danglings_on_mb_closure = NULL;
	this->stackings_on_branch = NULL;
	this->stackings_on_mb_closure = NULL;
	this->unpaired_forced = NULL;

	char line_buffer[MAX_HEADER_LENGTH];
	fgets(line_buffer, MAX_HEADER_LENGTH, seq_file);
	while(line_buffer[0] == ';')
	{
		fgets(line_buffer, MAX_HEADER_LENGTH, seq_file);
	}

	// Read label, contains new line character at the end of label.
	this->ctlabel = (char*)malloc(sizeof(char) * MAX_HEADER_LENGTH);
	strcpy(this->ctlabel, line_buffer);
	if(this->ctlabel[strlen(this->ctlabel) - 1] == '\n')
	{
		this->ctlabel[strlen(this->ctlabel) - 1] = 0;
	}
	this->check_set_label();
	
	//printf("seq label: %s\n", this->ctlabel);

	// Read and determine length of sequence.
	char cur_char = 0;
	this->numofbases = 0;

	// Start reading sequence data.
	while(1)
	{
		int ret = fscanf(seq_file, "%c", &cur_char);
		if(ret == EOF)
		{
			break;
		}

		if(cur_char == '1')
		{
			break;
		}

		if(cur_char != '\n' && cur_char != ' ')
		{
			this->numofbases++;
		}
	}

	//printf("Length of sequence is %d\n", this->numofbases);
	this->numseq = (int*)malloc(sizeof(int) * (this->numofbases + 1));
	this->nucs = (char*)malloc(sizeof(char) * (this->numofbases + 2));
	this->basepr = (int*)malloc(sizeof(int) * (this->numofbases + 1));
	this->unpaired_forced = (bool*)malloc(sizeof(bool) * (this->numofbases + 2));

	// Set file position to data position.
	// Cannot use fsetpos and fgetpos because for some reason they are messing up indices
	// when a linux text file is taken to a windows machine.
	fseek(seq_file, 0, SEEK_SET);

	// Read all information again before sequence data.
	fgets(line_buffer, MAX_HEADER_LENGTH, seq_file);
	while(line_buffer[0] == ';')
	{
		fgets(line_buffer, MAX_HEADER_LENGTH, seq_file);
	}

	this->nucs[0] = '#';
	int i = 1; // Sequence index, starts from 1.

	// Start reading sequence data.
	while(1)
	{
		// Read and validate input.
		int ret = fscanf(seq_file, "%c", &cur_char);
		if(ret == EOF)
		{
			break;
		}

		// Check end of sequence marker.
		if(cur_char == '1')
		{
			break;
		}

		// Process this nuc.
		if(cur_char != '\n' && cur_char != ' ')
		{
			map_nuc_IUPAC_code(cur_char, this->nucs[i], this->numseq[i], this->unpaired_forced[i]);			

			this->basepr[i] = 0; // No base pairing information.

			//printf("%c %d\n", this->nucs[i], this->numseq[i]);

			i++;
		}
	}

	// This is for ending sequences.
	this->nucs[i] = 0; 

	fclose(seq_file);
}

void t_structure::openfasta(const char* fasta_fp)
{
	// Very strict measure: Exit is sequence file is not verifiable.
	if(!this->verify_seq(fasta_fp))
	{
		printf("Could not verify sequence file %s @ %s(%d)\n", fasta_fp, __FILE__, __LINE__);
		exit(1);
	}

	FILE* fasta_file = open_f(fasta_fp, "r");
	if(fasta_file == NULL)
	{
		printf("fasta file %s does not exist @ %s(%d).\n", fasta_fp, __FILE__, __LINE__);
		exit(1);
	}

	this->numseq = NULL;
	this->nucs = NULL;
	this->basepr = NULL;
	this->danglings_on_branch = NULL;
	this->danglings_on_mb_closure = NULL;
	this->stackings_on_branch = NULL;
	this->stackings_on_mb_closure = NULL;

	char line_buffer[MAX_HEADER_LENGTH];
	fgets(line_buffer, MAX_HEADER_LENGTH, fasta_file);
	if(line_buffer[0] == '>')
	{
		// Copy label.
		this->ctlabel = (char*)malloc(sizeof(char) * MAX_HEADER_LENGTH);
		strcpy(this->ctlabel, &line_buffer[1]);
		if(this->ctlabel[strlen(this->ctlabel) - 1] == '\n')
		{
			this->ctlabel[strlen(this->ctlabel) - 1] = 0;
		}
	}
	this->check_set_label();

	// Read and determine length of sequence.
	char cur_char = 0;
	this->numofbases = 0;

	// Start reading sequence data.
	while(1)
	{
		int ret = fscanf(fasta_file, "%c", &cur_char);
		if(ret == EOF)
		{
			break;
		}

		// Found a new fasta sequence?
		if(cur_char == '>')
		{
			break;
		}

		if(cur_char != '\n' && cur_char != ' ')
		{
			this->numofbases++;
		}
	}

	//printf("Length of sequence is %d\n", this->numofbases);
	this->numseq = (int*)malloc(sizeof(int) * (this->numofbases + 1));
	this->nucs = (char*)malloc(sizeof(char) * (this->numofbases + 2));
	this->basepr = (int*)malloc(sizeof(int) * (this->numofbases + 1));
	this->unpaired_forced = (bool*)malloc(sizeof(bool) * (this->numofbases + 2));

	// Set file position to data position.
	// Cannot use fsetpos and fgetpos because for some reason they are messing up indices
	// when a linux text file is taken to a windows machine.
	fseek(fasta_file, 0, SEEK_SET);

	// Read captoin information.
	fgets(line_buffer, MAX_HEADER_LENGTH, fasta_file);

	int i = 1; // Sequence index, starts from 1.

	// Start reading sequence data.
	while(1)
	{
		// Read and validate input.
		int ret = fscanf(fasta_file, "%c", &cur_char);
		if(ret == EOF)
		{
			break;
		}

		// Check end of sequence marker.
		if(cur_char == '>')
		{
			break;
		}

		// Process this nuc.
		if(cur_char != '\n' && cur_char != ' ')
		{
			this->basepr[i] = 0; // No base pairing information.

			map_nuc_IUPAC_code(cur_char, this->nucs[i], this->numseq[i], this->unpaired_forced[i]);

			//printf("%c %d\n", this->nucs[i], this->numseq[i]);

			i++;
		}
	}

	// This is for ending sequences.
	this->nucs[i] = 0; 

	//printf("Read fasta file: %s (%d nucs)\n", this->nucs, this->numofbases);
	//getc(stdin);

	fclose(fasta_file);
}

bool t_structure::verify_ct(const char* ct_fp)
{
	return(true);
}

/*
Seq file should be like this:
;
[Empty line or comment or id ...]
[Empty line or sequence data]
*/
bool t_structure::verify_seq(const char* seq_fp)
{
	return(true);

	FILE* f_seq = open_f(seq_fp, "r");

	char line_buffer[MAX_HEADER_LENGTH];

	_fgets(line_buffer, MAX_HEADER_LENGTH, f_seq);

	// If the first character of first line is not semicolon, 
	// this is not a valid sequence file.
	if(line_buffer[0] != ';')
	{
		printf("Verification failed for sequence file %s @ %s(%d)\n", seq_fp, __FILE__, __LINE__);
		return(false);
	}

	int current_line_cnt = 2;
	int i_seq = 0;
	char seq_data[MAX_HEADER_LENGTH];

	// Read file and fill lines.
	while(1)
	{
		// Read next line starting with 2nd line.
		if(_fgets(line_buffer, MAX_HEADER_LENGTH, f_seq))
		{
			//printf("Current line_buffer: %s\n", line_buffer);

			// If currently read line is after 2nd line
			// the sequence data is being retrieved.
			if(current_line_cnt > 2)
			{
				for(int i = 0; i < (int)strlen(line_buffer); i++)
				{
					// If this is end of sequence, 
					if(seq_data[i_seq - 1] == '1') // Is sequence data already finished?
					{
						printf("Sequence data is ending before file ends, exiting at %s(%d)\n", __FILE__, __LINE__);
						return(false);
					}

					if(line_buffer[i] != '1' &&
						line_buffer[i] != 'A' &&
						line_buffer[i] != 'C' &&
						line_buffer[i] != 'G' &&
						line_buffer[i] != 'U' &&
						line_buffer[i] != 'T' &&
						line_buffer[i] != 'a' &&
						line_buffer[i] != 'c' &&
						line_buffer[i] != 'g' &&
						line_buffer[i] != 'u' &&
						line_buffer[i] != 't')
					{
						printf("Unknown nucleotide in sequence: %c, exiting at %s(%d)\n", line_buffer[i], __FILE__, __LINE__);
						return(false);
					}
					seq_data[i_seq++] = line_buffer[i];
				}
			}

			current_line_cnt++;
		}
		else
		{
			break;
		}
	}

	// If 2nd line is not read OR no sequence data is read,
	// return false.
	/*
	if(current_line_cnt < 3 ||		// Check if at least 3 lines are read.
		i_seq == 0 ||				// Check if sequence data is read.
		seq_data[i_seq - 1] != '1') // Check correct ending of seq_data
	{
		printf("Verification failed for sequence file %s @ %s(%d)\n", seq_fp, __FILE__, __LINE__);
		return(false);
	}
	*/

	if(current_line_cnt < 3)	// Check if at least 3 lines are read.
	{
		printf("Verification failed for sequence file %s @ %s(%d)\n", seq_fp, __FILE__, __LINE__);
		return(false);
	}

	if(i_seq == 0)				// Check if sequence data is read.
	{
		printf("Verification failed for sequence file %s @ %s(%d): No sequence data\n", seq_fp, __FILE__, __LINE__);
		return(false);
	}

	if(seq_data[i_seq - 1] != '1') // Check correct ending of seq_data
	{
		printf("Verification failed for sequence file %s @ %s(%d): %c\n", seq_fp, __FILE__, __LINE__, seq_data[i_seq - 1]);
		return(false);
	}

	fclose(f_seq);

	return(true);
}

char* _fgets(char* buffer, int n_max, FILE* f)
{
	char* ret = fgets(buffer, MAX_HEADER_LENGTH, f);

	//printf("buffer in _fgets: %s\n", buffer); 

	// Get rid of newline character.
	/* 
	012345
	AAAAA\n
	strlen is 6.
	*/
	if(buffer[strlen(buffer) - 1] == '\n')
	{
		buffer[strlen(buffer) - 1] = 0;
	}

	//printf("buffer in _fgets after fix: %s\n", buffer); 

	return(ret);
}




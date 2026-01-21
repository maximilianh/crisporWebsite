#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "structure.h"
#include <string.h>
#include <ctype.h>

bool _fgets(char* buffer, int n_max, FILE* f);

// Default constructor.
t_structure::t_structure()
{
    this->numofbases = 0;
	this->numseq = NULL;
	this->nucs = NULL;
	this->basepr = NULL;
	this->ctlabel = NULL;
	this->fp = NULL;
}

t_structure::t_structure(const char* _fp)
{
	char temp_fp[500];

	this->fp = (char*)malloc( (strlen(_fp) + 3) * sizeof(char) );
	strcpy(this->fp, _fp);
	strcpy(temp_fp, _fp);
	

	char* last_token;

	//printf("%s\n", _fp);

	// sense if fp is a seq file or a ct file.
	char* cur_token = strtok(temp_fp, "."); // Note that strtok does not work with pointers from stack.

	while(cur_token != NULL)
	{
		last_token = cur_token;
		cur_token = strtok(NULL, ".");
	}

	//printf("last token is: %s\n", last_token);

	char seq_str[] = "seq";
	char ct_str[] = "ct";
	if(strlen(last_token) == 3)
	{
		bool is_seq = true;
		for(int cnt = 0; cnt < 3; cnt++)
		{
			if(toupper(last_token[cnt]) != toupper(seq_str[cnt]))
			{
				is_seq = false;
			}
		}

		if(is_seq)
		{
			this->openseq(fp);
		}
	}
	else if(strlen(last_token) == 2)
	{
		bool is_ct = true;
		for(int cnt = 0; cnt < 2; cnt++)
		{
			if(toupper(last_token[cnt]) != toupper(ct_str[cnt]))
			{
				is_ct = false;
			}
		}

		if(is_ct)
		{
			this->openct(fp);
		}
	}
	else
	{
		printf("Could not determine file type of input @ %s(%d).\n", __FILE__, __LINE__);
		exit(0);
	}
}

void t_structure::openct(const char* ct_fp)
{
	FILE* ct_file = fopen(ct_fp, "r");
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
	fgets(this->ctlabel, MAX_HEADER_LENGTH, ct_file);

	//printf("ct label: %s\n", this->ctlabel);

	this->numseq = (int*)malloc(sizeof(int) * (this->numofbases + 1));
	this->nucs = (char*)malloc(sizeof(char) * (this->numofbases + 1));
	this->basepr = (int*)malloc(sizeof(int) * (this->numofbases + 1));

	// Read sequence data.
	for(int i = 1; i <= this->numofbases; i++)
	{
		int index;
		int some_val1;
		int some_val2;
		int some_val3;
		//                1  G 0  2  120 1
		fscanf(ct_file, "%d %c %d %d %d %d", &index, &this->nucs[i], &some_val1, &some_val2, &this->basepr[i], &some_val3);

		//printf("%c", this->nucs[i]);

		// Convert current base character into number value, from Dave's structure code.
		if (toupper(this->nucs[i]) == 'A') 
			this->numseq[i]=1;
		else if (toupper(this->nucs[i]) == 'C') 
			this->numseq[i]=2;
		else if (toupper(this->nucs[i]) == 'G') 
			this->numseq[i]=3;
		else if (toupper(this->nucs[i]) == 'U' || toupper(this->nucs[i]) == 'T') 
			this->numseq[i]=4;
		else if (toupper(this->nucs[i]) == 'I') 
			this->numseq[i]=5;
		else 
			this->numseq[i]=0;

		//printf("%d\n", this->basepr[i]);
	}

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

	FILE* seq_file = fopen(seq_fp, "r");
	if(seq_file == NULL)
	{
		printf("seq file %s does not exist @ %s(%d).\n", seq_fp, __FILE__, __LINE__);
		exit(1);
	}

	char line_buffer[MAX_HEADER_LENGTH];
	fgets(line_buffer, MAX_HEADER_LENGTH, seq_file);
	while(line_buffer[0] == ';')
	{
		fgets(line_buffer, MAX_HEADER_LENGTH, seq_file);
	}

	// Read label, contains new line character at the end of label.
	this->ctlabel = (char*)malloc(sizeof(char) * MAX_HEADER_LENGTH);
	strcpy(this->ctlabel, line_buffer);
	
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
			this->nucs[i] = cur_char;

			// Convert current base character into number value, from Dave's structure code.
			if (toupper(this->nucs[i]) == 'A') 
				this->numseq[i]=1;
			else if (toupper(this->nucs[i]) == 'C') 
				this->numseq[i]=2;
			else if (toupper(this->nucs[i]) == 'G') 
				this->numseq[i]=3;
			else if (toupper(this->nucs[i]) == 'U' || toupper(this->nucs[i]) == 'T') 
				this->numseq[i]=4;
			else if (toupper(this->nucs[i]) == 'I') 
				this->numseq[i]=5;
			else 
				this->numseq[i]=0;

			this->basepr[i] = 0; // No base pairing information.

			//printf("%c %d\n", this->nucs[i], this->numseq[i]);

			i++;
		}
	}

	// This is for ending sequences.
	this->nucs[i] = 0; 

	fclose(seq_file);
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

	FILE* f_seq = fopen(seq_fp, "r");

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
				for(int i = 0; i < strlen(line_buffer); i++)
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

bool _fgets(char* buffer, int n_max, FILE* f)
{
	bool ret = fgets(buffer, MAX_HEADER_LENGTH, f);

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


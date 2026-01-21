#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

bool check_file(char* fp)
{
	FILE* f_temp = fopen(fp, "r");
	if(f_temp == NULL)
	{		
		return(false);
	}

	fclose(f_temp);
	return(true);
}

// Check for valid CR LF's depending on OS,
// Should be run for all ASCII input files.
#define CR (0x0D)
#define LF (0x0A)
void validate_file(char* fp)
{
#ifdef _WIN32
	char cur_char;
	// Open file in binary.
	FILE* f_ip_bin = open_f(fp, "rb");
	while(fread(&cur_char, 1, 1, f_ip_bin) == 1)
	{
		if(cur_char == CR)
		{
			if(fread(&cur_char, 1, 1, f_ip_bin) == 1)
			{
				if(cur_char != LF)
				{
					// Just a warning here.
					printf("%s is not compatible with dos ascii files. CR+LF problem at %s(%d).\n", fp, __FILE__, __LINE__);
					//exit(0);
				}
			}
			else
			{
				// Just a warning here.
				printf("%s is not compatible with dos ascii files. CR+LF problem at %s(%d).\n", fp, __FILE__, __LINE__);
				//exit(0);
			}
		}
		else if(cur_char == LF) // If there is an immediate LF before seeing a CR, this is a linux file.
		{
			// Just a warning here.
			printf("%s is not compatible with dos ascii files. CR+LF problem at %s(%d).\n", fp, __FILE__, __LINE__);
			//exit(0);
		}

	}
	fclose(f_ip_bin);
#endif

#ifdef __unix__
	char cur_char;
	// Open file in binary.
	FILE* f_ip_bin = open_f(fp, "rb");
	while(fread(&cur_char, 1, 1, f_ip_bin) == 1)
	{
		// Linux files do not contain CR's.
		// They only contain LF's.
		if(cur_char == CR)
		{
			// Just a warning here.
			printf("%s is not compatible with Linux ascii files. CR+LF problem at %s(%d).\n", fp, __FILE__, __LINE__);
			//exit(0);
		}
	}
	fclose(f_ip_bin);
#endif

#ifdef __APPLE__
        char cur_char;
        // Open file in binary.
	FILE* f_ip_bin = open_f(fp, "rb");
        while(fread(&cur_char, 1, 1, f_ip_bin) == 1)
	  {
	    // Linux files do not contain CR's.
	    // They only contain LF's.
	    if(cur_char == CR)
	      {
		// Just a warning here.
		printf("%s is not compatible with Linux ascii files. CR+LF problem at %s(%d).\n", fp, __FILE__, __LINE__);                                                                                                                                   
	      }
	  }
        fclose(f_ip_bin);
#endif
}

FILE* open_f(const char* fp, const char* mode)
{
	if(fp == NULL || mode == NULL)
	{
		printf("Invalid arguments to open_f: %s.\n", fp);
		exit(0);
	}

	FILE* f = fopen(fp, mode);

	if(f == NULL)
	{
		if(mode[0] == 'r')
		{
			printf("Could not open %s for reading.\n", fp);
			exit(0);
		}
		else if(mode[0] == 'w')
		{
			printf("Could not open %s for writing.\n", fp);
			exit(0);
		}
		else
		{
			printf("Could not open %s for requested operation.\n", fp);
			exit(0);
		}
	}

	return(f);
}

const char* resolve_data_dir()
{
	// try to resolve the DATAPATH_ENV_VAR.
	const char* data_dir_from_env = getenv(DATAPATH_ENV_VAR);

	if(data_dir_from_env != NULL)
	{
		return(data_dir_from_env);
	}
	else
	{
		return(DATAPATH_DEFAULT);
	}

	printf("Could not resolve thermodynamics data directory.\n");
	exit(0);
}

char* x_fgets(char* buff, int size, FILE* file)
{
	if(fgets(buff, size, file) == NULL)
	{
		return(NULL);
	}

	if(buff[strlen(buff) - 1] == '\n')
	{
		int i_new_line_char = (int)strlen(buff) - 1;
		buff[i_new_line_char] = 0;
	}

	return(buff);
}


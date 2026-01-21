#include <string.h>
#include <limits.h>
#include "parts_compilation_directives.h"
#include "parts_paths.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ppf_cli.h"
#include "../../../src/phmm/structure/structure_object.h"
#include "../../../src/phmm/utils/file/utils.h"
#include "ppf_math.h"
//#include "ppf_math.h"

bool _DUMP_PPF_CLI_MESSAGES_ = false;

/*
pairwise_part_func [sequence 1 path] [sequence 2 path] -map/-pp -ap [aln prior file path] -fp1 [folding prior 1] -fp2 [folding prior 2] -ct1 [correct ct 1 path] -ct2 [correct ct 2 path]
*/

int t_ppf_cli::mode; 

t_ppf_cli::t_ppf_cli(char* seq1_fp, 
					 char* seq2_fp, 
					 int mode)
{
	this->init_paths_vars();

	// First two arguments should be paths to seq1 and seq2, which must exist.
	if(!check_file(seq1_fp))
	{
		//printf("sequence 1 could not be found @ %s, exiting @ %s(%d)\n", seq1_fp, __FILE__, __LINE__);
		//exit(0);
		this->last_error_code = ERR_SEQ1_NOT_EXISTENT;
		return;
	}
	else
	{
		// Do a validicity check on seq2 file.
		validate_file(seq1_fp);

		// Now safe to use this file.
		this->seq1_path = (char*)malloc(strlen(seq1_fp) + 3);
		strcpy(this->seq1_path, seq1_fp);

		if(_DUMP_PPF_CLI_MESSAGES_)
			printf("Seq1 path: %s\n", seq1_path);
	}

	// Test 2nd sequence path.
	if(!check_file(seq2_fp))
	{
		//printf("sequence 2 could not be found @ %s, exiting @ %s(%d)\n", seq2_fp, __FILE__, __LINE__);
		//exit(0);
		this->last_error_code = ERR_SEQ2_NOT_EXISTENT;
		return;
	}
	else
	{
		// Do a validicity check on seq2 file.
		validate_file(seq2_fp);

		// Now safe to use this file.
		this->seq2_path = (char*)malloc(strlen(seq2_fp) + 3);
		strcpy(this->seq2_path, seq2_fp);

		if(_DUMP_PPF_CLI_MESSAGES_)
			printf("Seq2 path: %s\n", seq2_path);
	}

	// Check if this is an map or pp calculation.
	if(mode == PARTS_RUN_MODE_MAP)
	{
		if(_DUMP_PPF_CLI_MESSAGES_)
			printf("Doing map calculation..\n");;

		//this->map = true;
		t_ppf_cli::mode = PARTS_RUN_MODE_MAP;
	}
	else if(mode == PARTS_RUN_MODE_PP)
	{
		if(_DUMP_PPF_CLI_MESSAGES_)
			printf("Doing pp calculation..\n");

		//this->map = false;
		t_ppf_cli::mode = PARTS_RUN_MODE_PP;
	}
	else if(mode == PARTS_RUN_MODE_STOCH_SAMPLE)
	{
		if(_DUMP_PPF_CLI_MESSAGES_)
			printf("Doing stochastic sampling..\n");

		t_ppf_cli::mode = PARTS_RUN_MODE_STOCH_SAMPLE; // Sampling based on stochastic tracebacking over pf arrays.
	}
	else if(mode == PARTS_RUN_MODE_INIT)
	{
		t_ppf_cli::mode = PARTS_RUN_MODE_INIT; // Sampling based on stochastic tracebacking over pf arrays.
	}
	else
	{
		this->last_error_code = ERR_MODE_NOT_VALID;
		return;
	}

	// Reset the prefixes.
	if(this->seq1_op_file_prefix == NULL)
	{
		this->seq1_op_file_prefix = new char[strlen("seq1") + 2];
		strcpy(this->seq1_op_file_prefix, "seq1");
	}

	if(this->seq2_op_file_prefix == NULL)
	{
		this->seq2_op_file_prefix = new char[strlen("seq2") + 2];
		strcpy(this->seq2_op_file_prefix, "seq2");
	}
}

void t_ppf_cli::set_output_prefixes(const char* seq1_prefix, const char* seq2_prefix){
  delete seq1_op_file_prefix;
  delete seq2_op_file_prefix;
  
  seq1_op_file_prefix = new char[strlen(seq1_prefix) + 2];
  strcpy(this->seq1_op_file_prefix, seq1_prefix);
  
  seq2_op_file_prefix = new char[strlen(seq2_prefix) + 2];
  strcpy(this->seq2_op_file_prefix, seq2_prefix);
  
  return;
}

t_ppf_cli::t_ppf_cli(char* conf_fp)
{
	this->init_paths_vars();

	this->load_conf_file(conf_fp);
}

t_ppf_cli::~t_ppf_cli()
{
	// Input sequence paths.
	if(seq1_path != NULL)
	{
		free(this->seq1_path);
	}

	if(seq2_path != NULL)
	{
		free(this->seq2_path);
	}

	if(seq1_op_file_prefix != NULL)
	{
		delete[] seq1_op_file_prefix;
	}

	if(seq2_op_file_prefix != NULL)
	{
		delete[] seq2_op_file_prefix;
	}

	if(seq1_map_ct_op != NULL)
	{
		free(this->seq1_map_ct_op);
	}
	if(seq2_map_ct_op != NULL)
	{
		free(this->seq2_map_ct_op);
	}
	if(map_aln_op != NULL)
	{
		free(this->map_aln_op);
	}

	if(seq1_pp_op != NULL)
	{
		free(this->seq1_pp_op);
	}

	if(seq2_pp_op != NULL)
	{
		free(this->seq2_pp_op);
	}

	if(seq1_sample_ct_op != NULL)
	{
		free(this->seq1_sample_ct_op);
	}

	if(seq2_sample_ct_op != NULL)
	{
		free(this->seq2_sample_ct_op);
	}

	if(sample_aln_op != NULL)
	{
		free(this->sample_aln_op);
	}

	if(seq1_SHAPE_path != NULL)
	{
		free(this->seq1_SHAPE_path);
	}
		
	if(seq2_SHAPE_path != NULL)
	{
		free(this->seq2_SHAPE_path);
	}

	// alignment and sequence ids for determining 
	// an identifier for the alignment.
	if(seq1_id != NULL)
	{
		free(this->seq1_id);
	}

	if(seq2_id != NULL)
	{
		free(this->seq2_id);
	}

	if(aln_id != NULL)
	{
		free(this->aln_id);
	}

	// aln prior file path, if existing.
	if(aln_prior_path != NULL)
	{
		free(this->aln_prior_path);
	}

	if(ins1_prior_path != NULL)
	{
		free(this->ins1_prior_path);
	}

	if(ins2_prior_path != NULL)
	{
		free(this->ins2_prior_path);
	}

	// fold prior 1 file path, if existing.
	if(fold_prior1_path != NULL)
	{
		free(this->fold_prior1_path);
	}

	if(str1_coinc_info_path != NULL)
	{
		free(this->str1_coinc_info_path);
	}

	// fold prior 2 file path, if existing.
	if(fold_prior2_path != NULL)
	{
		free(this->fold_prior2_path);
	}

	// unpairing paths.
	if(i_unpairing_prior1_path != NULL)
	{
		free(this->i_unpairing_prior1_path);
	}

	if(i_unpairing_prior2_path != NULL)
	{
		free(this->i_unpairing_prior2_path);
	}

	if(j_unpairing_prior1_path != NULL)
	{
		free(this->j_unpairing_prior1_path);
	}

	if(j_unpairing_prior2_path != NULL)
	{
		free(this->j_unpairing_prior2_path);
	}

	// correct ct1 path.
	if(correct_ct1_path != NULL)
	{
		free(this->correct_ct2_path);
	}

	// correct ct2 path.
	if(correct_ct2_path != NULL)
	{
		free(this->correct_ct2_path);
	}

	// The path of array dump file.
	if(array_dump_fp != NULL)
	{
		free(this->array_dump_fp);
	}

	if(loop_limits_fn = NULL)	
	{
		free(this->loop_limits_fn);
	}
}

t_ppf_cli::t_ppf_cli(int argc, char* argv[])
{
	if(argc < 4)
	{
		this->print_usage_msg();
		exit(0);
	}

	this->init_paths_vars();

	// First two arguments should be paths to seq1 and seq2, which must exist.
	if(!check_file(argv[1]))
	{
		//printf("sequence 1 could not be found @ %s, exiting @ %s(%d)\n", argv[1], __FILE__, __LINE__);
		//exit(0);
		this->last_error_code = ERR_SEQ1_NOT_EXISTENT;
		return;
	}
	else
	{
		// Do a validicity check on seq2 file.
		validate_file(argv[1]);

		// Now safe to use this file.
		this->seq1_path = (char*)malloc(strlen(argv[1]) + 3);
		strcpy(this->seq1_path, argv[1]);

		if(_DUMP_PPF_CLI_MESSAGES_)
			printf("Seq1 path: %s\n", seq1_path);
	}

	// Test 2nd sequence path.
	if(!check_file(argv[2]))
	{
		//printf("sequence 2 could not be found @ %s, exiting @ %s(%d)\n", argv[2], __FILE__, __LINE__);
		//exit(0);
		this->last_error_code = ERR_SEQ2_NOT_EXISTENT;
		return;
	}
	else
	{
		// Do a validicity check on seq2 file.
		validate_file(argv[2]);

		// Now safe to use this file.
		this->seq2_path = (char*)malloc(strlen(argv[2]) + 3);
		strcpy(this->seq2_path, argv[2]);

		if(_DUMP_PPF_CLI_MESSAGES_)
			printf("Seq2 path: %s\n", seq2_path);
	}

	// Check if this is an map or pp calculation.
	if(strcmp(argv[3], "-map") == 0)
	{
		if(_DUMP_PPF_CLI_MESSAGES_)
			printf("Doing map calculation..\n");;

		//this->map = true;
		t_ppf_cli::mode = PARTS_RUN_MODE_MAP;
	}
	else if(strcmp(argv[3], "-pp") == 0)
	{
		if(_DUMP_PPF_CLI_MESSAGES_)
			printf("Doing pp calculation..\n");

		//this->map = false;
		t_ppf_cli::mode = PARTS_RUN_MODE_PP;
	}
	else if(strcmp(argv[3], "-stochsample") == 0)
	{
		if(_DUMP_PPF_CLI_MESSAGES_)
			printf("Doing stochastic sampling..\n");

		t_ppf_cli::mode = PARTS_RUN_MODE_STOCH_SAMPLE; // Sampling based on stochastic tracebacking over pf arrays.
	}
	else
	{
		//printf("Could not determine what calculation to make from command line, unknown flag is %s\n", argv[3]); 
		//exit(0);
		this->last_error_code = ERR_MODE_NOT_VALID;
		return;
	}

	// Parsed mandatory flags, parse remaining information.
	int p_cnt = 4;
	while(p_cnt < argc)
	{
		//printf("%d\n", p_cnt);
		if(argv[p_cnt][0] == '-') // Is that a flag?
		{
			// Align prior?
			if(strcmp(argv[p_cnt], "-ap") == 0)
			{
				p_cnt++; // Alignment prior file path is in argv[p_cnt].
				
				// After -ap, have to define 3 alignment priors in aln ins1 ins2 order.

if(_DUMP_PPF_CLI_MESSAGES_)
				printf("Found alignment prior file: %s\n", argv[p_cnt]);

				if(!check_file(argv[p_cnt]) || !check_file(argv[p_cnt + 1]) || !check_file(argv[p_cnt + 2]))
				{
					printf("One or many of alignment prior files could not be found @ %s, %s, %s, exiting @ %s(%d)\n", argv[p_cnt], argv[p_cnt+1], argv[p_cnt+2], __FILE__, __LINE__);
					exit(0);
				}
				else
				{
					this->aln_prior_path = (char*)malloc(strlen(argv[p_cnt]) + 3);
					this->ins1_prior_path = (char*)malloc(strlen(argv[p_cnt + 1]) + 3);
					this->ins2_prior_path = (char*)malloc(strlen(argv[p_cnt + 2]) + 3);
					this->loop_limits_fn = (char*)malloc(sizeof(char) * (strlen(argv[p_cnt+3]) + 3));

					strcpy(this->aln_prior_path, argv[p_cnt]);
					strcpy(this->ins1_prior_path, argv[p_cnt+1]);
					strcpy(this->ins2_prior_path, argv[p_cnt+2]);
					strcpy(this->loop_limits_fn, argv[p_cnt+3]);
				}
	
				// Jump over next flag.
				p_cnt++; // jumps to ins1
				p_cnt++; // jumps to ins2
				p_cnt++; // jumps to ll fp.
				p_cnt++; // jumps to next one.
			}
			else if(strcmp(argv[p_cnt], "-fp1") == 0) // fold prior?
			{
				p_cnt++; // Alignment prior file path is in argv[p_cnt].

if(_DUMP_PPF_CLI_MESSAGES_)
				printf("Found fold prior1 file: %s\n", argv[p_cnt]);

				if(!check_file(argv[p_cnt]))
				{
					printf("Sequence 1 fold prior file could not be found @ %s, %s, %s, exiting @ %s(%d)\n", argv[p_cnt], argv[p_cnt+1], argv[p_cnt+2], __FILE__, __LINE__);
					exit(0);
				}
				else
				{
					this->fold_prior1_path = (char*)malloc(strlen(argv[p_cnt]) + 3);
					strcpy(this->fold_prior1_path, argv[p_cnt]);
				}
				// Jump over next flag.
				p_cnt++;
			}
			else if(strcmp(argv[p_cnt], "-fp2") == 0) // fold prior?
			{
				p_cnt++; // Alignment prior file path is in argv[p_cnt].

if(_DUMP_PPF_CLI_MESSAGES_)
				printf("Found fold prior2 file: %s\n", argv[p_cnt]);

				if(!check_file(argv[p_cnt]))
				{
					printf("Sequence 2 fold prior file could not be found @ %s, exiting @ %s(%d)\n", argv[p_cnt], __FILE__, __LINE__);
					exit(0);
				}
				else
				{
					this->fold_prior2_path = (char*)malloc(strlen(argv[p_cnt]) + 3);
					strcpy(this->fold_prior2_path, argv[p_cnt]);
				}
	
				// Jump over next flag.
				p_cnt++;				
			}
			else if(strcmp(argv[p_cnt], "-seq1_SHAPE_fp") == 0)
			{
				p_cnt++; // Alignment prior file path is in argv[p_cnt].

if(_DUMP_PPF_CLI_MESSAGES_)
				printf("seq1 SHAPE fp: %s\n", argv[p_cnt]);

				if(!check_file(argv[p_cnt]))
				{
					printf("Sequence 1 SHAPE data could not be found @ %s, exiting @ %s(%d)\n", argv[p_cnt], __FILE__, __LINE__);
					exit(0);
				}
				else
				{
					this->seq1_SHAPE_path = (char*)malloc(strlen(argv[p_cnt]) + 3);
					strcpy(this->seq1_SHAPE_path, argv[p_cnt]);
				}
	
				// Jump over next flag.
				p_cnt++;				
			}
			else if(strcmp(argv[p_cnt], "-seq2_SHAPE_fp") == 0)
			{
				p_cnt++; // Alignment prior file path is in argv[p_cnt].

if(_DUMP_PPF_CLI_MESSAGES_)
				printf("seq2 SHAPE fp: %s\n", argv[p_cnt]);

				if(!check_file(argv[p_cnt]))
				{
					printf("Sequence 2 SHAPE data could not be found @ %s, exiting @ %s(%d)\n", argv[p_cnt], __FILE__, __LINE__);
					exit(0);
				}
				else
				{
					this->seq2_SHAPE_path = (char*)malloc(strlen(argv[p_cnt]) + 3);
					strcpy(this->seq2_SHAPE_path, argv[p_cnt]);
				}
	
				// Jump over next flag.
				p_cnt++;				
			}
			else if(strcmp(argv[p_cnt], "-ct1") == 0) // Correct ct1 path?
			{
				p_cnt++; // Alignment prior file path is in argv[p_cnt].

if(_DUMP_PPF_CLI_MESSAGES_)
				printf("Found correct ct1 file: %s\n", argv[p_cnt]);

				if(!check_file(argv[p_cnt]))
				{
					printf("Sequence 1 correct ct file could not be found @ %s, exiting @ %s(%d)\n", argv[p_cnt], __FILE__, __LINE__);
					exit(0);
				}
				else
				{
					this->correct_ct1_path = (char*)malloc(strlen(argv[p_cnt]) + 3);
					strcpy(this->correct_ct1_path, argv[p_cnt]);

					// Load correct 1st structure.
					t_structure* ct1 = new t_structure(this->correct_ct1_path);


if(_DUMP_PPF_CLI_MESSAGES_)
{
					// Write ct1 dot plot.
					FILE* ct1_dot_plot_file = open_f("ct1_dp.txt", "w");
					for(int i1 = 1; i1 <= ct1->numofbases; i1++)
					{
						for(int i2 = 1; i2 <= ct1->numofbases; i2++)
						{
							if(ct1->basepr[i1] == i2)
							{
								fprintf(ct1_dot_plot_file, "1 ");
							}
							else
							{
								fprintf(ct1_dot_plot_file, "0 ");
							}
						}
						fprintf(ct1_dot_plot_file, "\n");
					}
					fclose(ct1_dot_plot_file);
}
				}	
	
				// Jump over next flag.
				p_cnt++;			
			}
			else if(strcmp(argv[p_cnt], "-ct2") == 0) // Correct ct2 path?
			{
				p_cnt++; // Alignment prior file path is in argv[p_cnt].

if(_DUMP_PPF_CLI_MESSAGES_)
				printf("Found correct ct2 file: %s\n", argv[p_cnt]);

				if(!check_file(argv[p_cnt]))
				{
					printf("Sequence 2 correct ct file could not be found @ %s, exiting @ %s(%d)\n", argv[p_cnt], __FILE__, __LINE__);
					exit(0);
				}
				else
				{
					this->correct_ct2_path = (char*)malloc(strlen(argv[p_cnt]) + 3);
					strcpy(this->correct_ct2_path, argv[p_cnt]);

					// Load correct 2nd structure.
					t_structure* ct2 = new t_structure(this->correct_ct2_path);

if(_DUMP_PPF_CLI_MESSAGES_)
{
					// Write ct2 dot plot.
					FILE* ct2_dot_plot_file = open_f("ct2_dp.txt", "w");
					for(int i1 = 1; i1 <= ct2->numofbases; i1++)
					{
						for(int i2 = 1; i2 <= ct2->numofbases; i2++)
						{
							if(ct2->basepr[i1] == i2)
							{
								fprintf(ct2_dot_plot_file, "1 ");
							}
							else
							{
								fprintf(ct2_dot_plot_file, "0 ");
							}
						}
						fprintf(ct2_dot_plot_file, "\n");
					}
					fclose(ct2_dot_plot_file);
}
				}
	
				// Jump over next flag.
				p_cnt++;				
			}
			else if(strcmp(argv[p_cnt], "-aw") == 0) // Set log alignment weight per nucleotide?
			{
				p_cnt++; // double log alignment weight per nucleotide is @ argv[p_cnt].

if(_DUMP_PPF_CLI_MESSAGES_)
				printf("Overriding log alignment weight per alignment position with %f\n", atof(argv[p_cnt]));

				this->log_aln_weight_per_aln_pos = atof(argv[p_cnt]);
					
				// Jump over next flag.
				p_cnt++;				
			}
            else if(strcmp(argv[p_cnt], "-rfpn") == 0) // Set log alignment weight per nucleotide?
            {
                p_cnt++; // double log alignment weight per nucleotide is @ argv[p_cnt].

if(_DUMP_PPF_CLI_MESSAGES_)
                printf("Overriding rescaling factor per nucleotide with %f\n", atof(argv[p_cnt]));

                this->cmd_rescaling_increment_factor_per_nucleotide = atof(argv[p_cnt]);

                // Jump over next flag.
                p_cnt++;
            }
            else if(strcmp(argv[p_cnt], "-seed") == 0) // Set log alignment weight per nucleotide?
            {
                    p_cnt++; // double log alignment weight per nucleotide is @ argv[p_cnt].

if(_DUMP_PPF_CLI_MESSAGES_)
                    printf("Setting stochastic sampling seed to %d.\n", atoi(argv[p_cnt]));

                    this->stoch_sampling_seed = atoi(argv[p_cnt]);

                    // Jump over next flag.
                    p_cnt++;
	    }
            else if(strcmp(argv[p_cnt], "-nsamp") == 0) // Set log alignment weight per nucleotide?
            {
                    p_cnt++; // double log alignment weight per nucleotide is @ argv[p_cnt].

if(_DUMP_PPF_CLI_MESSAGES_)
                    printf("Sample set size is %d.\n", atoi(argv[p_cnt]));

                    this->n_samples = atoi(argv[p_cnt]);

                    // Jump over next flag.
                    p_cnt++;
            }
			else if(strcmp(argv[p_cnt], "-seq1_op_prefix") == 0) // Set a file name for sampled cts: All sampled ct's are dumped into one file.
			{
				p_cnt++;

				this->seq1_op_file_prefix = (char*)malloc(strlen(argv[p_cnt]) + 3);
				strcpy(this->seq1_op_file_prefix, argv[p_cnt]);

                // Jump over next flag.
                p_cnt++;
			}
			else if(strcmp(argv[p_cnt], "-seq2_op_prefix") == 0) // Set a file name for sampled cts: All sampled ct's are dumped into one file.
			{
				p_cnt++;

				this->seq2_op_file_prefix = (char*)malloc(strlen(argv[p_cnt]) + 3);
				strcpy(this->seq2_op_file_prefix, argv[p_cnt]);

                // Jump over next flag.
                p_cnt++;
			}












			else if(strcmp(argv[p_cnt], "-seq1_map_ct_op") == 0) // Set a file name for sampled cts: All sampled ct's are dumped into one file.
			{
				p_cnt++;

				this->seq1_map_ct_op = (char*)malloc(strlen(argv[p_cnt]) + 3);
				strcpy(this->seq1_map_ct_op, argv[p_cnt]);

                // Jump over next flag.
                p_cnt++;
			}
			else if(strcmp(argv[p_cnt], "-seq2_map_ct_op") == 0) // Set a file name for sampled cts: All sampled ct's are dumped into one file.
			{
				p_cnt++;

				this->seq2_map_ct_op = (char*)malloc(strlen(argv[p_cnt]) + 3);
				strcpy(this->seq2_map_ct_op, argv[p_cnt]);

                // Jump over next flag.
                p_cnt++;
			}
			else if(strcmp(argv[p_cnt], "-map_aln_op") == 0) // Set a file name for sampled cts: All sampled ct's are dumped into one file.
			{
				p_cnt++;

				this->map_aln_op = (char*)malloc(strlen(argv[p_cnt]) + 3);
				strcpy(this->map_aln_op, argv[p_cnt]);

                // Jump over next flag.
                p_cnt++;
			}
			else if(strcmp(argv[p_cnt], "-seq1_pp_op") == 0) // Set a file name for sampled cts: All sampled ct's are dumped into one file.
			{
				p_cnt++;

				this->seq1_pp_op = (char*)malloc(strlen(argv[p_cnt]) + 3);
				strcpy(this->seq1_pp_op, argv[p_cnt]);

                // Jump over next flag.
                p_cnt++;
			}
			else if(strcmp(argv[p_cnt], "-seq2_pp_op") == 0) // Set a file name for sampled cts: All sampled ct's are dumped into one file.
			{
				p_cnt++;

				this->seq2_pp_op = (char*)malloc(strlen(argv[p_cnt]) + 3);
				strcpy(this->seq2_pp_op, argv[p_cnt]);

                // Jump over next flag.
                p_cnt++;
			}
			else if(strcmp(argv[p_cnt], "-seq1_sample_ct_op") == 0) // Set a file name for sampled cts: All sampled ct's are dumped into one file.
			{
				p_cnt++;

				this->seq1_sample_ct_op = (char*)malloc(strlen(argv[p_cnt]) + 3);
				strcpy(this->seq1_sample_ct_op, argv[p_cnt]);

                // Jump over next flag.
                p_cnt++;
			}
			else if(strcmp(argv[p_cnt], "-seq2_sample_ct_op") == 0) // Set a file name for sampled cts: All sampled ct's are dumped into one file.
			{
				p_cnt++;

				this->seq2_sample_ct_op = (char*)malloc(strlen(argv[p_cnt]) + 3);
				strcpy(this->seq2_sample_ct_op, argv[p_cnt]);

                // Jump over next flag.
                p_cnt++;
			}
			else if(strcmp(argv[p_cnt], "-sample_aln_op") == 0) // Set a file name for sampled cts: All sampled ct's are dumped into one file.
			{
				p_cnt++;

				this->sample_aln_op = (char*)malloc(strlen(argv[p_cnt]) + 3);
				strcpy(this->sample_aln_op, argv[p_cnt]);

                // Jump over next flag.
                p_cnt++;
			}







			else if(strcmp(argv[p_cnt], "-max_n_separation") == 0)
			{
				p_cnt++;
				this->max_n_separation_between_nucs = atoi(argv[p_cnt]);
				printf("Read maximum separation: %d\n", this->max_n_separation_between_nucs);
				p_cnt++;
			}
			else if(strcmp(argv[p_cnt], "-fold_env_prob_treshold") == 0)
			{
				p_cnt++;
				this->fold_env_prob_treshold = atof(argv[p_cnt]);
				p_cnt++;
			}
			else if(strcmp(argv[p_cnt], "-array_mem_limit_in_megs") == 0)
			{
				p_cnt++;
				this->array_mem_limit_in_megs = atof(argv[p_cnt]);
				p_cnt++;
			}
			else if(strcmp(argv[p_cnt], "-phmm_band_constraint_size") == 0)
			{
				p_cnt++;
				this->phmm_band_constraint_size = atoi(argv[p_cnt]);
				p_cnt++;
			}
			else if(strcmp(argv[p_cnt], "-str_coinc_env_prob_treshold") == 0)
			{
				p_cnt++;
				this->str_coinc_env_prob_treshold = atof(argv[p_cnt]);
				p_cnt++;
			}
			else if(strcmp(argv[p_cnt], "-use_array_files") == 0) // Use the files in array_dumps?
			{
				p_cnt++;
				this->use_array_files = true;
			}
			else if(strcmp(argv[p_cnt], "-save_array_files") == 0) // Use the files in array_dumps?
			{
				p_cnt++;
				this->save_array_files = true;
			}
			else
			{
				p_cnt++;
			}
		}
	}

	if(this->seq1_op_file_prefix == NULL)
	{
		this->seq1_op_file_prefix = (char*)malloc(sizeof(char) * (strlen("seq1") + 2));
		strcpy(this->seq1_op_file_prefix, "seq1");
	}

	if(this->seq2_op_file_prefix == NULL)
	{
		this->seq2_op_file_prefix = (char*)malloc(sizeof(char) * (strlen("seq2") + 2));
		strcpy(this->seq2_op_file_prefix, "seq2");
	}

	// Get seq and alignment id's.
	//this->get_seq_aln_ids();
}

void t_ppf_cli::init_paths_vars()
{
	// Initialize error code.
	this->last_error_code = NO_CLI_ERROR;

	this->seq1_id = NULL;
	this->seq2_id = NULL;
	this->aln_id = NULL;

	// Init paths.
	this->log_aln_weight_per_aln_pos = 1.0; // This is initially 1.0. If it is not read, will be set to 1.0 with same contributions for structure and alignment part.
	this->cmd_rescaling_increment_factor_per_nucleotide = 0.0; // This is initially 0.0, if it is not read, it is set to its default value defined in ppf_scale.h
	this->prior_to_posterior_prop_rate = 0.1; // By default, 10% of priors propagate into posteriors in sampling algorithm. In order to get pure information, set this paremeter to 0 in command line.
	this->n_samples = 1; // By default, number of samples is 1
	this->mode = PARTS_RUN_MODE_INIT; // Initially mode is nothing.
	this->array_mem_limit_in_megs = 0.0f;

	this->stoch_sampling_seed = 0xfffff;

	this->max_n_separation_between_nucs = 0x1fffffff;
	this->phmm_band_constraint_size = 0x1fffffff;

	this->seq1_path = NULL;
	this->seq2_path = NULL;

	this->seq1_map_ct_op = NULL;
	this->seq2_map_ct_op = NULL;
	this->map_aln_op = NULL;
	this->seq1_pp_op = NULL;
	this->seq2_pp_op = NULL;
	this->seq1_sample_ct_op = NULL;
	this->seq2_sample_ct_op = NULL;
	this->sample_aln_op = NULL;

	this->seq1_SHAPE_path = NULL;
	this->seq2_SHAPE_path = NULL;

	// A very low probability slightly greater than 0.0f 
	// to filter out base pairs with 0.0 pairing probability.
	this->fold_env_prob_treshold = 0.999f; 
	this->str_coinc_env_prob_treshold = 0.999f; 

	this->array_dump_fp = NULL;

	aln_prior_path = NULL;
	fold_prior1_path = NULL;
	fold_prior2_path = NULL;

	this->str1_coinc_info_path = NULL;

	i_unpairing_prior1_path = NULL;
	i_unpairing_prior2_path = NULL;
	j_unpairing_prior1_path = NULL;
	j_unpairing_prior2_path = NULL;
	correct_ct1_path = NULL;
	correct_ct2_path = NULL;
	aln_prior_path = NULL;
	ins1_prior_path = NULL;
	ins2_prior_path = NULL;

	loop_limits_fn = NULL;

	this->seq1_op_file_prefix = NULL;
	this->seq2_op_file_prefix = NULL;
}

// This checks for existence.
bool t_ppf_cli::check_file(char* file_path)
{
	FILE* _test = fopen(file_path, "r");
	if(_test == NULL)
	{
		return false;
	}
	else
	{
		fclose(_test); // Close opened sequence file.
		return true;
	}
}

int t_ppf_cli::GetErrorCode()
{
	return(this->last_error_code);
}

char* t_ppf_cli::GetErrorMessage(const int error_code)
{
	return(ppf_cli_error_msgs[error_code]);
}

//Usage depends on the scores used.
void t_ppf_cli::print_usage_msg()
{
	printf("USAGE: PARTS [sequence 1 path] [sequence 2 path] -map/-pp/-stochsample [OPTIONS] \n \
			OR\n\
			PARTS [configuration file path] \n \
		   OPTIONS: \n\
				-seq1_op_prefix [File name prefix for dumping sequence 1 related files]\n\
				-seq2_op_prefix [File name prefix for dumping sequence 2 related files] \n\
				-seq1_map_ct_op [File path for MAP structures of sequence 1] \n\
				-seq2_map_ct_op [File path for MAP structures of sequence 2] \n\
				-map_aln_op [File path for MAP alignment] \n\
				-seq1_pp_op [File path for posterior base pairing probabilities of sequence 1] \n\
				-seq2_pp_op [File path for posterior base pairing probabilities of sequence 2] \n\
				-seq1_sample_ct_op [File path for sampled structures of sequence 1] \n\
				-seq2_sample_ct_op [File path for sampled structures of sequence 2] \n\
				-sample_aln_op [File path for dumping sampled alignments] \n\
				-ap [aln_prior_file_path ins1_prior_file_path ins2_prior_file_path loop_limits_file_path] \n\
				-fp1 [folding prior probabilities for sequence 1] \n\
				-fp2 [folding prior probabilities for sequence 2] \n\
				-nsamp [number of stochastic traceback samples] \n\
				-seed [Pseudo-random number generator seed] \n\
				-ct1 [known ct file path for sequence 1] \n\
				-ct2 [known ct file path for sequence 2] \n\
				-aw [log alignment weight in pseudo-free energy computation] \n\
				-rfpn [rescaling factor per nucleotide] \n\
				-max_n_separation [maximum number separation between structurally co-incident nucleotides] \n\
				-fold_env_prob_treshold [minimum probability threshold for allocation of single pairing arrays] {ENTOR}\n\
				-array_mem_limit_in_megs [soft limit on memory to be allocated to arrays in Megabytes.\n\
				-phmm_band_constraint_size [Size of band constraint on all of HMM computations]\n\
				-str_coinc_env_prob_treshold [minimum probability threshold for allocation of co-incident positions]\n\n\
			Example of a configuration file: \n\
				seq1 RD0260.seq \n\
				seq2 RE2140.seq \n\
				seq1_op_prefix RD0260 \n\
				seq2_op_prefix RE2140 \n\
				mode stochsample \n\
				seed 1234\n\
				n_samp 1500\n\
			Example of a configuration file with explicit file names:\n\
				seq1 RD0260.seq \n\
				seq2 RE2140.seq \n\
				seq1_map_ct_op /home/my_map_ct1.ct \n\
				seq1_map_ct_op /home/my_map_ct2.ct \n\
				map_aln_op /home/map_aln.ct \n\
				mode map \n\
				n_samp 1500\n");
}

void t_ppf_cli::load_conf_file(char* conf_file)
{
	FILE* f_conf = fopen(conf_file, "r");

	if(f_conf == NULL)
	{
		this->last_error_code = ERR_CONF_FILE_NOT_EXISTENT;
		return;
	}

	//printf("loading configuration file %s\n", conf_file);
	//getc(stdin);

	// Read whole file.
	char cur_option[__MAX_PATH];
	char cur_val[__MAX_PATH];
	while(fscanf(f_conf, "%s %s", cur_option, cur_val) == 2)
	{
		//printf("Read %s %s\n", cur_option, cur_val);
		if(strcmp(cur_option, "seq1") == 0)
		{
			this->seq1_path = (char*)malloc(sizeof(char) * (strlen(cur_val) + 2));
			strcpy(this->seq1_path, cur_val);

			FILE* f_seq1 = fopen(this->seq1_path, "r");
			if(f_seq1 == NULL)
			{
				this->last_error_code = ERR_SEQ1_NOT_EXISTENT;
				fclose(f_conf);
				return;
			}
			else
			{
				fclose(f_seq1);
			}
		}
		else if(strcmp(cur_option, "seq2") == 0)
		{
			this->seq2_path = (char*)malloc(sizeof(char) * (strlen(cur_val) + 2));
			strcpy(this->seq2_path, cur_val);

			FILE* f_seq2 = fopen(this->seq2_path, "r");
			if(f_seq2 == NULL)
			{
				this->last_error_code = ERR_SEQ2_NOT_EXISTENT;
				fclose(f_conf);
				return;
			}
			else
			{
				fclose(f_seq2);
			}
		}
		else if(strcmp(cur_option, "mode") == 0)
		{
			if(strcmp(cur_val, "map") == 0)
			{
if(_DUMP_PPF_CLI_MESSAGES_)
				printf("Found mode map\n");
				this->mode = PARTS_RUN_MODE_MAP;
			}
			else if(strcmp(cur_val, "pp") == 0)
			{
				this->mode = PARTS_RUN_MODE_PP;
			}
			else if(strcmp(cur_val, "stochsample") == 0)
			{
				this->mode = PARTS_RUN_MODE_STOCH_SAMPLE;
			}
			else
			{
				this->last_error_code = ERR_MODE_NOT_VALID;
				fclose(f_conf);
				return;
			}
		}
		else if(strcmp(cur_option, "seq1_op_prefix") == 0)
		{
			this->seq1_op_file_prefix = (char*)malloc(sizeof(char) * (strlen(cur_val) + 2));
			strcpy(this->seq1_op_file_prefix, cur_val);
		}
		else if(strcmp(cur_option, "seq2_op_prefix") == 0)
		{
			this->seq2_op_file_prefix = (char*)malloc(sizeof(char) * (strlen(cur_val) + 2));
			strcpy(this->seq2_op_file_prefix, cur_val);
		}



		else if(strcmp(cur_option, "seq1_map_ct_op") == 0) // Set a file name for sampled cts: All sampled ct's are dumped into one file.
		{
			this->seq1_map_ct_op = (char*)malloc(strlen(cur_val) + 3);
			strcpy(this->seq1_map_ct_op, cur_val);
		}
		else if(strcmp(cur_option, "seq2_map_ct_op") == 0) // Set a file name for sampled cts: All sampled ct's are dumped into one file.
		{
			this->seq2_map_ct_op = (char*)malloc(strlen(cur_val) + 3);
			strcpy(this->seq2_map_ct_op, cur_val);
		}
		else if(strcmp(cur_option, "map_aln_op") == 0) // Set a file name for sampled cts: All sampled ct's are dumped into one file.
		{
			this->map_aln_op = (char*)malloc(strlen(cur_val) + 3);
			strcpy(this->map_aln_op, cur_val);
		}
		else if(strcmp(cur_option, "seq1_pp_op") == 0) // Set a file name for sampled cts: All sampled ct's are dumped into one file.
		{
			this->seq1_pp_op = (char*)malloc(strlen(cur_val) + 3);
			strcpy(this->seq1_pp_op, cur_val);
		}
		else if(strcmp(cur_option, "seq2_pp_op") == 0) // Set a file name for sampled cts: All sampled ct's are dumped into one file.
		{
			this->seq2_pp_op = (char*)malloc(strlen(cur_val) + 3);
			strcpy(this->seq2_pp_op, cur_val);
		}
		else if(strcmp(cur_option, "seq1_sample_ct_op") == 0) // Set a file name for sampled cts: All sampled ct's are dumped into one file.
		{
			this->seq1_sample_ct_op = (char*)malloc(strlen(cur_val) + 3);
			strcpy(this->seq1_sample_ct_op, cur_val);
		}
		else if(strcmp(cur_option, "seq2_sample_ct_op") == 0) // Set a file name for sampled cts: All sampled ct's are dumped into one file.
		{
			this->seq2_sample_ct_op = (char*)malloc(strlen(cur_val) + 3);
			strcpy(this->seq2_sample_ct_op, cur_val);
		}
		else if(strcmp(cur_option, "sample_aln_op") == 0) // Set a file name for sampled cts: All sampled ct's are dumped into one file.
		{
			this->sample_aln_op = (char*)malloc(strlen(cur_val) + 3);
			strcpy(this->sample_aln_op, cur_val);
		}
		else if(strcmp(cur_option, "nsamp") == 0)
		{
			this->n_samples = atoi(cur_val);
		}
                else if(strcmp(cur_option, "seed") == 0)
                {
                        this->stoch_sampling_seed = atoi(cur_val);
			//printf("Set stoch sampling seed to %d\n", this->stoch_sampling_seed);
                }
	}

	fclose(f_conf);

	//printf("Loaded conf file %s\n", conf_file);

	if(this->seq1_path == NULL)
	{
		this->last_error_code = ERR_SEQ1_NOT_EXISTENT;
		return;
	}

	if(this->seq2_path == NULL)
	{
		this->last_error_code = ERR_SEQ2_NOT_EXISTENT;
		return;
	}

	// Reset the prefixes.
	if(this->seq1_op_file_prefix == NULL)
	{
		this->seq1_op_file_prefix = (char*)malloc(sizeof(char) * (strlen("seq1") + 2));
		strcpy(this->seq1_op_file_prefix, "seq1");
	}

	if(this->seq2_op_file_prefix == NULL)
	{
		this->seq2_op_file_prefix = (char*)malloc(sizeof(char) * (strlen("seq2") + 2));
		strcpy(this->seq2_op_file_prefix, "seq2");
	}

}

//void t_ppf_cli::load_conf_file(char* conf_file)
//{
//	FILE* f_conf = open_f(conf_file, "r");
//	if(f_conf == NULL)
//	{
//		printf("Could not open configuration file %s @ %s(%d)\n", conf_file, __FILE__, __LINE__);
//		exit(0);
//	}
//
//	while(1)
//	{
//		char cur_option[1000];
//		char cur_value[1000];
//		if(fscanf("%s %s\n", cur_option, cur_value) != 2)
//		{
//			break;
//		}
//		else
//		{
//			if(strcmp(cur_option, "-proprate") == 0)
//			{
//				// Copy propagation rate constant from command line.
//				this->prior_to_posterior_prop_rate = atof(cur_value);
//			}
//			else if(strcmp(cur_option, "-fp1") == 0) // fold prior?
//			{
//if(_DUMP_PPF_CLI_MESSAGES_)
//				printf("Found fold prior1 file: %s\n", cur_value);
//
//				if(!check_file())
//				{
//					printf("Sequence 1 fold prior file could not be found @ %s exiting @ %s(%d)\n", cur_value, __FILE__, __LINE__);
//					exit(0);
//				}
//				else
//				{
//					this->fold_prior1_path = (char*)malloc(strlen(cur_value) + 3);
//					strcpy(this->fold_prior1_path, argv[p_cnt]);
//				}
//			}
//			else if(strcmp(cur_option, "-fp2") == 0) // fold prior?
//			{
//				p_cnt++; // Alignment prior file path is in argv[p_cnt].
//
//if(_DUMP_PPF_CLI_MESSAGES_)
//				printf("Found fold prior2 file: %s\n", cur_value);
//
//				if(!check_file(cur_value))
//				{
//					printf("Sequence 2 fold prior file could not be found @ %s, exiting @ %s(%d)\n", cur_value, __FILE__, __LINE__);
//					exit(0);
//				}
//				else
//				{
//					this->fold_prior2_path = (char*)malloc(strlen(cur_value) + 3);
//					strcpy(this->fold_prior2_path, argv[p_cnt]);
//				}
//	
//				// Jump over next flag.
//				p_cnt++;				
//			}
//			else if(strcmp(cur_option, "-str1_coinc_info") == 0) // fold prior?
//			{
//				p_cnt++; // Alignment prior file path is in argv[p_cnt].
//
//if(_DUMP_PPF_CLI_MESSAGES_)
//				printf("Found str co-inc info file for seq1: %s\n", argv[p_cnt]);
//
//				if(!check_file(argv[p_cnt]))
//				{
//					printf("Sequence 1 str co-inc info file could not be found @ %s, exiting @ %s(%d)\n", argv[p_cnt], __FILE__, __LINE__);
//					exit(0);
//				}
//				else
//				{
//					this->str1_coinc_info_path = (char*)malloc(strlen(argv[p_cnt]) + 3);
//					strcpy(this->str1_coinc_info_path, argv[p_cnt]);
//				}
//	
//				// Jump over next flag.
//				p_cnt++;				
//			}
//			else if(strcmp(cur_option, "-str2_coinc_info") == 0) // fold prior?
//			{
//				p_cnt++; // Alignment prior file path is in argv[p_cnt].
//
//if(_DUMP_PPF_CLI_MESSAGES_)
//				printf("Found str co-inc info file for seq2: %s\n", argv[p_cnt]);
//
//				if(!check_file(argv[p_cnt]))
//				{
//					printf("Sequence 2 str co-inc info file could not be found @ %s, exiting @ %s(%d)\n", argv[p_cnt], __FILE__, __LINE__);
//					exit(0);
//				}
//				else
//				{
//					this->str2_coinc_info_path = (char*)malloc(strlen(argv[p_cnt]) + 3);
//					strcpy(this->str2_coinc_info_path, argv[p_cnt]);
//				}
//	
//				// Jump over next flag.
//				p_cnt++;				
//			}
//			else if(strcmp(cur_option, "-ct1") == 0) // Correct ct1 path?
//			{
//				p_cnt++; // Alignment prior file path is in argv[p_cnt].
//
//if(_DUMP_PPF_CLI_MESSAGES_)
//				printf("Found correct ct1 file: %s\n", argv[p_cnt]);
//
//				if(!check_file(argv[p_cnt]))
//				{
//					printf("Sequence 1 correct ct file could not be found @ %s, exiting @ %s(%d)\n", argv[p_cnt], __FILE__, __LINE__);
//					exit(0);
//				}
//				else
//				{
//					this->correct_ct1_path = (char*)malloc(strlen(argv[p_cnt]) + 3);
//					strcpy(this->correct_ct1_path, argv[p_cnt]);
//
//					// Load correct 1st structure.
//					t_structure* ct1 = new t_structure(this->correct_ct1_path);
//
//
//if(_DUMP_PPF_CLI_MESSAGES_)
//{
//					// Write ct1 dot plot.
//					FILE* ct1_dot_plot_file = open_f("ct1_dp.txt", "w");
//					for(int i1 = 1; i1 <= ct1->numofbases; i1++)
//					{
//						for(int i2 = 1; i2 <= ct1->numofbases; i2++)
//						{
//							if(ct1->basepr[i1] == i2)
//							{
//								fprintf(ct1_dot_plot_file, "1 ");
//							}
//							else
//							{
//								fprintf(ct1_dot_plot_file, "0 ");
//							}
//						}
//						fprintf(ct1_dot_plot_file, "\n");
//					}
//					fclose(ct1_dot_plot_file);
//}
//				}	
//	
//				// Jump over next flag.
//				p_cnt++;			
//			}
//			else if(strcmp(cur_option, "-ct2") == 0) // Correct ct2 path?
//			{
//				p_cnt++; // Alignment prior file path is in argv[p_cnt].
//
//if(_DUMP_PPF_CLI_MESSAGES_)
//				printf("Found correct ct2 file: %s\n", argv[p_cnt]);
//
//				if(!check_file(argv[p_cnt]))
//				{
//					printf("Sequence 2 correct ct file could not be found @ %s, exiting @ %s(%d)\n", argv[p_cnt], __FILE__, __LINE__);
//					exit(0);
//				}
//				else
//				{
//					this->correct_ct2_path = (char*)malloc(strlen(argv[p_cnt]) + 3);
//					strcpy(this->correct_ct2_path, argv[p_cnt]);
//
//					// Load correct 2nd structure.
//					t_structure* ct2 = new t_structure(this->correct_ct2_path);
//
//if(_DUMP_PPF_CLI_MESSAGES_)
//{
//					// Write ct2 dot plot.
//					FILE* ct2_dot_plot_file = open_f("ct2_dp.txt", "w");
//					for(int i1 = 1; i1 <= ct2->numofbases; i1++)
//					{
//						for(int i2 = 1; i2 <= ct2->numofbases; i2++)
//						{
//							if(ct2->basepr[i1] == i2)
//							{
//								fprintf(ct2_dot_plot_file, "1 ");
//							}
//							else
//							{
//								fprintf(ct2_dot_plot_file, "0 ");
//							}
//						}
//						fprintf(ct2_dot_plot_file, "\n");
//					}
//					fclose(ct2_dot_plot_file);
//}
//				}
//	
//				// Jump over next flag.
//				p_cnt++;				
//			}
//			else if(strcmp(cur_option, "-pfdumpfp") == 0)
//			{
//				p_cnt++;
//
//				free(this->array_dump_fp);
//				this->array_dump_fp = (char*)malloc(sizeof(char) * (strlen(argv[p_cnt]) + 1));
//				strcpy(this->array_dump_fp, argv[p_cnt]);
//
//				// If the running mode is not dumping of pf arrays, 
//				// set the flag for loading this aray file.
//				if((this->mode & PARTS_RUN_MODE_DUMP_PF) != PARTS_RUN_MODE_DUMP_PF)
//				{
//					this->mode = this->mode | LOAD_PF_ARRAY_FLAG;
//				}
//
//				p_cnt++;
//			}
//			else if(strcmp(cur_option, "-aw") == 0) // Set log alignment weight per nucleotide?
//			{
//				p_cnt++; // double log alignment weight per nucleotide is @ argv[p_cnt].
//
//if(_DUMP_PPF_CLI_MESSAGES_)
//				printf("Overriding log alignment weight per alignment position with %f\n", atof(argv[p_cnt]));
//
//				this->log_aln_weight_per_aln_pos = atof(argv[p_cnt]);
//					
//				// Jump over next flag.
//				p_cnt++;				
//			}
//            else if(strcmp(cur_option, "-rfpn") == 0) // Set log alignment weight per nucleotide?
//            {
//                p_cnt++; // double log alignment weight per nucleotide is @ argv[p_cnt].
//if(_DUMP_PPF_CLI_MESSAGES_)
//                printf("Overriding rescaling factor per nucleotide with %f\n", atof(argv[p_cnt]));
//
//                this->cmd_rescaling_increment_factor_per_nucleotide = atof(argv[p_cnt]);
//
//                // Jump over next flag.
//                p_cnt++;
//            }
//            else if(strcmp(cur_option, "-nsamp") == 0) // Set log alignment weight per nucleotide?
//            {
//                    p_cnt++; // double log alignment weight per nucleotide is @ argv[p_cnt].
//if(_DUMP_PPF_CLI_MESSAGES_)
//                    printf("Sample set size is %d.\n", atoi(argv[p_cnt]));
//
//                    this->n_samples = atoi(argv[p_cnt]);
//
//                    // Jump over next flag.
//                    p_cnt++;
//            }
//			else if(strcmp(cur_option, "-sampled_cts1_op_fp") == 0) // Set a file name for sampled cts: All sampled ct's are dumped into one file.
//			{
//				p_cnt++;
//
//				this->sampled_cts1_op_fp = (char*)malloc(strlen(argv[p_cnt]) + 3);
//				strcpy(this->sampled_cts1_op_fp, argv[p_cnt]);
//
//                // Jump over next flag.
//                p_cnt++;
//			}
//			else if(strcmp(cur_option, "-sampled_cts2_op_fp") == 0) // Set a file name for sampled cts: All sampled ct's are dumped into one file.
//			{
//				p_cnt++;
//
//				this->sampled_cts2_op_fp = (char*)malloc(strlen(argv[p_cnt]) + 3);
//				strcpy(this->sampled_cts2_op_fp, argv[p_cnt]);
//
//                // Jump over next flag.
//                p_cnt++;
//			}
//			else if(strcmp(cur_option, "-sampled_cts_op_dir") == 0) // Set a directory for sampled cts: All structure samples are dumped into their own ct files.
//			{
//				p_cnt++;
//
//				this->sampled_cts_op_dir = (char*)malloc(strlen(argv[p_cnt]) + 3);
//				strcpy(this->sampled_cts_op_dir, argv[p_cnt]);
//
//				// Jump over next flag.
//                p_cnt++;
//			}
//			else if(strcmp(cur_option, "-sampled_alns_op_dir") == 0) // Set a directory for sampled alignments.
//			{
//				p_cnt++;
//
//				this->sampled_alns_op_dir = (char*)malloc(strlen(argv[p_cnt]) + 3);
//				strcpy(this->sampled_alns_op_dir, argv[p_cnt]);
//
//                // Jump over next flag.
//                p_cnt++;
//			}
//			else if(strcmp(cur_option, "-map_ct1_op_fp") == 0) // Set a file path for MAP structure of sequence 1.
//			{
//				p_cnt++;
//
//				this->map_ct_1_op_fp = (char*)malloc(strlen(argv[p_cnt]) + 3);
//				strcpy(this->map_ct_1_op_fp, argv[p_cnt]);
//
//				// Jump over next flag.
//                p_cnt++;
//			}
//			else if(strcmp(cur_option, "-map_ct2_op_fp") == 0) // Set a file path for MAP structure of sequence 2.
//			{
//				p_cnt++;
//
//				this->map_ct_2_op_fp = (char*)malloc(strlen(argv[p_cnt]) + 3);
//				strcpy(this->map_ct_2_op_fp, argv[p_cnt]);
//
//                // Jump over next flag.
//                p_cnt++;
//			}
//			else if(strcmp(cur_option, "-map_aln_op_fp") == 0) //  // Set a file path for MAP alignment.
//			{
//				p_cnt++;
//
//				this->map_alns_op_fp = (char*)malloc(strlen(argv[p_cnt]) + 3);
//				strcpy(this->map_alns_op_fp, argv[p_cnt]);
//
//                // Jump over next flag.
//                p_cnt++;
//			}
//			else if(strcmp(cur_option, "-max_n_separation") == 0)
//			{
//				p_cnt++;
//				this->max_n_separation_between_nucs = atoi(argv[p_cnt]);
//				p_cnt++;
//			}
//			else if(strcmp(cur_option, "-fold_env_prob_treshold") == 0)
//			{
//				p_cnt++;
//				this->fold_env_prob_treshold = atof(argv[p_cnt]);
//				p_cnt++;
//			}
//			else if(strcmp(cur_option, "-array_mem_limit_in_megs") == 0)
//			{
//				p_cnt++;
//				this->array_mem_limit_in_megs = atof(argv[p_cnt]);
//				p_cnt++;
//			}
//			else if(strcmp(cur_option, "-phmm_band_constraint_size") == 0)
//			{
//				p_cnt++;
//				this->phmm_band_constraint_size = atoi(argv[p_cnt]);
//				p_cnt++;
//			}
//			else if(strcmp(cur_option, "-str_coinc_env_prob_treshold") == 0)
//			{
//				p_cnt++;
//				this->str_coinc_env_prob_treshold = atof(argv[p_cnt]);
//				p_cnt++;
//			}
//			else if(strcmp(cur_option, "-str_coinc_est_sample_size") == 0)
//			{
//				p_cnt++;
//				this->str_coinc_estimator_sample_size = atoi(argv[p_cnt]);
//				p_cnt++;
//			}
//			else
//			{
//				p_cnt++;
//			}
//		}
//	};
//
//	fclose(f_conf);
//}

bool t_ppf_cli::have_correct_ct1()
{
	if(this->correct_ct1_path == NULL)
	{
		return(false);
	}
	else
	{
		return(true);
	}
}

bool t_ppf_cli::have_correct_ct2()
{
	if(this->correct_ct2_path == NULL)
	{
		return(false);
	}
	else
	{
		return(true);
	}
}

bool t_ppf_cli::have_fold_prior1()
{
	if(this->fold_prior1_path == NULL)
	{
		return(false);
	}
	else
	{
		return(true);
	}
}

bool t_ppf_cli::have_fold_prior2()
{
	if(this->fold_prior2_path == NULL)
	{
		return(false);
	}
	else
	{
		return(true);
	}
}

bool t_ppf_cli::have_aln_prior()
{
	if(this->aln_prior_path == NULL)
	{
		return(false);
	}
	else
	{
		return(true);
	}
}

// seq ids: Change all .'s to _'s and 
// tokenize string using / or \ as delimiters,
// get last token.
void t_ppf_cli::get_seq_aln_ids()
{
	char* cur_token = NULL;
	char previous_token[300];
	char _seq1_path[500];
	char _seq2_path[500];

	// Copy temporary strings.
	strcpy(_seq1_path, this->seq1_path);
	strcpy(_seq2_path, this->seq2_path);

	// Get seq1 id.
	cur_token = strtok(_seq1_path, "/\\");

	if(cur_token != NULL)
	{
		strcpy(previous_token, cur_token);
	}

	while(cur_token != NULL)
	{
		cur_token = strtok(NULL, "/\\");

		if(cur_token != NULL)
		{
			strcpy(previous_token, cur_token);
		}
	}

	this->seq1_id = (char*)malloc(sizeof(char) * (strlen(previous_token) + 3));
	// Have seq name in previous token.
	strcpy(this->seq1_id, previous_token);

	// Replace dots with underscores.
	for(int cnt = 0; cnt < strlen(this->seq1_id); cnt++)
	{
		if(this->seq1_id[cnt] == '.')
		{
			this->seq1_id[cnt] = '_';
		}
	}

if(_DUMP_PPF_CLI_MESSAGES_)
	printf("Sequence 1 ID: %s\n", this->seq1_id);

	// Get seq2 id.
	cur_token = strtok(_seq2_path, "/\\");

	if(cur_token != NULL)
	{
		strcpy(previous_token, cur_token);
	}

	while(cur_token != NULL)
	{
		cur_token = strtok(NULL, "/\\");

		if(cur_token != NULL)
		{
			strcpy(previous_token, cur_token);
		}
	}

	this->seq2_id = (char*)malloc(sizeof(char) * (strlen(previous_token) + 3));
	// Have seq name in previous token.
	strcpy(this->seq2_id, previous_token);

	// Replace dots with underscores.
	for(int cnt = 0; cnt < strlen(this->seq2_id); cnt++)
	{
		if(this->seq2_id[cnt] == '.')
		{
			this->seq2_id[cnt] = '_';
		}
	}

if(_DUMP_PPF_CLI_MESSAGES_)
	printf("Sequence 2 ID: %s\n", this->seq2_id);

	// Form alignment id.
	this->aln_id = (char*)malloc( sizeof(char) * (strlen(this->seq1_id) + strlen(this->seq2_id) + 3) );

	sprintf(this->aln_id, "%s_%s", this->seq1_id, this->seq2_id);

if(_DUMP_PPF_CLI_MESSAGES_)
	printf("Alignment ID: %s\n", this->aln_id);
}

#include <string.h>
#include <limits.h>
#include <time.h>
#include <stdio.h>
#include "ppf_timer.h"
#include "ppf_cli.h"
#include "parts_paths.h"
#include "../../../src/phmm/utils/file/utils.h"

t_ppf_timer::t_ppf_timer(t_ppf_cli* ppf_cli)
{
	this->start = (double)clock();
	this->ppf_cli = ppf_cli;
}

double t_ppf_timer::get_elapsed_time()
{
	double current = (double)clock();
	double elapsed_time =  (current - this->start) / CLOCKS_PER_SEC;
	return(elapsed_time);
}

void t_ppf_timer::log_elapsed_time(char* event_log)
{
	double elapsed_time = this->get_elapsed_time();

	FILE* timer_file = open_f(TIME_OP_FP, "a");
	// Log this elapsed time 
	fprintf(timer_file, "%s %s %10f\n", ppf_cli->aln_id, event_log, elapsed_time);

	fclose(timer_file);
}


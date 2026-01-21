#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "ppf_progress_bar.h"
#include "ppf_cli.h"
#include <string.h>

/*
class t_ppf_progress_bar
{
public:
	bool up_count;
	int range;
	char progress_char;
	t_ppf_cli* ppf_cli;
	t_ppf_progress_bar(t_ppf_cli* _ppf_cli, char _progress_char, bool _up_count, int _range);
	~t_ppf_progress_bar();

	// update_bar: redraws the bar.
	void update_bar(int progress);
};
*/

t_ppf_progress_bar::t_ppf_progress_bar(t_ppf_cli* _ppf_cli, char _progress_char, bool _up_count, int _range)
{
	this->ppf_cli = _ppf_cli;
	this->progress_char = _progress_char;
	this->up_count = _up_count;
	this->range = _range; // This is the counting range of progress bar, this value can be any int.
}

t_ppf_progress_bar::~t_ppf_progress_bar()
{
	// Destruct.
}

void t_ppf_progress_bar::update_bar(int progress)
{
	// Calculate # of progress_char(s) to print.
	int n_progress_chars;

	if(this->up_count)
	{
		n_progress_chars = (int) (((double)N_FULL_PROGRESS_CHAR / (double)range) * (double)progress);
	}
	else
	{
		n_progress_chars = (int) (((double)N_FULL_PROGRESS_CHAR / (double)range) * (double)(range - progress));
	}

	int perc_progress = 0;
	if(this->up_count)
	{
		perc_progress = (int)( ((double)progress / (double)range) * 100);
	}
	else
	{
		perc_progress = (int)( ((double)(range - progress) / (double)range) * 100);
	}

	printf(   "\r%d%%[", perc_progress);

	for(int cnt = 0; cnt < n_progress_chars - 1; cnt++)
	{
		printf("%c", this->progress_char);
	}

	// The n_proress_chars th character is the exciter_char.
	char exciter_chars[] = "\\|/-";

	if(this->up_count)
	{
		printf("%c", exciter_chars[progress % strlen(exciter_chars)]);
	}
	else
	{
		printf("%c", exciter_chars[(range - progress) % strlen(exciter_chars)]);
	}

	for(int cnt = n_progress_chars; cnt < N_FULL_PROGRESS_CHAR; cnt++)
	{
		printf(" ", this->progress_char);
	}

	//printf("]%c      ", exciter_chars[progress % strlen(exciter_chars)]);
	printf("]         ");
	fflush(stdout);
}


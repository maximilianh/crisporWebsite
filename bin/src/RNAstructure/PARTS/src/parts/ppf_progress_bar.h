#ifndef _PPF_PROGRESS_BAR_
#define _PPF_PROGRESS_BAR_

#define N_FULL_PROGRESS_CHAR (60) // Put 60 progress chars for completion.

class t_ppf_cli;

// A simple progress bar on command line.
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

#endif // _PPF_PROGRESS_BAR_


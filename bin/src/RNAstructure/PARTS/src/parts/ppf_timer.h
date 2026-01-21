#ifndef _PPF_TIMER_
#define _PPF_TIMER_

class t_ppf_cli;

class t_ppf_timer
{
public:
	double start;
	double end;
	t_ppf_cli* ppf_cli;

	t_ppf_timer(t_ppf_cli* ppf_cli);
	~t_ppf_timer();

	double get_elapsed_time();

	void log_elapsed_time(char* event_log);
};

#endif //_PPF_TIMER_


#ifndef _PP_RESULTS_
#define _PP_RESULTS_

class t_ppf_cli;

class t_pp_results
{
public:
	t_pp_results(t_ppf_cli* _ppf_cli, double** _posterior_bp1_probabilities, double** _posterior_bp2_probabilities);

	t_ppf_cli* ppf_cli;
	double** posterior_bp1_probabilities;
	double** posterior_bp2_probabilities;
};

#endif // _PP_RESULTS_


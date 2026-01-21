#include "pp_results.h"
#include <stdio.h>
#include <stdlib.h>

/*
#ifndef _PP_RESULTS_
#define _PP_RESULTS_

class t_ppf_cli;

class t_pp_results
{
public:
	t_pp_results(t_ppf_cli* _ppf_cli, double** _posterior_bp_probabilities);

	ppf_cli* ppf_cli;
	double** posterior_bp_probabilities;
};

#endif // _PP_RESULTS_*/

t_pp_results::t_pp_results(t_ppf_cli* _ppf_cli, double** _posterior_bp1_probabilities, double** _posterior_bp2_probabilities)
{
	this->ppf_cli = _ppf_cli;
	this->posterior_bp1_probabilities = _posterior_bp1_probabilities;
	this->posterior_bp2_probabilities = _posterior_bp2_probabilities;
}


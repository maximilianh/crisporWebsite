#include <math.h>
#include "ppf_operators.h"
#include "ppf_math.h"
#include "../../../src/phmm/phmm.h"

double ppf_max_callback(double x, double y)
{
	return(MAX(x,y));
}

double ppf_sum_callback(double x, double y)
{
	return(SUM(x,y));
}

bool ppf_compare_callback(double x, double y)
{
	return(COMPARE(x,y));
}

bool ppf_geq_callback(double x, double y)
{
	return(GEQ(x,y));
}

double ppf_sub_callback(double x, double y)
{
	return(SUB(x,y));
}

double ppf_exclude_callback(double x, double y)
{
	if(COMPARE(x, y))
	{
		return(x);
	}
	else
	{
		return(MAX(x,y));
	}
}


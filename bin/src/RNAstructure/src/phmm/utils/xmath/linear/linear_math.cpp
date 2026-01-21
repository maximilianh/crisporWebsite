#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include "nnm_math.h"
#include "../log/xlog_math.h"
#include "linear_math.h"

double lin_sum(double num1, double num2)
{
	return(num1 + num2);
}

double lin_mul(double num1, double num2)
{
	return(num1 * num2);
}

double lin_sub(double num1, double num2)
{
	return(num1 - num2);
}

double lin_div(double num1, double num2)
{
	return(num1 / num2);
}

double lin_pow(double num1, double num2)
{
	return(pow(num1, num2));
}

double lin_max(double num1, double num2)
{
	if(lin_compare(num1, num2))
	{
		return(num1);
	}
	else if(num1 > num2)
	{
		return(num1);
	}
	else
	{
		return(num2);
	}
}

bool lin_compare(double val1, double val2)
{
	// Perfect match returns true immediately.
	if(val1 == val2)
	{
		return(true);
	}

	// Do not do logarithmic check if any of the values are 0.
	if(val1 == 0 || val2 == 0)
	{
		return(false);
	}

	// Do a logarithmic type of epsilon-comparison on values.
	double log1 = xlog(val1);
	double log2 = xlog(val2);

	// If xlogs can be equated directly, return true.
	if(log1 == log2)
	{
		return(true);
	}

	//if(log1 <= log2 + XLOG_EPSILON && log1 >= log2 - XLOG_EPSILON)
	if(fabs(log1 - log2) <= XLOG_EPSILON)
	{
		return(true);
	}

	//printf("fabs: %.15f\nXLOG_EPSILON: %.15f\n", fabs(log1 - log2), XLOG_EPSILON);
	return(false);
}

// >= operator for double values.
bool lin_geq(double val1, double val2)
{
	if(lin_compare(val1, val2))
	{
		return(true);
	}

	if(val1 > val2)
	{
		return(true);
	}

	return(false);
}

bool lin_gt(double val1, double val2)
{
	// Do not do logarithmic check if any of the values are 0.
	if(val1 == 0 || val2 == 0)
	{
		if(lin_compare(val1, val2))
		{
			return(false);
		}

		if(val1 == 0)
		{
			return(false);
		}
		else
		{
			return(true);
		}
	}

	// Do a logarithmic type of epsilon-comparison on values.
	double log1 = xlog(val1);
	double log2 = xlog(val2);

	// If xlogs can be equated directly, return true.
	if(log1 == log2)
	{
		return(false);
	}

	if(log1 > log2 + XLOG_EPSILON)
	{
		return(true);
	}
	else
	{
		return(false);
	}
}

double get_linear_comp_prec()
{
	double prec = 1.0f;
	double one = 1.0f;
	while(1)
	{
		if(lin_compare(one, (one + prec)))
		{
			printf("%lf = %lf + %G\n", one, one, prec);
			break;
		}

		prec = prec / 2.0f;
	}

	return(prec);
}


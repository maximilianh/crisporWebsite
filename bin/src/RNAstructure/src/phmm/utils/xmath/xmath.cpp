#include "xmath.h"
#include <stdio.h>
#include <stdlib.h>
#include "log/xlog_math.h"
#include "linear/linear_math.h"

t_xmath::t_xmath(int _math_mode)
{
	this->math_mode = _math_mode;

	// Set the callbacks.
	if(this->math_mode == MATH_MOD_LIN)
	{
		this->sum = lin_sum;
		this->sub = lin_sub;
		this->mul = lin_mul;
		this->div = lin_div;
		this->pow = lin_pow;
		this->max = lin_max;

		this->compare = lin_compare;
		this->geq = lin_geq;
		this->gt = lin_gt;

		this->precision = get_linear_comp_prec;
	}
	else if(this->math_mode == MATH_MOD_LOG)
	{
		this->sum = xlog_sum;
		this->sub = xlog_sub;
		this->mul = xlog_mul;
		this->div = xlog_div;
		this->pow = xlog_pow;
		this->max = xlog_max;

		this->compare = xlog_comp;
		this->geq = xlog_geq;
		this->gt = xlog_gt;

		this->precision = get_xlog_comp_prec;
	}
	else
	{
		printf("Could not determine math type %d @ %s(%d)\n", this->math_mode, __FILE__, __LINE__);
		exit(0);
	}
}

t_xmath::~t_xmath()
{}
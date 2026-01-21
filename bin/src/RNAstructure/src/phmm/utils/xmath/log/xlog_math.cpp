#include <math.h>
#include "xlog_math.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

// MOVED TO HEADER: xlog_sum, xlog_sub, xlog_mul, xlog_div, xlog_pow 

// Determine whether the first xlog value represents a number that is greater 
//    than or equal to the number represented by the second xlog value 
//    or if the xlog values are equal within an absolute tolerance interval.
bool xlog_geq(const double& log1, const double& log2) {
	return log1 > log2 || xlog_comp(log1,log2);
}

// Determine whether the first xlog value represents a number that is greater 
//    than or equal to the number represented by the second xlog value 
//    AND the two xlog values are NOT equal to each other within an absolute
//    tolerance interval.
bool xlog_gt(const double& log1, const double& log2) {
	return log1 > LOG_OF_ZERO && log1+XLOG_EPSILON > log2;
}

// Compare two xlog values checking for intervals based on epsilon value.
// The interpretation of epsilon based comparison is: (x/y) < ~1 & (y/x) < ~1
bool xlog_comp(const double& log1, const double& log2)
{
	return log1 == log2 || 
		   (log1 <= LOG_OF_ZERO && log2 <= LOG_OF_ZERO) ||
		   (log1 <= log2 + XLOG_EPSILON && log1 >= log2 - XLOG_EPSILON);
}

// Get the precision of xlog comparison.
double get_xlog_comp_prec()
{
	double log_prec = log(1.0);
	const double ln1 = log(1);
	const double ln2 = log(2);
	while(1)
	{
		if(xlog_comp(ln1, xlog_sum(ln1, log_prec))) {
			printf("%lf = %lf + %G\n", ln1, ln1, log_prec);
			break;
		}
		log_prec = xlog_div(log_prec, ln2);
	}
	return log_prec;
}

// Return the xlog value that represents the maximum of the two numbers 
// represented by the xlog value arugments.
double xlog_max(const double& log_val1, const double& log_val2) {
	return std::max(log_val1, log_val2);
}

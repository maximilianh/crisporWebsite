// Determine the definitions of linear/log calculations.

#ifndef _PPF_MATH_
#define _PPF_MATH_

#include "parts_compilation_directives.h"

#ifdef _LINEAR_COMPUTATIONS_
#include "../../../src/phmm/utils/xmath/linear/linear_math.h"
	//#include <math.h>
	#define ZERO (0.0)
	#define CONVERT_FROM_LIN(x) (x) // Convert a linear value into current space.
	#define CONVERT_FROM_LOG(x) (exp(x)) // Convert a log value into current space.
	#define CONVERT_TO_LIN(x) (x) // Convert a linear value into current space.
	#define CONVERT_TO_LOG(x) (xlog(x)) // Convert a log value into current space.
	#define EPS LINEAR_EPSILON
	// Define four basic operations as four macros to maximize speed.
	#define SUM(log1, log2) ((log1) + (log2))
	#define MUL(log1, log2) ((log1) * (log2))
	#define DIV(log1, log2) ((log1) / (log2))
	#define SUB(log1, log2) ((log1) - (log2))
	#define MUL3(log1, log2, log3) (log1 * log2 * log3)
	#define MUL4(log1, log2, log3, log4) (log1 * log2 * log3 * log4)
	#define SUM3(log1, log2, log3) (log1 + log2 + log3)
	#define SUM4(log1, log2, log3, log4) (log1 + log2 + log3 + log4)

	#define COMPARE(val1, val2) ( lin_compare((val1), (val2)) )

	#define GEQ(x,y) (lin_geq(x,y))
	#define GT(x,y) (lin_gt(x,y))
#endif //LINEAR_CALCULATIONS

#ifdef _LOG_COMPUTATIONS_
#include "../../../src/phmm/utils/xmath/log/xlog__math.h"
	#define ZERO (xlog(0.0))
	#define CONVERT_FROM_LIN(x) (xlog(x)) // Convert a linear value into current space.
	#define CONVERT_FROM_LOG(x) (x) // Convert a log value into current space.
	#define CONVERT_TO_LIN(x) (xexp(x)) // Convert a linear value into current space.
	#define CONVERT_TO_LOG(x) (x) // Convert a log value into current space.
	#define EPS XLOG_EPSILON
	#define SUM(log1, log2) ( xlog_sum((log1), (log2)) )
	#define MUL(log1, log2) ( xlog_mul((log1), (log2)) )
	#define DIV(log1, log2) ( xlog_div((log1), (log2)) )
	#define SUB(log1, log2) ( xlog_sub((log1), (log2)) )
	#define COMPARE(log1, log2) ( xlog_comp((log1), (log2)) )
	#define MUL3(log1, log2, log3) (xlog_mul(log1, xlog_mul(log2, log3)))
	#define MUL4(log1, log2, log3, log4) (xlog_mul(log1, xlog_mul(log2, xlog_mul(log3, log4))))
	#define SUM3(log1, log2, log3) (xlog_sum(log1, xlog_sum(log2, log3)))
	#define SUM4(log1, log2, log3, log4) (xlog_sum(log1, xlog_sum(log2, xlog_sum(log3, log4))))

	#define GEQ(x,y) (xlog_geq(x,y))
	#define GT(x,y) (xlog_gt(x,y))
#endif // LOG_CALCULATIONS

#endif // _PPF_MATH_


#ifndef _XLOG_MATH_
#define _XLOG_MATH_

//#include "math.h"
#include <cmath>
#include <limits>
#include <iostream>
#include <cfloat>
#include <cassert>
#include "../../../../defines.h"
#include <stdexcept>

// This header contains extended logarithmic/exponential functions derived from 
// methods described in TP Mann's 2006 paper:
//   "Numerically stable hidden Markov model implementation." (Tobias P. Mann, 2006)
//   http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf
// Our functions are more optimized than those Mann described, and they use a
// different value to represent LOG_OF_ZERO.
// The functions here are used to implement PPF_Math, which is used in the PHMM forward-backward algorithm.
// They were later adapted to implement the log-scale calculations used in the partition function (pfunction).

#define LOG_OF_ZERO_LEGACY (-1.0 * pow(2.0, 30))

#ifdef USE_XLOG_ZERO
	static const double LOG_OF_ZERO = -log(DBL_MAX)*1000;
	#define IS_LOG_ZERO(x) (x<=LOG_OF_ZERO)
#else
	// Log-scale values less than or equal to LOG_OF_ZERO are assumed to represent zero.
	static const double LOG_OF_ZERO = -std::numeric_limits<double>::infinity();
	#define IS_LOG_ZERO(x) false
	//#define LOG_OF_ZERO -1073741824 // -2^(30), can be represented exactly.
#endif
static const double LOG_OF_ONE = 0.0;

#if defined USE_LOGP1_LOOKUP_SUM || defined USE_LOGP1_LOOKUP_DIFF
#include "log_lookup.h"
#endif

#ifdef USE_LOGP1_LOOKUP_SUM
static log_lookup_sum LogTable = log_lookup_sum(-50, LOG_LOOKUP_ENTRIES); // use default number of steps.
#endif
#ifdef USE_LOGP1_LOOKUP_DIFF
static log_lookup_sub LogTableDiff = log_lookup_sub(50,50000); // no longer needed, because subtractions have been removed from pfunction
#endif

// Calculates log1p(exp(x))  -- used for calculating xlog_sum.
inline double log1pexp(const double& x) noexcept {
	if (IS_LOG_ZERO(x)) return 0.0; 	// if x <= LOG_OF_ZERO, then exp(x) is 0, and log(1+0)=log(1)=0, so return 0.
	#ifdef USE_LOGP1_LOOKUP_SUM
	return LogTable.exp_lu(x);
	#else
	return log1p(exp(x));
	#endif
}

// calculates log1p(-exp(x))  -- used for calculating xlog_sub.
inline double log1pNexp(const double& x) noexcept {
	if (IS_LOG_ZERO(x)) return 0.0; 	// if x <= LOG_OF_ZERO, then exp(x) is 0, and log(1-0)=log(1)=0, so return 0.
	#ifdef USE_LOGP1_LOOKUP_DIFF
	return LogTableDiff.exp_lu(x);
	#else
	return log1p(-exp(x));
	#endif
}

// Convert a number (in standard, linear scale) to a log-scale value (i.e. stored as the natural logarithm of the original number).
double xlog(const double& linear_value) noexcept;
// Convert a log-scale value back to its original linear-scale number (by taking the exponent of it, base e).
double xexp(const double& xlog) noexcept;
// Calculates the product of two log-scale values.
double xlog_mul(const double& xlog1, const double& xlog2);
// Calculates the quotient of two log-scale values.
double xlog_div(const double& xlog1, const double& xlog2);
// Calculates the sum of two log-scale values.
double xlog_sum(const double& xlog1, const double& xlog2);
// Calculates the difference between two log-scale values.
double xlog_sub(const double& xlog1, const double& xlog2);
// Calculates the result of raising a log-scale value to a power. (The power should be in linear scale.)
double xlog_pow(const double& xlog, const double& power);

inline double xlog(const double& value) noexcept {
	#ifdef USE_XLOG_ZERO
	if (value == 0.0) return LOG_OF_ZERO; // the log function itself will throw an exception if value < 0.
	#endif
	return log(value); // the log function itself will throw an exception if value < 0.
}

inline double xexp(const double& xlog) noexcept {
	#ifdef USE_XLOG_ZERO
	if (xlog <= LOG_OF_ZERO) return 0.0;
	#endif
	return exp(xlog);
}

// Computes log(exp(a)+exp(b))
inline double xlog_sum(const double& a, const double& b) {
	// Derivation:   log(exp(a)-exp(b))  
	//             = log(exp(a) * (1-exp(b-a)) )
	//             = a+log(1-exp(b-a))
	//             = a+log1p(-exp(b-a))
    #ifdef USE_XLOG_ZERO
	if(IS_LOG_ZERO(a)) return b;
    if(IS_LOG_ZERO(b)) return a;
	#endif
	// Note: The test of a>b is important when A or B is greater than MAX_DOUBLE.
	// e.g. A=1E+500 and B=1 --- exp(a-b) will overflow, while exp(b-a) will not.
	// As long as a>=b, the value of exp(b-a) will always be in the range (0,1], 
	// so it will not overflow (even when one or both numbers are greater than 
	// MAX_DOUBLE in the linear scale).
	// In fact, the result will always be max(a,b)+Q where Q is in the range [0..ln(2)].
	return a>b ? a+log1p(exp(b-a)): b+log1p(exp(a-b));
	/* 
	// -----Possible alternatives (todo: test for speed)-----
    // -----alt.1----- (temporary vars for min,max)
	const double mx = std::max(a,b), mn=std::min(a,b);
    return mx+log1pexp(mn-mx);
	// -----alt.2----- (implicit compiler-defined temporary vars for min,max)
    return std::max(a,b)+log1pexp(std::min(a,b)-std::max(a,b));
	// -----alt.3----- (max and -abs)
	return std::max(a,b)+log1pexp(-std::abs(a-b));
	*/
}

inline double xlog_sum2(const double& a, const double& b) {
	// Derivation:   log(exp(a)-exp(b))  
	//             = log(exp(a) * (1-exp(b-a)) )
	//             = a+log(1-exp(b-a))
	//             = a+log1p(-exp(b-a))
    #ifdef USE_XLOG_ZERO
	if(IS_LOG_ZERO(a)) return b;
    if(IS_LOG_ZERO(b)) return a;
	#endif
	// Note: The test of a>b is important when A or B is greater than MAX_DOUBLE.
	// e.g. A=1E+500 and B=1 --- exp(a-b) will overflow, while exp(b-a) will not.
	// As long as a>=b, the value of exp(b-a) will always be in the range (0,1], 
	// so it will not overflow (even when one or both numbers are greater than 
	// MAX_DOUBLE in the linear scale).
	// In fact, the result will always be max(a,b)+Q where Q is in the range [0..ln(2)].
	return a>b ? a+log1pexp(b-a): b+log1pexp(a-b);
	/* 
	// -----Possible alternatives (todo: test for speed)-----
    // -----alt.1----- (temporary vars for min,max)
	const double mx = std::max(a,b), mn=std::min(a,b);
    return mx+log1pexp(mn-mx);
	// -----alt.2----- (implicit compiler-defined temporary vars for min,max)
    return std::max(a,b)+log1pexp(std::min(a,b)-std::max(a,b));
	// -----alt.3----- (max and -abs)
	return std::max(a,b)+log1pexp(-std::abs(a-b));
	*/
}
inline double xlog_sum3(const double& a, const double& b) {
	double c, d, x;
	if (a<b){
		c = b;
		d = a;
	}
	else{
		c = a;
		d = b;
	}
	x = c-d;

	if (x > 11.8624794162)
		return d + x;
	if (x < double(3.3792499610))
	{
		if (x < double(1.6320158198))
		{
			if (x < double(0.6615367791))
				return d + ((double(-0.0065591595)*x+double(0.1276442762))*x+double(0.4996554598))*x+double(0.6931542306);
			return d + ((double(-0.0155157557)*x+double(0.1446775699))*x+double(0.4882939746))*x+double(0.6958092989);
		}
		if (x < double(2.4912588184))
			return d + ((double(-0.0128909247)*x+double(0.1301028251))*x+double(0.5150398748))*x+double(0.6795585882);
		return d + ((double(-0.0072142647)*x+double(0.0877540853))*x+double(0.6208708362))*x+double(0.5909675829);
	}
	if (x < double(5.7890710412))
	{
		if (x < double(4.4261691294))
			return d + ((double(-0.0031455354)*x+double(0.0467229449))*x+double(0.7592532310))*x+double(0.4348794399);
		return d + ((double(-0.0010110698)*x+double(0.0185943421))*x+double(0.8831730747))*x+double(0.2523695427);
	}
	if (x < double(7.8162726752))
		return d + ((double(-0.0001962780)*x+double(0.0046084408))*x+double(0.9634431978))*x+double(0.0983148903);
	return d + ((double(-0.0000113994)*x+double(0.0003734731))*x+double(0.9959107193))*x+double(0.0149855051);
}

// Computes log(exp(a)-exp(b))
inline double xlog_sub(const double& a, const double& b) { // (previously xlog_sub)
	// Derivation:   log(exp(a)-exp(b))  
	//             = log(exp(a) * (1-exp(b-a)) )
	//             = a+log(1-exp(b-a))
	//             = a+log1p(-exp(b-a))
	#ifdef USE_XLOG_ZERO
	// `a` must be >= `b` so that exp(a) >= exp(b), otherwise the value inside the log would be negative.
	if (b<=LOG_OF_ZERO) return a; // Important to do this test, because (a<b) is NOT an error if (b <= LOG_OF_ZERO)
	if (a<b) throw std::runtime_error("Subtraction of xlog values resulted in an unrepresentable negative number. (in " __FILE__ ")");
	return a==b ? LOG_OF_ZERO : a+log1pNexp(b-a);  // note that (b-a) is always < 0, which is required because we are computing log(1-exp(b-a))
	#else
	return a+log1pNexp(b-a);  // note that (b-a) is always < 0, which is required because we are computing log(1-exp(b-a))
	#endif
}

inline double xlog_mul(const double& log1, const double& log2) {
	#ifdef USE_XLOG_ZERO
	return (log1 <= LOG_OF_ZERO || log2 <= LOG_OF_ZERO) ? LOG_OF_ZERO : log1+log2;
	#else
	return log1+log2;
	#endif
}

// Returns 0 if log1 is 0 no matter what log2 is.
inline double xlog_div(const double& log1, const double& log2) {
	#ifdef USE_XLOG_ZERO
	if(log1 <= LOG_OF_ZERO) return LOG_OF_ZERO;
	if(log2 <= LOG_OF_ZERO) throw std::runtime_error("Division by xlog zero-value (in " __FILE__ ")" );
	#endif
	return log1-log2;
}

inline double xlog_pow(const double& log_value, const double& pow) {
	#ifdef USE_XLOG_ZERO
	return log_value <= LOG_OF_ZERO ? LOG_OF_ZERO : log_value * pow;
	#else
	return log_value * pow;
	#endif
}



// ###### PPF_Math ###### //
// The following are legacy functions from PPF_Math (used by phmm, e.g. src/phmm/utils/xmath/xmath.cpp)

//The epsilon value that defines the absolute tolerance interval used to compare 
// xlog values in the PPF_Math functions xlog_comp, xlog_gt, and xlog_geq.
// This was previously set to 1E-8, which assumed that because exp(1E-8) is a
// small number, that it made for a reliable comparison.
#define XLOG_EPSILON 1E-10

// Determine if two xlog values are equal, within an absolute tolerance interval given by XLOG_EPSILON.
bool xlog_comp(const double& xlog1, const double& xlog2);

// Determine whether the first xlog value represents a number that is greater 
//    than or equal to the number represented by the second xlog value 
//    or if the xlog values are equal within an absolute tolerance interval.
bool xlog_geq(const double& xlog1, const double& xlog2);

// Determine whether the first xlog value represents a number that is greater 
//    than or equal to the number represented by the second xlog value 
//    AND the two xlog values are NOT equal to each other within an absolute
//    tolerance interval.
bool xlog_gt(const double& xlog1, const double& xlog2);

// Get the precision of xlog comparison.
double get_xlog_comp_prec();

// Return the xlog value that represents the maximum of the two numbers 
// represented by the xlog value arugments.
double xlog_max(const double& xlog1, const double& xlog2);

// ###### ^ PPF_Math ^ ###### //

#endif // _XLOG_MATH_

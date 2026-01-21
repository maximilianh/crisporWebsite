#ifndef PFUNCTION_MATH_H
#define PFUNCTION_MATH_H

#include "defines.h"
#include "phmm/utils/xmath/log/xlog_math.h"
#include "log_double.h"
#include <type_traits>

//PFPRECISION is the datatype used by the partition function. Options are:  double, log_double, float, or EXTENDED_DOUBLE
#ifdef PF_LOG_CLASS
	#define PFPRECISION log_double
	//typedef log_double PFPRECISION;
#else
	#define PFPRECISION double
	//typedef double PFPRECISION;
	#ifdef PF_LINEAR
		#define PF_NUC_SCALING
	#else
		#define PF_LOG_CALC // If PFPRECISION is defined as double, also define PF_LOG_CALC to use log-based values instead of linear scaling.
	#endif
#endif

// #define LOG_LOOKUP_ENTRIES 10000 // default number of lookup entries (aka steps) in log_lookup_sum (xlog_math,log_lookup)
 
 // binary_operation is the type of a pointer to a function that takes 
 // two arguments of the same type and returns a value of the same type.
template<typename T> using binary_operation = T(*)(const T&,const T&);

// Perform left-associative aggregate (aka "left fold") using a binary function.
// Operation should be a binary_operation.
template<typename R, binary_operation<R> Operation, typename T>  inline T binary_reduce(const T& t) { return t; }
template<typename R, binary_operation<R> Operation, typename T, typename U, typename... Rest> 
inline R binary_reduce(const T& t, const U& u, Rest... rest) {
	return binary_reduce<R,Operation>(Operation(t,u), rest...);
}
// template<typename T, T func(const T&,const T&)>
// struct make_operation { 
// 	inline T operator()(const T& t, const T& u)  { return func(t,u); }
// };

// template<typename Operation, typename T>  inline T binary_reduce(T t) { return t; }
// template<typename Operation, typename T, typename U, typename... Rest> 
// inline auto binary_reduce(const T& t, const U& u, Rest... rest) 
// 	-> typename std::common_type<T,U>::type {
// 	return binary_reduce<Operation>(Operation()(t,u), rest...);
// }


#ifdef PF_LOG_CALC
	// Calculate the partition function using explicit log values (typed as doubles)

	// // Terminal template pattern for log_product
    // template <typename T>
	// inline T log_product(T t) {	
    //     return t;
    // }

	// // Calculate the product of two log-scale values (recursive template)
	// template <typename T, typename... Rest>
	// inline T log_product(T t, Rest... rest) {
	// 	return xlog_mul(t, log_product(rest...));
	// }

    // // Terminal template pattern for log_sum
	// template <typename T>
	// inline T log_sum(T t) { return t; }

	// // Calculate the sum of two log-scale values (recursive template)
	// template <typename T, typename... Rest>
	// inline T log_sum(T t, Rest... rest) {
	// 	return xlog_sum(t,log_sum(rest...));
	// }
	#define ZERO LOG_OF_ZERO
	#define ONE 0.0
	#define SUM(...)  binary_reduce<PFPRECISION,xlog_sum2>(__VA_ARGS__)
	#define DIFF(...)  xlog_sub(__VA_ARGS__)
	#define PROD(...) binary_reduce<PFPRECISION,xlog_mul>(__VA_ARGS__)
	#define DIV(...) xlog_div(__VA_ARGS__)
	#define POWER(...) xlog_pow(__VA_ARGS__)

	// Convert a linear-scale value to a log-scale value (only used when PF_LOG_CALC is defined -- otherwise returns arg as-is)
	#define TO_XLOG(X) xlog(X)
    // Convert a log-scale value to linear scale (only used when PF_LOG_CALC is defined -- otherwise returns arg as-is)
	#define TO_LINEAR(X) xexp(X)
	// Convert a log-scale value to PFPRECISION. Used only when reading log-scale values
	// from CUDA into CPU code. If PF_LOG_CALC is defined, then RNA also uses log-scale
	// values, so no conversion is necessary.
	#define XLOG_TO_PF(X) X
	
    // Take the log of a PFPRECISION value.
    // If PF_LOG_CALC is defined, this would be xlog(TO_LINEAR(log_value)) = xlog(xepx(log_value)) 
    //   = log_value   .. i.e. return the argument as-is
	#define PF_LOG_PFPRECISION(X) X
	
    // PF_EXP_LINEAR -- Take the exp of a *LINEAR-scale* PFPRECISION value 
	// (i.e. one that has NOT been converted into log-scale).
	// Therefore we can consider the number *already* implicitly the exp of the
	// log-scale representation of the argument -- i.e. exp(log(X))
    // So if PF_LOG_CALC is defined, just return the argument as-is.
	#define PF_EXP_LINEAR(X) X
#else
	// Calculate the partition function using either linear scale or implicit log values (typed as log_double)

    // // Terminal template pattern for linear_sum
	// template <typename T>
	// inline PFPRECISION linear_sum(const T& t) { 
	// 	return t;
	// }

	// // Linear-scale counterpart to log_sum -- simply calculates the sum of all arguments. (recursive template)
	// template <typename T, typename... Rest>
	// inline PFPRECISION linear_sum(const T& t, Rest... rest) {
	// 	return t + linear_sum(rest...);
	// }

    // // Terminal template pattern for linear_product
	// template <typename T>
	// inline PFPRECISION linear_product(const T& t) {
	// 	return t;
	// }

	// // Linear-scale counterpart to log_product -- simply calculates the product of all arguments. (recursive template)
	// template <typename T, typename... Rest>
	// inline PFPRECISION linear_product(T t, Rest... rest) {
	// 	return t * linear_product(rest...);
	// }

	//template<typename T, typename T>
	inline PFPRECISION linear_sum(const PFPRECISION& left, const PFPRECISION& right) { return left + right; }
	inline PFPRECISION linear_mul(const PFPRECISION& left, const PFPRECISION& right) { return left * right; }

    // Linear counterpart of log_sub -- calculates the difference between left and right (i.e.: left - right)
	inline PFPRECISION linear_sub(const PFPRECISION& left, const PFPRECISION& right) {
		return left-right;
	}
    // Linear counterpart of log_div -- calculates the quotient of left and right (i.e.: left / right)
	inline PFPRECISION linear_div(const PFPRECISION& left, const PFPRECISION& right) {
		return left / right;
	}

	#define SUM(...)  binary_reduce<PFPRECISION,linear_sum>(__VA_ARGS__)
	#define DIFF(...) linear_sub(__VA_ARGS__)
	#define PROD(...) binary_reduce<PFPRECISION,linear_mul>(__VA_ARGS__)
	#define DIV(...)  linear_div(__VA_ARGS__)
	#define POWER(...) pow(__VA_ARGS__)
	#define PF_EXP_LINEAR(...) exp(__VA_ARGS__)      // exp of a PFPRECISION value that has NOT been converted to log scale.
	#define PF_LOG_PFPRECISION(...) log(__VA_ARGS__) // log of a PFPRECISION value

	// Explicitly convert PFPRECISION (e.g. log_double) to double
	#define TO_LINEAR(X) (double)(X) // static_cast<double>(X)

	// XLOG_TO_PF -- Convert a log-scale value to PFPRECISION. Used only when reading log-scale values
	//     from CUDA into CPU code. 
	//     Behavior depends on what PFPRECISION is defined as (e.g. double vs log_double)
	#ifdef PF_LOG_CLASS
		// Convert log value to PFPRECISION (i.e. log_double) -- Initialize a new log_double from the log-scale double.
		#define XLOG_TO_PF(X) ld_from_xlog(X)
		#define TO_XLOG(X)    ld_from_linear(X) // explicitly convert double to PFPRECISION (e.g. log_double)
		#define ZERO log_double::zero()
		#define ONE log_double::one()
	#else
		// Convert log value to PFPRECISION (i.e.  linear double) -- So take exp.
		#define XLOG_TO_PF(X) xexp(X)
		// For initializing PFPRECISION with a linear-scale value
		#define TO_XLOG(X)  (PFPRECISION)(X) //static_cast<PFPRECISION>(X)
		#define ZERO 0.0
		#define ONE 1.0
	#endif
#endif //PF_LOG_CALC

#ifdef PF_NUC_SCALING
	// Using linear scale (PF_LINEAR) 
	// So per-nucleotide scaling is requried in the partition function to prevent overflow.
	#define PFSCALE(...) pfscale(__VA_ARGS__)
	// Scales the value by scaling_factor raised to the power nucleotide_count. I.e.: value * pow(scaling_factor, nucleotide_count)
	inline PFPRECISION pfscale(const PFPRECISION &value, const double& scaling_factor, const int nucleotide_count) {
		return value*pow(scaling_factor, nucleotide_count);
	}
	#define SCALING , data->scaling
	#define SCALING2 , scaling
	#define TWOSCALING , twoscaling
	#define POWSCALING(X) , pow(data->scaling,X)
	#define POWSCALING2(X) , pow(scaling,X)
#else
	// Using log scale (PF_LOG_CLASS or PF_LOG_CALC) 
	// So do not need per-nucleotide scaling in the partition function.
	#define PFSCALE(X, ...) X
	#define SCALING
	#define SCALING2
	#define TWOSCALING
	#define POWSCALING(X)
	#define POWSCALING2(X)
#endif

// Additional scaling macros
#define PFUNSCALE(V,S,P) PFSCALE(V,S,-P)
#define PFSCALE2(V,S) PFSCALE(V,S,2)
#define PFUNSCALE2(V,S) PFSCALE(V,S,-2)
#define PFSCALE_DS(V,P) PFSCALE(V,data->scaling,-P)
#define PFSCALE_S(V,P) PFSCALE(V,scaling,-P)

// Explicit conversion of double to PFPRECISION (used for troubleshooting log_double)
#define TO_PFP(X) TO_XLOG(X)

#endif //PFUNCTION_MATH_H

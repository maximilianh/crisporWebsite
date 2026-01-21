#if !defined(BOLTZMANN_H)
#define BOLTZMANN_H
#include "../src/phmm/utils/xmath/log/xlog_math.h"

/*======================================================================
// define boltzman() *once* in the code -- this is the place right now. 
// is partition.h is a better place?
//
// also removed a version of boltzman that accepts ints -- this
//  causes weird results, e.g., in calculating partition function
//  and base pair probabilities with SHAPE or 'EX' data. (rhiju das, 2011)
//
========================================================================*/
#ifdef PF_LOG_CALC
inline double boltzman(double i, double temp) {
	if (i>=INFINITE_ENERGY) return LOG_OF_ZERO;
	else return -i/conversionfactor/(RKC*temp);
}
#else
inline double boltzman(double i, double temp) {
	if (i>=INFINITE_ENERGY) return 0;
	else return exp(-i/conversionfactor/(RKC*temp));
}
#endif
// inline log_double boltzman(double i, double temp) {
// 	log_double out;
// 	if (i>=INFINITE_ENERGY) out._value = LOG_OF_ZERO;
// 	else out._value = -i/conversionfactor/(RKC*temp);
// 	return out;
// }

// inline PFPRECISION boltzman(double i, PFPRECISION temp) {

// 	if (i==INFINITE_ENERGY) return 0;
// 	else return (PFPRECISION) (exp((-(i)/((double)conversionfactor))/((double) RKC * (double) temp)));

// }

#endif



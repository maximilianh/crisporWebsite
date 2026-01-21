#ifndef ZERO_H
#define ZERO_H
//#include "logdouble.h"
#include "../src/phmm/utils/xmath/log/xlog_math.h"

template<typename T>
inline T zero() {}

template <> inline int zero<int>(){
	return 10000000;
}

template <> inline log_double zero<log_double>(){
	return 0.0;
}

template <> inline double zero<double>(){
	return 0.0;
}

template<typename T>
inline T one() {}

template <> inline int one<int>(){
	return 0;
}

template <> inline log_double one<log_double>(){
	return 1.0;
}

template <> inline double one<double>(){
	return 1.0;
}

#define IDENTITY zero<T>()
#define ONE one<T>()

//template <> real_t zero<real_t>(){
//	return Z;
//}

#endif

#if defined(WIN32) && !defined(__GNUC__)
	#include <math.h>
#else
	#include <mathimf.h>
#endif

double gammafunction(double arg) {
	return tgamma(arg);
}


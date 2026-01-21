#ifndef INT_H
#define INT_H

#include <math.h>
#include <stdlib.h>
#if defined SHORT

typedef short int_t;
#define conversion_factor 10.0f
//#define round(x) (floor(x>0? x+0.5 : x-0.5))
#define RF "%d"
#define RCONST(x) x
#define ALMOST_ZERO RCONST(1e-7)
#define MATHFN(x) x##f

#elif defined INT

typedef int int_t;
#define conversion_factor 100.0f
//#define round(x) (floor(x>0? x+0.5 : x-0.5))
#define RF "%d"
#define RCONST(x) x
#define ALMOST_ZERO RCONST(1e-15)
#define MATHFN(x) x

#elif defined LONG

typedef long int int_t;
#define conversion_factor 100.0f
//#define round(x) (floor(x>0? x+0.5 : x-0.5))
#define RF "%d"
#define RCONST(x) x##L
#define ALMOST_ZERO RCONST(1e-19)
#define MATHFN(x) x##l

#else

#error Must define one of SHORT, INT, or LONG

#endif

//inline int_t STR_TO_INT(const char* x,char** y) ((int_t) round(strtof(x,y)*conversion_factor))
#define round(x) (floor( x+0.5f )
//#define STR_TO_INT(x,y) (round(strtof(x,y)*conversion_factor)
#define STR_TO_INT(x,y) floor(strtof(x,y)*conversion_factor + 0.5f)
#define LOG MATHFN(log)
#define EXP MATHFN(exp)
#define LOG1P MATHFN(log1p)
#define COS MATHFN(cos)
#define SIN MATHFN(sin)
#define SQRT MATHFN(sqrt)
#define RINT MATHFN(rint)
#define FABS MATHFN(fabs)
#define HALF RCONST(0.5)
#define PI RCONST(3.1415926535897932384626433832795029)
#define INF (14000) //the GPU doesn't like INT_MAX or SHORT_MAX
#define NOT_A_NUMBER 14000 //an arbitrary large number to represent infinity

//static inline real_t sq(real_t x) { return x*x; }

#endif /* INT_H */

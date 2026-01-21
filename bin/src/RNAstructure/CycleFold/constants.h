#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__
#include <map>
#include "../src/log_double.h"

/*
enum NCM_type {single3,single4,single5,single6,
	double21,double12,double13,double22,double23,
	double25,double31,double32,double33,double34,
	double35,double42,double43,double44,double52,double53};

enum nucs {AA,AU,AG,AC,UA,UU,UG,UC,GA,GU,GG,GC,CA,CU,CG,CC};
enum pairtype {cWW,cWH,cHW,cWS,cSW,cHH,cHS,cSS,cSH,
				tWW,tWH,tHW,tWS,tSW,tHH,tHS,tSS,tSH};
*/
typedef std::string NCM_type;
typedef std::string pairtype;
typedef std::string nucs;
typedef log_double real_t;

enum orientation {INTERIOR, EXTERIOR};
const int num_pairtypes = 18;
const int num_NCMs = 21;
const int num_orientations = 2;
const int BIGGEST_NCM = 12;


#endif //__CONSTANTS_H__

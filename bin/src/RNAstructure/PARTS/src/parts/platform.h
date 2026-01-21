#if !defined _PLATFORM_ 
#define _PLATFORM_

#include "structure.h"
#include <math.h>

#define sgifix strcat(line," ");

//#define CHAR_BIT (8)

inline int min(int x, int y)
{
	int r = y + ((x - y) & ((x - y) >> (sizeof(int) * CHAR_BIT - 1))); // min(x, y)

	return(r);
}

inline int max (int x,int y)
{
	int r = x - ((x - y) & ((x - y) >> (sizeof(int) * CHAR_BIT - 1))); // max(x, y)

	return(r);
}

int pow10(int i);

#endif

/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 */

#ifndef PLATFORM_H
#define PLATFORM_H

#define sgifix strcat(line," ");

inline int min(int i, int j) {
  return i < j ? i : j;
}

inline int max (int i, int j) {
  return i < j ? j : i;
}

inline int pow10(int n) {
	int j = 1;

	for (int i = 0; i < n; i++) {
		j *= 10;
	}

	return j;
}



#endif

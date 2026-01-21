/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005,2006
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 */

#ifndef DYNALIGNARRAY_H
#define DYNALIGNARRAY_H

#include "defines.h"

// This class encapsulates the large 4-dimensional arrays used
// by Dynalign

//Note that highlimit[i] and lowlimit[i] are the spans of allowed alignments in seq2 for nucletide 
//i from sequence 1.  IMPORTANT:  These arrays must persist until after the call of ~dynalignarray.
class dynalignarray {
	private:
  short *Lowlimit,*Highlimit;          
  short N,N2,Ndiff; // length of sequence one and two
  bool optimalonly;
  short infinite;

	public:
		integersize ****array;
		dynalignarray(short n1, short n2, short *lowlimit, short *highlimit, bool Optimalonly=false);
		dynalignarray();
		~dynalignarray();
		void allocate(short n1, short n2, short *lowlimit, short *highlimit, bool Optimalonly=false);
		integersize &f(short i, short j,short k, short l);
};

inline integersize &dynalignarray::f(short i, short j,short k, short l) {
  if (i>N&&j>N) {
    i -= N;
    j -= N;
    k -= N2;
    l -= N2;
  }
  return array[i][j][k][l];
}

#endif

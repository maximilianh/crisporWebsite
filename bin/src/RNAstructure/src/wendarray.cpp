/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005,2006
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 */

//The wendarray class provides functionality for encapsulating the w3 and w5 arrays.

#include "wendarray.h"

#include "defines.h"


wendarray::wendarray() {
	

}

wendarray::wendarray(short n1, short n2, short *lowlimit, short *highlimit) {
	allocate(n1, n2, lowlimit, highlimit);

}

// Functions for wendarray.  Used for w3 and w5
void wendarray::allocate(short n1, short n2, short *lowlimit, short *highlimit) {
  N1 = n1;
  N2 = n2;
  Lowlimit = lowlimit;
  


  short Low,High;

  array = new integersize *[n1+2];
  for (int i = 0; i <= n1 + 1; ++i) {
    
    
	//Changed by DHM 11/11/2010
	//No longer limit the bounds with w3 and w5.
	//This is required because w3 and w5 don't require the ends to align.
	//if (i<=N1) {
	//	Low=lowlimit[i];
	//	High=highlimit[i];
	//}
	
	//else {
		Low=0;
		High=N2+2;
	//}
	array[i] = new integersize[High-Low+2/*2*m+2*/];
	for (int j = 0; j <= High-Low; ++j) {
      array[i][j]=DYNALIGN_INFINITY;
    }

	//No shift is needed now.
	//array[i]=array[i]-Low;//shift pointer for access
  
  }

  
}

wendarray::~wendarray() {
  int i;
	
  for (i=0;i<=N1+1;++i) {
	  //restore pointer:
	  //if(i<=N1) array[i]=array[i]+Lowlimit[i];
	  //else array[i]=array[i];
	  delete[] array[i];
  }
  delete[] array;
  

}

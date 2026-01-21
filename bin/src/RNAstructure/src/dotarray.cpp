/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 */

#include "dotarray.h"

#include "defines.h"

dotarray::dotarray(int size) {
	short int i,j;

	//initialize the array
  array = new integersize *[size+1];

  for (i=0;i<=(size);i++)  {
   	array[i] = new integersize [i+1];
  }

  for (i=0;i<=size;i++) {
    for (j=0;j<=i;j++) {
      array[i][j] = DYNALIGN_INFINITY;
    }
  }

  store = size;
}

dotarray::~dotarray() {
 	short int i;

  for (i=0;i<=store;i++) {
    delete[] array[i];
  }
  delete[] array;
}

integersize &dotarray::dot(int i, int j) {
  return array[j][i];
}

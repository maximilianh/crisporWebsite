/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 */

#ifndef DOTARRAY_H
#define DOTARRAY_H

#include "defines.h"

// dotarray encapsulates the array needed to store dot plot information

class dotarray {
private:
  integersize **array;
  short int store;

public:
  dotarray(int size);
  ~dotarray();

  integersize &dot(int i, int j);
};

#endif

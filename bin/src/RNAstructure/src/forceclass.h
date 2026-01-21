
#ifndef FORCECLASS_H
#define FORCECLASS_H

#include "structure.h"

// forceclass encapsulates a large 2-d arrays used by the dynamic
// algorithm to enforce folding constraints
class forceclass {
private:
  int Size;

public:
  int k;
  char **dg;
  
  // the constructor allocates the space needed by the arrays
  forceclass(int size);

  // the destructor deallocates the space used
  ~forceclass();
  
  // f is an integer function that references the correct element of
  // the array
  char &f(int i, int j);
};

inline char &forceclass::f(int i, int j) 
{
  if (i > j) std::swap(i,j);
  if (i > Size) 
  {
    i -= Size;
    j -= Size;
  }
  return dg[i][j-i];
}

//! force is used to prepare arrays for the function dynamic, used during the
//!	       fill routines - it coordinates the force...() functions above
void force(structure *ct,forceclass *fce, bool *lfce);
//the force... functions are used to initialize the arrays used to apply
//constraints during the folding process
void forcepair(int x,int y,structure *ct,forceclass *v);
void forcedomain(int x,int y,structure *ct,forceclass *v);
void forcesingle(int x,structure* ct,forceclass *v);
void forcedbl(int dbl,structure* ct,forceclass *w,bool *v);
void forceinter(int dbl,structure* ct,forceclass *w);
void forceinterefn(int dbl,structure* ct,forceclass *w);

//returns true if pair i and j is not a GU pair and the most adjacent pairs are not GU
//	used with chemical modification data to make sure that some contributions are not over-counted
bool notgu(int i, int j, structure *ct);

#endif

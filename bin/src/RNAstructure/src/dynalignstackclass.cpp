/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 */

#include "dynalignstackclass.h"

#include "defines.h"

void dynalignstackclass::allocate_stack() {
  stackenergy = new integersize[max];
  stack=new short int *[max];
  for (short i = 0; i < max; i++) {
    stack[i] = new short int[4];
  }
  openness = new bool [max];
#ifdef DYNALIGN_II
    vmodness = new bool [max];
#else
#endif
}

dynalignstackclass::dynalignstackclass(short int stacksize) {
  max = stacksize;
  size = 0;
  allocate_stack();
}
#ifdef DYNALIGN_II
bool dynalignstackclass::pull(short int *i,short int *j, short int *a, short int *b, 
                              integersize *energy, bool *open,bool *if_vmod)
#else
bool dynalignstackclass::pull(short int *i,short int *j, short int *a, short int *b, 
                              integersize *energy, bool *open)
#endif
 {
  if (size==0) {
    return false;
  }

  size--;
  *i = stack[size][0];
  *j = stack[size][1];
  *a = stack[size][2];
  *energy = stackenergy[size];
  *b = stack[size][3];
  *open = openness[size];
#ifdef DYNALIGN_II
    *if_vmod = vmodness[size];
#else
#endif
  return true;
}
#ifdef DYNALIGN_II
void dynalignstackclass::push(short int i,short int j, short int k, short int l, 
                              integersize energy,bool open, bool if_vmod)
#else
void dynalignstackclass::push(short int i,short int j, short int k, short int l, 
                              integersize energy, bool open)
#endif
{
  short index;

  if (size == max) {
    //allocate more space:
    dynalignstackclass *temp;
    temp = new dynalignstackclass(max);
    for (index=0;index<max;index++) {
#ifdef DYNALIGN_II
         temp->push(stack[index][0],stack[index][1],stack[index][2],stack[index][3],stackenergy[index],vmodness[size],openness[index]);
#else
      temp->push(stack[index][0],stack[index][1],stack[index][2],stack[index][3],stackenergy[index],openness[index]);
#endif
    }
    delete_array();
    max = 2*max;

    allocate_stack();
    for (index=0;index<(max/2);index++) {
#ifdef DYNALIGN_II
        temp->pull(&stack[index][0],&stack[index][1],&stack[index][2],&stack[index][3],&stackenergy[index],&vmodness[size],&openness[index]);
#else
      temp->pull(&stack[index][0],&stack[index][1],&stack[index][2],&stack[index][3],&stackenergy[index],&openness[index]);
#endif
    }
  }
  stack[size][0] = i;
  stack[size][1] = j;
  stack[size][2] = k;
  stackenergy[size] = energy;
  stack[size][3] = l;
#ifdef DYNALIGN_II
    vmodness[size] = if_vmod;
#else
#endif
  openness[size] = open;
  size++;
}
	
void dynalignstackclass::delete_array() {
  for (short i = 0; i < max; i++) {
    delete[] stack[i];
  }
  delete[] stack;

  delete[] stackenergy;
  delete[] openness;
#ifdef DYNALIGN_II
  delete[] vmodness;
#else
#endif
}

dynalignstackclass::~dynalignstackclass() {
  delete_array();
}

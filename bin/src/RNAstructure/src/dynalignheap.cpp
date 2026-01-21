/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 */

#include "dynalignheap.h"

#include "defines.h"

dynalignheap::dynalignheap(int allocate) {
  size=0;
  heapi = new short [allocate];
  heapj = new short [allocate];
  heapa = new short [allocate];
  heapb = new short [allocate];
  max=allocate;
  heape = new integersize [allocate];
}

dynalignheap::~dynalignheap() {
  delete[] heapi;
  delete[] heapj;
  delete[] heapa;
  delete[] heapb;
  delete[] heape;
}

void dynalignheap::push(short i, short j, short a, short b, integersize e) {
  dynalignheap *temp;
  int index;

  if (size==max) {
    temp=new dynalignheap(max);
    for (index=0;index<max;index++) temp->push(heapi[index],heapj[index],heapa[index],heapb[index],heape[index]);
    delete[] heapi;
    delete[] heapj;
    delete[] heapa;
    delete[] heapb;
    delete[] heape;

    max=10*max;
    heapi = new short [max];
    heapj = new short [max];
    heapa = new short [max];
    heapb = new short [max];
    heape = new integersize [max];

    for (index=0;index<size;index++) temp->read(index,&(heapi[index]),&(heapj[index]),&(heapa[index]),&(heapb[index]),
                                                &(heape[index]));
    delete temp;
  }
  heapi[size]=i;
  heapj[size]=j;
  heapa[size]=a;
  heapb[size]=b;
  heape[size]=e;
  size++;
}
void dynalignheap::read(int position, short *i,short *j, short *a, short *b, integersize *e) {
  *i = heapi[position];
  *j = heapj[position];
  *a = heapa[position];
  *b = heapb[position];
  *e = heape[position];
}

integersize dynalignheap::peak(int position) {
  return heape[position];
}

void dynalignheap::swap(int p1, int p2) {
  short i,j,a,b,e;
  i = heapi[p1];
  j = heapj[p1];
  a = heapa[p1];
  b = heapb[p1];
  e = heape[p1];

  heapi[p1] = heapi[p2];
  heapj[p1] = heapj[p2];
  heapa[p1] = heapa[p2];
  heapb[p1] = heapb[p2];
  heape[p1] = heape[p2];

  heapi[p2] = i;
  heapj[p2] = j;
  heapa[p2] = a;
  heapb[p2] = b;
  heape[p2] = e;
}

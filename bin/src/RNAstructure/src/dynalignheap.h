/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 */

#ifndef DYNALIGNHEAP_H
#define DYNALIGNHEAP_H

#include "defines.h"

class dynalignheap {
private:
	int max;
	short *heapi,*heapj, *heapa,*heapb;
	integersize *heape;

public:
	int size;
	dynalignheap(int allocate=10000);
	~dynalignheap();
	void push(short i, short j, short a, short b, integersize e);
	void read(int position, short *i,short *j, short *a, short *b, integersize *e);
	integersize peak(int position);
	void swap(int p1, int p2);
};

#endif

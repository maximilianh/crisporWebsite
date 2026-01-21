/*
 * Code modified from: 
 * 	Dynalign: by David Mathews, copyright 2002, 2003,2004,2005
 *	Contributors: Chris Connett and Andrew Yohn, 2006
 *	Contributors: Jason Gloor, 2007
 */

#ifndef EXPECTMAXSTACK_H
#define EXPECTMAXSTACK_H

#include "defines.h"

class expectMaxStack {
private:
	void allocate_stack();

public:
	int size, **stack, maximum;

	expectMaxStack(int stacksize = 50);
	~expectMaxStack();

	bool pull(int *i, int *j);
	void push(int i, int j);
	
	void delete_array();
};

#endif

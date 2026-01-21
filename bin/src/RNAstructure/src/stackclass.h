/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 */

#ifndef STACKCLASS_H
#define STACKCLASS_H

#include "defines.h"

class stackclass 
{
private:
	void allocate_stack();

public:
	short size, **stack, maximum;
	integersize *stackenergy;

	stackclass(short int stacksize = 50);
	~stackclass();

	bool pull(short int *i, short int *j, short int *open, integersize *energy, short int *pair);
	void push(short int i, short int j, short int open, integersize energy, short int pair);
	
	void delete_array();
};

#endif

/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 */

#ifndef DYNALIGNSTACKCLASS_H
#define DYNALIGNSTACKCLASS_H

#include "defines.h"

class dynalignstackclass {
	short **stack;
	int size, max;
	integersize *stackenergy;
	bool *openness;
        void allocate_stack();
#ifdef DYNALIGN_II
        bool *vmodness;
#else
#endif
public:
	dynalignstackclass(short int stacksize = 50);
#ifdef DYNALIGN_II
    	bool pull(short int *i,short int *j, short int *a, short int *b, 
		  integersize *energy, bool *open,bool *if_vmod);
	void push(short int i,short int j, short int a, 
		  short int b, integersize energy,bool open=false, bool if_vmod=false);
#else
	bool pull(short int *i,short int *j, short int *a, short int *b, 
            integersize *energy, bool *open);
	void push(short int i,short int j, short int a, 
            short int b, integersize energy, bool open=false);
#endif	
	void delete_array();
	~dynalignstackclass();
};

#endif

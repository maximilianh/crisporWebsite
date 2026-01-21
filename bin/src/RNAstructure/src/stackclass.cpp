/*
 *  Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 *  Contributors: Chris Connett and Andrew Yohn, 2006
 */

#include "stackclass.h"
#include "defines.h"

void stackclass::allocate_stack() 
{
	stackenergy = new integersize [maximum];
	stack = new short int *[maximum];
	for (short i = 0; i < maximum; i++)
	{
		stack[i] = new short int[4];
	}
}

stackclass::stackclass(short int stacksize) 
{
	maximum = stacksize;
	size = 0;
	allocate_stack();
}

bool stackclass::pull(	short int *i,
						short int *j, 
						short int *open, 
						integersize *energy, 
						short int *pair) 
{		
	if (size == 0) return false;
	else 
	{
		size--;
		*i = stack[size][0];
		*j = stack[size][1];
		*open = stack[size][2];
		*pair = stack[size][3];
		*energy = stackenergy[size];
		return true;
	}
}
	
void stackclass::push(	short int i, 
						short int j, 
						short int open, 
						integersize energy, 
						short int pair	)
{
	
	if (size == maximum) 
	{
		//allocate more space:
		stackclass *temp;
		temp = new stackclass(maximum);
		
		for (short k = 0; k < maximum; k++) 
		{
			temp->push(stack[k][0], stack[k][1], stack[k][2], stackenergy[k], stack[k][3]);
		}

		delete_array();
		
		maximum = 2 * maximum;

		allocate_stack();
		
		for (short k = 0; k < (maximum/2); k++) 
		{
			temp->pull(&(stack[k][0]), &(stack[k][1]), &(stack[k][2]), &(stackenergy[k]), &(stack[k][3]));
		}

		delete temp;
	}
		
	stack[size][0] = i;
	stack[size][1] = j;
	stack[size][2] = open;
	stack[size][3] = pair;
	stackenergy[size] = energy;
	size++;
}
	
void stackclass::delete_array() 
{
	for (short i = 0; i < maximum; i++) 
	{
		delete[] stack[i];
	}
	delete[] stack;
	delete[] stackenergy;
}

stackclass::~stackclass() 
{
	delete_array();
}

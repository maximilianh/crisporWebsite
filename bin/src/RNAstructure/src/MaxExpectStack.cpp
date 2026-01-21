/*
 * This class was derived from:
 * 	Dynalign:  by David Mathews, copyright 2002, 2003,2004,2005
 *	Contributors: Chris Connett and Andrew Yohn, 2006
 *	Contributors: Jason Gloor, 2007
 *
 *	The stack takes i and j integers
 */

#include "MaxExpectStack.h"

#include "defines.h"

void expectMaxStack::allocate_stack() {
	int i;
  
	stack = new int *[maximum];
	for (i=0;i<maximum;i++) stack[i] = new int [2];
}

expectMaxStack::expectMaxStack(int stacksize) {
	maximum = stacksize;
	size = 0;
	allocate_stack();
}

bool expectMaxStack::pull(int *i,int *j) {
		
	if (size==0) return false;
	else {
		size--;
		*i = stack[size][0];
		*j = stack[size][1];
		return true;
	}
}
	
void expectMaxStack::push(int i,int j){
	int k;

	if (size == maximum) {
		//allocate more space:
		expectMaxStack *temp;
		temp = new expectMaxStack(maximum);
		for (k=0;k<maximum;k++) {
			temp->push(stack[k][0],stack[k][1]);
		}
		delete_array();
		maximum = 2*maximum;

		allocate_stack();
		for (k=0;k<(maximum/2);k++) {
			temp->pull(&(stack[k][0]),&(stack[k][1]));
		}

		delete temp;
	}
		
	stack[size][0] = i;
	stack[size][1] = j;
	size++;
}

void expectMaxStack::delete_array() {
	for (int i = 0; i < maximum; i++) {
 		delete[] stack[i];
  	}
	delete[] stack;
}

expectMaxStack::~expectMaxStack() {
	delete_array();
}

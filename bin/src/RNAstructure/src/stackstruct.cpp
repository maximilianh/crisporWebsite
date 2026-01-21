

#include "stackstruct.h"

void push(stackstruct *stack,int a,int b,int c,int d)
{
  (stack->sp)++;
  stack->stk[stack->sp][0]= a;
  stack->stk[stack->sp][1]= b;
  stack->stk[stack->sp][2]= c;
  stack->stk[stack->sp][3]= d;
}

void pull(stackstruct *stack,int *i,int *j,int *open,int *null,int *stz)
{
	if (stack->sp==0) {
		*stz = 1;
		return;
	}
	else {
    *stz = 0;
    *i = stack->stk[stack->sp][0];
    *j = stack->stk[stack->sp][1];
    *open= stack->stk[stack->sp][2];
    *null= stack->stk[stack->sp][3];
    stack->sp--;
	}
}

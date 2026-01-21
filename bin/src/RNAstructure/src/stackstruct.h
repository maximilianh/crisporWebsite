

#ifndef STACKSTRUCT_H
#define STACKSTRUCT_H

// this structure contains a stack of data, used by functions that
// analyze a structure piecewise
struct stackstruct {
	int stk[101][4],sp;
};

//push info onto the stack
void push(stackstruct *stack,int a,int b,int c,int d);

//pull info from the stack
void pull(stackstruct *stack,int *i,int *j,int *open,int *null,int *stz);

#endif

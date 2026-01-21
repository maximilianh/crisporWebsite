#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "ppf_tb_stack.h"

bool _DUMP_TB_STACK_MESSAGES_ = false;

// Constructor, do not do anything?
t_ppf_tb_stack::t_ppf_tb_stack()
{
	this->n_nodes = 0; // There is nothing on the stack.
	this->tos = new t_stack_node(0, 0, 0, 0, TRACE_stack_empty); // Special empty stack node, which is popped once to finish traceback.
}

t_ppf_tb_stack::~t_ppf_tb_stack()
{
	delete(this->tos);
}

// Allocate a node and set structure to first node.
void t_ppf_tb_stack::push_str(int i, int j, int k, int l, char array_id)
{
	// Allocate a new node.	
	t_stack_node* new_node = new t_stack_node(i, j, k, l, array_id);

	// Put new node to top of stack, tos.
	t_stack_node* last_tos = this->tos; // Copy stack.
	this->tos = new_node; // Set top of stack.
	this->tos->next = last_tos; // Set next.

	if(_DUMP_TB_STACK_MESSAGES_)
	{
		printf("Pushed %s(%d, %d, %d, %d)\n", trace_names[array_id], i,j,k,l);
	}
}

// Get first node and delete it from stack.
// Return indices set to -1 on empty stack.
void t_ppf_tb_stack::pop_str(int& i, int& j, int& k, int& l, char& array_id)
{
	// Copy tos to arguments.
	i = this->tos->i;
	j = this->tos->j;
	k = this->tos->k;
	l = this->tos->l;
	array_id = this->tos->array_id;

	// Now have to delete tos and set tos to next one.
	t_stack_node* to_be_deleted_tos = this->tos;
	this->tos = to_be_deleted_tos->next; // Set new tos.

	delete(to_be_deleted_tos);
}

t_stack_node::t_stack_node(int i, int j, int k, int l, char array_id)
{
	this->i = i;
	this->j = j;
	this->k = k;
	this->l = l;

	this->array_id = array_id;

	this->next = NULL; // Set next of any newly created stack node to NULL.
}

// Destruct.
t_stack_node::~t_stack_node()
{}

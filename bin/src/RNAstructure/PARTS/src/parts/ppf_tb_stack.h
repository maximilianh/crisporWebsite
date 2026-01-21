#ifndef _TB_STACK_
#define _TB_STACK_

// Forwrd declare stack node.
class t_stack_node; 

// Array identifiers for traceback.
enum
{
	TRACE_SS,
	TRACE_WL,
	TRACE_W,
	TRACE_WMB,
	TRACE_WMBL,
	TRACE_V_mhe,
	TRACE_W_ij_mhi,
	TRACE_W_kl_mhi,
	TRACE_W_ext,
	TRACE_stack_empty // When traceback is finished.
};

static char trace_names[29][100] = 
{
	"SS",
	"WL",
	"W",
	"WMB",
	"WMBL",
	"V_mhe",
	"W_ij_mhi",
	"W_kl_mhi",
	"W_ext",
	"EMPTY_STACK"
};

// t_ppf_tb_stack: 
// Optimal traceback stack for handling of 
class t_ppf_tb_stack
{
public:
	int n_nodes;

	t_stack_node* tos;

	t_ppf_tb_stack();
	~t_ppf_tb_stack();

	// Allocate a node and set structure to first node.
	void push_str(int i, int j, int k, int l, char array_id);

	// Get first node and delete it from stack.
	void pop_str(int& i, int& j, int& k, int& l, char& array_id);
};

// A stack node.
class t_stack_node
{
public:
	int i;
	int j;
	int k;
	int l;
	char array_id;
	t_stack_node* next;

	t_stack_node(int i, int j, int k, int l, char array_id);
	~t_stack_node();
};

#endif // _TB_STACK_

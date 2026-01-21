#ifndef _ARRAY_FILE_MANAGER_
#define _ARRAY_FILE_MANAGER_

#include <vector> 

using namespace std;

enum{V_FILE, WL_FILE, W_FILE, WMBL_FILE, WMB_FILE, VMHE_FILE, N_ARRAYS_LOGGED};
static char array_file_names[N_ARRAYS_LOGGED][20] = {"V", "WL", "W", "WMBL", "WMB", "VMHE"};

#define ARRAY_IP_DIR "array_dumps/in"
#define ARRAY_OP_DIR "array_dumps/out"

class t_ppf_cli;
class t_seq_man;

class t_array_file_manager
{
public:
	t_array_file_manager(t_seq_man* _seq_man);
	~t_array_file_manager();

	t_seq_man* seq_man;

	vector<FILE*>* in_array_files;
	vector<FILE*>* out_array_files;

	bool write_array_value(int out_array, double val);
	bool read_array_value(int in_array, double& val);
};

#endif // _ARRAY_FILE_MANAGER_




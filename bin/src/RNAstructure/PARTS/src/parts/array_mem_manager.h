#ifndef _ARRAY_MEM_MANAGER_
#define _ARRAY_MEM_MANAGER_

#define N_ARRAYS (9)

#define LOG_THRESHOLD_UPDATE (0.001)

#define MEG_BYTES (1024 * 1024)
#define GIG_BYTES (1024 * 1024 * 1024)

class t_seq_man;
class t_ppf_cli;
struct t_aln_env_result;

class t_array_mem_manager
{
public:
	t_array_mem_manager(t_seq_man* _seq_man, t_ppf_cli* _ppf_cli);
	~t_array_mem_manager();

	t_seq_man* seq_man;
	t_ppf_cli* ppf_cli;

	t_aln_env_result* get_env_threshold_per_mem_limit(double array_mem_limit);
	double get_mem_per_array();
	double get_total_mem_size();
};

#endif // _ARRAY_MEM_MANAGER_

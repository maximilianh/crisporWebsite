#ifndef _MAP_RESULTS_
#define _MAP_RESULTS_

class t_structure;
class t_MAP_structures;
class t_MAP_alignment;
class t_MAP_mhr_info;
class t_ppf_cli;
class t_seq_man;

// This structure represents a joint of all map results.
class t_map_results
{
public:
	// Construct the results from map structures and map alignment.
	t_map_results(t_seq_man* _seq_man, 
			 t_ppf_cli* _ppf_cli, 
			 t_MAP_structures* _map_strs, 
			 t_MAP_alignment* _map_alignment, 
			 t_MAP_mhr_info* _map_mhr_info);

	t_MAP_structures* map_strs;
	t_MAP_alignment* map_alignment;
	t_MAP_mhr_info* map_mhr_info;
	t_ppf_cli* ppf_cli;
	t_seq_man* seq_man;

	~t_map_results();

	// Mutate map_str's into two structure classes.
	t_structure* str1;
	t_structure* str2;
	void mutate_map_strs();
};

#endif // _MAP_RESULTS_


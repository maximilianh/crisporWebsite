#ifndef _SEED_MANAGER_
#define _SEED_MANAGER_

#include <vector>
using namespace std;

// Global seed manager for RNG's. 

class t_seed_manager
{
public:
	t_seed_manager();
	t_seed_manager(char* seeds_fp);
	~t_seed_manager();

	static vector<int>* seeds;
	static int i_next_seed;

	// Generte a unique seed.
	static int seed_me();

	static void dump_seeds();
};

#endif // _SEED_MANAGER_


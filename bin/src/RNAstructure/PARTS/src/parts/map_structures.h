// This header file encapsulates various decoded structures, like map,
// map and wrt other criterion.
#ifndef _MAP_STRUCTURES_
#define _MAP_STRUCTURES_

class t_seq_man;
class t_ppf_cli;

class t_MAP_structures
{
public:
	int N1;
	int N2;
	int* seq1_map_ct_bps;
	int* seq2_map_ct_bps;
	t_seq_man* seq_man;
	t_ppf_cli* ppf_cli;

	// Constructors and destructors.
	t_MAP_structures(t_seq_man* seq_man, t_ppf_cli* ppf_cli);
	~t_MAP_structures();

	// Add base pairs to structures while doing MAP decoding.
	void add_bp_to_MAP_ct1(int i, int j);
	void add_bp_to_MAP_ct2(int k, int l);

	// Dump MAP decoded ct's in ct format.
	void dump_map_cts();
};

#endif // _MAP_STRUCTURES_

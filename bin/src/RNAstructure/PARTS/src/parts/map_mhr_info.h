#ifndef _MHR_INFO_
#define _MHR_INFO_

enum{FP_SIDE, TP_SIDE};

class t_MAP_structures;
class t_MAP_alignment;
class t_ppf_cli;

#define GAP_CHAR ('.')

// This class represents the MHR's that are inside the MAP structural alignment.
// This information can be traced from MAP structures and MAP alignment, after PARTS MAP traceback.
class t_MAP_mhr_info
{
public:
	t_MAP_structures* map_structures;
	t_MAP_alignment* map_alignment;
	t_ppf_cli* ppf_cli;
	int n_mhrs; // Number of MHR's in the structural alignment.
	int* mhrs1; // mhrs[i] indicates which mhr ith nucleotide in 1st sequence is in.
	int* mhrs2; // mhrs[k] indicates which mhr kth nucleotide in 2nd sequence is in.
	int* mhr1_sides;
	int* mhr2_sides;

	t_MAP_mhr_info(t_ppf_cli* ppf_cli, t_MAP_structures* map_structures, t_MAP_alignment* map_alignment);
	~t_MAP_mhr_info();



	// Trace and determine mhr of each nucleotide in both sequences.
	void trace_MAP_mhrs();
	void _trace_MAP_mhrs();
	void __trace_MAP_mhrs();
};

#endif // _MHR_INFO_



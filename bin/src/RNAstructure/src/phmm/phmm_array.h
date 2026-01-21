#ifndef _PHMM_ARRAY_
#define _PHMM_ARRAY_

/*
t_phmm_array is a 3D array that contains state and index information, for storing 
forward/backward/ML array computation results.
*/

class t_phmm_array
{
public:
	t_phmm_array(int _n1, int _n2, int _phmm_band_constraint_size, bool mallocate);
	~t_phmm_array();

	int n1;
	int n2;

	double n_bytes_alloced;

	static int low_phmm_limit(int i, int n1, int n2, int phmm_band_constraint_size);
	static int high_phmm_limit(int i, int n1, int n2, int phmm_band_constraint_size);
	static bool check_phmm_boundary(int i, int k, int n1, int n2, int phmm_band_constraint_size);
	bool check_phmm_boundary(int i, int k);

	void set_hmm_array_banded_limits();

	int* low_phmm_array_limits;
	int* high_phmm_array_limits;

	int phmm_band_constraint_size;

	// Indexed by i1, k, i_state
	double*** array;

	double& x(int i, int k, int state);
};

#endif // _PHMM_ARRAY_



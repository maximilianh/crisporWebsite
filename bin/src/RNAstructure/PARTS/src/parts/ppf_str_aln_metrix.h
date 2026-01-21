#ifndef _PPF_STR_ALN_METRIX_
#define _PPF_STR_ALN_METRIX_

class t_ppf_str_aln_metrix
{
public:
	int n_samples;

	t_ppf_str_aln_metrix(int _n_samples);
	~t_ppf_str_aln_metrix();

	// bds: Base Pair Difference: Summate number of base pair
	// between structures of 1st sequence and 2nd sequence.
	double get_bds_distance(int i1, int i2, int* bps_i1_1, int* bps_i2_1, int* bps_i1_2, int* bps_i2_2);

	void dump_distances();

};


#endif // _PPF_STR_ALN_METRIX_
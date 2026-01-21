#ifndef _PPF_SCALE_
#define _PPF_SCALE_

// PPF scaling handler for linear ppf calculation.

class t_ppf_V;
class t_ppf_WL;
class t_ppf_WR;
class t_ppf_W;
class t_ppf_WMBL;
class t_ppf_WMBR;
class t_ppf_WMB;
class t_ppf_W_mhi;
class t_ppf_V_mhe;

class t_spf_array; // For determining scaling factors for each nucleotide.
class t_seq_man; // Sequence manager to interface with sequence data.
class t_frag_aln_enum_array;
class t_ppf_cli;
class t_ppf_loops;

// http://www.cprogramming.com/tutorial/floating_point/understanding_floating_point_representation.html
#define MAX_ARRAY_VALUE_LIMIT (pow(10.0, 280)) // This max limit is the limit of rescaling, when values reach this limit, a rescaling is needed to avoid overflow.
#define MIN_ARRAY_VALUE_LIMIT (pow(10.0, -280)) // This min limit is the limit of rescaling, when values reach this limit, a rescaling is needed to avoid overflow.
#define DEFAULT_PER_NUCLEOTIDE_RESCALE_INCREMENT (1.2) // Rescaling factor of each nucleotide. 

class t_ppf_scale
{
	// Those should be inaccessible to other classes.
public:
	t_seq_man* seq_man;
	t_ppf_cli* ppf_cli;
	t_ppf_loops* ppf_loops;

	//double rescale_factor; // Initial rescale factor.
	double rescaling_increment_factor_per_nucleotide;

	double* seq1_scales; // Scaling factors of each nucleotide in 1st sequence.
	double* seq2_scales; // Scaling factors of each nucleotide in 2nd sequence.

	t_ppf_scale(t_ppf_loops* ppf_loops);
	t_ppf_scale(t_ppf_loops* ppf_loops, double init_rescale_factor);

	// Calculate scaling factor for each nucleotide using spf interface.
	void scaling_factor_by_nuc();

	// Rescale check function, calls each array to rescale if going out of precision for
	// one of the arrays.
	bool check_pp_ppf_array_rescale(int i, int j, int k, int l);
	bool check_pp_ppf_external_array_rescale(int i, int j, int k, int l);

	bool check_back_frag_array_rescale(int i, int j, int k, int l,
								t_frag_aln_enum_array* back_frag_array);

	// Read scaler values by nucleotide index and 
	// i1 is nucleotide index of sequence 1 and i2 is ...
	double get_scale_value(int i1, int i2); 

	double get_partial_scale_factor(int s1, int e1, int s2, int e2);
	double get_log_partial_scale_factor(int s1, int e1, int s2, int e2);

	double cumulative_rescale_external_factor(int i, int j, int k, int l);
	
	// Following function returns log value of an array value, which is calculated in 
	// linear space. Since calculated in linear space, array value contains a partial scaling factor.
	// This function is good for dumping values of arrays.
	double get_unscaled_log_value_by_indices(int i, int j, int k, int l, double scaled_value);

	void copy_scales_to_spf_by_id(int id, double* spf_scales);

	double cumulative_rescale_factor(int i, int j, int k, int l);
	double cumulative_rescale_factor_exc_ik(int i, int j, int k, int l);
	void increment_scales_by_rescale(bool up_scale);

	~t_ppf_scale();
};

#endif // _PPF_SCALE_

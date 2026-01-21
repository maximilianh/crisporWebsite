#ifndef _FOLDING_CONSTRAINTS_
#define _FOLDING_CONSTRAINTS_

#define POS_MEM_NOT_EXIST (SHRT_MAX)

class t_structure;

/*
Manipulate and interface with the structural constraints.
*/
class t_folding_constraints
{
public:
	t_folding_constraints(t_structure* _rna_seq);
	t_folding_constraints(char* _seq_fp);
	t_folding_constraints(t_structure* _rna_seq, double** linear_pp, double pp_threshold);
	t_folding_constraints(char* _seq_fp, double** linear_pp, double pp_threshold);
	t_folding_constraints(t_folding_constraints* _folding_constraint);
	~t_folding_constraints();

	static int pairable[5][5];

	void alloc_init_maps();
	void free_maps();

	void mallocate_ptr_reloc_maps();
	void free_ptr_reloc_maps();

	//t_energy_loops* energy_loops;
	t_structure* rna_seq;

	bool** same_loop_map;
	bool** str_coinc_map;
	bool** pairing_map;
	bool* forbid_non_v;

	// Pointer relocation maps, derived from pairing probabilities.
	//static short** compute_paired_ptr_reloc_map(t_structure* rna_seq, double pp_threshold, double** pairing_probs);
	//static short** compute_coinc_ptr_reloc_map(t_structure* rna_seq, double pp_threshold, double** pairing_probs);
	//static short** compute_paired_ptr_reloc_map(char* seq_fp, double pp_threshold, double** pairing_probs);
	//static short** compute_coinc_ptr_reloc_map(char* seq_fp, double pp_threshold, double** pairing_probs);
	short** coinc_pointer_relocation_map; // Use to allocate unpaired (only co-incident locations)
	short** paired_pointer_relocation_map; // Use to allocate paired arrays.

	// Force pairing of two nucleotides.
	bool force_pairing(int i, int j);

	bool forbid_non_v_emission(int i);
	bool forbid_non_v_emission(int i, int j);

	// Check the closability of following loops.
	bool check_internal_loop(int i, int j, int inner_i, int inner_j);
	bool check_hairpin_loop(int i, int j);

	void dump_constraints();

	// Compute the saturated base pair for the pairing probabilities.
	void compute_saturated_structure(double** pairing_probs);
	int* saturated_str_bps;

private:
	void compute_ptr_reloc_maps(double pp_threshold, double** pairing_probs);
};

#endif // _FOLDING_CONSTRAINTS_


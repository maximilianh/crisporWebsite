#ifndef _MAP_ALIGNMENT_
#define _MAP_ALIGNMENT_

// Handles alignment detemrination for map alignment.

class t_seq_man;
class t_ppf_cli;

class t_MAP_alignment
{
public:
	t_seq_man* seq_man;
	t_ppf_cli* ppf_cli;
	int** seq1_alns; // an array which holds alignment info for sequence 1
	int** seq2_alns; // an array which holds alignment info for sequence 2

	t_MAP_alignment(t_seq_man* seq_man, t_ppf_cli* ppf_cli);
	~t_MAP_alignment();

	// Following function are for setting alignment values.
	void set_seq1_ins(int seq1_index, int seq2_index);
	void set_seq2_ins(int seq1_index, int seq2_index);
	void set_aln(int seq1_index, int seq2_index);

	// Get the length of the sequence alignment.
	int get_l_aln();

	// Dump map alignment.
	void dump_map_alignment();

	// Alignment strings.
	char* aln_str1;
	char* aln_str2;

	int l_aln;

	int* aln_index_line1;
	int* aln_index_line2;
};

#endif // _MAP_ALIGNMENT_


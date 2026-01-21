#ifndef _PFE_SAMPLED_ALIGNMENT_
#define _PFE_SAMPLED_ALIGNMENT_

#include <vector>

using namespace std;

// Handles alignment element of one sampled structural alignment.
class t_seq_man;
class t_ppf_cli;

struct t_aln_site
{
	int i1;
	int i2;
};

class t_stoch_sampled_alignment
{
public:
	t_seq_man* seq_man;
	t_ppf_cli* ppf_cli;

	vector<t_aln_site*>* alignment_path;

	t_stoch_sampled_alignment(t_seq_man* seq_man, t_ppf_cli* ppf_cli);
	~t_stoch_sampled_alignment();

	// Following function are for setting alignment values.
	void set_seq1_ins(int seq1_index, int seq2_index);
	void set_seq2_ins(int seq1_index, int seq2_index);
	void set_aln(int seq1_index, int seq2_index);

	// Dump map alignment.
	void dump_sampled_alignment(int _sample_id);

	void reset_aln();

	// Alignment strings.
	char* aln_str1;
	char* aln_str2;	
};

bool cmp_aln_site(t_aln_site* aln_site1, t_aln_site* aln_site2);

#endif // _MAP_ALIGNMENT_




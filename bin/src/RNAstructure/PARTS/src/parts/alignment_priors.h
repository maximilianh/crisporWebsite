#ifndef _ALIGNMENT_PRIORS_
#define _ALIGNMENT_PRIORS_

class t_seq_man;
class t_ppf_cli;

// This class manages alignment priors and what to return as a prior (i.e. LOG scores of alignment position states.)
// Basically encapsulates a 3 dimensional matrix just as forward-backward array.
// It is to be used extensively in PPF and also PHMM. If serialization is implemented then
// it might be helpful 
class t_aln_priors
{
public:
	int N1;
	int N2;
	double*** priors; // Accession is the same as forward-backward arrays: [seq1 index][seq2 index][state].

	t_ppf_cli* ppf_cli;

	double n_bytes_alloced;

	// Constructor can allocate and initialize prior array.
	// Allocate empty, UNIFORM priors for each INDIVIDUAL ALIGNMENT: Contribution ofd any alignment is same as any other one without depending on length of alignment.
	// Nonuniformity (i.e. length dependency) comes when priors are set from hmm model.
	t_aln_priors(t_seq_man* seq_man, bool mallocate); 
	t_aln_priors(t_seq_man* seq_man, t_ppf_cli* ppf_cli, bool mallocate);
	t_aln_priors(t_seq_man* seq_man, t_ppf_cli* ppf_cli, double** aln_probs, double** ins1_probs, double** ins2_probs, bool mallocate);

	~t_aln_priors();

	// Array allocator used by both ctors.
	void alloc_array(bool mallocate);
	void free_array();

	double x(int i1, int i2, int state);
};

#endif // _ALIGNMENT_PRIORS_

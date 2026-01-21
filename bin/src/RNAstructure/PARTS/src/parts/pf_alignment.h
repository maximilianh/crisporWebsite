/* 
t_pf_alignment class: encapsulates sequence alignment inference from structural alignment model.
This class infers the sequence alignment information that will be prior to actual pairwise HMM model 
from structrual alignmet model, namely the pairwise partition function. 

The sequence alignment probabilities (ALN, INS1, INS2) are encapsulated in 3 dimensional arrays as implemented
in pairwise HMM in order to satisfy COMPAREatibility.

Inference is acCOMPARElished by changing inside-outside nature of pairwise partition function into a forward-backward
calculation type of a relation. However it seems impossible to infer exact forward and backward arrays from 
partition function since there is long reange relations in pairwise partition function calculation which
breaks markovian assumptions in forward-backward algorithm.

Instead of exact forward and backward arrays, sequence alignment probabilities for 3 different states are calculated,
more precisely, the enumeration scores for 3 states are calculated. So that at the end I should have sth like this:

state_probs_array[cnt1][cnt2][ALN] refers to SUM of scores of all possible structures which include cnt1 aligned to cnt2.
state_probs_array[cnt1][cnt2][INS1] refers to SUM of scores of all possible structures where there is an INS1 at cnt1, cnt2
alignment position. 
INS2 is similar.

assuming everything works ok until now, one last problem is inferring normalized probabilities from above scores. It 
should be noted that we can get exact likelihood ratios from above scores as in turbo decoding however meanings of 
likelihoods are not clear to me now.
*/
#ifndef _SEQ_ALN_ARRAY_
#define _SEQ_ALN_ARRAY_

// Following are for different states to access sequence alignment array.
#define STATE_ALN (0)
#define STATE_INS1 (1)
#define STATE_INS2 (2)

class t_seq_man; // Sequence manager to interface with sequence data.

// Need following two classes in order to calculate sequence alignment probabilities.
class t_ppf_W; //
class t_ppf_WMB; //

class t_pf_alignment
{
public:
	int N1;
	int N2;

	t_pf_alignment(t_seq_man* seq_man); // Constructor, allocates and inits probability array.
	~t_pf_alignment(); // Destructor.

	// Accession to proabilities array is as following:
	// 1st index is seq1 index, 2nd index is seq2 index, 3rd index is state index.
	double*** state_probs_array; // Actual array which holds sequence alignment probabilities.

	// This function calls following 3 function to infer 
	void infer_seq_aln_probs(t_ppf_W* W, t_ppf_WMB* WMB); // Just two of these are enough.

	// Following 3 functions calculate different state probabilities at alignment positions.
	void infer_seq_ALN_probs(t_ppf_W* W, t_ppf_WMB* WMB); // 
	void infer_seq_INS1_probs(t_ppf_W* W, t_ppf_WMB* WMB); // 
	void infer_seq_INS2_probs(t_ppf_W* W, t_ppf_WMB* WMB); // 
};

#endif // _SEQ_ALN_ARRAY_

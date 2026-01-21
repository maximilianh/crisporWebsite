#ifndef _P_ALIGNMENT_
#define _P_ALIGNMENT_

class t_p_alignment
{
public:
	t_p_alignment(char* _seq1_aln_line, char* _seq2_aln_line);
	~t_p_alignment();

	char* seq1_aln_line;
	char* seq2_aln_line;

	double get_aln_similarity(char gap_symbol);
};

#endif // _P_ALIGNMENT_


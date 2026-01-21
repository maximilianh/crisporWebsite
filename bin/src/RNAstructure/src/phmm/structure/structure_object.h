#ifndef _STRUCTURE_
#define _STRUCTURE_

#include <vector>
using namespace std;

#define MIN_LOOP (3)
#define MAX_INT_LOOP_SIZE (30)

#define MAX_HEADER_LENGTH (1000)

// Define the pairability array.

/*
Extended implementation of t_structure class to read stacking/danglings from the file.
*/
class t_structure
{
public:
    int numofbases;
	int* numseq;
	char* nucs;
	int* basepr;
	char* ctlabel;
	//char* fp;

	bool* unpaired_forced;

	int* stackings_on_branch;
	int* danglings_on_branch; 

	int* stackings_on_mb_closure;
	int* danglings_on_mb_closure; 

	t_structure(const char* fp); // Can take a seq or a ct file.
	t_structure(t_structure* structure); // Copy constructor.
	t_structure(const char* ct_label, vector<char>* nucs, bool fix_seq_chars=true);
	t_structure(); // Defaults constructor.
	~t_structure();

	// Compare two structures and return true if they are equal.
	static bool cmp_seq(t_structure* str1, t_structure* str2);
	static vector<t_structure*>* read_multi_seq(const char* const multi_seq_fp, bool fix_seq_chars=true);

	void openct(const char* ct_fp);
	void openseq(const char* seq_fp);
	void openfasta(const char* seq_fp);

	void check_set_label();

	bool verify_ct(const char* ct_fp);
	bool verify_seq(const char* seq_fp);

	static void map_nuc_IUPAC_code(char raw_nuc, 
									char& trans_nuc, 
									int& num, 
									bool& force_unpaired);
};

#endif // _STRUCTURE_



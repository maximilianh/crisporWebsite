#ifndef PROBSCAN_CLASS
#define PROBSCAN_CLASS

#include "RNA.h"
#include "../src/defines.h"
#include <vector>
#include "../src/phmm/utils/xmath/log/xlog_math.h"




//types which contain nucleotide indices which specify a structure

//! A ProbScan hairpin closed by i and j
typedef struct hp {double probability;
                                     int i;
                                     int j;} hairpin_t;

//! A ProbScan internal loop closed by i and j on the exterior and k and l in the interior, where i<k<l<j
typedef struct il {double probability;
                                     int i;
                                     int j;
                                     int k;
                                     int l;} internal_loop_t;
//! A ProbScan base pair stack
typedef struct bp {double probability;
                                     int i;
                                     int j;
                                     int k;
                                     int l;} basestack_t;

//! A ProbScan multibranch loop

//! A multibranch loop, represented by a vector pairs.
//! branches[0] must close the multibranch loop
typedef struct mb{double probability;
                std::vector<std::pair<int,int> > branches;} multibranch_loop_t;

//functions to make loops
//return a loop with probability and closing nuc indices provided in arguments
hairpin_t hairpin(double p,int i, int j);
internal_loop_t internal_loop(double p,int i, int j, int k, int l);
basestack_t basestack(double p,int i, int j, int k, int l);
//to build a multibranch loop, call multibranch_loop with the closing base pair
//then call add_branch with the branches, from 5' to 3'
multibranch_loop_t multibranch_loop(int i, int j);
void add_branch(multibranch_loop_t& mb,int k, int l);

//display loops and their probabilities to stdout, fixed = true is fixed width/fixed = false is scientific
void show_hairpins(vector<hairpin_t>,bool fixed);
void show_stacks(vector<basestack_t>, bool fixed);
void show_internal_loops(vector<internal_loop_t>, bool fixed);
void show_bulge_loops(vector<internal_loop_t>, bool fixed);
void show_mbl(multibranch_loop_t mbl,bool fixed);

//! A ProbScan multibranch loop element.

//! Represents an element of a multibranch loop
//! a pair or an unpaired nucleotide,used in multibranch loop probability calculation
//! either a helix closed by i and j or an unpaired nucleotide at i
class mb_element{
 public:
    int i;
    int j;
    bool is_a_pair;
    mb_element(std::pair<int,int> h) : i(h.first),j(h.second),is_a_pair(true) {}
    mb_element(int nuc) : i(nuc),j(0),is_a_pair(false) {}
};

class ProbScan : public RNA
{
 public:
    //!Constructor - user provides a sequence as a c string.

    //! The partition function will be calculated. If the sequence is long, this may take some time.
    //!	Capitalization makes no difference.
    //! Note that sequences will subsequently be indexed starting at 1 (like a biologist), so that the 0th position in the sequence array will be nucleotide 1.
    //!	\param sequence is a NULL terminated c string containing the nucleotide sequence.
	//!	\param alphabetName is a cstring that indicates the folding libary, where the default, "rna", is the standard RNA parameters.
    ProbScan(std::string sequence, const char* const alphabetName="rna");
    //!Constructor - user provides a filename for existing file as a c string.

    //!	The existing file, specified by filename, can either be a a sequence file or an RNAstructure partition function save file.
    //!	Therefore, the user provides a flag for the file:from_sequence_file
    //!		true => .seq or .fasta file, false => partition function save (.pfs) file.
    //! If the input file is a sequence the partition function will be calculated. If the sequence is long, this may take some time.
    //!	This constructor generates internal error codes that can be accessed by GetErrorCode() after the constructor is called.  0 = no error.
    //! The errorcode can be resolved to a c string using GetErrorMessage.
    //!	If a sequence filename is provided, the contructor needs to be explicitly told the folding alphabet, via alphabetName, because files do not store this information.
    //! Note also that save files explicitly store the thermodynamic parameters, therefore changing the backbone type as compared to the original calculation will not change structure predictions.
    //! \param filename is null terminated c string containing the path to the input file.
    //! \param from_sequence_file is a bool which tells the constructor whether we are initializing from a sequence file, in which case the partition function must be calculated
    //!	\param alphabetName is a cstring that indicates the folding libary, where the default, "rna", is the standard RNA parameters.
    ProbScan(const char filename[],bool from_sequence_file, const char* const alphabetName="rna");

    //!	Returns probability of a hairpin closed at a specific position

    //!\param i The 5' nucleotide closing the hairpin
    //!\param j The 3' nucleotide closing the hairpin
    //!\return A double containing the probability of the hairpin
    double probability_of_hairpin(int i,int j);

	//!\param i The 5' nucleotide closing the hairpin
	//!\param j The 3' nucleotide closing the hairpin
	//!\param ip The 5'-most nucleotide closing the stem
	//!\param jp The 3'-most nucleotide closing the stem
	//!\return A double containing the probability of the hairpin stem-loop
	double probability_of_stemloop(int i, int j, int ip, int jp);
//search over all possible hairpins, return a vector of hairpin_t
//for all hairpins with probability>threshold

    //!	Calculates the probabilities of all possible hairpins in this sequence

    //!\param min The minimum size of a hairpin
    //!\param max The maximum size of a hairpin
    //!\param threshold The minimum probability for candidate hairpins
    //!\return A vector of hairpin objects, containing the positions of the hairpins and their probabilities
    std::vector<hairpin_t> probability_of_all_hairpins(int min,int max,double threshold);

    //!	Returns probability of an internal loop or bulge loop closed at a specific position

    //!\param i The 5' nucleotide closing the loop on the exterior
    //!\param j The 3' nucleotide closing the loop on the exterior
    //!\param k The 5' nucleotide closing the loop on the interior
    //!\param l The 3' nucleotide closing the loop on the interior
    //!\return A double containing the probability of the internal loop
    double probability_of_internal_loop(int i,int j,int k,int l);


    //!	Calculates the probabilities of all possible internal loops and/or bulge loops in this sequence

    //!\param threshold the minimum probability of candidate loops
    //!\param mode a string which indicates what type of loops should be searched for. Allowed values are "internal", "bulge", and "both"
    //!\return A vector of internal loop objects, containing the positions of the loops and their probabilities
    std::vector<internal_loop_t> probability_of_all_internal_loops(double threshold,
                                    std::string mode=std::string("both"));

    //! Calculates probability of a base pair stack closed at a specific position
    //! Note that this is a special case of probability_of_helix where the size is set to 1

    //!\param i The 5' nucleotide closing the stack
    //!\param j The 3' nucleotide closing the stack
    //!\return A double containing the probability of the stack
    double probability_of_stack(int i,int j);

    //!	Calculates probability of an helix at a specific position

    //!\param i The 5' nucleotide closing the helix on the exterior
    //!\param j The 3' nucleotide closing the helix on the exterior
    //!\param how_many_stacks The number of base pair STACKS in the helix (this is the number of pairs minus 1)
    //!\return A double containing the probability of the helix
    double probability_of_helix(const int i, const int j, const int how_many_stacks);


    //!	Calculates the probabilities of all possible helices in this sequence of a specific length

    //!\param threshold the minimum probability of candidate helices
    //!\param length the number of base pair stacks to search for
    //!\return A vector of helix objects, containing the positions of the helices and their probabilities
    std::vector<basestack_t> probability_of_all_helices(double threshold,int length);

    //!	Calculates probability of a multibranch loop at a specific position

    //!\param mb A multibranch loop object, containing a vector of pairs describing the multibranch loop. These can be created with the multibranch_loop function. See the text interface for the ProbScan program for an example of usage.
    //!\return A double containing the probability of the multibranch loop
    double probability_of_multibranch_loop(const multibranch_loop_t& mb);
 private:
//calculate equilibrium constant for a multibranch loop defined by
//a multibranch_loop_t for use in probability calculation
    PFPRECISION equilibrium_constant_for_multibranch_loop(const multibranch_loop_t&);
//helper functions for Kmb calculation
    std::vector<mb_element> construct_mb_element_array(const multibranch_loop_t&);
};
//print element array for debugging multibranch calculation
void show_mb_element_array(vector<mb_element>);
#endif

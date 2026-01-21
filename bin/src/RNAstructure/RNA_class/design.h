


#if !defined(DESIGN_H)
#define DESIGN_H

#include "RNA.h"
#include "../src/random.h"

#include <cstdlib>
#include <vector>

//! design Class.
/*!
	The design class provides an entry point for all the sequence design routines in RNAstructure.
*/

//Note the stylized comments provide facility for automatic documentation via doxygen.


using namespace std;


class design: public RNA {


	public:

		//******************************
		//Constructors:
		//******************************

		//!Constructor.  Programmer provides a filename for a ct file.

		//! Read a ct with name filename.  This is a template for the structure.  The sequence in the ct file will be ignored.
		//! If IsRNA is true (default), this is RNA, and it is DNA otherwise.

		//! \param filename is a null terminated cstring that provides the filename and path.
		//! \param alphabet is c-string that specifies the alphabet (e.g. rna, dna, or custom)
		design(const char filename[], const char* const alphabet);

		//! Design the sequence. 
		//! This function will choose a sequence with low ensemble defect to fold to the structure read in the constructor.
		//! The sequence is stored in the underlying RNA class.
		//! \param pernucdefect is the maximum allowed ensemble defect per nucleotide.
		//! \param random determines the method for choosing sequence.  If false, use pre-selected sequences, and otherwise (true), randomly choose the sequence.
		//! \param maxdepth is the maximum extent to which the structure will be sub-divided in the binary decomposition.  The default is 5.
		//! \param heuristic if true, avoids performing a partition function on the complete sequence.  Default is false.
		//! \param MaxRedesignC specifies the maximum number of redesigns.
		//! \param MaxMutateC is the maximum number of nucleotide mutations to make.
		//! \param MaxLeafRedesignC is the maximum number of time a leaf will be redesigned when building the tree.
		//! \param randomSeed is the seed for the random number generator.
		//! \return An int error code that can be resolved to a message using RNA::GetErrorMessage.  A return of 0 means no errors.
		int design_sequence(double &pernucdefect, const bool random, const bool bias, const char prob_info[], const int maxdepth = 5, bool heuristic=false, int MaxRedesignC=10,int MaxMutateC=4, int MaxLeafRedesignC=3, long randomSeed = 1L, bool allow_isolated=false);

		//! Change the defaults for maximum redesigns.
		//! This function allows the programmer to change the maximum number
		//! of redesigns at three levels from the defaults.
		//! The defaults are to redesign leaves (leaf) up to 3 times, parent nodes
		//! (parent) up to 10 times, and to mutate nucleotides for defect-weighted 
		//! mutation (mutate) up to 4 times.
		//! \param leaf is an int that specifies the maximum number of times the leaf is re-optimized at random.
		//! \param parent is an int that specifies the maximum number of redesigns per parent node.
		//! \param mutate is an int that specifies the maximum number of times a nucleotide will be mutated during defect-weighted reoptimization.
		void SpecifyRedesignLimits(int leaf, int parent, int mutate);

		//! Change the default for whether weighted-defect optimization should occur.
		//! This function allows the programmer to change whether the ascent back up the tree (if check of weighted defect fails)
		//! should only redesign a fraction of the sequence.
		//! The default is for this to occur.
		//! \param DefectWeighted is a bool that specifies whether defect-weighted re-optimization should occur (true=yes, false=no).
		void SpecifyWeightedDefect(bool DefectWeighted);


	private:
		//binary decomposition
		//This function is called to decompose a structure from nucleotide nucstart to nucend, where nucleotides missingstart to missingend (inclusive) are cut out.
		//The missing region must be contiguous.
		//The tree is stored in tree, indexed by level (first) and then position (second).
		//maxdepth is the maximum number of binary splits to be performed, if possible.  tree stores 0 at depths and positions that could not be decomposed.
		//currentdepth is the depth in the decomposition to be performed.
		//The function calls itself recursively from currentdepth to maxdepth.
		//numbering is int that labels the fragment at currentdepth-1.  It determines the labels that will be used for each fragment of the split at currentdepth. 
		void decompose(int nucstart,int nucend,int currentdepth, int maxdepth, int **tree, int missingstart=0, int missingend=0);
		
		//check if a sequence fragment is close enough to half the remaining sequence to perform decomposition
		bool closeenoughtocut(int i,int j,int nucstart,int nucend,int missingstart,int missingend,double CLOSENESS);

		//mark the binary decomposition tree
		void marktree(int beststart,int bestend,int nucstart,int nucend,int missingstart,int missingend,int currentdepth,int **tree);

		//decide if the fragment from currentstart to currentend is better than a fragment from missingstart to missingend in that it is closer to half the sequence.
		void bestdecomposition(int nuctstart, int nucend, int currentstart, int currentend, int *beststart, int *bestend, int missingstart, int missingend);

		//Once the structure is decomposed and the decomposition is in int **tree, this function selects the sequence that will fold
		//Returns the ensemble defect or the reulting design
		double SelectSequence(int **tree, bool random, bool bias,int depth,const double pernucdefect, long seed=1, bool allowisolated=false);

		//once the structure is decomposed and the decomposition is in int **tree, this function selects the sequence that will fold
		//Returns the ensemble defect or the reulting design
		double SelectSequenceHeuristic(int **tree, bool random, bool bias, int depth,const double pernucdefect, long seed=1, bool allowisolated=false);

		//Using the tree, find all fragments at level level and place on the stack
		void FindFragments(int **tree,int level, int start, int stop, int missingstart, int missingstop, vector<int> *stackstart, vector<int> *stackend, vector<int> *stackmissingstart, vector<int> *stackmissingend, vector<int> *stackfragmentdepth);


		//Fill in the sequence for fragment from start to end, excluding nucs from missingstart to missingend
		void FillSequence(int start, int end, int missingstart, int missingend, bool random, bool bias, randomnumber *dice, vector<vector<string> > &Helices, vector<vector<string> > &Loops); 

		//For j in the complete sequence, map j to the index in a fragment starting at start, with a piece from missingstart to missingend absent, where the absent fragment is filled with XXXXXX 
		int MapNuctoFragment(int j,int start,int missingstart,int missingend);
		int MapFragmenttoNuc(int j,int start,int missingstart,int missingend);
		void Mutation(int maxDefPos,int start,int missingstart,int missingend, char* sequence,vector<int>& Mutated);
		void GetDefect(int start, int end, int missingstart, int missingend, vector<double>& def, double& defect, RNA* fragment);
		void Debug1(int start, int end, int missingstart, int missingend, char* sequence, RNA* fragment);
		void PlaceSeqOnStack(vector<int> *stackstart, vector<int> *stackend, vector<int> *stackmissingstart, vector<int> *stackmissingend, vector<int> *stackfragmentdepth);
		void StoreMutation(int start, int end, int missingstart, int missingend, char* sequence);
		void StoreBestSequence(int start, int end, int missingstart, int missingend, char** sequence, int fragmentdepth);
		char tonuc(int i);
		int toint(char i);
		double leafdesign(int start,int end,int missingstart,int missingend, bool random, bool bias, randomnumber *dice, vector<vector<string> > *Helices,vector<vector<string> > *Loops, double pernucdefect,Thermodynamics *thermo, bool allowisolated);

		//Structure that does Defect Weighted Leaf Optimization
		void LeafOptimize(const double pernucdefect, randomnumber& dice, double& defect, vector<double>& def, RNA* fragment, int start, int end, int missingstart, int missingend, char* sequence, Thermodynamics *thermo, bool allowisolated);


		int randomnuc(randomnumber *dice);//return the integer for a random nucleotide choice
		void randompair(int *i, int *j, randomnumber *dice);//set integer sequence representations for a random pair choice

		int MaxRedesign;
		int MaxMutate;
		int MaxLeafRedesign;
		bool defectweighted;//indicate if defect-weighted choice of fragment should occur, default is true
		//vector<string> helices;
		//vector<string> loops;
		int numbering;




		
		//These need to be read from a data file, but for now, these are hard-wired (in the constructor):
		//First, these are pairs that are allowed, and the bias towards choosing them:
		//vector<vector<int> > pairs; //first dimension is the pair index, second dimension is nucleotide IDs.
		vector<vector<double>> pairbias; //GIve the probablity of choosing the pair with index i

		//Now indicate single stranded nucleotides allowed, and the bias:
		//vector<int> singles;//A list of nucleotides allowed for being single stranded
		vector<double> singlebias;//For each entry in singles, this array gives the probability for sampling

		bool bias_in_leaf_refinement;//Track if a sequence bias should be used in leaf refinement.

};

inline int pow(const int i, const int j) {
	int total,count;

	if (j==0) return 1;
	total = i;

	for (count=2;count<=j;count++) {
		total = total * i;
	}

	return total;

}

#endif //!defigned (DESIGN_H)


/*
 * A program that designs a nucleic acid sequence, with low ensemble defect, 
 * that is expected to fold into the structure specified by the user.
 * The designed strand can be composed of either DNA or RNA.
 *
 * (c) 2015 Mathews Lab, University of Rochester Medical Center.
 * Written by Richard M. Watson
 */

#ifndef DESIGN_INTERFACE_H
#define DESIGN_INTERFACE_H

#include <vector>

#include "../RNA_class/RNA.h"
#include "../RNA_class/design.h"


class designInterface {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	designInterface();

	/*
	 * Name:        parse
	 * Description: Parses command line arguments to determine what options are required for a particular calculation.
	 * Arguments:
	 *     1.   The number of command line arguments.
	 *     2.   The command line arguments themselves.
	 * Returns:
	 *     True if parsing completed without errors, false if not.
	 */
	bool parse( int argc, char** argv );

	/*
	 * Name:        run
	 * Description: Run calculations.
	 */
	int run();

	/**
	 * Shows interface parameters (for debugging)
	 */
	void showParams();

 private:
	// Private variables.

	// Description of the calculation type.
	string calcType;

	// The input ct file which describes the structure to design. The sequence in the ct file will be ignored.
	string ctFile;           

	// Determines the method for choosing sequence -- use pre-selected sequences (true) or randomly choose the sequence (false)
	bool preselectSequences;

	bool biasDesign;


	//bool to track whether isolated base pairs are allowed
	bool allow_isolated;

	string biasProbFile;
	
	// Name of the alphabet (e.g. "rna" or "dna" etc).
	string alphabet;	

	//// if true, avoids performing a partition function on the complete sequence.  Default is false.
	//bool useHeuristicMode;  

	// The maximum allowed ensemble defect per nucleotide.
	double defectThreshold;

	// The maximum extent to which the structure will be sub-divided in the binary decomposition.  The default is 5.
	int maxDepth;
	
	// The output ct structure file.
	string outFile;

	//// Turn pseudoknot design on or off
	//bool designPseudoknot; 

	// The maximum number of redesigns per parent node. The default is 10.
	int maxRedesign;

	// The maximum number of times a nucleotide will be mutated during defect-weighted reoptimization. The default is 4.
	int maxMutate;
	
	// The maximum number of times the leaf is re-optimized at random. The default is 3.
	int maxLeafRedesign;

	// Print sequence to standard output.
	//bool stdPrint;      Not used. Instead, see docs for --output flag

	// Whether or not to time the design process.
	bool useTimer;

	// Random seed
	long randSeed;

	// The program displays the random seed, unless it was specified by the user. This keeps track of whether or not they did.
	bool setRandomSeed;

	//void DesignPseudoknot(bool IsRNA, RNA* rna, RNA* rna1, vector<vector<string> > &Helices, vector<vector<string> > &Loops, vector<int> &pairsBroken, int &sequenceLength, design* des, long randomSeed = 1L);
};

#endif /* DESIGN_INTERFACE_H */

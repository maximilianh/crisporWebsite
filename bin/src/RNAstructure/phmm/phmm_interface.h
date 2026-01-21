/*
 * A program that folds a strand of nucleic acids.
 * This strand of nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */


#include "../RNA_class/RNA.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"
#include "../src/phmm/phmm_aln.h"

class phmm_interface {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	phmm_interface();

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
	void run();

 private:
	// Private variables.

	// Input and output file names.
	string seq1;          // The input sequence files
	string seq2;
	string outfile;		//file to write output

	//  Flag signifying whether to use maximum likelihood alignment
	bool ML;
	//  Flag signifying whether to print log probabilities (default no)
	bool LP;

};


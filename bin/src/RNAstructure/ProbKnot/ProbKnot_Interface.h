/*
 * A program that finds pseudoknots in a strand of nucleic acids.
 * These nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#ifndef PROBKNOT_H
#define PROBKNOT_H

#include "../RNA_class/RNA.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

class ProbKnot {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	ProbKnot();

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

	// Description of the calculation type.
	string calcType;

	// Input and output file names.
	string inFile;           // The input file.
	string ctFile;           // The output ct file.

	// Flag signifying if calculation handles RNA (true) or DNA (false).
	bool isRNA;
    
    // Boolean flag signifying if input is an ensemble structure file (true) or not (false).
	bool isEnsemble;

	// Boolean flag signifying if input is a sequence file (true) or not (false).
	bool isSequence;

	// The number of iterations for the algorithm to do.
	int iterations;

	// The minimum helix length allowed.
	int minHelixLength;

	//probability threshold needed to accept pair, defaults to zero 
	double threshold; 
};

#endif /* PROBKNOT_H */

/*
 * A program that folds a strand of nucleic acids from a previously folded save file.
 * This strand of nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#ifndef REFOLD_H
#define REFOLD_H

#include "../RNA_class/RNA.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

class refold {
 public:
 	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	refold();

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
	string saveFile;         // The input save file.
	string ctFile;           // The output ct file.

	// The maximum number of structures to generate.
	int maxStructures;

	// The maximum percent energy difference.
	double percent;

	// The window size for calculation.
	int windowSize;
};

#endif /* REFOLD_H */

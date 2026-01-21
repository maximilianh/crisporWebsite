/*
 * A program that creates structures based on levels of probable pairs.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#ifndef PROBABLE_PAIRS_H
#define PROBABLE_PAIRS_H

#include "../RNA_class/RNA.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

class ProbablePair {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	ProbablePair();

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
	string input;            // The input plot data file.
	string ctFile;           // The output ct file.

	// Flag signifying if calculation handles RNA (true) or DNA (false).
	bool isRNA;

	// Boolean flag signifying if input is a sequence file (true) or not (false).
	bool isSequence;

	// Float designating what threshold should be used.
	float threshold;
};

#endif /* PROBABLE_PAIR_H */

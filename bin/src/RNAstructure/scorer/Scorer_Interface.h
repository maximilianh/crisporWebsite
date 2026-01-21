/*
 * A program that scores two structures and outputs their sensitivity and PPV.
 * These structures can be composed of either DNA or RNA.
 *
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#ifndef SCORER_INTERFACE_H
#define SCORER_INTERFACE_H

#include <iomanip>

#include "../RNA_class/RNA.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"
#include "../src/score.h"

class Scorer_Interface {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	Scorer_Interface();

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
	string predicted;        // The input predicted ct file.
	string accepted;         // The input accepted ct file.
	string output;           // The output scoring file.

	// Boolean determining if exact comparisons are mandated (true) or if slippage is allowed (false).
	bool exact;

	// The number of the predicted structure to compare with the accepted structure, if there is more than one structure and a particular one is specified.
	int number;

	// Boolean determining if printing out to the screen is activated (true) or not (false).
	bool print;
};

#endif /* SCORER_INTERFACE_H */

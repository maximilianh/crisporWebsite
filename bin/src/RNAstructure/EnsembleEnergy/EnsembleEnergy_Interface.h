/*
 * A program that prints out an ensemble energy of nucleic acid input.
 *
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#ifndef ENSEMBLE_ENERGY_H
#define ENSEMBLE_ENERGY_H

#include "../RNA_class/RNA.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

class EnsembleEnergy_Interface {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	EnsembleEnergy_Interface();

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

	// Input file name.
	string input;

	// Flag signifying if calculation handles RNA (true) or DNA (false).
	bool isRNA;

	// Boolean flag signifying if input is a sequence file (true) or not (false).
	bool isSequence;

	// Boolean flag signifying if output progress messages will be suppressed.
	bool isSilent;
};

#endif /* ENSEMBLE_ENERGY_H */

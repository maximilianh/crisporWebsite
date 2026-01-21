/*
 * A program that calculates thermodynamic parameters from a set of oligonucleotides.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#ifndef OLIGOSCREEN_H
#define OLIGOSCREEN_H

#include "../RNA_class/Oligowalk_object.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

class oligoscreen {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	oligoscreen();

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
	string listFile;         // The input oligo list file.
	string reportFile;       // The output report file.

	// Flag signifying if calculation handles RNA (true) or DNA (false).
	bool isRNA;

	// The temperature at which calculation occurs.
	double temperature;
};

#endif /* OLIGOSCREEN_H */

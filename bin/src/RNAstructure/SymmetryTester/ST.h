#pragma once

#include "../RNA_class/RNA.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

class ST_Interface {
public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	ST_Interface();

	/*
	 * Name:        parse
	 * Description: Parses command line arguments to determine what options are required for a particular calculation.
	 * Arguments:
	 *     1.   The number of command line arguments.
	 *     2.   The command line arguments themselves.
	 * Returns:
	 *     True if parsing completed without errors, false if not.
	 */
	bool parse(int argc, char** argv);

	/*
	 * Name:        run
	 * Description: Run calculations.
	 */
	bool run();

private:
	// Private variables.

	


	

	// Nucleic acid type (rna, dna or custom)
	string alphabet;

	
};
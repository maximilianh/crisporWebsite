/*
 * A program that converts a dot bracket file to a CT file.
 *
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#ifndef DOT2CT_INTERFACE_H
#define DOT2CT_INTERFACE_H

#include "../RNA_class/RNA.h"
#include "../src/ParseCommandLine.h"

class dot2ct_Interface {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	dot2ct_Interface();

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
	bool run();

 private:
	// Private variables.

	// Description of the calculation type.
	string calcType;

	// Input and output file names.
	string bracketFile;      // The input dot bracket file.
	string ctFile;           // The output ct file.
	bool quiet;
};

#endif /* DOT2CT_INTERFACE_H */

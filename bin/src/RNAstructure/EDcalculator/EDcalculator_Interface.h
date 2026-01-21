/*
 * A program that designs a nucleic acid sequence, with low ensemble defect, 
 * that is expected to fold into the structure specified by the user.
 * The designed strand can be composed of either DNA or RNA.
 *
 * (c) 2015 Mathews Lab, University of Rochester Medical Center.
 * Written by Richard M. Watson
 */

#ifndef EDCALCULATOR_INTERFACE_H
#define EDCALCULATOR_INTERFACE_H

#include <vector>
#include <string>
#include "../RNA_class/RNA.h"

class EDcalculator_Interface {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	EDcalculator_Interface();

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

 private:
	// Private variables.

	// The input ct file which describes the structure to design. The sequence in the ct file will be ignored.
	string ctFile;           

	// Name of the alphabet (e.g. "rna" or "dna" etc).
	string alphabet,nucfilename;

	int structurenumber,nuc_start,nuc_end;

	bool raw;

	string constraintFile;

	string outputFile;

	//bool to track whether isolated base pairs are allowed
	bool allow_isolated;
};

#endif /* EDCALCULATOR_INTERFACE_H */

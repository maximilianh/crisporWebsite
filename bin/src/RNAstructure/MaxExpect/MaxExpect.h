/*
 * A program that predicts structures composed of probable base pairs and single-stranded nucleotides.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter.
 */

#ifndef MAX_ACCURACY_H
#define MAX_ACCURACY_H

#include "../RNA_class/RNA.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

class MaxExpect {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	MaxExpect();

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

	// The weight given to base pairs.
	double gamma;

	// Nucleic acid type (rna, dna or custom)
	string alphabet;

	// Boolean flag signifying if input is a sequence file (true) or not (false).
	bool isSequence;

	// The maximum number of structures to generate.
	int maxStructures;

	// The maximum percent energy difference.
	double percent;

	// The window size for calculation.
	int windowSize;
};

#endif /* MAX_ACCURACY_H */

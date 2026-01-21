/*
 * A program that folds two strands of nucleic acids.
 * These strands of nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#ifndef BIFOLD_H
#define BIFOLD_H

#include "../RNA_class/HybridRNA.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

class bifold {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	bifold();

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
	string seqFile1;         // The first input sequence file.
	string seqFile2;         // The second input sequence file.
	string ctFile;           // The output ct file.
	string saveFile;         // The optional output save file.
	
	// Nucleic acid type (rna, dna or custom)
	string alphabet;

	// Flag signifying if intramolecular pairs are forbidden (true) or not (false).
	bool forbidIntramolecular;

	// The maximum internal bulge loop size.
	int maxLoop;
	
	// The maximum number of structures to generate.
	int maxStructures;

	// The maximum percent energy difference.
	double percent;

	// The temperature at which calculation occurs.
	double temperature;
	
	// The window size for calculation.
	int windowSize;

	//Indicate if the input is a list of sequences, rather than two sequences
	bool IsList;
};

#endif /* BIFOLD_H */

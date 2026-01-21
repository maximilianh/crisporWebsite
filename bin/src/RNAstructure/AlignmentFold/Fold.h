/*
 * A program that folds a strand of nucleic acids.
 * This strand of nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#ifndef FOLD_INTERFACE_H
#define FOLD_INTERFACE_H

#include "../RNA_class/RNA.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

class Fold_Interface {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	Fold_Interface();

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
	string seqFile;          // The input sequence file.
	string ctFile;           // The output ct file.
	string saveFile;         // The optional output save file.

	string constraintFile;   // The optional folding constraints file.
	string experimentalFile; // The optional input bonus file.
	string SHAPEFile;        // The optional SHAPE constraints file.
	string DSHAPEFile;		//The optional differential SHAPE file.
	string DMSFile;        // The optional DMS constraints file.
    string DMSNTFile;
	string CMCTFile;        // The optional CMCT constraints file.

	string singleOffsetFile; // The optional single strand offset file.
	string doubleOffsetFile; // The optional double strand offset file.

	// The input offset.
	double experimentalOffset;

	// The input scaling.
	double experimentalScaling;

	// The intercept for SHAPE constraints.
	double intercept;

	// The intercept for single-stranded SHAPE constraints.
	double interceptSingle;

	// Nucleic acid type (rna, dna or custom)
	string alphabet;

	//  Flag signifying whether only the mfe structure is needed
	bool quickfold;

	//  Flag signifying whether O(N^4) internal loop search should be used
	bool simple_iloops;

	//  Flag signifying whether coaxial stacking recursions should be used
	bool disablecoax;

	// The maximum pairing distance.
	int maxDistance;

	// The maximum internal bulge loop size.
	int maxLoop;

	// The maximum number of structures to generate.
	int maxStructures;

	// The maximum percent energy difference.
	double percent;

	// The slope for SHAPE constraints.
	double slope;

	// The slope for differential SHAPE data.
	double Dslope;

	// The slope for single-stranded SHAPE constraints.
	double slopeSingle;

	// The temperature at which calculation occurs.
	double temperature;

	// The number of bootstraping iterations to be done.
	double bootstrap;

	// The window size for calculation.
	int windowSize;

	// Allow the user to override the sequence name.
	string sequenceName;

	// suppress unnecessary output
	bool quiet;

	// write structure in dot-bracket notation.
	bool useBracketNotation;

	//Auxiliary function used to sample a SHAPE, DMS, CMCT file.
	string sample_file(string shapefile, int numnuc, int iter);
};

#endif /* FOLD_INTERFACE_H */

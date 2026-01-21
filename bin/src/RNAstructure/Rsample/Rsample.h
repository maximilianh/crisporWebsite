/*
 * A program that calculates the partition function for a strand of nucleic acids.
 * This strand of nucleic acids can be composed of either DNA or RNA.
 *
 */

#ifndef RSAMPLE_H
#define RSAMPLE_H

#include <string>
#include "../RNA_class/RNA.h"
#include "../RNA_class/RsampleData.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

class RsampleInterface {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	RsampleInterface();

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

	// Input and output file names.
	string seqFile;          // The input sequence file.
	string pfsFile;          // The output partition function save file.
	string SHAPEFile;        // The SHAPE constraints file.
	
	// Names of paired-end, paired-middle and unpaired files
	string peFile;
	string pmFile;
	string upFile;


	// The Cparam for SHAPE constraints.
	double Cparam;
	
	// The Offset for SHAPE constraints.
	double Offset;

	// Flag signifying if calculation handles RNA (true) or DNA (false).
	bool isRNA;

	// The maximum pairing distance.
	int maxDistance;

    // number of samples for stochastic sampling
	int numsamples;

	// The temperature at which calculation occurs.
	double temperature;

	//the random number seed 
	int seed;


	//Flag to use the DMS distributions instead of SHAPE
	bool isDMS;


	//Flag to use the DMS distributions, capped at a maximum of 4.0, instead of SHAPE
	double max;

	// Nucleic acid type (rna, dna or custom)
	string alphabet;


	std::string constraintFile;   // The optional folding constraints file.
};

#endif /* RSAMPLE_H */

/*
 * A program that calculates the probscan function for a strand of nucleic acids.
 * This strand of nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#ifndef PROBSCANINTERFACE
#define PROBSCANINTERFACE

#include "../RNA_class/RNA.h"
#include "../RNA_class/ProbScan.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

class probscanInterface {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	probscanInterface();

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
    void processLine(std::string input);

	/*
	 * Name:        run
	 * Description: Run calculations.
	 */
	void run();

 private:
	// Private variables.

	// Description of the calculation type.
	std::string calcType;

    //nuc indices for scan
    std::vector<int>mb;

	// Input and output file names.
	std::string seqFile;          // The input sequence file.
	std::string pfsFile;          // The output probscan function save file.
	std::string inputFile;
    std::string loop_file;   //input file with multibranch loops
    std::string pairSpecification;
	std::string pairstemloopSpecification;
	bool fromSequence;

    // prediction mode options
    bool hairpin;
    bool helix;
    int numstacks;
    bool internal;
    bool bulge;
    bool multibranch;
    bool pairs;
	bool stemloop;

	double threshold;

	// The temperature at which calculation occurs.
	double temperature;

	// Nucleic acid type (rna, dna or custom)
	string alphabet;

	//fixed digits or scientific
	bool fixed;
};

#endif /* PROBSCANINTERFACE */

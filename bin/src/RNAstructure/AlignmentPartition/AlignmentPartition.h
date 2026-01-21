/*
 * A program that performs AlignmentPartition on multiple sequence alignment of homologous RNA.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Abhinav Mittal
 */

#ifndef PARTITION_H
#define PARTITION_H

#include "../RNA_class/RNA.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

#ifdef SMP
	#define PARTITION_EXE_NAME "AlignmentPartition-smp"
#else
	#define PARTITION_EXE_NAME "AlignmentPartition"
#endif

class AlignmentPartitionInterface {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	 AlignmentPartitionInterface();

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
	string seqFile;          // The input sequence file.
	string pfsFile;          // The output partition function save file.

	string constraintFile;   // The constraints file.
	string doubleOffsetFile; // The optional double strand offset file.
	string experimentalFile; // The experimental pair bonus file.
	string SHAPEFile;        // The SHAPE constraints file.

	// The experimental pair bonus offset.
	double experimentalOffset;

	// The experimental pair bonus scaling.
	double experimentalScaling;

	// The intercept for SHAPE constraints.
	double intercept;

	// Nucleic acid type
	string alphabet;

	// The maximum pairing distance.
	int maxDistance;

	// The slope for SHAPE constraints.
	double slope;

	// The temperature at which calculation occurs.
	double temperature;

	// Flag indicating whether coaxial stacking recursions are used
	bool disablecoax;

	// Suppress progress and other unnecessary output.
	bool quiet;

	// base pairing cutoff
	float bp_cutoff;

	// maximum size of internal loop
	int maxinter;
};

#endif /* PARTITION_H */

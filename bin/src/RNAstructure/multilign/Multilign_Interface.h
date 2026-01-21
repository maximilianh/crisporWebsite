/*
 * A program that predicts structures from multiple sequences using progressive iterations of Dynalign.
 *
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#ifndef MULTILIGN_INTERFACE_H
#define MULTILIGN_INTERFACE_H

#include <sstream>

#include "../RNA_class/Multilign_object.h"
#include "../src/configfile.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"
#include "../src/TProgressDialog.h"

class Multilign_Interface {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	Multilign_Interface();

	/*
	 * Name:        parse
	 * Description: Parse command line arguments to determine what options are
	 *              required for a partitcular calculation.
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

	/**
	 * Name:        usage
	 * Description: Print out a special usage message for this interface.
	 * Arguments:
	 *     1. parser
	 *        The command line parser.
	 */
	void usage( ParseCommandLine* parser );

 private:
	// Private variables important for multiple sequence calculations.

	// The calculation type description.
	string calcType;

	// 2D vector of input file names.
	vector< vector<string> > files;

	// Boolean flag, true if SHAPE data is being used, false if not.
	bool hasSHAPE;

	// The number of processors the calculation is run on.
	// By default, this is 1 (serial mode), but it may be more than that if the calculation is run in parallel.
	int processors;

	// The number of sequences used as input.
	int seqNumber;

 private:
	// Private variables.

	// The multiple alignment file.
	string alignmentFile;

	// The alignment window size.
	int alignmentWindow;

	// The base pairing window size.
	int basepairWindow;

	// The gap penalty.
	double gap;

	// The sequence index that will be considered the "main" sequence.
	int indexSeq;

	// The intercept for SHAPE constraints.
	double intercept;

	// Boolean flag specifying whether base pair inserts should be allowed (true) or not (false).
	bool inserts;

	// Boolean flag specifying whether the nucleic acid type is RNA (true) or DNA (false).
	bool isRNA;

	// The number of iteraions over which Multilign runs.
	int iterations;

	// Boolean flag specifying whether to keep intermediate files (dsv and aout) or not after calulation.
	bool keepIntermediate;

	// Boolean flag specifying whether to do local folding (true) or global folding (false).
	bool local;

	// The maximum DSV change allowed.
	double maxDsvChange;

	// The maximum number of pairs allowed.
	int maxPairs;

	// The default M (separation) parameter.
	double maxSeparation;

	// The maximum number of structures.
	int maxStructures;

	// The maximum percent energy difference.
	double percent;

	// The maximum percent difference in folding free energy change from single sequence folding for pairs allowed in a subsequent calculation.
	double percentSingle;

	// Boolean flag specifying whether sequences should be randomized (true) or not (false).
	bool random;

	// The slope for SHAPE constraints.
	double slope;

	// The termperature at which calculation occurs.
	double temperature;
};

#endif /* MULTILIGN_INTERFACE_H */

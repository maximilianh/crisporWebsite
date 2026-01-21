/*
 * A program that predicts structures using the TurboFold algorithm.
 *
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#ifndef TURBOFOLD_INTERFACE_H
#define TURBOFOLD_INTERFACE_H

#include <sstream>

#include "../src/TurboFold_object.h"
#include "../src/configfile.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"
#include "../src/TProgressDialog.h"
#include "../src/phmm/structure/structure_object.h"

class TurboFold_Interface {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	TurboFold_Interface();
    ~TurboFold_Interface();

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
	int run();

	/**
	 * Name:        usage
	 * Description: Print out a special usage message for this interface.
	 * Arguments:
	 *     1. parser
	 *        The command line parser.
	 */
	void usage( ParseCommandLine &parser );

 private:
	// Private variables important for multiple sequence calculations.

	// The calculation type description.
	string calcType;

	// vectors of file names.
	vector<string> sequenceFiles, outputCtFiles, outputPfsFiles, shapeFiles, outputStartPfsFiles; //, rsampleFiles;

    string conf_file;
	
	// Boolean flag, true if SHAPE data is being used, false if not.
	bool hasSHAPE;

	// The number of processors the calculation is run on.
	// By default, this is 1 (serial mode), but it may be more than that if the calculation is run in parallel.
	int processors;

	//vector<double> parameters;

	// The maximum pairing distance.
	int distance;

	// The intercept for SHAPE constraints.
	double intercept;

	// The MEA mode maximum number of structures.
	int maxStructures;

	// The MEA mode gamma value.
	double meaGamma;

	// The ProbKnot mode minimum helix length.
	int minHelixLength;

	// The mode that TurboFold runs under.
	string mode;

	// The MEA mode maximum percent energy difference.
	double percent;

	// The ProbKnot mode number of iterations.
	int pkIterations;

	// The slope for SHAPE constraints.
	double slope;

	// The termperature at which calculation occurs.
	double temperature;

	// The Threshold mode probable pair threshold.
	double threshold;

	// The gamma value for TurboFold.
	double turboGamma;

	// The iterations value for TurboFold.
	double turboIterations;

	// The MEA mode window size.
	int windowSize;

	// Initialize the default alignment output filename.
    string OutAln;

    // Initialize the default alignment output format.
    string AlnFormat;

    // Initialize the default alignment output max column number.
    int ColumnNumber;

    // Rsample parameters.
    // The Cparam for SHAPE constraints.
    double Cparam;
    
    // The Offset for SHAPE constraints.
    double Offset;

    // number of samples for stochastic sampling
    int rsample_numsamples;

    // Names of paired-end, paired-middle and unpaired files
    string peFile;
    string pmFile;
    string upFile;

    // Name of partition function output file.
    //string pfsFile;

    bool is_RSample_mode;
   
    int rsample_seed;

	vector<vector<double> > shape_reactivities;

	vector<t_structure*> *fasta_sequences;
};

#endif /* TURBOFOLD_INTERFACE_H */

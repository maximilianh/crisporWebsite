/*
 * ShapeKnots, a program that predicts RNA secondary structures with pseudoknots.
 *
 * (c) 2013 
 * Mathews Lab, University of Rochester Medical Center
 * Weeks Lab, The University at North Carolina at Chapel Hill
 * Code contributors: Wayne Higgins, Stanislav Bellaousov, David H. Mathews
 */

#ifndef ShapeKnots_Interface_H
#define ShapeKnots_Interface_H

#include "../src/ShapeKnots.h"

#include <iomanip>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>

#include "../src/rna_library.h"
#include "../src/structure.h"
#include "../src/algorithm.h"
#include "../RNA_class/RNA.h"
#include "../src/ParseCommandLine.h"

//using namespace std;

class ShapeKnots_Interface {
 public:
	// Public constructor and methods.
	
	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	ShapeKnots_Interface();

	/*
	 * Name:        parse
	 * Description: Parses command line arguments to determine what options are required for a particular calculation.
	 * Arguments:
	 *     1.   The number of command line arguments.
	 *     2.   The command line arguments themselves.
	 * Returns:
	 *     True if parsing completed without errors, false if not.
	 */
	bool Parse( int argc, char** argv );

	/*
	 * Name:        run
	 * Description: Run calculations.
	 */
	bool run();

 private:
	//Define private variables.

	string seqFile;//Name of a sequence file containing input data.

	string ctFile;//Name of a CT file to which output will be written.

	string doubleOffsetFile;//Name of double-stranded offset file.

	string constraintFile;// The optional folding constraints file.

	string DMSFile;// The optional DMS constraints file.
    string DMSNTFile;// The optional DMSNT constraints file.

	double Dslope;//Differential slope parameters used with differential SHAPE restraints(kcal/mol).

    string DSHAPEFile;//Name of a differential SHAPE restraint file.

	int InMaxStructures;//Maximum number of internally generated structures for each call of the dynamic programming algorithm.

	int InPercent;//Maximum percent difference in folding free energy change for internally generated suboptimal structures for each call of the dynamic programming algorithm.

	int InWindowSize;//Window size for the internally generated suboptimal structures for each call of the dynamic programming algorithm.

	int OutMaxStructures;//Maximum number of structures to be outputted.

	int OutPercent;//Maximum percent difference in folding free energy change for generating suboptimal structures in the output.

	double P1, P2;//Pseudoknot energy model parameters (kcal/mol).

	int finallistSize; //Maximum number of helices to be processed.

	string SHAPEFile;//Name of a SHAPE restraints file.

	double slope, intercept;//Slope and intercept parameters used with SHAPE restraints(kcal/mol).

	string singleOffsetFile;//Name of single-stranded offset file.

	int OutWindowSize;//Window size for outputted suboptimal structures.

	bool ifWindowOptions;//Tracks if the OutWindowSize flag was specified by the user.

    string experimentalFile; // The optional input bonus file.
	// The input offset.
	double experimentalOffset;
	// The input scaling.
	double experimentalScaling;

};
#endif /* ShapeKnots_Interface_H */

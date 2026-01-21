/*
 * A program that calculates a probability dot plot.
 * This class can write image output to Postscript or SVG.
 * It can also write to a dot plot text file.
 *
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter.
 */

#ifndef PROBABILITY_PLOT_H
#define PROBABILITY_PLOT_H

#include "../RNA_class/RNA.h"
#include "../src/DotPlotHandler.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

class ProbabilityPlot {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	ProbabilityPlot();

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
	string inputFile;        // The input plot data file.
	string outputFile;       // The output file.

	// The number of legend entries allowed in the plot.
	int entries;

	// Boolean flag, true if SVG is to be written, false if not.
	bool isSVG;

	// Boolean flag signifying whether the input file should be a standard text file.
	bool logPlot;

	// Boolean flag signifying whether the input file should be a matrix.
	bool matrixPlot;

	// The maximum bound of the plot.
	double maxBound;

	// The minimum bound of the plot.
	double minBound;

	// Boolean flag signifying whether a text output file should be written.
	bool writeText;

	// Option configuring the output of descriptions. (see LegendDescriptionSettings for details)
	LegendDescriptionSettings descriptionSettings;
};

#endif /* PROBABILITY_PLOT_H */

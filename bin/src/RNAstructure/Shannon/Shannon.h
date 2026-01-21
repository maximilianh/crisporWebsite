/*
  * A program that calculates the Shannon entropy of each nucleotide in a sequence.
 
 *
 * (c) 2021 Mathews Lab, University of Rochester Medical Center.
 * Written by David H. Mathews
 */

#ifndef PROBABILITY_PLOT_H
#define PROBABILITY_PLOT_H

#include "../RNA_class/RNA.h"
#include "../src/DotPlotHandler.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

class Shannon {
public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	Shannon();

	/*
	 * Name:        parse
	 * Description: Parses command line arguments to determine what options are required for a particular calculation.
	 * Arguments:
	 *     1.   The number of command line arguments.
	 *     2.   The command line arguments themselves.
	 * Returns:
	 *     True if parsing completed without errors, false if not.
	 */
	bool parse(int argc, char** argv);

	/*
	 * Name:        run
	 * Description: Run calculations.
	 */
	void run();

private:
	// Private variables.


	// Input and output file names.
	string inputFile;        // The input plot data file.
	string outputFile;       // The output file.


	// Option configuring the output of descriptions. (see LegendDescriptionSettings for details)
	LegendDescriptionSettings descriptionSettings;
};

#endif /* PROBABILITY_PLOT_H */

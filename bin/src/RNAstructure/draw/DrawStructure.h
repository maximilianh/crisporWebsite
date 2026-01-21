/*
 * A program that draws a structure and writes image output.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#ifndef DRAWSTRUCTURE_H
#define DRAWSTRUCTURE_H

#include "../src/ParseCommandLine.h"
#include "../src/StructureImageHandler.h"

class DrawStructure {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	DrawStructure();

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
	string inputFile;        // The input structure ct file.
	string outputFile;       // The output image file.

	string probabilityFile;  // The optional probability annotation file.
	string loopProbabilityFile;  // The optional loop probability annotation file.
	string SHAPEFile;        // The optional SHAPE annotation file.

	// Boolean flag, true if a circular structure should be drawn, false if not.
	bool circular;

	// Boolean flag, true if circles should be around nucleotides, false if not.
	bool encircle;

	// Boolean flag, true if a linear structure should be drawn, false if not.
	bool flat;

	// Boolean flag, true if SVG is to be written, false if not.
	bool isSVG;

	// Boolean flag, true if the structures should be drawn counterclockwise, false if not.
	bool levorotatory;

	// The number of a specific structure to output, or the beginning number (if outputting a range of structures with -nn).
	int number;

	// The ending number of a range of structures.
	int endnumber;

	// The maximum number of structures to draw (used for RNAstructureWeb to limit responses.)
	int maxStructures;

	// Boolean flag, true if the probability annotation file should be text, false if not.
	bool textAnnotation;

	// Option configuring the output of descriptions. (see LegendDescriptionSettings for details)
	LegendDescriptionSettings descriptionSettings;
};

#endif /* DRAWSTRUCTURE_H */

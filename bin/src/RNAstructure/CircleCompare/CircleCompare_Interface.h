/*
 * A program that compares two structures in a circular layout.
 *
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#ifndef CIRCLECOMPARE_INTERFACE_H
#define CIRCLECOMPARE_INTERFACE_H

#include "../src/ParseCommandLine.h"
#include "../src/StructureComparedImageHandler.h"

class CircleCompare_Interface {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	CircleCompare_Interface();

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
	string predicted;        // The input predicted ct file.
	string accepted;         // The input accepted ct file.
	string output;           // The output Postscript image file.

	string probabilityFile;  // The optional probability annotation file.
	string SHAPEFile;        // The optional SHAPE annotation file.

	// Boolean flag signifying if the alternative color scheme is used (true) or not (false).
	bool alternative;

	// Boolean flag, true if circles should be around nucleotides, false if not.
	bool encircle;

	// Boolean flag signifying exact pairs in scoring (true) or not (false).
	bool exact;

	// Boolean flag, true if file names are to be shown in addition to descriptions, false if not.
	bool filenames;

	// Boolean flag, true if SVG is to be written, false if not.
	bool isSVG;

	// Boolean flag, true if the structures should be drawn counterclockwise, false if not.
	bool levorotatory;

	// The number of the predicted structure to compare with the accepted structure, if there is more than one structure and a particular one is specified.
	int number;

	// The number of predicted structures to compare against.
	int predictedCount;

	// Boolean flag, true if probability annotation is read for the predicted structure, false if for the accepted structure.
	bool probabilityForPredicted;

	// Boolean flag, true if the probability annotation file should be text, false if not.
	bool textAnnotation;
};

#endif /* CIRCLECOMPARE_INTERFACE_H */

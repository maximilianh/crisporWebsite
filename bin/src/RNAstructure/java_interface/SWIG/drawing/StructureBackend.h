/*
 * A header file for a class that holds structure methods for Java drawing.
 *
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter.
 */

#ifndef STRUCTURE_BACKEND_H
#define STRUCTURE_BACKEND_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//Stored in java_interface/SWIG/drawing
#include "../../../RNA_class/RNA.h"
#include "../../../src/ErrorChecker.h"
#include "../../../src/StructureImageHandler.h"

using namespace std;

class StructureBackend {
 public:
	// Public methods.

	/**
	 * Name:        Constructor.
	 * Description: Initializes the class.
	 */
	StructureBackend();

	/**
	 * Name:        addAnnotationProbability
	 * Description: Read probability annotation into structures.
	 * Arguments:
	 *     1. file
	 *        The annotation file to read from.
	 * Returns:
	 *     True if annotation was read correctly, false if not.
	 */
	bool addAnnotationProbability( string file );

	/**
	 * Name:        addAnnotationSHAPE
	 * Description: Read SHAPE annotation into structures.
	 * Arguments:
	 *     1. file
	 *        The annotation file to read from.
	 * Returns:
	 *     True if annotation was read correctly, false if not.
	 */
	bool addAnnotationSHAPE( string file );

	/**
	 * Name:        flip
	 * Description: Flip all structures horizontally.
	 */
	void flip();

    /**
	 * Name:  setFlipped
	 * Description: Flip structures horizontally whose 'isFlipped' property 
	 * differs from the flipped argument.
	 */
	void setFlipped(bool flipped);

	/**
	 * Name:        getStructureData
	 * Description: Get the data for a particular structure.
	 * Arguments:
	 *     1. number
	 *        The structure to get data for.
	 */
	string getStructureData( int number );

	/**
	 * Name:        readStructureData
	 * Description: Read structures into this class.
	 * Arguments:
	 *     1. file
	 *        The structure file to read from.
	 * Returns:
	 *     The empty string if structures were read correctly. An error message otherwise.
	 */
	string readStructureData( string file );

	/**
	 * Name:        removeAnnotation
	 * Description: Remove annotation from all structures.
	 */
	void removeAnnotation();

	/*
	 * Name:        setNucleotidesCircled
	 * Description: Set nucleotides circled or uncircled on all structures.
	 * Arguments:
	 *     1. circled
	 *        True if nucleotides should be circled, false if not.
	 */ 
	void setNucleotidesCircled( bool circled );

	/*
	 * Name:        setBackboneStyle
	 * Description: Draw the RNA backbone in a specific style --
	 *              Normal (Radial), Circular, or Flat (Linear).
	 *
	 * Arguments:
	 *     1. style an integer that indicates the backbone style:
	 *              0: Normal (Radial)
	 *              1: Flat (linear)
	 *              2: Circular
	 */ 
	void setBackboneStyle( int style );

	/*
	 * Name:        writeDotBracketFile
	 * Description: Write a dot bracket file.
	 * Arguments:
	 *     1. inFile
	 *        The file to get dot bracket structures from.
	 *     2. outFile
	 *        The file to write dot bracket structures to.
	 * Returns:
	 *     A string showing the completion status.
	 */ 
	string writeDotBracketFile( string inFile, string outFile );

	/*
	 * Name:        writeHelixFile
	 * Description: Write a helix file from a structure.
	 * Arguments:
	 *     1. inFile
	 *        The file to get dot bracket structures from.
	 *     2. outFile
	 *        The file to write dot bracket structures to.
	 *     3. number
	 *        The structure number.
	 * Returns:
	 *     A string showing the completion status.
	 */ 
	string writeHelixFile( string inFile, string outFile, int number );

	/*
	 * Name:        writePostscriptFile
	 * Description: Write a Postscript file from a structure.
	 * Arguments:
	 *     1. outFile
	 *        The file to write Postscript structures to.
	 *     2. number
	 *        The structure number.
	 */ 
	void writePostscriptFile( string outFile, int number );

	/*
	 * Name:        writeSVGFile
	 * Description: Write a SVG file from a structure.
	 * Arguments:
	 *     1. outFile
	 *        The file to write Postscript structures to.
	 *     2. number
	 *        The structure number.
	 */ 
	void writeSVGFile( string outFile, int number );

	/*
	 * Name:        getStructureCount
	 * Description: Return the number of structures (i.e. RNA::
	 */ 
	int getStructureCount();

 private:
	// Private variables.

	// The file holding structure data.
	string structureDataFile;

	// The vector of StructureImageHandlers.
	vector<StructureImageHandler> structureHandlers;

	// The number of structure handlers held in this class.
	int structures;
};

#endif /* STRUCTURE_BACKEND_H */

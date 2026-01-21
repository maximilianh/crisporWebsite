/*
 * A class that holds structure image data and can write those images to files if necessary.
 *
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter.
 */

#ifndef STRUCTURE_IMAGE_HANDLER_H
#define STRUCTURE_IMAGE_HANDLER_H

#include <limits>
#include <map>
#include <vector>

#include "../RNA_class/RNA.h"
#include "DrawingConstants.h"
#include "ErrorChecker.h"

using namespace std;

class StructureImageHandler {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes private variables.
         */
	StructureImageHandler();

	/*
	 * Name:        addAnnotationProbability
	 * Description: Add probability annotation to the structure.
	 * Arguments:
	 *     1. file
	 *        The file that holds probability annotation information.
	 *     2. text
	 *        True if the annotation file is a dot plot text file, false if not.
	 * Returns:
	 *     A string showing the completion status.
	 */ 
	string addAnnotationProbability( string file, bool text = false );

	/*
	 * Name:        addAnnotationLoopProbability
	 * Description: Add loop probability annotation to the structure.
	 * Arguments:
	 *     1. file
	 *        The file that holds probability annotation information.
     *     2. r
     *        An RNA object with the relevant structure information
	 * Returns:
	 *     A string showing the completion status.
	 */ 
	string addAnnotationLoopProbability( string file, RNA& r , int structurenumber=1);


	/*
	 * Name:        addAnnotationSHAPE
	 * Description: Add SHAPE annotation to the structure.
	 * Arguments:
	 *     1. file
	 *        The file that holds SHAPE annotation information.
	 * Returns:
	 *     A string showing the completion status.
	 */ 
	string addAnnotationSHAPE( string file );

	/*
	 * Name:        flipHorizontally
	 * Description: Flip the structure image horizontally.
	 */ 
	void flipHorizontally();
	void setFlipped(bool flipped);

	/*
	 * Name:        readCircular
	 * Description: Read a structure and determine its coordinates for a
	 *              circular layout.
	 * Arguments:
	 *     1. file
	 *        The file that holds structure information.
	 *     2. number
	 *        The structure to get information about.
	 * Returns:
	 *     A string showing the completion status.
	 */ 
	string readCircular( string file, int number );

	/*
	 * Name:        readLinear
	 * Description: Read a structure and determine its coordinates for a
	 *              linear layout.
	 * Arguments:
	 *     1. file
	 *        The file that holds structure information.
	 *     2. number
	 *        The structure to get information about.
	 * Returns:
	 *     A string showing the completion status.
	 */ 
	string readLinear( string file, int number );

	/*
	 * Name:        readRadial
	 * Description: Read a structure and determine its coordinates for a
	 *              radial layout.
	 * Arguments:
	 *     1. file
	 *        The file that holds structure information.
	 *     2. number
	 *        The structure to get information about.
	 * Returns:
	 *     A string showing the completion status.
	 */ 
	string readRadial( string file, int number );

	/*
	 * Name:        removeAnnotation
	 * Description: Remove annotation from nucleotides.
	 */
	void removeAnnotation();

	/*
	 * Name:        setNucleotidesCircled
	 * Description: Set nucleotides to be surrounded by circles when drawn.
	 * Arguments:
	 *     1. encircle
	 *        True if nucleotides should be circled, false if not.
	 */
	void setNucleotidesCircled( bool encircle );

	/*
	 * Name:        toString
	 * Description: Return a string representation of the structure handled by
	 *              this object. This is usually a Java convention, not a C++
	 *              one, but it is useful to get a snapshot of this structure
	 *              and access all information about it with a single easy call.
	 */
	string toString();

	/*
	 * Name:        writePostscript
	 * Description: Write a structure as a Postscript image.
	 * Arguments:
	 *     1. file
	 *        The file name to write.
	 *     2. append
	 *        True if structures should be appended on to this file, false
	 *        if a new structure file should be created for writing.
	 *     3. pagenumber
	 *        If multiple structures are being written, the pagenumber of the structure to be written	 
	 *     3. totalpages
	 *        The total number of structures to be written
	 */
	void writePostscript( string file, bool append, int pagenumber=1, int totalpages=1);

	/*
	 * Name:        writeSVG
	 * Description: Write a structure as an SVG image.
	 * Arguments:
	 *     1. file
	 *        The file name to write.
	 */
	void writeSVG( string file );

	// Returns the description (which was read in from the structure file)
	string getDescription();
	
	// Set the description that is written in the legend of the output figure.
	void setDescription(const string desc);

 private:
	// Private image writing workhorse function.

	/*
	 * Name:        writeImageFile
	 * Description: Write a structure as an image.
	 *              This method encapsulates the common strategy for
	 *              writing structure images in different formats.
	 * Arguments:
	 *     1. file
	 *        The file name to write.
	 *     2. append
	 *        True if structures should be appended on to this file, false
	 *        if a new structure file should be created for writing.
	 *     3. isSVG
	 *        True if the file should be SVG, false if not.
	 *     3. pagenumber
	 *        If multiple structures are being written, the pagenumber of the structure to be written	 
	 *     3. totalpages
	 *        The total number of structures to be written
	 */
	void writeImageFile( string file, bool append, bool isSVG, int pagenumber=1, int totalpages=1);
    std::map<int,int> getPairs();
	/**
	* Clear pairs, annotations, coordinates, legend, extras, etc 
	* in preparation for reading in a new structure.
	*/
	void reset();

 protected:
	// Protected data structures.

	// Annotation colors.
	vector<string> annotations;

	// Backbone, nucleotide, and label coordinates.
	vector<string> coordinates;

	// Legend text.
	vector<string> legend;

	// Legend text colors.
	vector<string> legendColors;

	// Any extra data that is necessary to place on the image.
	// Note that this data is assumed to be formatted properly for a given
	// image type.
	vector<string> extras;

	// Pairing curve coordinates.
	vector<string> pairs;

 protected:
	// Protected variables.

	// A boolean flag specifying if the structure is bimolecular (true) or
	// not (false).
	bool bimolecular;

	// A boolean flag specifying if the nucleotides are circled (true) or
	// not (false).
	bool circleNucs;

	// The description of the image.
	string description;

	// The maximum X bound of the image.
	double maxX;

	// The maximum Y bound of the image.
	double maxY;

	// whether or not the structure has been flipped
	bool isFlipped;
};

#endif /* STRUCTURE_IMAGE_HANDLER_H */

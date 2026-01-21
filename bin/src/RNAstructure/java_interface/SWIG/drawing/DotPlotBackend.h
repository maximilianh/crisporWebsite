/*
 * A header file for a class that holds dot plot methods for Java drawing.
 *
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter.
 */

#ifndef DOT_PLOT_BACKEND_H
#define DOT_PLOT_BACKEND_H

#include <math.h>

//This file is in java_interface/SWIG/drawing, so the path to the root RNAstructure directory is ../../../
#include "../../../RNA_class/Dynalign_object.h"
#include "../../../src/DotPlotHandler.h"
#include "../../../src/ErrorChecker.h"

using namespace std;

class DotPlotBackend {
 public:
	// Public methods.

	/**
	 * Name:        Constructor.
	 * Description: Initializes the class.
	 */
	DotPlotBackend();

	/**
	 * Name:        Destructor.
	 * Description: Frees dynamic variables where necessary.
	 */
	~DotPlotBackend();

	/*
	 * Name:        getBounds
	 * Description: Get the maximum bounds of the plot.
	 * Returns:
	 *     The maximum bounds of the plot, as a string.
	 */
	string getBounds();

	/*
	 * Name:        getColors
	 * Description: Get the current, minimum, and maximum colors allowed in
	 *              this plot.
	 * Returns:
	 *     The colors, as a string.
	 */
	string getColors();

	/*
	 * Name:        getCurrentMax
	 * Description: Get the current maximum plot value.
	 * Returns:
	 *     The value.
	 */
	double getCurrentMax();

	/*
	 * Name:        getDefaultMin
	 * Description: Get the default minimum plot value.
	 * Returns:
	 *     The value.
	 */
	double getDefaultMin();

	/*
	 * Name:        getDefaultMax
	 * Description: Get the default maximum plot value.
	 * Returns:
	 *     The value.
	 */
	double getDefaultMax();

	/*
	 * Name:        getCurrentMin
	 * Description: Get the current minimum plot value.
	 * Returns:
	 *     The value.
	 */
	double getCurrentMin();

	/**
	 * Name:        getDotData
	 * Description: Get the data for a particular dot.
	 * Arguments:
	 *     1. i
	 *        The dot index in sequence 1, one-indexed.
	 *     2. j
	 *        The dot index in sequence 2, one-indexed.
	 * Returns:
	 *     The dot's data.
	 */
	string getDotData( int i, int j );

	/**
	 * Name:        getDotDataByLocation
	 * Description: Get the data for a particular dot, based on where it was
	 *              clicked and the plot's scale.
	 * Arguments:
	 *     1. x
	 *        The clicked location's X coordinate.
	 *     2. y
	 *        The clicked location's Y coordinate.
	 *     3. scale
	 *        The scale of the plot.
	 * Returns:
	 *     The dot's data.
	 */
	string getDotDataByLocation( int x, int y, double scale );

	/**
	 * Name:        getGridLine
	 * Description: Get the data for a particular grid line.
	 * Arguments:
	 *     1. number
	 *        The index of the line, one-indexed.
	 * Returns:
	 *     The grid line data.
	 */
	string getGridLine( int number );

	/**
	 * Name:        getLegendData
	 * Description: Get the dot plot's legend data, with entries separated by
	 *              a new line character.
	 * Returns:
	 *     The dot plot's legend data.
	 */
	string getLegendData();

	/**
	 * Name:        getLength
	 * Description: Get the dot plot's length.
	 * Returns:
	 *     The dot plot's length.
	 */
	int getLength();

	/**
	 * Name:        readData
	 * Description: Read the data for a particular dot plot.
	 * Arguments:
	 *     1. file
	 *        The dot plot file to read from.
	 *     2. strand
	 *        The strand read.
	 * Returns:
	 *     True if data was read correctly, false if not.
	 */
	bool readData( string file, int strand );

	/*
	 * Name:        setColors
	 * Description: Set the current number of colors allowed in this plot.
	 * Arguments:
	 *     1. colors
	 *        The number of colors.
	 */
	void setColors( int colors );

	/*
	 * Name:        setRange
	 * Description: Set the range of the dot plot.
	 * Arguments:
	 *     1. minimum
	 *        The lowest allowable value.
	 *     2. maximum
	 *        The highest allowable value.
	 */
	void setRange( double minimum, double maximum );

	/*
	 * Name:        writePostscriptFile
	 * Description: Write a Postscript file from a plot.
	 * Arguments:
	 *     1. outFile
	 *        The file to write a Postscript plot to.
	 */ 
	void writePostscriptFile( string outFile );

	/*
	 * Name:        writeStructuresFile
	 * Description: Write a probable structures file from a plot.
	 * Arguments:
	 *     1. file
	 *        The file the plot was drawn from.
	 *     2. outFile
	 *        The file to write an SVG plot to.
	 */ 
	string writeStructuresFile( string file, string outFile );

	/*
	 * Name:        writeSVGFile
	 * Description: Write a SVG file from a plot.
	 * Arguments:
	 *     1. outFile
	 *        The file to write an SVG plot to.
	 */ 
	void writeSVGFile( string outFile );

	/**
	 * Name:        writeTextFile
	 * Description: Write a text file from a dot plot.
	 * Arguments:
	 *     1. file
	 *        The dot plot file to write.
	 */
	void writeTextFile( string file );

 private:
	// Private variables.

	// The current and default range values for the dot plot.
	double currentMin, currentMax;
	double defaultMin, defaultMax;

	// The dot plot handler that holds dot data.
	DotPlotHandler* handler;

	// The dot plot length.
	int length;

	// The number of colors in the plot.
	int plotColors;

	// The nucleotide sequence.
	vector<string> sequence;

	// The map of values.
	map<string,double> values;
};

#endif /* DOT_PLOT_BACKEND_H */

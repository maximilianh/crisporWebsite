/*
 * A header file for a class that holds dot plot data and can write those images to files if necessary.
 *
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 * Redone in 2012.
 * Written by Jessica S. Reuter.
 */

#include <fstream>
#include <iomanip>
#include <limits>
#include <math.h>
#include <map>
#include <ostream>
#include <stdlib.h>
#include <utility>
#include <vector>
#include <algorithm>
#include "DrawingConstants.h"

/*
 * A namespace that holds dot plot constants.
 */
namespace DotPlotConstants {

	// Dividers.
	const string DIVIDER_ENERGY = "Free Energy (kcal/mol)";
	const string DIVIDER_PROBABILITY = "-log10(Probability)";

	// Accepted numbers of legend entries.
	const int ENTRIES_MINIMUM = 3;
	const string ENTRIES_MINIMUM_STRING = "3";
	const int ENTRIES_DEFAULT = 5;
	const string ENTRIES_DEFAULT_STRING = "5";
	const int ENTRIES_MAXIMUM = 15;
	const string ENTRIES_MAXIMUM_STRING = "15";

	// Postscript legend text constants.
	const int LEGEND_TEXT_SIZE_PS = 15;
	const string LEGEND_TEXT_OPEN_PS = "gsave /sfm2 { findfont exch makefont setfont } bind def [15 0 0 -15 0 0] /Courier-Bold sfm2";
	const string LEGEND_TEXT_CLOSE_PS = "grestore";

	// SVG legend text constants.
	const int LEGEND_TEXT_SIZE_SVG = 19;
	const string LEGEND_TEXT_OPEN_SVG = "<g font-family=\"monospace\" font-size=\"19\">";
	const string LEGEND_TEXT_CLOSE_SVG = "</g>";

	// Miscellaneous properties.
	const int LABEL_LINE_LENGTH = 10;
	const int PLOT_STRETCH = 4;
};

// Namespace usage declarations.
using namespace std;
using namespace DotPlotConstants;

/*
 * The DotPlotHandler class.
 */
class DotPlotHandler {
 public:
	// Public methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes private variables.
	 * Arguments:
	 *     1. text
	 *        The plot description.
	 *     2. size
	 *        The plot size.
	 *     3. triangular
	 *        True if the plot is triangular, false if not.
	 */
	DotPlotHandler( string text, int size, bool triangular = true );

	/*
	 * Name:        addDotValue
	 * Description: Add a dot to this handler.
	 * Arguments:
	 *     1. i
	 *        The i index of the dot.
	 *     2. j
	 *        The j index of the dot.
	 *     3. value
	 *        The dot value.
	 */
	void addDotValue( int i, int j, double value );

	/*
	 * Name:        getDotData
	 * Description: Get a particular dot's data: its location and color.
	 * Arguments:
	 *     1. i
	 *        The i index of the dot.
	 *     2. j
	 *        The j index of the dot.
	 * Returns:
	 *     The drawing data for the dot.
	 */
	string getDotData( int i, int j );

	/*
	 * Name:        setDescription
	 * Description: Set the plot description.
	 * Arguments:
	 *     1. newDescription
	 *        The description to set for the plot.
	 */
	void setDescription( string newDescription );

	/*
	 * Name:        setLegend
	 * Description: Set the legend properties.
	 * Arguments:
	 *     1. entries
	 *        The number of entries (colors) in the legend.
	 */
	void setLegend( int entries );

	/*
	 * Name:        setLegendDivider
	 * Description: Set the legend divider.
	 * Arguments:
	 *     1. text
	 *        The divider inside each legend entry.
	 */
	void setLegendDivider( string text );

	/*
	 * Name:        setLegendMaximum
	 * Description: Set the legend maximum value.
	 * Arguments:
	 *     1. maximum
	 *        The new legend maximum.
	 */
	void setLegendMaximum( double maximum );

	/*
	 * Name:        setLegendMinimum
	 * Description: Set the legend minimum value.
	 * Arguments:
	 *     1. minimum
	 *        The new legend minimum.
	 */
	void setLegendMinimum( double minimum );

	/*
	 * Name:        toString
	 * Description: Return a string representation of the dot plot handled by
	 *              this object. This is usually a Java convention, not a C++
	 *              one, but it is useful to get a snapshot of this plot
	 *              and access all information about it with a single easy call.
	 * Arguments:
	 *     1. includeDots
	 *        True if this string should include dot data, false if not.
	 *        In general, dot data isn't necessary because plots can get big.
	 *        Default is false.
	 * Returns:
	 *     A representation of the image as a string.
	 */
	string toString( bool includeDots = false );

	/*
	 * Name:        writePostscriptImage
	 * Description: Write the plot as a Postscript image.
	 * Arguments:
	 *     1. file
	 *        The file to write.
	 */
	void writePostscriptImage( string file );

	/*
	 * Name:        writeSVGImage
	 * Description: Write the plot as a SVG image.
	 * Arguments:
	 *     1. file
	 *        The file to write.
	 */
	void writeSVGImage( string file );

	/*
	 * Name:        writeTextFile
	 * Description: Write the plot as a text file.
	 * Arguments:
	 *     1. file
	 *        The file to write.
	 */
	void writeTextFile( string file );

 private:
	// Private image writing workhorse function.

	/*
	 * Name:        writeImageFile
	 * Description: Write the plot as some type of image file.
	 * Arguments:
	 *     1. file
	 *        The file to write.
	 *     2. isSVG
	 *        True if the image is SVG, false if not.
	 */
	void writeImageFile( string file, bool isSVG );

 private:
	// Private variables.

	// The current maximum and minimum plot range.
	double currentMax;
	double currentMin;

	// The default maximum and minimum plot range.
	double defaultMax;
	double defaultMin;

	// The plot description.
	string description;

	// The vector of dot values.
	vector< vector<double> > dots;

	// The vector of grid line data.
	vector<string> gridData;

	// The legend divider and ranges.
	string legendDivider;
	vector<string> legendRanges;
	vector< pair<double,double> > legendRangesMap;

	// The maximum bounds of the plot.
	int maxX;
	int maxY;
};

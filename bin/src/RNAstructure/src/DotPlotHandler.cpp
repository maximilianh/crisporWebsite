/*
 * An implementation file for a class that holds dot plot data and can write those images to files if necessary.
 *
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 * Redone in 2012.
 * Written by Jessica S. Reuter.
 */

#include "DotPlotHandler.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
DotPlotHandler::DotPlotHandler( string text, int size, bool triangular ) {

	// Set the description.
	description = text;

	// Initialize the default minimum and maximum.
	defaultMin = numeric_limits<double>::infinity();
	defaultMax = defaultMin * -1;

	// Initialize the current minimum and maximum.
	currentMin = defaultMin;
	currentMax = defaultMax;

	// Initialize the dot plot data structure.
	for( int i = 1; i <= size; i++ ) {
		vector<double> row;
		for( int j = 1; j <= size; j++ ) {
			row.push_back( defaultMin );
		}
		dots.push_back( row );
	}

	// Write the grid lines.
	// Go through each possible index on the adjusted length.
	// Only write a grid line if the index marks the first nucleotide, last nucleotide, or every tenth nucleotide.
	// In this case, every 10th nucleotide is actually every 10th index multiplied by the plot stretch, to make room for dots.
	// Also, if the raw plot length is not divisible by 10, the last grid line divisible by 10 is not labeled.
	// This is to make sure that the final grid line, which shows the plot line, is clearly seen.
	int adjustedPlotLength = size * PLOT_STRETCH;
	int block = 10 * PLOT_STRETCH;
	int currentGridLine = 0;
	int labelGap = 5;
	int adjust = PLOT_STRETCH;
	for( int i = 1; i <= adjustedPlotLength; i++ ) {
		if( ( i == 1 ) || ( i == adjustedPlotLength ) || ( i % block == 0 ) ) {

			// Determine what the next grid line label number should be.
			// A label number of 0 means no label should be written.
			int label = currentGridLine * 10;
			if( i == 1 ) { label = 1; }
			else if( i == adjustedPlotLength ) { label = dots.size(); }
			else if( adjustedPlotLength - i < block ) { label = 0; }

			// Determine the adjustment away from the grid border for the grid line, if it needs to be something other than the default.
			if( i != 1 ) { adjust = i; }

			// Determine the coordinates of the horizontal grid line.
			int x1 = ( triangular ) ? BORDER + adjust : BORDER;
			int x2 = BORDER + adjustedPlotLength + LABEL_LINE_LENGTH;
			int y = BORDER + TEXTSIZE + LABEL_LINE_LENGTH + adjust;
			int hx = x2 + labelGap;
			int hy = y + ( GLYPHSIZE / 2 );

			// Build the horizontal grid line and add it to the grid data.
			stringstream horizontal( stringstream::in | stringstream::out );
			horizontal << x1 << " " << y << " " << x2 << " " << y;
			if( label != 0 ) { horizontal << " " << label << " " << hx << " " << hy; }
			gridData.push_back( horizontal.str() );

			// Determine the coordinates of the vertical grid line.
			int places =
				( label >= 10000 ) ? 5 :
				( label >= 1000 ) ? 4 :
				( label >= 100 ) ? 3 :
				( label >= 10 ) ? 2 : 1;
			int x = BORDER + adjust;
			int y1 = BORDER + TEXTSIZE + labelGap;
			int y2 = ( triangular ) ? y1 + adjust : y1 + adjustedPlotLength;
			y2 += ( LABEL_LINE_LENGTH - labelGap );
			int vx = x - ( (GLYPHSIZE * places) / 2 );
			int vy = BORDER + TEXTSIZE;

			// Build the vertical grid line and add it to the grid data.
			stringstream vertical( stringstream::in | stringstream::out );
			vertical << x << " " << y1 << " " << x << " " << y2;
			if( label != 0 ) { vertical << " " << label << " " << vx << " " << vy; }
			gridData.push_back( vertical.str() );

			// Increment the grid line.
			currentGridLine++;
		}

		// If the last grid line was written and the plot is triangular, write a diagonal line to complete the plot border.
		if( triangular && ( i == adjustedPlotLength ) ) {
			adjust = PLOT_STRETCH;
			int x1 = BORDER + adjust;
			int y1 = BORDER + TEXTSIZE + LABEL_LINE_LENGTH + adjust;
			int x2 = ( x1 + adjustedPlotLength ) - adjust;
			int y2 = ( y1 + adjustedPlotLength ) - adjust;
			stringstream diagonal( stringstream::in | stringstream::out );
			diagonal << x1 << " " << y1 << " " << x2 << " " << y2;
			gridData.push_back( diagonal.str() );
		}
	}

	// Determine the maximum bounds of the image.
	int places =
		( size >= 10000 ) ? 5 :
		( size >= 1000 ) ? 4 :
		( size >= 100 ) ? 3 :
		( size >= 10 ) ? 2 : 1;
	maxX = BORDER + adjustedPlotLength + LABEL_LINE_LENGTH + labelGap + ( places * GLYPHSIZE ) + BORDER;
	maxY = BORDER + TEXTSIZE + labelGap + LABEL_LINE_LENGTH + adjustedPlotLength + BORDER;
}

///////////////////////////////////////////////////////////////////////////////
// Add a dot plot value at a specific place.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::addDotValue( int i, int j, double value ) {

	// Set the dot value.
	dots[j-1][i-1] = value;

	// Adjust bounds using this value if necessary.
	if( value != numeric_limits<double>::infinity() ) {
		if( value < defaultMin ) {
			defaultMin = value;
			currentMin = defaultMin;
		}

		if( value > defaultMax ) {
			defaultMax = value;
			currentMax = defaultMax;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
// Get information about a dot from a specific place.
///////////////////////////////////////////////////////////////////////////////
string DotPlotHandler::getDotData( int i, int j ) {

	// Get the dot value.
	double value = dots[j-1][i-1];

	// If the value is outside the current range, return.
	double epsilon = numeric_limits<double>::epsilon();
	bool min1 = ( currentMin <= value );
	bool min2 = ( fabs( currentMin - value ) < epsilon );
	bool max1 = ( currentMax >= value );
	bool max2 = ( fabs( currentMax - value ) < epsilon );
	if( (( min1 || min2 ) && ( max1 || max2 )) == false ) { return ""; }

	// If the value is valid, determine its data.
	if( value != numeric_limits<double>::infinity() ) {

		// Using the dot value, determine its color.
		// Set the dot color to red as the default.
		// If the legend has multiple entries, determine which area the dot falls in.
		double red = 1.0, green = 0.0, blue = 0.0;
		if( legendRangesMap.size() > 1 ) {
			for( int k = 1; k <= legendRangesMap.size(); k++ ) {

				// Get the range string and the range pair data for this entry.
				string entry = legendRanges[k-1];
				double low = legendRangesMap[k-1].first;
				double high = legendRangesMap[k-1].second;

				// If the dot value is in this entry's range, set the color.
				if( ( low <= value ) && ( value <= high ) ) {
					stringstream entryStream( legendRanges[k-1] );
					string temp;
					for( int i = 1; i <= 2; i++ ) { entryStream >> temp;}
					entryStream >> red;
					entryStream >> green;
					entryStream >> blue;
				}
			}
		}

		// Determine the X and Y coordinates of the dot.
		double x = BORDER + (j*4)-1;
		double y = BORDER + TEXTSIZE + LABEL_LINE_LENGTH + (i*4)-1;

		// Create a string that holds the dot data, and return it.
		stringstream valueStream( stringstream::in | stringstream::out );
		valueStream << fixed << setprecision( numeric_limits<double>::digits10 ) << value;
		stringstream dotStream( stringstream::in | stringstream::out );
		dotStream << x << " " << y << " " << red << " " << green << " " << blue;
		return dotStream.str();
	}

	// Otherwise, return an empty string.
	return "";
}

///////////////////////////////////////////////////////////////////////////////
// Set a special plot description.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::setDescription( string newDescription ) {
	description = newDescription;
}

///////////////////////////////////////////////////////////////////////////////
// Set values in the dot plot legend.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::setLegend( int entries ) {

	// If the number of entries is out of the acceptable range, return.
	if( entries < ENTRIES_MINIMUM ) { return; }
	if( entries > ENTRIES_MAXIMUM ) { return; }

	// Clear the legend vector and map.
	legendRanges.clear();
	legendRangesMap.clear();

	// Calculate the range each legend entry covers.
	// Also, determine if there is an odd or even number of entries.
	double increment = ( currentMax - currentMin ) / (double)entries;
	bool odd = ( entries % 2 == 1 );
	double dynamics = ( odd ) ? ( entries - 3 ) / 2 : ( entries - 2 ) / 2;
	double percent = 1 / ( dynamics + 1 );

	// Add the legend entries.
	for( int i = 1; i <= entries; i++ ) {

		// Determine the amount of red in this entry's color, from 1 (complete color) to 0 (no color).
		double red = 1.0 - ( (i-1) * percent );
		if( red < 0.0 ) { red = 0.0; }

		// Determine the amount of blue in this entry's color, from 1 (complete color) to 0 (no color).
		double blue = 1.0 - ( (entries - i) * percent );
		if( blue < 0.0 ) { blue = 0.0; }

		// Determine the amount of green in this entry's color, from 1 (complete color) to 0 (no color).
		double green = 1.0;
		if( red > 0.0 ) { green -= red; }
		if( blue > 0.0 ) { green -= blue; }

		// Determine the lower bound of this entry.
		double low = currentMin;
		low += ( increment * ( i - 1 ) );

		// Determine the upper bound of this entry.
		double high = currentMin;
		high += ( increment * i );
		if( i == entries ) { high = currentMax; }

		// Add the legend entry to the legend vector.
		stringstream entryStream( stringstream::in | stringstream::out );
		if( legendDivider == DIVIDER_ENERGY ) {
			entryStream << fixed << setprecision( 1 );
		}
		if (legendDivider == DIVIDER_PROBABILITY) {
			// LOGP_ZERO rounds numbers to 0 if they are within  +/-LOGP_EPSILON (e.g. 1E-12) of zero. This helps avoid differences in output due to small precision differences between platforms.
			entryStream //<< setprecision(4) 
			            << LOGP_ZERO(low) << " " << LOGP_ZERO(high);
		} else
			entryStream << low << " " << high;
		entryStream << fixed << setprecision( 2 ) << " " << red << " " << green << " " << blue;
		legendRanges.push_back( entryStream.str() );

		// Add the legend entry's range to the legend map.
		legendRangesMap.push_back( make_pair( low, high ) );
	}
}

///////////////////////////////////////////////////////////////////////////////
// Set the legend divider.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::setLegendDivider( string text ) {
	legendDivider = text;
}

///////////////////////////////////////////////////////////////////////////////
// Set the maximum value in the dot plot.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::setLegendMaximum( double newMaximum ) {

	// If the current maximum is out of range, return.
	if( newMaximum < defaultMin ) { return; }
	if( newMaximum > defaultMax ) { return; }
	if( newMaximum < currentMin ) { return; }

	// Set the maximum.
	currentMax = newMaximum;
}

///////////////////////////////////////////////////////////////////////////////
// Set the minimum value in the dot plot.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::setLegendMinimum( double newMinimum ) {

	// If the current minimum is out of range, return.
	if( newMinimum < defaultMin ) { return; }
	if( newMinimum > defaultMax ) { return; }
	if( newMinimum > currentMax ) { return; }

	// Set the minimum.
	currentMin = newMinimum;
}

///////////////////////////////////////////////////////////////////////////////
// Get a string representation of the dot plot.
///////////////////////////////////////////////////////////////////////////////
string DotPlotHandler::toString( bool includeDots ) {

	// Create the string stream.
	stringstream stream( stringstream::in | stringstream::out );

	// Show the maximum plot bounds.
	stream << "Max Bounds: " << "(" << maxX << "," << maxY << ")" << endl;

	// Show the description.
	stream << "Description: " << description << endl;

	// Show the legend.
	stream << "Legend: " << legendDivider << endl;
	for( unsigned int i = 1; i <= legendRanges.size(); i++ ) {

		// Get the legend data.
		double low, high, red, green, blue;
		stringstream legendStream( legendRanges[i-1] );
		legendStream >> low;
		legendStream >> high;
		legendStream >> red;
		legendStream >> green;
		legendStream >> blue;

		// Convert the colors to RGB.
		red = floor( red * 255.0 );
		green = floor( green * 255.0 );
		blue = floor( blue * 255.0 );

		// Build the legend entry and add it to the string.
		stringstream buildStream( stringstream::in | stringstream::out );
		string divider2 = ( i == legendRanges.size() ) ? " <= " : " <  ";
		buildStream << low << " <= " << legendDivider
		            << divider2 << high << " "
		            << red << " " << green << " " << blue;
		stream << buildStream.str() << endl;
	}

	// Show the grid lines.
	stream << "Grid lines:" << endl;
	for( unsigned int i = 1; i <= gridData.size(); i++ ) {
		stream << gridData[i-1] << endl;
	}

	// If dots should be included in this string, do so.
	if( includeDots == true ) {
		stream << "Dots:" << endl;
		int size = dots.size();
		for( int i = 1; i <= size; i++ ) {
			for( int j = 1; j <= size; j++ ) {
				string dotData = getDotData( i, j );
				if( dotData != "" ) {
					stream << "(" << i << "," << j << "): " << dotData << endl;
				}
			}
		}
	}

	// Return the string from the stream.
	return stream.str();
}

///////////////////////////////////////////////////////////////////////////////
// Write a dot plot image file.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::writeImageFile( string file, bool isSVG ) {

	// Initialize an index used to handle conversion template indices later.
	size_t index = string::npos;

	// Set black as the grid and label color.
	string gridColor = getColorString( BLACK, isSVG );

	// Determine the dot size.
	stringstream dotSizeStream( stringstream::in | stringstream::out );
	dotSizeStream << PLOT_STRETCH - 2;
	string dotSizeString = dotSizeStream.str();

	// Determine the dimension by which the image should be scaled.
	double absoluteXBound = ( !isSVG ) ? XBOUND_PS : XBOUND_SVG;
	double absoluteYBound = ( !isSVG ) ? YBOUND_PS : YBOUND_SVG;
	double scaleX = absoluteXBound / maxX;
	double scaleY = absoluteYBound / ( maxY - ( 2 * TEXTSIZE_BIGGER ) );
	double scale = min( scaleX, scaleY );

	// Open the output file and write the start of the image.
	ofstream out( file.c_str() );
	string startMarker = ( !isSVG ) ? createStartPS() : createStartSVG();
	out << startMarker << endl;

	// Open the scaled area for the plot.
	string scaleStart = ( !isSVG ) ? SCALE_OPEN_PS : SCALE_OPEN_SVG;
	string scaleString;
	stringstream scaleStream( stringstream::in | stringstream::out );
	scaleStream << scale;
	scaleStream >> scaleString;
	while( ( index = scaleStart.find( SCALEFACTOR ) ) != string::npos ) { scaleStart = scaleStart.replace( index, SCALEFACTOR.size(), scaleString ); }
	out << scaleStart << endl;

	// Write the grid lines.
	// Size the grid line text the same as the legend text.
	string gridResizeOpen = ( !isSVG ) ? LEGEND_TEXT_OPEN_PS : LEGEND_TEXT_OPEN_SVG;
	out << gridResizeOpen << endl;
	for( int i = 1; i <= gridData.size(); i++ ) {

		// Get the data from the next grid line.
		string x1, y1, x2, y2, number = "", numberX, numberY;
		stringstream gridStream( gridData[i-1] );
		gridStream >> x1;
		gridStream >> y1;
		gridStream >> x2;
		gridStream >> y2;
		if( gridStream >> number ) {
			gridStream >> numberX;
			gridStream >> numberY;
		}

		// Write the next line which makes up the grid, and its number if applicable.
		string line = ( !isSVG ) ? LINE_PS : LINE_SVG;
		while( ( index = line.find( COLOR ) ) != string::npos ) { line = line.replace( index, COLOR.size(), gridColor ); }
		while( ( index = line.find( LINEWEIGHT ) ) != string::npos ) { line = line.replace( index, LINEWEIGHT.size(), "1" ); }
		while( ( index = line.find( STARTX ) ) != string::npos ) { line = line.replace( index, STARTX.size(), x1 ); }
		while( ( index = line.find( STARTY ) ) != string::npos ) { line = line.replace( index, STARTY.size(), y1 ); }
		while( ( index = line.find( ENDX ) ) != string::npos ) { line = line.replace( index, ENDX.size(), x2 ); }
		while( ( index = line.find( ENDY ) ) != string::npos ) { line = line.replace( index, ENDY.size(), y2 ); }
		out << line << endl;
		if( number != "" ) {
			string text = ( !isSVG ) ? TEXT_PS : TEXT_SVG;
			while( ( index = text.find( COLOR ) ) != string::npos ) { text = text.replace( index, COLOR.size(), gridColor ); }
			while( ( index = text.find( LOCX ) ) != string::npos ) { text = text.replace( index, LOCX.size(), numberX ); }
			while( ( index = text.find( LOCY ) ) != string::npos ) { text = text.replace( index, LOCY.size(), numberY ); }
			while( ( index = text.find( TEXTSTRING ) ) != string::npos ) { text = text.replace( index, TEXTSTRING.size(), number ); }
			out << text << endl;
		}
	}
	string gridResizeClose = ( !isSVG ) ? LEGEND_TEXT_CLOSE_PS : LEGEND_TEXT_CLOSE_SVG;
	out << gridResizeClose << endl;

	// Write the dots.
	int size = dots.size();
	for( int i = 1; i <= size; i++ ) {
		for( int j = 1; j <= size; j++ ) {

			// Get the data for the next dot.
			string dotData = getDotData( i, j );
			if( dotData != "" ) {
				stringstream dotDataStream( dotData );

				// Get the X and Y coordinates of the dot.
				string x, y;
				dotDataStream >> x;
				dotDataStream >> y;

				// Check if the dot is in the proper range to draw.
				double epsilon = numeric_limits<double>::epsilon();
				double value = dots[j-1][i-1];
				bool min1 = ( currentMin <= value );
				bool min2 = ( fabs( currentMin - value ) < epsilon );
				bool max1 = ( currentMax >= value );
				bool max2 = ( fabs( currentMax - value ) < epsilon );

				// If the value is in the proper range, do the dot drawing.
				if( ( min1 || min2 ) && ( max1 || max2 ) ) {

					// Determine the color of the dot.
					string red, green, blue;
					dotDataStream >> red;
					if( isSVG ) {
						double redVal;
						stringstream redStream( stringstream::in | stringstream::out );
						redStream << red;
						redStream >> redVal;
						redVal *= 255;
						stringstream redStream2( stringstream::in | stringstream::out );
						redStream2 << fixed << setprecision( 0 ) << redVal;
						red = redStream2.str();
					}
					dotDataStream >> green;
					if( isSVG ) {
						double greenVal;
						stringstream greenStream( stringstream::in | stringstream::out );
						greenStream << green;
						greenStream >> greenVal;
						greenVal *= 255;
						stringstream greenStream2( stringstream::in | stringstream::out );
						greenStream2 << fixed << setprecision( 0 ) << greenVal;
						green = greenStream2.str();
					}
					dotDataStream >> blue;
					if( isSVG ) {
						double blueVal;
						stringstream blueStream( stringstream::in | stringstream::out );
						blueStream << blue;
						blueStream >> blueVal;
						blueVal *= 255;
						stringstream blueStream2( stringstream::in | stringstream::out );
						blueStream2 << fixed << setprecision( 0 ) << blueVal;
						blue = blueStream2.str();
					}

					// Determine the next dot color.
					string color = ( !isSVG ) ? COLOR_TEMPLATE_PS : COLOR_TEMPLATE_SVG;
					color = color.replace( color.find( "RED" ), 3, red );
					color = color.replace( color.find( "GREEN" ), 5, green );
					color = color.replace( color.find( "BLUE" ), 4, blue );

					// Draw the next dot.
					string rectData = ( !isSVG ) ? RECTANGLE_PS : RECTANGLE_SVG;
					while( ( index = rectData.find( COLOR ) ) != string::npos ) { rectData = rectData.replace( index, COLOR.size(), color ); }
					while( ( index = rectData.find( LOCX ) ) != string::npos ) { rectData = rectData.replace( index, LOCX.size(), x ); }
					while( ( index = rectData.find( LOCY ) ) != string::npos ) { rectData = rectData.replace( index, LOCY.size(), y ); }
					while( ( index = rectData.find( WIDTH ) ) != string::npos ) { rectData = rectData.replace( index, WIDTH.size(), dotSizeString ); }
					while( ( index = rectData.find( HEIGHT ) ) != string::npos ) { rectData = rectData.replace( index, HEIGHT.size(), dotSizeString ); }
					out << rectData << endl;
				}
			}
		}
	}

	// Close the scaled area for the plot.
	string scaleClose = ( !isSVG ) ? SCALE_CLOSE_PS : SCALE_CLOSE_SVG;
	out << scaleClose << endl;

	// Open the resizing group for text.
	// Also determine what the size of the resized text is.
	int resizedSize = ( !isSVG ) ? LEGEND_TEXT_SIZE_PS : LEGEND_TEXT_SIZE_SVG;
	string resizeOpen = ( !isSVG ) ? LEGEND_TEXT_OPEN_PS : LEGEND_TEXT_OPEN_SVG;
	out << resizeOpen << endl;

	// Write the legend into the image.
	unsigned int legendPlusOne = legendRanges.size() + 1;
	for( unsigned int i = 1; i <= legendRanges.size(); i++ ) {

		// Determine the Y location for this entry string.
		unsigned int next = i - 1;
		if( currentMin == currentMax ) { next = legendRanges.size() - 1; }
		double location = absoluteYBound - ( resizedSize * (legendPlusOne - next) );
		stringstream entryStream( stringstream::in | stringstream::out );
		entryStream << location;
		string entryYString = entryStream.str();

		// Pull out the information for this particular entry.
		// Note that colors are specified on a scale of 0 (no color) to 1 (complete color).
		// While this is OK for Postscript, SVG requires some conversion, and other formats may require it in the future too.
		string low, high, red, green, blue;
		stringstream info( legendRanges[i-1] );
		info >> low;
		info >> high;
		info >> red;
		if( isSVG ) {
			double redVal;
			stringstream redStream( stringstream::in | stringstream::out );
			redStream << red;
			redStream >> redVal;
			redVal *= 255;
			stringstream redStream2( stringstream::in | stringstream::out );
			redStream2 << fixed << setprecision( 0 ) << redVal;
			red = redStream2.str();
		}
		info >> green;
		if( isSVG ) {
			double greenVal;
			stringstream greenStream( stringstream::in | stringstream::out );
			greenStream << green;
			greenStream >> greenVal;
			greenVal *= 255;
			stringstream greenStream2( stringstream::in | stringstream::out );
			greenStream2 << fixed << setprecision( 0 ) << greenVal;
			green = greenStream2.str();
		}
		info >> blue;
		if( isSVG ) {
			double blueVal;
			stringstream blueStream( stringstream::in | stringstream::out );
			blueStream << blue;
			blueStream >> blueVal;
			blueVal *= 255;
			stringstream blueStream2( stringstream::in | stringstream::out );
			blueStream2 << fixed << setprecision( 0 ) << blueVal;
			blue = blueStream2.str();
		}

		// Determine the next entry string.
		string entryDividerStart = ( !isSVG ) ? " <= " : " &#60;= ";
		string entryDividerEnd = ( !isSVG ) ? " <  " : " &#60;  ";
		if( i == legendRanges.size() ) { entryDividerEnd = entryDividerStart; }
		stringstream textStream( stringstream::in | stringstream::out );
		if( currentMin < currentMax ) { textStream << low << entryDividerStart << legendDivider << entryDividerEnd << high; }
		else { textStream << legendDivider << " = " << currentMin; }
		string text = textStream.str();

		// Determine the next entry color.
		string color = ( !isSVG ) ? COLOR_TEMPLATE_PS : COLOR_TEMPLATE_SVG;
		color = color.replace( color.find( "RED" ), 3, red );
		color = color.replace( color.find( "GREEN" ), 5, green );
		color = color.replace( color.find( "BLUE" ), 4, blue );

		// Build the next legend entry string and write it to the file.
		string entry = ( !isSVG ) ? TEXT_PS : TEXT_SVG;
		while( ( index = entry.find( COLOR ) ) != string::npos ) { entry = entry.replace( index, COLOR.size(), color ); }
		while( ( index = entry.find( LOCX ) ) != string::npos ) { entry = entry.replace( index, LOCX.size(), BORDER_STRING ); }
		while( ( index = entry.find( LOCY ) ) != string::npos ) { entry = entry.replace( index, LOCY.size(), entryYString ); }
		while( ( index = entry.find( TEXTSTRING ) ) != string::npos ) { entry = entry.replace( index, TEXTSTRING.size(), text ); }
		out << entry << endl;

		// If only one specific entry should be written, break out of the loop.
		if( currentMin == currentMax ) { break; }
	}

	// Write the image description.
	if( !description.empty() ) {

		// Create the description Y location string.
		string descYString;
		stringstream descStream( stringstream::in | stringstream::out );
		descStream << ( absoluteYBound - resizedSize + 2 );
		descStream >> descYString;

		// Create the description string and write it in the file.
		string descData = ( !isSVG ) ? TEXT_PS : TEXT_SVG;
		string descText = ( isSVG ) ? description : LegendDescriptionSettings::escapeBackSlashes(description); // for PS output, we have to escape backslashes (e.g. in windows file paths)
		while( ( index = descData.find( COLOR ) ) != string::npos ) { descData = descData.replace( index, COLOR.size(), getColorString( BLACK, isSVG ) ); }
		while( ( index = descData.find( LOCX ) ) != string::npos ) { descData = descData.replace( index, LOCX.size(), BORDER_STRING ); }
		while( ( index = descData.find( LOCY ) ) != string::npos ) { descData = descData.replace( index, LOCY.size(), descYString ); }
		while( ( index = descData.find( TEXTSTRING ) ) != string::npos ) { descData = descData.replace( index, TEXTSTRING.size(), descText ); }
		out << descData << endl;
	}

	// Close the resizing group for text.
	string resizeClose = ( !isSVG ) ? LEGEND_TEXT_CLOSE_PS : LEGEND_TEXT_CLOSE_SVG;
	out << resizeClose << endl;

	// Finish the file.
	string endMarker = ( !isSVG ) ? END_MARKER_PS : END_MARKER_SVG;
	out << endMarker << endl;

	// Close the output stream.
	out.close();
}

///////////////////////////////////////////////////////////////////////////////
// Write a Postscript dot plot image.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::writePostscriptImage( string file ) {
	writeImageFile( file, false );
}

///////////////////////////////////////////////////////////////////////////////
// Write an SVG dot plot image.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::writeSVGImage( string file ) {
	writeImageFile( file, true );
}

///////////////////////////////////////////////////////////////////////////////
// Write a dot plot text file.
///////////////////////////////////////////////////////////////////////////////
void DotPlotHandler::writeTextFile( string file ) {

	// Open the output stream to the text file and write the header.
	ofstream out( file.c_str() );
	int length = dots.size();
	out << length << endl << "i\tj\t" << legendDivider << endl;

	// Loop through each possible pair/dot.
	for( int i = 1; i <= length; i++ ) {
		for( int j = 1; j <= length; j++ ) {

			// If the value is within range, write the dot.
			double value = dots[j-1][i-1];
			double epsilon = numeric_limits<double>::epsilon();
			bool min1 = ( currentMin <= value );
			bool min2 = ( fabs( currentMin - value ) < epsilon );
			bool max1 = ( currentMax >= value );
			bool max2 = ( fabs( currentMax - value ) < epsilon );
			if( ( min1 || min2 ) && ( max1 || max2 ) ) {
				out << i << "\t" << j << "\t" << value << endl;
			}
		}
	}

	// Close the written text file.
	out.close();

}

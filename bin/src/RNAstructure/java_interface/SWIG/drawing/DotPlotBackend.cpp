/*
 * An implementation file for a class that holds dot plot methods for
 * Java drawing.
 *
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter.
 */

#include "DotPlotBackend.h"
#include <iostream>

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
DotPlotBackend::DotPlotBackend() {

	// Initialize the default minimum and maximum.
	defaultMin = numeric_limits<double>::infinity();
	defaultMax = defaultMin * -1;

	// Initialize the current minimum and maximum.
	currentMin = defaultMin;
	currentMax = defaultMax;

	// Initialize the other data members.
	handler = 0;
	length = 0;
	plotColors = ENTRIES_DEFAULT;
}

///////////////////////////////////////////////////////////////////////////////
// Destructor.
///////////////////////////////////////////////////////////////////////////////
DotPlotBackend::~DotPlotBackend() {

	if( handler != 0 ) { delete handler; }
}

///////////////////////////////////////////////////////////////////////////////
// Get the dot plot bounds.
///////////////////////////////////////////////////////////////////////////////
string DotPlotBackend::getBounds() {

	// Get the handler's data as a string.
	string data = handler->toString();

	// Determine the x bounds and shave off the borders, then add a 5 pixel
	// border on each side.
	size_t xBound = data.find_first_of( "(" ) + 1;
	size_t xLength = data.find_first_of( "," ) - xBound;
	int x = 0;
	stringstream xStream( data.substr( xBound, xLength ) );
	xStream >> x;
	x -= ( BORDER * 2 );
	x += 10;

	// Determine the y bounds and shave off the borders, then add a 5 pixel
	// border on each side.
	size_t yBound = data.find_first_of( "," ) + 1;
	size_t yLength = data.find_first_of( ")" ) - yBound;
	int y = 0;
	stringstream yStream( data.substr( yBound, yLength ) );
	yStream >> y;
	y -= ( BORDER * 2 );
	y += 10;

	// Build and return the bounds string.
	stringstream boundsStream( stringstream::in | stringstream::out );
	boundsStream << x << " " << y;
	return boundsStream.str();
}

///////////////////////////////////////////////////////////////////////////////
// Get the colors allowed in the plot.
///////////////////////////////////////////////////////////////////////////////
string DotPlotBackend::getColors() {

	stringstream stream( stringstream::in | stringstream::out );
	stream << plotColors << " " << ENTRIES_MINIMUM << " " << ENTRIES_MAXIMUM;
	return stream.str();
}

///////////////////////////////////////////////////////////////////////////////
// Get the current maximum plot value.
///////////////////////////////////////////////////////////////////////////////
double DotPlotBackend::getCurrentMax() {

	return currentMax;
}

///////////////////////////////////////////////////////////////////////////////
// Get the current minimum plot value.
///////////////////////////////////////////////////////////////////////////////
double DotPlotBackend::getCurrentMin() {

	return currentMin;
}

///////////////////////////////////////////////////////////////////////////////
// Get the default minimum plot value.
///////////////////////////////////////////////////////////////////////////////
double DotPlotBackend::getDefaultMin() {

	return defaultMin;
}

///////////////////////////////////////////////////////////////////////////////
// Get the default maximum plot value.
///////////////////////////////////////////////////////////////////////////////
double DotPlotBackend::getDefaultMax() {

	return defaultMax;
}

///////////////////////////////////////////////////////////////////////////////
// Get information about a dot from a specific place.
///////////////////////////////////////////////////////////////////////////////
string DotPlotBackend::getDotData( int i, int j ) {

	// Get the data directly from the handler.
	// If it's empty, return.
	string dataString = handler->getDotData( i, j );
	if( dataString == "" ) { return ""; }

	// Split the data into its pieces.
	double x, y, red, green, blue;
	stringstream data( dataString );
	data >> x;
	data >> y;
	data >> red;
	data >> green;
	data >> blue;

	// Shave off the border, then add a 10 pixel border.
	x = ( x - BORDER ) + 10;
	y = ( y - BORDER ) + 10;

	// Convert the colors to rgb.
	red *= 255;
	red = floor( red );
	green *= 255;
	green = floor( green );
	blue *= 255;
	blue = floor( blue );

	// Rebuild the data string and return it.
	stringstream rebuilt( stringstream::in | stringstream::out );
	rebuilt << x << " " << y << " "
	        << red << " " << green << " " << blue;
	return rebuilt.str();
}

///////////////////////////////////////////////////////////////////////////////
// Get information about a dot from a specific place.
///////////////////////////////////////////////////////////////////////////////
string DotPlotBackend::getDotDataByLocation( int x, int y, double scale ) {

	// Determine where the dot would be if scaled to 100%.
	double xScaled = (double)x / scale;
	double yScaled = (double)y / scale;

	// Go through each dot to see if it could be where the user clicked.
	// If a dot is identified, return its data.
	for( int i = 1; i <= length; i++ ) {
		for( int j = 1; j <= length; j++ ) {

			// Determine the X and Y coordinates of the dot.
			double dotX = BORDER + (j*4)-1;
			double dotY = BORDER + TEXTSIZE + LABEL_LINE_LENGTH + (i*4)-1;

			// Shave off the border, then add a 10 pixel border.
			dotX = ( dotX - BORDER ) + 10;
			dotY = ( dotY - BORDER ) + 10;

			// Check to make sure the point is inside the dot; if it is, return
			// the dot's data.
			if( ( dotX <= xScaled ) && ( xScaled <= ( dotX + 3 ) ) ) {
				if( ( dotY <= yScaled ) && ( yScaled <= ( dotY + 3 ) ) ) {

					// Build the dot pair to get the dot value.
					stringstream pair( stringstream::in | stringstream::out );
					pair << i << "-" << j;

					// If the value is outside the default range, record this fact.
					double value = values[pair.str()];
					double epsilon = numeric_limits<double>::epsilon();
					bool min1 = ( defaultMin <= value );
					bool min2 = ( fabs( defaultMin - value ) < epsilon );
					bool max1 = ( defaultMax >= value );
					bool max2 = ( fabs( defaultMax - value ) < epsilon );
					bool outOfRange = 
						( (( min1 || min2 ) && ( max1 || max2 )) == false );

					// Build the data string.
					stringstream dotOut( stringstream::in | stringstream::out );
					dotOut << i << "(" << sequence[i-1] << ") -- "
					       << j << "(" << sequence[j-1] << "): ";
					if( outOfRange == true ) { dotOut << "infinity"; }
					else { dotOut << values[pair.str()]; }

					// Return the data string.
					return dotOut.str();
				}
			}
		}
	}

	// Return an empty string if no dot was identified.
	return "";
}

///////////////////////////////////////////////////////////////////////////////
// Get a grid line.
///////////////////////////////////////////////////////////////////////////////
string DotPlotBackend::getGridLine( int number ) {

	// Get the handler's data as a string.
	string fullData = handler->toString();
	size_t gridStart = fullData.find( "Grid lines:" );
	size_t gridEnd = fullData.find( "Legend:" );
	size_t substringLength = gridEnd - gridStart;
	string data = fullData.substr( gridStart, substringLength );

	// Pull out the grid line.
	int counter = 0;
	istringstream stream( data );
	string line;
	while( getline( stream, line ) ) {
		if( counter == number ) {

			// Get raw data from the grid line.
			stringstream lineStream( line );
			double x1, y1, x2, y2, labelX, labelY;
			int label = 0;
			lineStream >> x1;
			lineStream >> y1;
			lineStream >> x2;
			lineStream >> y2;
			if( lineStream >> label ) {
				lineStream >> labelX;
				lineStream >> labelY;
			}

			// Shave off the border, then add a 10 pixel border.
			x1 = ( x1 - BORDER ) + 10;
			y1 = ( y1 - BORDER ) + 10;
			x2 = ( x2 - BORDER ) + 10;
			y2 = ( y2 - BORDER ) + 10;
			if( label > 0 ) {
				labelX = ( labelX - BORDER ) + 10;
				labelY = ( labelY - BORDER ) + 10;
			}

			// Return the rebuilt grid line.
			stringstream rebuilt( stringstream::in | stringstream::out );
			rebuilt << x1 << " " << y1 << " " << x2 << " " << y2 << " ";
			if( label > 0 ) {
				rebuilt << label << " " << labelX << " " << labelY;
			}
			return rebuilt.str();
		}
		counter++;
	}

	// Return an empty string if no gridline was found.
	return "";
}

///////////////////////////////////////////////////////////////////////////////
// Get the plot's legend data.
///////////////////////////////////////////////////////////////////////////////
string DotPlotBackend::getLegendData() {

	string fullData = handler->toString();
	size_t gridStart = fullData.find( "Legend:" );
	size_t gridEnd = fullData.find( "Grid lines:" );
	size_t substringLength = gridEnd - gridStart;
	string data = fullData.substr( gridStart, substringLength );
	data = data.substr( data.find_first_of( "\n" ) + 1 );
	data = data.substr( 0, data.find_last_of( "\n" ) );
	return data;
}

///////////////////////////////////////////////////////////////////////////////
// Get the dot plot length.
///////////////////////////////////////////////////////////////////////////////
int DotPlotBackend::getLength() {

	return length;
}

///////////////////////////////////////////////////////////////////////////////
// Read dot plot data.
///////////////////////////////////////////////////////////////////////////////
bool DotPlotBackend::readData( string file, int strand ) {
	// Check to make sure the file type is a valid dot plot type.
	// If not, return false.
	string extension = file.substr( file.find_last_of( "." ) + 1 );
	bool text = ( extension == "dp" );
	bool dynalignSave = ( extension == "dsv" );
	bool partitionSave = ( extension == "pfs" );
	bool foldingSave = ( extension == "sav" );
	if( text || dynalignSave || partitionSave || foldingSave ) { /* Valid */ }
	else { return false; }

	try {
	// If the data file is a dot plot text file, read it.
	if( text ) {

		// Open the text file stream to read data.
		// If an error occurred, close the file and return false.
		ifstream inFile( file.c_str() );
		if( !( inFile.is_open() ) ) {
			inFile.close();
			return false;
		}

		// Pull the plot length from the first line.
		// If an error occurred, close the file and return false.
		string firstLine;
		getline( inFile, firstLine );
		stringstream lengthStream( firstLine );
		lengthStream >> length;
		if( length <= 0 ) {
			inFile.close();
			return false;
		}

		// Pull the plot divider from the end of the second line.
		// This is the third column of that line.
		string divider = "";
		string dividerLine;
		getline( inFile, dividerLine );
		stringstream dividerStream( dividerLine );
		for( int i = 1; i <= 3; i++ ) { dividerStream >> divider; }

		// Create the dot plot handler.
		handler = new DotPlotHandler( file.c_str(), length );
		handler->setLegendDivider( divider );

		// For all the rest of the lines, pull out the dot data.
		string dataLine = "";
		while( true ) {

			// Get the next line.
			// If it's empty, break out of the loop.
			// Otherwise, build it into a stringstream to get data from.
			getline( inFile, dataLine );
			if( dataLine == "" ) { break; }
			stringstream dataStream( dataLine );

			// Read the i index from the line.
			// If an error occurred, close the file and return false.
			// The dot plot handler doesn't need to be deleted because that's
			// taken care of in the destructor.
			int i = 0;
			dataStream >> i;
			if( i <= 0 ) {
				inFile.close();
				return false;
			}

			// Read the j index from the line.
			// If an error occurred, close the file and return false.
			// The dot plot handler doesn't need to be deleted because that's
			// taken care of in the destructor.
			int j = 0;
			dataStream >> j;
			if( j <= 0 ) {
				inFile.close();
				return false;
			}

			// Read the dot value from the line.
			// If an error occurred, close the file and return false.
			// The dot plot handler doesn't need to be deleted because that's
			// taken care of in the destructor.
			double value = numeric_limits<double>::infinity();
			dataStream >> value;

			// Add the value to the handler.
			handler->addDotValue( i, j, value );
			stringstream val( stringstream::in | stringstream::out );
			val << i << "-" << j;
			values[val.str()] = value;
			if( value != numeric_limits<double>::infinity() ) {
				if( value < defaultMin ) {
					defaultMin = value;
					currentMin = defaultMin;
				} else if( value > defaultMax ) {
					defaultMax = value;
					currentMax = defaultMax;
				}
			}
		}

		// Since there's no way to know the sequence this data came from, fill
		// the sequence vector with "?".
		for( int i = 1; i <= length; i++ ) { sequence.push_back( "?" ); }

		// Close the file stream after data was read correctly.
		inFile.close();
	}

	// If the data file is a Dynalign save file, read it.
	else if( dynalignSave ) {

		// Create the Dynalign object that holds the dot plot data.
		Dynalign_object* object = new Dynalign_object( file.c_str() );
		ErrorChecker<Dynalign_object>* checker =
			new ErrorChecker<Dynalign_object>( object );
		string error = checker->returnError();

		// If an error occurred creating the Dynalign object, delete it and its
		// error checker, then return false.
		if( error != "" ) {
			delete checker;
			delete object;
			return false;
		}

		// Build the dot plot handler using the strand's sequence length.
		length =
			( strand == 1 ) ? object->GetRNA1()->GetSequenceLength() :
			object->GetRNA2()->GetSequenceLength();
		handler = new DotPlotHandler( file, length );

		// Add all possible dots to the dot plot. 
		for( int i = 1; i < length; i++ ) {
			for( int j = i+1; j <= length; j++ ) { //j must always be greater than i!   //i+1
				double energy = object->GetBestPairEnergy( strand, i, j );
				if( energy > 0.0 ) {
					energy = numeric_limits<double>::infinity();
				} else {
					stringstream stream( stringstream::in | stringstream::out );
					stream << fixed << setprecision( 1 ) << energy;
					stream >> energy;
				}
				handler->addDotValue( i, j, energy );
				stringstream val( stringstream::in | stringstream::out );
				val << i << "-" << j;
				values[val.str()] = energy;
				if( energy != numeric_limits<double>::infinity() ) {
					if( energy < defaultMin ) {
						defaultMin = energy;
						currentMin = defaultMin;
					} else if( energy > defaultMax ) {
						defaultMax = energy;
						currentMax = defaultMax;
					}
				}
			}
		}

		// Set the plot divider and legend.
		handler->setLegendDivider( DIVIDER_ENERGY );
		handler->setLegend( ENTRIES_DEFAULT );

		// Get the plotted strand and determine its sequence.
		RNA* rna = ( strand == 1 ) ? object->GetRNA1() : object->GetRNA2();
		for( int i = 1; i <= length; i++ ) {
			stringstream seqStream( stringstream::in | stringstream::out );
			seqStream << rna->GetNucleotide( i );
			sequence.push_back( seqStream.str() );
		}

		// Set the description based on the strand number.
		string seqDescription = rna->GetCommentString();
		string seqString = ( strand == 1 ) ? "1" : "2";
		string fullDescription =
			file + " -- Sequence " + seqString + ": " + seqDescription;
		handler->setDescription( fullDescription );

		// Delete the error checker and Dynalign object.
		delete checker;
		delete object;
	}

	// If the data file is a partition function save file, read it.
	else if( partitionSave ) {
		// Create the RNA object that holds the dot plot data.
		RNA* object = new RNA( file.c_str(), FILE_PFS );
		ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( object );
		string error = checker->returnError();
		// If an error occurred creating the RNA object, delete it and its
		// error checker, then return false.
		if( error != "" ) {
			delete checker;
			delete object;
			return false;
		}

		// Build the dot plot handler using the strand's sequence length.
		length = object->GetSequenceLength();
		handler = new DotPlotHandler( file, length );
		cout << "DotPlotBackend-Before Loop\n";
		// Add all possible dots to the dot plot.
		try {
		for( int i = 1; i < length; i++ ) {
			//for( int j = 1; j <= length; j++ ) {        //RMW: j should never be less than i.
			for( int j = i+1; j <= length; j++ ) {
			
				double probability = object->GetPairProbability( i, j );
				if( ( probability == 0.0 ) || ( probability == -0.0 ) ) {
					probability = numeric_limits<double>::infinity();
				} else {
					probability = -log10( probability );
				}
				if (i==length/2 && j==length/2)
					cout << "Dot plot mid-point\n";

				if (i==3*length/4 && j==3*length/4)
					cout << "Dot plot 3/4 point\n";

				if (i==length-2 && j==length-2)
					cout << "Dot plot -2-point\n";

				if (i==length)
					cout << "Dot plot i=length point\n";

				if (i==length-1 && j==length-1)
					cout << "Dot plot near-end-point\n";

				handler->addDotValue( i, j, probability );
				stringstream val( stringstream::in | stringstream::out );
				val << i << "-" << j;
				values[val.str()] = probability;
				if( probability != numeric_limits<double>::infinity() ) {
					if( probability < defaultMin ) {
						defaultMin = probability;
						currentMin = defaultMin;
					} else if( probability > defaultMax ) {
						defaultMax = probability;
						currentMax = defaultMax;
					}
				}
			}
		}
		} catch (...) {
			cout << "DOT PLOT ERROR" << endl;
		}
		cout << "DotPlotBackend-After Loop\n";
		// Set the plot divider and legend.
		handler->setLegendDivider( DIVIDER_PROBABILITY );
		handler->setLegend( ENTRIES_DEFAULT );

		// Determine the plotted strand's sequence.
		for( int i = 1; i <= length; i++ ) {
			stringstream seqStream( stringstream::in | stringstream::out );
			seqStream << object->GetNucleotide( i );
			sequence.push_back( seqStream.str() );
		}
		// Delete the error checker and RNA object.
		delete checker;
		delete object;
	}

	// If the data file is a folding save file, read it.
	else if( foldingSave ) {

		// Create the RNA object that holds the dot plot data.
		RNA* object = new RNA( file.c_str(), FILE_SAV );
		ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( object );
		string error = checker->returnError();

		// If an error occurred creating the RNA object, delete it and its
		// error checker, then return false.
		if( error != "" ) {
			delete checker;
			delete object;
			return false;
		}

		// Build the dot plot handler using the strand's sequence length.
		length = object->GetSequenceLength();
		handler = new DotPlotHandler( file, length );

		// Add all possible dots to the dot plot.
		for( int i = 1; i < length; i++ ) {
			for( int j = i+1; j <= length; j++ ) { //j must always be greater than i.
				double energy = object->GetPairEnergy( i, j );
				if( energy > 0.0 ) {
					energy = numeric_limits<double>::infinity();
				} else {
					stringstream stream( stringstream::in | stringstream::out );
					stream << fixed << setprecision( 1 ) << energy;
					stream >> energy;
				}
				handler->addDotValue( i, j, energy );
				stringstream val( stringstream::in | stringstream::out );
				val << i << "-" << j;
				values[val.str()] = energy;
				if( energy != numeric_limits<double>::infinity() ) {
					if( energy < defaultMin ) {
						defaultMin = energy;
						currentMin = defaultMin;
					} else if( energy > defaultMax ) {
						defaultMax = energy;
						currentMax = defaultMax;
					}
				}
			}
		}

		// Set the plot divider and legend.
		handler->setLegendDivider( DIVIDER_ENERGY );
		handler->setLegend( ENTRIES_DEFAULT );

		// Determine the plotted strand's sequence.
		for( int i = 1; i <= length; i++ ) {
			stringstream seqStream( stringstream::in | stringstream::out );
			seqStream << object->GetNucleotide( i );
			sequence.push_back( seqStream.str() );
		}

		// Delete the error checker and RNA object.
		delete checker;
		delete object;
	}

	} catch (exception& ex) {
		cout << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
		cout << "DotPlotBackend readData Error: " << ex.what() << endl;
		cout << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	} catch (int ex) {
		cout << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
		cout << "DotPlotBackend readData Error: (int) " << ex << endl;
		cout << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	} catch (string& ex) {
		cout << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
		cout << "DotPlotBackend readData Error: (str) " << ex << endl;
		cout << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	} catch (char* ex) {
		cout << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
		cout << "DotPlotBackend readData Error:  chr) " << ex << endl;
		cout << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	} catch (...) {
		cout << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
		cout << "DotPlotBackend readData Error (Unknown)";
		cout << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	}

	// If the method got to this point, it means that data was read correctly.
	// Return true.
	return true;
}

///////////////////////////////////////////////////////////////////////////////
// Set the colors allowed in the plot.
///////////////////////////////////////////////////////////////////////////////
void DotPlotBackend::setColors( int colors ) {

	plotColors = colors;
	handler->setLegend( colors );
}

///////////////////////////////////////////////////////////////////////////////
// Set the range for the plot.
///////////////////////////////////////////////////////////////////////////////
void DotPlotBackend::setRange( double minimum, double maximum ) {

	// Set the minimum.
	handler->setLegendMinimum( minimum );
	if( maximum < defaultMin ) { return; }
	if( maximum > defaultMax ) { return; }
	if( maximum < currentMin ) { return; }
	currentMax = maximum;

	// Set the maximum.
	handler->setLegendMaximum( maximum );
	if( minimum < defaultMin ) { return; }
	if( minimum > defaultMax ) { return; }
	if( minimum > currentMax ) { return; }
	currentMin = minimum;

	// Redo the legend.
	handler->setLegend( plotColors );
}

///////////////////////////////////////////////////////////////////////////////
// Write a Postscript file.
///////////////////////////////////////////////////////////////////////////////
void DotPlotBackend::writePostscriptFile( string outFile ) {

	handler->writePostscriptImage( outFile );
}

///////////////////////////////////////////////////////////////////////////////
// Write a probable structures file.
///////////////////////////////////////////////////////////////////////////////
string DotPlotBackend::writeStructuresFile( string file, string outFile ) {

	// Create the RNA object and its error checker.
	RNA* strand = new RNA( file.c_str(), FILE_PFS );
	ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );
	string error = checker->returnError();

	// If no error occurred, calculate structures.
	if( error == "" ) {
		int mainCalcError = strand->PredictProbablePairs( 0 );
		error = checker->returnError( mainCalcError );
	}

	// If no error occurred, write a CT file.
	if( error == "" ) {
		int writeError = strand->WriteCt( outFile.c_str() );
		error = checker->returnError( writeError );
	}

	// Delete the error checker and data structure, and return any error.
	delete checker;
	delete strand;
	return error;
}

///////////////////////////////////////////////////////////////////////////////
// Write an SVG file.
///////////////////////////////////////////////////////////////////////////////
void DotPlotBackend::writeSVGFile( string outFile ) {

	handler->writeSVGImage( outFile );
}

///////////////////////////////////////////////////////////////////////////////
// Write a dot plot text file.
///////////////////////////////////////////////////////////////////////////////
void DotPlotBackend::writeTextFile( string file ) {

	handler->writeTextFile( file );
}

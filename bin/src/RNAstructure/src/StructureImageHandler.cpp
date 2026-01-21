/*
 * An implementation file for a class that holds structure image data and can write those images to files if necessary.
 *
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter.
 */

#include "StructureImageHandler.h"
#include "loop_utils.h"
#include "../RNA_class/RNA.h"
#include "../RNA_class/ProbScan.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
StructureImageHandler::StructureImageHandler() {
	bimolecular = false;
	description = "";
	maxY = maxX = numeric_limits<double>::infinity() * -1.0;
}

///////////////////////////////////////////////////////////////////////////////
// Add probability annotation to the structure.
///////////////////////////////////////////////////////////////////////////////
string StructureImageHandler::addAnnotationProbability(
	string file, bool text ) {

	// Pull out the pairs from the structure.
	map<int, int> pairingMap = getPairs();

	// Create the vector of annotation probabilities.
	vector<double> probabilities;

	// If the annotation file is a partition function save file, use an RNA
	// object to read in probabilities.
	if( text == false ) {

		// Create the data strand and its error checker.
		RNA* strand = new RNA( file.c_str(),FILE_PFS );
		ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );

		// Check the data strand for errors, and if one occurred, delete the data structures and return the error.
		string result = checker->returnError();
		if( result != "" ) {
			delete strand;
			delete checker;
			return result;
		}

		// If the annotation sequence length is too long for the given coordinates, return an error.
		int length = strand->GetSequenceLength();
		if( length > coordinates.size() ) {
			delete strand;
			delete checker;
			return "Cannot apply this annotation file to the given sequence.";
		}

		// For each nucleotide, get its probability and set it in the probabilities vector.
		for( int i = 1; i <= length; i++ ) {

			// Initialize the probability and initialize the nucleotide's pair, if any.
			double probability = 1.0;
			int pair = 0;
			if( pairingMap.find( i ) != pairingMap.end() ) { pair = pairingMap[i]; }

			// If the nucleotide is paired, set its pairing probability.
			if( pair != 0 ) {
				if( pair > i ) { probability = strand->GetPairProbability( i, pair ); }
				else { probability = strand->GetPairProbability( pair, i ); }
			}

			// If the nucleotide is not paired, calculate its probability of being single-stranded.
			else {

				// For each nucleotide up to but not including the current one, subtract the probability that the current nucleotide will be paired to any of them.
				for( int j = 1; j < i; j++ ) {
					probability = probability - strand->GetPairProbability( j, i );
				}

				// For each nucleotide after but not including the current one, subtract the probability that the current nucleotide will be paired to any of them.
				for( int j = i + 1; j <= length; j++ ) {
					probability = probability - strand->GetPairProbability( i, j );
				}
			}

			// Add the probability to the vector.
			probabilities.push_back( probability );
		}

		// Delete the data strand and error checker.
		delete strand;
		delete checker;
	}

	// If the annotation file is a text file, read annotation data from that.
	else {

		// Declare a variable for line reading.
		string line;

		// Attempt to open the annotation file.
		// If it couldn't be open, return an error.
		ifstream in( file.c_str() );
		if( !in.is_open() ) { return "Error opening annotation file."; }

		// Read in the first line to get the length of the annotation sequence.
		// If the sequence length cannot be read, return an error.
		// If the annotation sequence length is too long for the given coordinates, return an error.
		int length;
		getline( in, line );
		stringstream lengthStream( line );
		if( lengthStream >> length ) {
			if( length > coordinates.size() ) {
				return "Cannot apply this annotation file to the given sequence.";
			}
		} else { return "Error reading annotation file."; }

		// Skip the second line, which has column headers.
		getline( in, line );

		// Initialize the probabilities as -1, so they'll all be black later.
		for( int i = 1; i <= length; i++ ) { probabilities.push_back( -1 ); }

		// For all the other lines, read them in.
		int i, j;
		double value;
		while( !in.eof() ) {

			// Get the next line, and put it into a stream.
			getline( in, line );
			stringstream lineStream( line );
			if( !( lineStream.str() == "" ) ) {

				// Read in the first index and return an error if it didn't happen correctly.
				if( !( lineStream >> i ) ) {
					in.close();
					return "Error reading annotation file.";
				}

				// Read in the second index and return an error if it didn't happen correctly.
				if( !( lineStream >> j ) ) {
					in.close();
					return "Error reading annotation file.";
				}

				// Read in the probability value and return an error if it didn't happen correctly.
				if( !( lineStream >> value ) ) {
					in.close();
					return "Error reading annotation file.";
				}

				// If the pairing map contains the read pair, save the probability.
				if( pairingMap.find( i ) != pairingMap.end() ) {
					if( pairingMap[i] == j ) {
						double converted = pow( 10, -value );
						probabilities[i-1] = converted;
						probabilities[j-1] = converted;
					}
				}
			}
		}

		// Close the annotation file.
		in.close();
	}

	// Clear any existing annotation data.
	annotations.clear();
	legend.clear();
	legendColors.clear();

	// For each probability, add its appropriate color to the annotations vector.
	unsigned int numProbs = probabilities.size();
	for( int i = 1; i <= numProbs; i++ ) {
		double probability = probabilities[i-1];
		string color =
			( probability >= 0.99 ) ? RED :
			( probability >= 0.95 ) ? ORANGE :
			( probability >= 0.90 ) ? YELLOW :
			( probability >= 0.80 ) ? GREEN :
			( probability >= 0.70 ) ? LIGHT_GREEN :
			( probability >= 0.60 ) ? LIGHT_BLUE :
			( probability >= 0.50 ) ? BLUE :
			( probability > 0 ) ? PINK :
			BLACK;
		annotations.push_back( color );
	}

	// Set the probability annotation legend text.
	legend.push_back( "      Probability >= 99%" );
	legend.push_back( "99% > Probability >= 95%" );
	legend.push_back( "95% > Probability >= 90%" );
	legend.push_back( "90% > Probability >= 80%" );
	legend.push_back( "80% > Probability >= 70%" );
	legend.push_back( "70% > Probability >= 60%" );
	legend.push_back( "60% > Probability >= 50%" );
	legend.push_back( "50% > Probability" );

	// Set the probability annotation colors.
	legendColors.push_back( RED );
	legendColors.push_back( ORANGE );
	legendColors.push_back( YELLOW );
	legendColors.push_back( GREEN );
	legendColors.push_back( LIGHT_GREEN );
	legendColors.push_back( LIGHT_BLUE );
	legendColors.push_back( BLUE );
	legendColors.push_back( PINK );

	// Return an empty string to show that no error occurred.
	return "";
}

template<typename T>
void fill_probabilities(const vector<T>& loops,
                        map<int,double>& probabilities,
                        ProbScan& p)
{
    for(size_t i=0;i<loops.size();i++){
        const vector<int> n = loops[i].nucs();
        for(size_t j=0;j<n.size();j++){
            probabilities[n[j]-1] = loops[i].getProbability(p);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// Add loop probability annotation to the structure.
///////////////////////////////////////////////////////////////////////////////
string StructureImageHandler::addAnnotationLoopProbability(
	string file,RNA& r, int structurenumber ) {

	// Pull out the pairs from the structure.
	map<int, int> pairingMap = getPairs();

	// Create the table of annotation probabilities.
	map<int,double> probabilities;

	// use a probscan object to read in probabilities.
    ProbScan p = ProbScan( file.c_str(), false);
    ErrorChecker<ProbScan> checker = ErrorChecker<ProbScan>( &p );

    // Check the data strand for errors, and if one occurred, delete the data structures and return the error.
    string result = checker.returnError();
    if( result != "" ) {
        return result;
    }

    // If the annotation sequence length is too long for the given coordinates, return an error.
    int length = p.GetSequenceLength();
    if( length > (int) coordinates.size() ) {
        return "Cannot apply this annotation file to the given sequence.";
    }

    for(int i=0;i<r.GetSequenceLength();i++){
        //initialize all probabilities to a negative value
        //so that loops with no probability information get colored black
        probabilities[i] = -1.0;
    }
    //fill the probabilities using the probscan object
    fill_probabilities(loop::find_hairpins(r,structurenumber),probabilities,p);
    fill_probabilities(loop::find_internals(r,structurenumber),probabilities,p);
    fill_probabilities(loop::find_multibranch(r,structurenumber),probabilities,p);
    fill_probabilities(loop::find_stems(r,structurenumber),probabilities,p);

	// Clear any existing annotation data.
	annotations.clear();
	legend.clear();
	legendColors.clear();

	// For each probability, add its appropriate color to the annotations vector.
	for( int i = 0; i < r.GetSequenceLength(); i++ ) {
		double probability = probabilities.at(i);
		string color =
			( probability >= 0.90 ) ? RED :
			( probability >= 0.80 ) ? ORANGE :
			( probability >= 0.70 ) ? YELLOW :
			( probability >= 0.60 ) ? GREEN :
			( probability >= 0.50 ) ? LIGHT_GREEN :
			( probability >= 0.30 ) ? LIGHT_BLUE :
			( probability >= 0.10 ) ? BLUE :
			( probability > 0 ) ? PINK :
			BLACK;
		annotations.push_back( color );
	}

	// Set the probability annotation legend text.
	legend.push_back( "      Probability >= 90%" );
	legend.push_back( "90% > Probability >= 80%" );
	legend.push_back( "80% > Probability >= 70%" );
	legend.push_back( "70% > Probability >= 60%" );
	legend.push_back( "60% > Probability >= 50%" );
	legend.push_back( "50% > Probability >= 30%" );
	legend.push_back( "30% > Probability >= 10%" );
	legend.push_back( "10% > Probability" );

	// Set the probability annotation colors.
	legendColors.push_back( RED );
	legendColors.push_back( ORANGE );
	legendColors.push_back( YELLOW );
	legendColors.push_back( GREEN );
	legendColors.push_back( LIGHT_GREEN );
	legendColors.push_back( LIGHT_BLUE );
	legendColors.push_back( BLUE );
	legendColors.push_back( PINK );

	// Return an empty string to show that no error occurred.
	return "";
}

///////////////////////////////////////////////////////////////////////////////
// Add SHAPE annotation to the structure.
///////////////////////////////////////////////////////////////////////////////
string StructureImageHandler::addAnnotationSHAPE( string file ) {

	annotations.clear(); // Clear any existing annotation data.
	annotations.resize(coordinates.size(), GRAY); // fill entire annotations array with the default value.
	legend.clear();
	legendColors.clear();

	// Attempt to open the annotation file.
	// If it couldn't be open, return an error.
	ifstream in( file.c_str() );
	if( !in.good() ) { return "Error opening annotation file."; }

	// Read in the annotation data.
	unsigned long linenumber = 0;
	string line;
	while( getline( in, line ) ) {
		linenumber++;
		// Skip blank lines.
		trim(line); if(line.empty()) continue;
		// Put the line contents into a string stream.
		stringstream lineStream( stringstream::in | stringstream::out );
		lineStream.str(line);

		// Read the index affected by SHAPE data from the line.
		// If the index is valid, increase the size of the annotation array up to and including to the index.
		// If the next index isn't a proper number, return a formatting error.
		unsigned int index = 0;
		if( !(lineStream >> index) ) return sfmt("Invalid (non-numeric) index in SHAPE file \"%s\" on line %li", file.c_str(), linenumber);
		if (index > coordinates.size())
				return sfmt("SHAPE data index %u is greater than sequence length (%u bases) in SHAPE file \"%s\" on line %li.", index, coordinates.size(), file.c_str(), linenumber);

		// If the value is valid set the annotation color, otherwise return a formatting error.
		string &color = annotations[index-1];
		double shape;
		if( !(lineStream >> shape) ) return sfmt("Invalid (non-numeric) data value in SHAPE file \"%s\" on line %li.", file.c_str(), linenumber);

		if (!lineStream.eof())  return sfmt("Invalid SHAPE data format (or extra text) in SHAPE file \"%s\" on line %li.", file.c_str(), linenumber);

		if (shape > 0.85)      color = RED;
		else if (shape >= 0.4) color = ORANGE;
		else if (shape >= -500) color = BLACK;
		// otherwise, leave it default (GRAY)
	}

	// Once file reading is done, close the file.
	in.close();

	// Set the SHAPE annotation legend text.
	legend.push_back( "      SHAPE >= 0.85" );
	legend.push_back( "0.85 > SHAPE >= 0.4" );
	legend.push_back( "0.4 > SHAPE" );
	legend.push_back( "No Data" );

	// Set the SHAPE annotation colors.
	legendColors.push_back( RED );
	legendColors.push_back( ORANGE );
	legendColors.push_back( BLACK );
	legendColors.push_back( GRAY );

	// Return an empty string to show that no error occurred.
	return "";
}

///////////////////////////////////////////////////////////////////////////////
// Flip the structure image horizontally.
///////////////////////////////////////////////////////////////////////////////
void StructureImageHandler::flipHorizontally() {
	isFlipped=!isFlipped;
	//cout << "Flipping. Now: " << isFlipped << endl;
	// Flip the X components of the coordinates vector.
	for( unsigned int i = 1; i <= coordinates.size(); i++ ) {

		// If a coordinate is blank, skip to the next non-blank one.
		while( coordinates[i-1] == "" ) { i++; }

		// Pull the data out of the next pair string.
		stringstream outStream( coordinates[i-1] );
		string nuc;
		double nucX, nucY, backX, backY, labelX, labelY, num;
		outStream >> nuc;
		outStream >> nucX;
		outStream >> nucY;
		outStream >> backX;
		outStream >> backY;
		outStream >> labelX;
		outStream >> labelY;
		outStream >> num;

		// Flip coordinates where necessary.
		double convertedNucX = maxX - nucX;
		double convertedBackX = ( backX != -1 ) ? maxX - backX : backX;
		double convertedLabelX = ( labelX != -1 ) ? maxX - labelX : labelX;

		// Rebuild the coordinate string with converted data.
		stringstream inStream( stringstream::in | stringstream::out );
		inStream << nuc << " "
		         << convertedNucX << " " << nucY << " "
		         << convertedBackX << " " << backY << " "
		         << convertedLabelX << " " << labelY << " " << num;
		coordinates[i-1] = inStream.str();
	}

	// Flip the X components of the pairs vector.
	for( unsigned int i = 1; i <= pairs.size(); i++ ) {

		// Pull the data out of the next pair string.
		stringstream outStream( pairs[i-1] );
		string pairing, colorData;
		double x1, y1, centerX, centerY, x2, y2;
		outStream >> pairing;
		outStream >> x1;
		outStream >> y1;
		outStream >> centerX;
		outStream >> centerY;
		outStream >> x2;
		outStream >> y2;
		getline( outStream, colorData );

		// Rebuild the pair string with converted data.
		stringstream inStream( stringstream::in | stringstream::out );
		inStream << pairing << " "
		         << ( maxX - x1 ) << " " << y1 << " "
		         << ( maxX - centerX ) << " " << centerY << " "
		         << ( maxX - x2 ) << " " << y2;
		if( colorData != "" ) { inStream << " " << colorData; }
		pairs[i-1] = inStream.str();
	}
}
void StructureImageHandler::setFlipped(bool flipped) {
	if (isFlipped != flipped) flipHorizontally();
}

void StructureImageHandler::reset() {
	// Clear the coordinates and pairs.
	coordinates.clear();
	pairs.clear();
	maxY = maxX = numeric_limits<double>::infinity() * -1.0;
	isFlipped=false;
}

///////////////////////////////////////////////////////////////////////////////
// Read circular structure data.
///////////////////////////////////////////////////////////////////////////////
string StructureImageHandler::readCircular( string file, int number ) {
	reset();

	// Create the data strand and its error checker.
	RNA* strand = new RNA( file.c_str(), FILE_CT, DT_RNA, true,true);
	ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );

	// Check the data strand for errors, and if one occurred, delete the data structures and return the error.
	string result = checker->returnError();
	if( result != "" ) {
		delete strand;
		delete checker;
		return result;
	}

	// If the desired structure doesn't exist, return an error.
	int numStructures = strand->GetStructureNumber();
	if( ( number < 1 ) || ( number > numStructures ) ) {
		delete strand;
		delete checker;
		return "Structure number out of range.";
	}

	// Initialize the sequence length and description.
	int length = strand->GetSequenceLength();
	description = strand->GetCommentString( number );

	// Determine how many places the sequence length takes up.
	int places =
		( length >= 10000 ) ? 5 :
		( length >= 1000 ) ? 4 :
		( length >= 100 ) ? 3 :
		( length >= 10 ) ? 2 : 1;
	places++;

	// Determine whether the structure is bimolecular.
	for( int i = 1; i <= length; i++ ) {
		if( strand->GetNucleotide( i ) == 'I' ) {
			bimolecular = true;
			i += length;
		}
	}

	// Calculate the layout data for the circle.
	int radius = ( TEXTSIZE * (length + 2) ) / 4;
	double slice = 360.0 / ( (double)(length+2) );
	int lineLength = ( !bimolecular ) ? 3 * TEXTSIZE : 0;
	double translateToCenter = (double)radius + (double)lineLength + (double)( places * GLYPHSIZE );
	double centerCircleX = translateToCenter;
	double centerCircleY = translateToCenter;
	maxX = 2 * translateToCenter;
	maxY = 2 * translateToCenter;

	// Calculate the nucleotide coordinates, and that of their labels, if necessary.
	for( int i = 1; i <= length; i++ ) {

		// Create a string stream that holds nucleotide data.
		stringstream stream( stringstream::in | stringstream::out );
		stream.str( "" );

		// Get the next nucleotide, and only find its coordinates if it's not part of the intramolecular linker.
		char nuc = strand->GetNucleotide( i );
		if( nuc != 'I' ) {

			// Add the nucleotide to the stream.
			stream << nuc << " ";

			// Determine what the angle and X and Y coordinates of this nucleotide should be.
			double radians = (slice * (double)i) * (PI / 180.0) - (PI / 2.0); 
			double x = ( cos( radians ) * (double)radius ) + translateToCenter;
			double y = ( sin( radians ) * (double)radius ) + translateToCenter;
			stream << " " << x << " " << y << " ";

			// Add "-1 -1" as the backbone coordinates, since circular structures aren't drawn with a backbone.
			stream << "-1 -1 ";

			// If the structure isn't bimolecular, put the nucleotide number label coordinates and number in the string stream.
			// If not, just add "-1 -1 -1" as placeholders.
			if( ( !bimolecular ) && ( i % 10 == 0 ) ) {

				// Determine the coordinates for the end of the label line.
				double xLabel = x + ( lineLength * cos( radians ) );
				double yLabel = y + ( lineLength * sin( radians ) );
				stream << xLabel << " " << yLabel << " " << i;
			} else { stream << "-1 -1 -1"; }
		}

		// Pull the string out of the string stream and add it to the coordinates vector.
		coordinates.push_back( stream.str() );
	}

	// Calculate the location of the nucleotide pairings.
	// Only calculate a nucleotide pairing location if the nucleotide is paired to begin with.
	for( int i = 1; i <= length; i++ ) {
		int pair = strand->GetPair( i, number );
		if( pair > i ) {

			// Initialize the variables needed to find pairing data.
			double x1, y1, x2, y2;
			string temp;
			istringstream stream( coordinates[i-1] );
			istringstream stream2( coordinates[pair-1] );

			// Get the X and Y coordinates of the first nucleotide.
			stream >> temp;
			stream >> x1;
			stream >> y1;

			// Get the X and Y coordinates of the second nucleotide.
			stream2 >> temp;
			stream2 >> x2;
			stream2 >> y2;

			// Calculate the midpoint between the nucleotides.
			double midX = ( x1 + x2 ) / 2.0;
			double midY = ( y1 + y2 ) / 2.0;

			// Calculate quantities necessary to determine the altitude (adjusted distance between midpoint and circle center) of a pair.
			double distanceBetweenCenterAndMidpoint = sqrt( pow( centerCircleX - midX, 2 ) + pow( centerCircleY - midY, 2 ) );
			double distanceBetweenNucleotides = sqrt( pow( x2 - x1, 2 ) + pow( y2 - y1, 2 ) );
			double angleBetweenMidpointAndCircleCenter = atan2( midY - centerCircleY, midX - centerCircleX );

			// Calculate the altitude, and its X/Y components.
			double altitude = sqrt( 0.75 * pow ( distanceBetweenNucleotides, 2 ) );
			double altitudeX = cos( angleBetweenMidpointAndCircleCenter ) * altitude * ( distanceBetweenCenterAndMidpoint / (double)radius );
			double altitudeY = sin( angleBetweenMidpointAndCircleCenter ) * altitude * ( distanceBetweenCenterAndMidpoint / (double)radius );

			// Calculate the X and Y coordinates of the control point at half the altitude.
			double centerX = midX - ( altitudeX / 2.0 );
			double centerY = midY - ( altitudeY / 2.0 );

			// Create a string stream to hold the pairing points in order, then extract the data as one string and add it to the pairings vector.
			stringstream pairData( stringstream::in | stringstream::out );
			pairData << i << "-" << pair << " " << x1 << " " << y1 << " " << centerX << " " << centerY << " " << x2 << " " << y2;
			pairs.push_back( pairData.str() );
		}
	}

	// Delete the data strand and error checker, then return an empty string to show that no error occurred.
	delete strand;
	delete checker;
	return "";
}

///////////////////////////////////////////////////////////////////////////////
// Read linear structure data.
///////////////////////////////////////////////////////////////////////////////
string StructureImageHandler::readLinear( string file, int number ) {

	// Clear the coordinates and pairs.
	reset();

	// Create the data strand and its error checker.
	RNA* strand = new RNA( file.c_str(), FILE_CT, DT_RNA, true, true);
	ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );

	// Check the data strand for errors, and if one occurred, delete the data structures and return the error.
	string result = checker->returnError();
	if( result != "" ) {
		delete strand;
		delete checker;
		return result;
	}

	// If the desired structure doesn't exist, return an error.
	int numStructures = strand->GetStructureNumber();
	if( ( number < 1 ) || ( number > numStructures ) ) {
		delete strand;
		delete checker;
		return "Structure number out of range.";
	}

	// Initialize the sequence length and description.
	int length = strand->GetSequenceLength();
	description = strand->GetCommentString( number );

	// Determine whether the structure is bimolecular.
	for( int i = 1; i <= length; i++ ) {
		if( strand->GetNucleotide( i ) == 'I' ) {
			bimolecular = true;
			i += length;
		}
	}

	// Calculate the maximum Y bound.
	double backbonePieceLength = (double)( GLYPHSIZE * 2 );
	double stretchedLineLength = length * backbonePieceLength;
	int lineLength = ( !bimolecular ) ? 3 * TEXTSIZE : 0;
	maxY = stretchedLineLength + ( 2 * BORDER );

	// Create constants that determine the X coordinates of the nucleotides and their labels, which are always the same.
	double nucX = lineLength + (double)( 2 * GLYPHSIZE );
	double nucLabelX = nucX - lineLength;

	// Calculate the nucleotide coordinates, and that of their labels, if necessary.
	for( int i = 1; i <= length; i++ ) {

		// Create a string stream that holds nucleotide data.
		stringstream stream( stringstream::in | stringstream::out );
		stream.str( "" );

		// Get the next nucleotide, and only find its coordinates if it's not part of the intramolecular linker.
		char nuc = strand->GetNucleotide( i );
		if( nuc != 'I' ) {

			// Add the nucleotide to the stream.
			stream << nuc << " ";

			// Determine what the Y coordinate of this nucleotide should be and put it in the stream with the standard X coordinate.
			double y = BORDER + ( backbonePieceLength * ( i - 1 ) );
			stream << " " << nucX << " " << y << " ";

			// Determine what the Y coordinate of this nucleotide's backbone should be and put it in the stream with the standard X coordinate.
			// If the nucleotide is the last one, just add "-1 -1" as placeholders.
			double backY = y + backbonePieceLength;
			if( i != length ) { stream << nucX << " " << backY << " "; }
			else { stream << "-1 -1 "; }

			// If the structure isn't bimolecular, put the nucleotide number label coordinates and number in the string stream.
			// If not, just add "-1 -1 -1" as placeholders.
			if( ( !bimolecular ) && ( i % 10 == 0 ) ) { stream << nucLabelX << " " << y << " " << i; }
			else { stream << "-1 -1 -1"; }
		}

		// Pull the string out of the string stream and add it to the coordinates vector.
		coordinates.push_back( stream.str() );
	}

	// Calculate the location of the nucleotide pairings.
	// Only calculate a nucleotide pairing location if the nucleotide is paired to begin with.
	for( int i = 1; i <= length; i++ ) {
		int pair = strand->GetPair( i, number );
		if( pair > i ) {

			// Initialize the variables needed to find pairing data.
			double y1, y2;
			string temp;
			istringstream stream( coordinates[i-1] );
			istringstream stream2( coordinates[pair-1] );

			// Get the Y coordinate of the first nucleotide.
			for( int j = 1; j <= 2; j++ ) { stream >> temp; }
			stream >> y1;

			// Get the Y coordinate of the second nucleotide.
			for( int j = 1; j <= 2; j++ ) { stream2 >> temp; }
			stream2 >> y2;

			// Calculate 70% of the distance between nucleotides plus the X baseline (X coordinate of control point).
			// Calculate the midpoint between the nucleotides (Y coordinate of control point).
			double controlX = ( ( y2 - y1 ) * 0.7 ) + nucX;
			double controlY = ( y1 + y2 ) / 2.0;

			// Calculate the maximum X bound using the control point.
			maxX = max( maxX, controlX + lineLength + ( 2 * GLYPHSIZE ) );

			// Create a string stream to hold the pairing points in order, then extract the data as one string and add it to the pairings vector.
			stringstream pairData( stringstream::in | stringstream::out );
			pairData << i << "-" << pair << " " << nucX << " " << y1 << " " << controlX << " " << controlY << " " << nucX << " " << y2;
			pairs.push_back( pairData.str() );
		}
	}

	// Delete the data strand and error checker, then return an empty string to show that no error occurred.
	delete strand;
	delete checker;
	return "";
}

///////////////////////////////////////////////////////////////////////////////
// Read radial structure data.
///////////////////////////////////////////////////////////////////////////////
string StructureImageHandler::readRadial( string file, int number ) {

	// Clear the coordinates and pairs.
	reset();

	// Create the data strand and its error checker.
	RNA* strand = new RNA( file.c_str(), FILE_CT, DT_RNA, true,true);
	ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );

	// Check the data strand for errors, and if one occurred, delete the data structures and return the error.
	string result = checker->returnError();
	if( result != "" ) {
		delete strand;
		delete checker;
		return result;
	}

	// If the desired structure doesn't exist, return an error.
	int numStructures = strand->GetStructureNumber();
	if( ( number < 1 ) || ( number > numStructures ) ) {
		delete strand;
		delete checker;
		return "Structure number out of range.";
	}

	// Initialize the sequence length and description.
	int length = strand->GetSequenceLength();
	description = strand->GetCommentString( number );

	// If the structure has pseudoknots in it, it can't be drawn radial, so clean up data for this method, call circular drawing, and return that routine's result.
	if( strand->ContainsPseudoknot( number ) ) {
		delete strand;
		delete checker;
		return readCircular( file, number );
	}

	// Calculate all the drawing coordinates and determine whether the structure is bimolecular.
	int drawSize = TEXTSIZE + 4;
	strand->DetermineDrawingCoordinates( drawSize, drawSize, number );
	for( int i = 1; i <= length; i++ ) {
		if( strand->GetNucleotide( i ) == 'I' ) {
			bimolecular = true;
			i += length;
		}
	}

	// Go through all the calculated nucleotide and label locations to calculate how nucleotides should be translated to fit on the image.
	int minX = numeric_limits<int>::infinity(), minY = minX;
	for( int i = 1; i <= length; i++ ) {

		// Get the nucleotide location coordinates.
		int x = strand->GetNucleotideXCoordinate( i );
		int y = strand->GetNucleotideYCoordinate( i );

		// Check whether the next nucleotide's location should update the offset or max bounds.
		if( strand->GetNucleotide( i ) != 'I' ) {
			if( x < minX ) { minX = x; }
			if( x > maxX ) { maxX = x; }

			if( y < minY ) { minY = y; }
			if( y > maxY ) { maxY = y; }
		}

		// Check whether the next nucleotide label's location should update the offset or max bounds.
		if( !bimolecular ) {

			// Get the nucleotide label coordinates.
			int xLabel = strand->GetLabelXCoordinate( i );
			int yLabel = strand->GetLabelYCoordinate( i );

			// If the nucleotide label is valid, check its location.
			if( !( (xLabel == 0) && (yLabel == 0) ) ) {

				// Determine how many places the label takes up, based on the number.
				int places =
					( i >= 10000 ) ? 5 :
					( i >= 1000 ) ? 4 :
					( i >= 100 ) ? 3 :
					( i >= 10 ) ? 2 : 1;
				places++;

				// Calculate the number width.
				double boxWidth = 0.0;
				for( int i = 1; i <= places; i++ ) { boxWidth += GLYPHSIZE; }

				// If the X coordinate of the label line is out of bounds, adjust the bounds accordingly.
				if( (xLabel - boxWidth) < minX ) { minX = xLabel - boxWidth; }
				if( (xLabel + boxWidth) > maxX ) { maxX = xLabel + boxWidth; }

				// If the Y coordinate of the label line is out of bounds, adjust the bounds accordingly.
				if( (yLabel - TEXTSIZE) < minY ) { minY = yLabel - TEXTSIZE; }
				if( (yLabel + TEXTSIZE) > maxY ) { maxY = yLabel + TEXTSIZE; }
			}
		}
	}

	// Create X and Y offsets to make sure the lowest X and/or Y coordinate is 0.
	// Then, edit the minimum and maximum bounds.
	int offsetX = minX * -1, offsetY = minY * -1;
	minX += offsetX;
	maxX += offsetX;
	minY += offsetY;
	maxY += offsetY;
	maxX += ( BORDER * 2 );
	maxY += ( BORDER * 2 );

	// Calculate the nucleotide coordinates, and that of their labels, if necessary.
	for( int i = 1; i <= length; i++ ) {

		// Create a string stream that holds nucleotide data.
		stringstream stream( stringstream::in | stringstream::out );
		stream.str( "" );

		// Get the next nucleotide, and only find its coordinates if it's not part of the intramolecular linker.
		char nuc = strand->GetNucleotide( i );
		if( nuc != 'I' ) {

			// Add the nucleotide to the stream.
			stream << nuc << " ";

			// Determine what the raw X and Y coordinates of this nucleotide should be.
			int x = strand->GetNucleotideXCoordinate( i ) + offsetX + BORDER;
			int y = strand->GetNucleotideYCoordinate( i ) + offsetY + BORDER;
			stream << x << " " << y << " ";

			// Calculate the backbone line segment's coordinates if a line segment should be drawn.
			// Otherwise, just add "-1 -1" as placeholders.
			if( ( i != length ) && ( strand->GetNucleotide( i ) != 'I' ) && ( strand->GetNucleotide( i + 1 ) != 'I' ) ) {
				int backX = strand->GetNucleotideXCoordinate( i + 1 ) + offsetX + BORDER;
				int backY = strand->GetNucleotideYCoordinate( i + 1 ) + offsetY + BORDER;
				stream << backX << " " << backY << " ";
			} else { stream << "-1 -1 "; }

			// If the structure isn't bimolecular, put the nucleotide number label coordinates and number in the string stream.
			// If not, just add "-1 -1 -1" as placeholders.
			if( !bimolecular ) {

				// Get the nucleotide label coordinates.
				int xLabel = strand->GetLabelXCoordinate( i );
				int yLabel = strand->GetLabelYCoordinate( i );

				// If the nucleotide label is valid, calculate its data.
				if( !( (xLabel == 0) && (yLabel == 0) ) ) {
					xLabel = strand->GetLabelXCoordinate( i ) + offsetX + BORDER;
					yLabel = strand->GetLabelYCoordinate( i ) + offsetY + BORDER;
					stream << xLabel << " " << yLabel << " " << i;
				} else { stream << "-1 -1 -1"; }
			} else { stream << "-1 -1 -1"; }
		}

		// Pull the string out of the string stream and add it to the coordinates vector.
		coordinates.push_back( stream.str() );
	}

	// Calculate the location of the nucleotide pairings.
	// Only calculate a nucleotide pairing location if the nucleotide is paired to begin with.
	for( int i = 1; i <= length; i++ ) {
		int pair = strand->GetPair( i, number );
		if( pair > i ) {

			// Initialize the variables needed to find pairing data.
			double x1, y1, x2, y2;
			istringstream stream( coordinates[i-1] );
			istringstream stream2( coordinates[pair-1] );
			string temp;

			// Get the X and Y coordinates of the first nucleotide.
			stream >> temp;
			stream >> x1;
			stream >> y1;

			// Get the X and Y coordinates of the second nucleotide.
			stream2 >> temp;
			stream2 >> x2;
			stream2 >> y2;

			// Calculate the midpoint between the nucleotides.
			double midX = ( x1 + x2 ) / 2.0;
			double midY = ( y1 + y2 ) / 2.0;

			// Create a string stream to hold the pairing points in order, then extract the data as one string and add it to the pairings vector.
			stringstream pairData( stringstream::in | stringstream::out );
			pairData << i << "-" << pair << " " << x1 << " " << y1 << " " << midX << " " << midY << " " << x2 << " " << y2;
			pairs.push_back( pairData.str() );
		}
	}

	// Flip the coordinates horizontally so the structure is rendered clockwise by default.
	flipHorizontally();

	// Delete the data strand and error checker, then return an empty string to show that no error occurred.
	delete strand;
	delete checker;
	return "";
}

///////////////////////////////////////////////////////////////////////////////
// Remove annotation from this structure.
///////////////////////////////////////////////////////////////////////////////
void StructureImageHandler::removeAnnotation() {
	annotations.clear();
	legend.clear();
	legendColors.clear();
}

///////////////////////////////////////////////////////////////////////////////
// Set nucleotides circled.
///////////////////////////////////////////////////////////////////////////////
void StructureImageHandler::setNucleotidesCircled( bool encircle ) {
	circleNucs = encircle;
}

///////////////////////////////////////////////////////////////////////////////
// Get a string representation of the structure.
///////////////////////////////////////////////////////////////////////////////
string StructureImageHandler::toString() {

	// If no structure data has been read, return an empty string.
	if( coordinates.size() == 0 ) { return ""; }

	// Create the string stream.
	stringstream stream( stringstream::in | stringstream::out );

	// Initialize the string to show absent data.
	string absent = "N/A";

	// Get the number of nucleotides.
	unsigned int nucleotides = coordinates.size();

	// Add each nucleotide's data.
	for( unsigned int i = 1; i <= coordinates.size(); i++ ) {

		// Get the nucleotide's data.
		string nucChar, nucX, nucY, backX, backY, labelX, labelY, labelNum;
		stringstream dataStream( coordinates[i-1] );
		dataStream >> nucChar;
		dataStream >> nucX;
		dataStream >> nucY;
		dataStream >> backX;
		dataStream >> backY;
		dataStream >> labelX;
		dataStream >> labelY;
		dataStream >> labelNum;

		// Add the nucleotide number, character and placement data.
		stream << "Nucleotide " << i << " --" << endl
		       << "\tCharacter: " << nucChar << endl
		       << "\tPlaced at: " << "(" << nucX << "," << nucY << ")" << endl;

		// Add the color.
		string colorString = getColorString( BLACK, true );
		if( annotations.size() != 0 ) { colorString = getColorString( annotations[i-1], true ); }
		colorString = colorString.substr( colorString.find_first_of( "(" ) + 1 );
		colorString = colorString.substr( 0, colorString.find_first_of( ")" ) );
		stream << "\tColor (RGB): " << colorString << endl;

		// Add the backbone location.
		string locationString = absent;
		if( backX != "-1" ) { locationString = "(" + backX + "," + backY + ")"; }
		stream << "\tBackbone stretches to: " << locationString << endl;

		// Add the label string.
		string labelString = absent;
		if( labelX != "-1" ) { labelString = "(" + labelX + "," + labelY + ")"; }
		stream << "\tLabel placed at: " << labelString << endl;

		// Show the coordinates of the pair and its control point.
		string pairString = absent;
		string controlString = absent;
		for( unsigned int j = 1; j <= pairs.size(); j++ ) {
			size_t dash = pairs[j-1].find_first_of( "-" );
			stringstream testStream( pairs[j-1].substr( 0, dash ).c_str() );
			unsigned int test;
			testStream >> test;
			if( test == i ) {
				string temp, pair;
				string controlX, controlY;
				stringstream dataStream2( pairs[j-1].substr( dash+1 ) );
				dataStream2 >> pair;
				for( int k = 1; k <= 2; k++ ) { dataStream2 >> temp; }
				dataStream2 >> controlX;
				dataStream2 >> controlY;
				pairString = pair;
				controlString = "(" + controlX + "," + controlY + ")";
				break;
			}
		}
		stream << "\tPaired to: " << pairString << endl
		       << "\tControl point: " << controlString << endl;
	}

	// Show the maximum structure bounds.
	stream << "Max Bounds: " << "(" << maxX << "," << maxY << ")" << endl;

	// Show the description.
	stream << "Description: " << description << endl;

	// If a legend exists, show it as well.
	if( legend.size() > 0 ) {
		stream << "Legend:" << endl;
		for( unsigned int j = 1; j <= legend.size(); j++ ) {
			stream << legend[j-1] << " -- " << getColorString( legendColors[j-1], true ) << endl;
		}
	}

	// Return the string from the stream.
	return stream.str();

}

///////////////////////////////////////////////////////////////////////////////
// Write an image file.
///////////////////////////////////////////////////////////////////////////////
void StructureImageHandler::writeImageFile( string file, bool append, bool isSVG, int pagenumber, int pages) {

	// Initialize an index used to handle conversion template indices later.
	size_t index = string::npos;

	// Open the output file with the proper mode.
	ios_base::openmode mode = ( append ) ? ios_base::app : ios_base::trunc;
	ofstream out( file.c_str(), mode );

	// Write the start of the image.
	string startMarker = ( !isSVG ) ? createStartPS(pagenumber,pages) : createStartSVG();
	out << startMarker << endl;

	// If no pairs exist in this structure, set a description to say this.
	if( pairs.size() == 0 ) {
		string noStructureString = "The given structure contains no pairs.";
		string textData = ( !isSVG ) ? TEXT_PS : TEXT_SVG;
		while( ( index = textData.find( COLOR ) ) != string::npos ) { textData = textData.replace( index, COLOR.size(), getColorString( BLACK, isSVG ) ); }
		while( ( index = textData.find( LOCX ) ) != string::npos ) { textData = textData.replace( index, LOCX.size(), TEXTSIZE_BIGGER_STRING ); }
		while( ( index = textData.find( LOCY ) ) != string::npos ) { textData = textData.replace( index, LOCY.size(), TEXTSIZE_BIGGER_STRING ); }
		while( ( index = textData.find( TEXTSTRING ) ) != string::npos ) { textData = textData.replace( index, TEXTSTRING.size(), noStructureString ); }
		out << textData << endl;
	}

	// If pairs do exist in the file, write the structure image.
	else {

		// Convert each piece of coordinates data to a string that holds the nucleotide, its color, and its label, where applicable.
		vector<string> coordinatesTransformed;
		for( unsigned int i = 1; i <= coordinates.size(); i++ ) {

			// If a coordinate is blank, skip to the next non-blank one.
			while( coordinates[i-1] == "" ) { i++; }

			// Declare variables for conversion data.
			string nuc;
			double backX, backY, nucX, nucY, labelX, labelY;
			string backXString, backYString, nucXString, nucYString, labelXString, labelYString;
			int num;

			// Read in the conversion data.
			istringstream dataStream( coordinates[i-1] );
			dataStream >> nuc;
			dataStream >> nucX;
			dataStream >> nucY;
			dataStream >> backX;
			dataStream >> backY;
			dataStream >> labelX;
			dataStream >> labelY;
			dataStream >> num;

			// Convert the coordinate data to strings for later. 
			stringstream dataStringStream( stringstream::in | stringstream::out );
			dataStringStream << nucX << " " << nucY << " " << backX << " " << backY << " " << labelX << " " << labelY;
			dataStringStream >> nucXString;
			dataStringStream >> nucYString;
			dataStringStream >> backXString;
			dataStringStream >> backYString;
			dataStringStream >> labelXString;
			dataStringStream >> labelYString;

			// Get the raw color of the nucleotide.
			string color = BLACK;
			if( annotations.size() > 0 ) { color = annotations[i-1]; }

			// Determine what the X and Y coordinates of the nucleotide character label should be, as numbers and strings.
			int charLabelX = nucX - ( TEXTSIZE / 3 );
			int charLabelY = nucY + ( TEXTSIZE / 3 );
			if( nuc == "g" ) { charLabelY -= 3; }
			string charLabelXString, charLabelYString;
			stringstream charConvertStream( stringstream::in | stringstream::out );
			charConvertStream << charLabelX << " " << charLabelY;
			charConvertStream >> charLabelXString;
			charConvertStream >> charLabelYString;

			// Determine where the nucleotide number label should be, if necessary.
			string boxXString, boxYString, boxWidthString;
			string numLabelXString, numLabelYString, numString;
			if( (labelX != -1.0) && (labelY != -1.0) ) {

				// Determine how many places the label takes up, based on the number.
				int places =
					( num >= 10000 ) ? 5 :
					( num >= 1000 ) ? 4 :
					( num >= 100 ) ? 3 :
					( num >= 10 ) ? 2 : 1;

				// Calculate the number width.
				double boxWidth = (double)places * (double)GLYPHSIZE;

				// Calculate the position of the label box and its text.
				double boxX = labelX - ( boxWidth / 2.0 );
				double boxY = labelY - ( (double)NUC_LABEL_HEIGHT / 2.0 );
				double numLabelX = boxX;
				double numLabelY = boxY + (double)TEXTSIZE - (double)NUC_POS_ADJUSTMENT;

				// Convert the label box data into strings.
				stringstream labelBoxStream( stringstream::in | stringstream::out );
				labelBoxStream << boxX << " " << boxY << " " << boxWidth;
				labelBoxStream >> boxXString;
				labelBoxStream >> boxYString;
				labelBoxStream >> boxWidthString;

				// Convert the nucleotide label data into strings.
				stringstream nucStream( stringstream::in | stringstream::out );
				nucStream << numLabelX << " " << numLabelY << " " << num;
				nucStream >> numLabelXString;
				nucStream >> numLabelYString;
				nucStream >> numString;
			}

			// Open the conversion output stream.
			stringstream out( stringstream::in | stringstream::out );

			// Draw a backbone segment, if necessary.
			if( (backX != -1) && (backY != -1) ) {
				string backData = ( !isSVG ) ? LINE_PS : LINE_SVG;
				while( ( index = backData.find( COLOR ) ) != string::npos ) { backData = backData.replace( index, COLOR.size(), getColorString( BLACK, isSVG ) ); }
				while( ( index = backData.find( LINEWEIGHT ) ) != string::npos ) { backData = backData.replace( index, LINEWEIGHT.size(), "1" ); }
				while( ( index = backData.find( STARTX ) ) != string::npos ) { backData = backData.replace( index, STARTX.size(), nucXString ); }
				while( ( index = backData.find( STARTY ) ) != string::npos ) { backData = backData.replace( index, STARTY.size(), nucYString ); }
				while( ( index = backData.find( ENDX ) ) != string::npos ) { backData = backData.replace( index, ENDX.size(), backXString ); }
				while( ( index = backData.find( ENDY ) ) != string::npos ) { backData = backData.replace( index, ENDY.size(), backYString ); }
				out << backData << " ";
			}

			// Draw a number label for the nucleotide, if necessary.
			if( (labelX != -1) && (labelY != -1) ) {
				string lineData = ( !isSVG ) ? LINE_PS : LINE_SVG;
				while( ( index = lineData.find( COLOR ) ) != string::npos ) { lineData = lineData.replace( index, COLOR.size(), getColorString( BLACK, isSVG ) ); }
				while( ( index = lineData.find( LINEWEIGHT ) ) != string::npos ) { lineData = lineData.replace( index, LINEWEIGHT.size(), "1" ); }
				while( ( index = lineData.find( STARTX ) ) != string::npos ) { lineData = lineData.replace( index, STARTX.size(), nucXString ); }
				while( ( index = lineData.find( STARTY ) ) != string::npos ) { lineData = lineData.replace( index, STARTY.size(), nucYString ); }
				while( ( index = lineData.find( ENDX ) ) != string::npos ) { lineData = lineData.replace( index, ENDX.size(), labelXString ); }
				while( ( index = lineData.find( ENDY ) ) != string::npos ) { lineData = lineData.replace( index, ENDY.size(), labelYString ); }
				string rectData = ( !isSVG ) ? RECTANGLE_PS : RECTANGLE_SVG;
				while( ( index = rectData.find( COLOR ) ) != string::npos ) { rectData = rectData.replace( index, COLOR.size(), getColorString( WHITE, isSVG ) ); }
				while( ( index = rectData.find( LOCX ) ) != string::npos ) { rectData = rectData.replace( index, LOCX.size(), boxXString ); }
				while( ( index = rectData.find( LOCY ) ) != string::npos ) { rectData = rectData.replace( index, LOCY.size(), boxYString ); }
				while( ( index = rectData.find( WIDTH ) ) != string::npos ) { rectData = rectData.replace( index, WIDTH.size(), boxWidthString ); }
				while( ( index = rectData.find( HEIGHT ) ) != string::npos ) { rectData = rectData.replace( index, HEIGHT.size(), NUC_LABEL_HEIGHT_STRING ); }
				string textData = ( !isSVG ) ? TEXT_PS : TEXT_SVG;
				while( ( index = textData.find( COLOR ) ) != string::npos ) { textData = textData.replace( index, COLOR.size(), getColorString( BLACK, isSVG ) ); }
				while( ( index = textData.find( LOCX ) ) != string::npos ) { textData = textData.replace( index, LOCX.size(), numLabelXString ); }
				while( ( index = textData.find( LOCY ) ) != string::npos ) { textData = textData.replace( index, LOCY.size(), numLabelYString ); }
				while( ( index = textData.find( TEXTSTRING ) ) != string::npos ) { textData = textData.replace( index, TEXTSTRING.size(), numString ); }
				out << lineData << " " << rectData << " " << textData << " ";
			}

			// Draw the nucleotide letter label.
			string loopData = ( !isSVG ) ? CIRCLE_PS : CIRCLE_SVG;
			string circleOutline = ( circleNucs ) ? getColorString( BLACK, isSVG ) : getColorString( WHITE, isSVG );
			while( ( index = loopData.find( LINEWEIGHT ) ) != string::npos ) { loopData = loopData.replace( index, LINEWEIGHT.size(), "1" ); }
			while( ( index = loopData.find( LOCX ) ) != string::npos ) { loopData = loopData.replace( index, LOCX.size(), nucXString ); }
			while( ( index = loopData.find( LOCY ) ) != string::npos ) { loopData = loopData.replace( index, LOCY.size(), nucYString ); }
			while( ( index = loopData.find( BACKGROUND ) ) != string::npos ) { loopData = loopData.replace( index, BACKGROUND.size(), getColorString( WHITE, isSVG ) ); }
			while( ( index = loopData.find( OUTLINE ) ) != string::npos ) { loopData = loopData.replace( index, OUTLINE.size(), circleOutline ); }
			while( ( index = loopData.find( RADIUS ) ) != string::npos ) { loopData = loopData.replace( index, RADIUS.size(), CIRCLE_RADIUS_STRING ); }
			string nucsData = ( !isSVG ) ? TEXT_PS : TEXT_SVG;
			while( ( index = nucsData.find( COLOR ) ) != string::npos ) { nucsData = nucsData.replace( index, COLOR.size(), getColorString( color, isSVG ) ); }
			while( ( index = nucsData.find( LOCX ) ) != string::npos ) { nucsData = nucsData.replace( index, LOCX.size(), charLabelXString ); }
			while( ( index = nucsData.find( LOCY ) ) != string::npos ) { nucsData = nucsData.replace( index, LOCY.size(), charLabelYString ); }
			while( ( index = nucsData.find( TEXTSTRING ) ) != string::npos ) { nucsData = nucsData.replace( index, TEXTSTRING.size(), nuc ); }
			out << loopData << " " << nucsData;

			// Set the coordinates in the coordinates vector.
			coordinatesTransformed.push_back( out.str() );
		}

		// Convert each piece of pairings data to a string that holds the path which draws the pair, as well as its color.
		vector<string> pairsTransformed;
		for( unsigned int i = 1; i <= pairs.size(); i++ ) {

			// Create the variables necessary to convert data.
			string x1, y1, centerX, centerY, x2, y2;
			istringstream stream( pairs[i-1] );
			string temp, color;

			// Initialize the data.
			stream >> temp;
			stream >> x1;
			stream >> y1;
			stream >> centerX;
			stream >> centerY;
			stream >> x2;
			stream >> y2;
			if( !( stream >> color ) ) { color = BLACK; }

			// Convert the data.
			string data = ( !isSVG ) ? CURVE_PS : CURVE_SVG;
			while( ( index = data.find( CURVEWEIGHT ) ) != string::npos ) { data = data.replace( index, CURVEWEIGHT.size(), "3" ); }
			while( ( index = data.find( X1 ) ) != string::npos ) { data = data.replace( index, X1.size(), x1 ); }
			while( ( index = data.find( Y1 ) ) != string::npos ) { data = data.replace( index, Y1.size(), y1 ); }
			while( ( index = data.find( CONTROLX ) ) != string::npos ) { data = data.replace( index, CONTROLX.size(), centerX ); }
			while( ( index = data.find( CONTROLY ) ) != string::npos ) { data = data.replace( index, CONTROLY.size(), centerY ); }
			while( ( index = data.find( X2 ) ) != string::npos ) { data = data.replace( index, X2.size(), x2 ); }
			while( ( index = data.find( Y2 ) ) != string::npos ) { data = data.replace( index, Y2.size(), y2 ); }
			while( ( index = data.find( COLOR ) ) != string::npos ) { data = data.replace( index, COLOR.size(), getColorString( color, isSVG ) ); }
			pairsTransformed.push_back( data );
		}

		// Determine the dimension by which the image should be scaled.
		double absoluteXBound = ( !isSVG ) ? XBOUND_PS : XBOUND_SVG;
		double absoluteYBound = ( !isSVG ) ? YBOUND_PS : YBOUND_SVG;
		double adjustedYBound = absoluteYBound - ( TEXTSIZE_BIGGER * ( legend.size() + 2 ) );
		double scaleX = absoluteXBound / maxX;
		double scaleY = adjustedYBound / maxY;
		double scale = min( scaleX, scaleY );

		// Open the scaled area.
		string scaleStart = ( !isSVG ) ? SCALE_OPEN_PS : SCALE_OPEN_SVG;
		string scaleString;
		stringstream scaleStream( stringstream::in | stringstream::out );
		scaleStream << scale;
		scaleStream >> scaleString;
		while( ( index = scaleStart.find( SCALEFACTOR ) ) != string::npos ) { scaleStart = scaleStart.replace( index, SCALEFACTOR.size(), scaleString ); }
		out << scaleStart << endl;

		// Draw the pairs and coordinates into the image.
		for( unsigned int i = 1; i <= pairsTransformed.size(); i++ ) { out << pairsTransformed[i-1] << endl; }
		for( unsigned int i = 1; i <= coordinatesTransformed.size(); i++ ) { out << coordinatesTransformed[i-1] << endl; }

		// Close the scaled area.
		string scaleClose = ( !isSVG ) ? SCALE_CLOSE_PS : SCALE_CLOSE_SVG;
		out << scaleClose << endl;

		// If any extra data is to be written into the file, do it.
		// Note that the data is assumed to be formatted properly for the given input file type, and is not scaled.
		int extraLines = extras.size();
		if( extraLines != 0 ) {
			for( int i = 1; i <= extraLines; i++ ) {
				out << extras[i-1] << endl;
			}
		}

		// Write out the beginning of the legend and/or description text resizing group.
		string resizeStart = ( !isSVG ) ? LEGEND_RESIZE_START_PS : LEGEND_RESIZE_START_SVG;
		out << resizeStart << endl;

		// If a legend should be written into the image, do so.
		if( legend.size() > 0 ) {

			// Write each entry into the legend.
			unsigned int legendPlusOne = legend.size() + 1;
			for( unsigned int i = 1; i <= legend.size(); i++ ) {

				// Determine the Y location for this entry string.
				unsigned int next = i - 1;
				double location = absoluteYBound - ( TEXTSIZE_LEGEND_BIGGER * (legendPlusOne - next) );
				stringstream entryStream( stringstream::in | stringstream::out );
				entryStream << location;
				string entryYString = entryStream.str();

				// Build the next legend entry string and write it to the file.
				string entry = ( !isSVG ) ? TEXT_PS : TEXT_SVG;
				while( ( index = entry.find( COLOR ) ) != string::npos ) { entry = entry.replace( index, COLOR.size(), getColorString( legendColors[next], isSVG ) ); }
				while( ( index = entry.find( LOCX ) ) != string::npos ) { entry = entry.replace( index, LOCX.size(), BORDER_STRING ); }
				while( ( index = entry.find( LOCY ) ) != string::npos ) { entry = entry.replace( index, LOCY.size(), entryYString ); }
				while( ( index = entry.find( TEXTSTRING ) ) != string::npos ) { entry = entry.replace( index, TEXTSTRING.size(), legend[next] ); }
				out << entry << endl;
			}
		}

		// Write the image description, if necessary.
		if( description != "" ) {

			// Trim the description and remove any control characters.
			const size_t whiteIndex1 = description.find_first_not_of( " \n\r\t" );
			if( whiteIndex1 != string::npos ) { description = description.substr( whiteIndex1 ); }
			const size_t whiteIndex2 = description.find_last_not_of( " \n\r\t" );
			if( whiteIndex2 != string::npos ) { description = description.substr( 0, whiteIndex2 + 1 ); }
			while( (index = description.find( '(', 0 )) != string::npos ) { description = description.erase( index, 1 ); }
			while( (index = description.find( ')', 0 )) != string::npos ) { description = description.erase( index, 1 ); }
			while( (index = description.find( '<', 0 )) != string::npos ) { description = description.erase( index, 1 ); }
			while( (index = description.find( '>', 0 )) != string::npos ) { description = description.erase( index, 1 ); }

			// If the sanitized description is still a nonzero length, write it.
			if( description != "" ) {

				// Create the description Y location string.
				string descYString;
				stringstream descStream( stringstream::in | stringstream::out );
				descStream << ( absoluteYBound - TEXTSIZE_LEGEND_BIGGER + 5 );
				descStream >> descYString;

				// Trim the description if necessary so it fits on one line.
				int maxDescriptionSize = ( !isSVG ) ? DESC_PS : DESC_SVG;
				if( (int)description.size() > maxDescriptionSize ) {
					description = description.substr( 0, maxDescriptionSize - 3 );
					description += "...";
				}

				// Create the description string and write it in the file.
				string descData = ( !isSVG ) ? TEXT_PS : TEXT_SVG;
				while( ( index = descData.find( COLOR ) ) != string::npos ) { descData = descData.replace( index, COLOR.size(), getColorString( BLACK, isSVG ) ); }
				while( ( index = descData.find( LOCX ) ) != string::npos ) { descData = descData.replace( index, LOCX.size(), BORDER_STRING ); }
				while( ( index = descData.find( LOCY ) ) != string::npos ) { descData = descData.replace( index, LOCY.size(), descYString ); }
				while( ( index = descData.find( TEXTSTRING ) ) != string::npos ) { descData = descData.replace( index, TEXTSTRING.size(), description ); }
				out << descData << endl;
			}
		}

		// Write out the end of the legend and/or description text resizing group.
		string resizeEnd = ( !isSVG ) ? LEGEND_RESIZE_END_PS : LEGEND_RESIZE_END_SVG;
		out << resizeEnd << endl;
	}

	// Finish the file.
	string endMarker = ( !isSVG ) ? END_MARKER_PS : END_MARKER_SVG;
	out << endMarker << endl;

	// Close the output stream.
	out.close();
}

map<int,int> StructureImageHandler::getPairs(){
	map<int, int> pairingMap;
	for( unsigned int i = 1; i <= pairs.size(); i++ ) {

		// Split off the string at the beginning of the pairing data that shows pairing indices.
		string pairDataString;
		stringstream pairSplitStream( pairs[i-1] );
		pairSplitStream >> pairDataString;

		// Get the data from the pairing string.
		string indexString, pairString;
		stringstream indexSplitStream( pairDataString );
		getline( indexSplitStream, indexString, '-' );
		getline( indexSplitStream, pairString );

		// Get the nucleotides that make the next pair.
		int index, pair;
		stringstream indexStream( indexString );
		indexStream >> index;
		stringstream pairStream( pairString );
		pairStream >> pair;

		// Save the pair in both directions.
		pairingMap[index] = pair;
		pairingMap[pair] = index;
	}
    return pairingMap;
}

///////////////////////////////////////////////////////////////////////////////
// Write a Postscript image.
///////////////////////////////////////////////////////////////////////////////
void StructureImageHandler::writePostscript( string file, bool append ,int pagenumber, int pages) {
	writeImageFile( file, append, false , pages);
}

///////////////////////////////////////////////////////////////////////////////
// Write an SVG image.
///////////////////////////////////////////////////////////////////////////////
void StructureImageHandler::writeSVG( string file ) {
	writeImageFile( file, false, true );
}

// Returns the description (which was read in from the structure file)
string StructureImageHandler::getDescription() { return description; }

// Set the description that is written in the legend of the output figure.
void StructureImageHandler::setDescription(const string desc) { description = desc; }

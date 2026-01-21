/*
 * A program that calculates a probability dot plot.
 * This class can write image output to Postscript or SVG.
 * It can also write to a dot plot text file.
 *
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter.
 */

#include "ProbabilityPlot.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
ProbabilityPlot::ProbabilityPlot() {

	// Initialize the calculation type description.
	calcType = "Probability dot plot";

	// Initialize all alternate input file format flags to false.
	logPlot = false;
	matrixPlot = false;

	// Initialize SVG image writing, and text file writing to false.
	isSVG = false;
	writeText = false;

	// Initialize the number of legend entries and the dot plot bounds.
	entries = ENTRIES_DEFAULT;
	minBound = numeric_limits<double>::infinity() * -1;
	maxBound = numeric_limits<double>::infinity();
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool ProbabilityPlot::parse( int argc, char* argv[] ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "ProbabilityPlot" );
	parser->addParameterDescription( "input file", "The name of the input file that holds base pairing probabilities. This file may be one of the following file types. 1) Partition function save file (binary file). 2) Matrix file (plain text). Note that in order to use a matrix file, the \"--matrix\" flag must be specified. 3) Dot plot file (plain text). This file is in the standard format exported by all dot plot interfaces when the \"text\" option is used. Note that in order to use a dot plot file, the \"--log10\" flag must be specified." );
	parser->addParameterDescription( "output file", "The name of a file to which output will be written. Depending on the options selected, this may be one of the following file types. 1) A Postscript image file. 2) An SVG image file. 3) A plain text file." );

	// Add the entries option.
	vector<string> entriesOptions;
	entriesOptions.push_back( "-e" );
	entriesOptions.push_back( "-E" );
	entriesOptions.push_back( "--entries" );
	parser->addOptionFlagsWithParameters( entriesOptions, "Specifies the number of colors in the dot plot. Default is " + ENTRIES_DEFAULT_STRING + " colors. Minimum is " + ENTRIES_MINIMUM_STRING + " colors. Maximum is " + ENTRIES_MAXIMUM_STRING + " colors." );

	// Add the log probabilities option.
	vector<string> logOptions;
	logOptions.push_back( "--log10" );
	parser->addOptionFlagsNoParameters( logOptions, "Specifies that the input file format is a dot plot text file of log10 base pair probabilities. Giving this flag with one of the text options would give a file identical to the input file." );

	// Add the matrix option.
	vector<string> matrixOptions;
	matrixOptions.push_back( "--matrix" );
	parser->addOptionFlagsNoParameters( matrixOptions, "Specifies that the input file format is a plain text matrix of base pair probabilities." );

	// Add the maximum option.
	vector<string> maximumOptions;
	maximumOptions.push_back( "-max" );
	maximumOptions.push_back( "-MAX" );
	maximumOptions.push_back( "--maximum" );
	parser->addOptionFlagsWithParameters( maximumOptions, "Specifies the maximum value that is viewable in the plot. Default is the largest allowable point in a given data file. If the given value is greater than the default, it is ignored." );

	// Add the minimum option.
	vector<string> minimumOptions;
	minimumOptions.push_back( "-min" );
	minimumOptions.push_back( "-MIN" );
	minimumOptions.push_back( "--minimum" );
	parser->addOptionFlagsWithParameters( minimumOptions, "Specifies the minimum value that is viewable in the plot. Default is the smallest allowable point in a given data file. If the given value is less than the default, it is ignored." );

	// Add the SVG option.
	vector<string> svgOptions;
	svgOptions.push_back( "--svg" );
	parser->addOptionFlagsNoParameters( svgOptions, "Specify that the output file should be an SVG image file, rather than a Postscript image file." );

	// Add the text option.
	vector<string> textOptions;
	textOptions.push_back( "-t" );
	textOptions.push_back( "-T" );
	textOptions.push_back( "--text" );
	parser->addOptionFlagsNoParameters( textOptions, "Specifies that output should be a dot plot (text) file." );

	// Add the description option.
	vector<string> descOptions;
	descOptions.push_back( "--desc" );
	parser->addOptionFlagsWithParameters( descOptions, "Configure the output of descriptions. Valid values are: (1) \"\" or \"~none\" -- Do not write a description (2) \"~file\" -- If the default description corresponds to a file or path, use only the base name of the path (i.e. no directory or file extension). (3) \"~~\" or \"~default\" -- Use the default description (this is the same as not specifying the flag) (4) \"~list|DESC1|DESC2|DESC3\" -- use this syntax when the output is expected to have more than one plot. It specifies a list of descriptions will be applied in the order given. The character immediately after \"~list\" will be used as the separator (i.e. it need not be the bar (|) character. (5) Any other value is assumed to be the literal description you want to have displayed in the plot legend.");

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		inputFile = parser->getParameter( 1 );
		outputFile = parser->getParameter( 2 );
	}

	// Get the entries option.
	if( !parser->isError() ) {
		parser->setOptionInteger( entriesOptions, entries );
		if( ( entries < ENTRIES_MINIMUM ) || ( entries > ENTRIES_MAXIMUM ) ) {
			if( entries < ENTRIES_MINIMUM ) { cerr << "Too few plot entries given." << endl; }
			else { cerr << "Too many plot entries given." << endl; }
			parser->setError();
		}
	}

	// Get the log probabilities option.
	if( !parser->isError() ) { logPlot = parser->contains( logOptions ); }

	// Get the matrix option.
	if( !parser->isError() ) { matrixPlot = parser->contains( matrixOptions ); }

	// Get the plot bounds options.
	if( !parser->isError() ) {
		parser->setOptionDouble( minimumOptions, minBound );
		parser->setOptionDouble( maximumOptions, maxBound );
		if( minBound > maxBound ) { parser->setErrorSpecialized( "Minimum plot value cannot be greater than maximum plot value." ); }
	}

	// Get the SVG option.
	if( !parser->isError() ) { isSVG = parser->contains( svgOptions ); }

	// Get the text option.
	if( !parser->isError() ) { writeText = parser->contains( textOptions ); }

	// Get the description option.
	if( !parser->isError() && parser->contains( descOptions ) )
		descriptionSettings.parse(parser->getOptionString(descOptions, false)); 

	// If both the log probabilities and the text option are specified, show an error.
	// This is because the output file would be identical to the input file.
	if( logPlot && writeText ) {
		cerr << "Cannot specify log10 probabilities input when writing a text file." << endl
		     << "The output file would be identical to the input file." << endl;
		parser->setError();
	}

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}


///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void ProbabilityPlot::run() {

	/*
	 * Read in dot plot data.
	 */

	// Print a message that says the dot plot data is being read.
	cout << "Reading dot plot data..." << flush;

	// Initialize the dot plot handler and an error string.
	DotPlotHandler* plotHandler = 0;
	string error = "";
	string description = descriptionSettings.getDescription(inputFile, true);

	// If the input file is a matrix file, read raw dot plot probabilities from a square matrix.
	if( matrixPlot ) {

		// Open the text file stream to read data.
		// If the file was not opened correctly, set an error.
		string line;
		ifstream inFile( inputFile.c_str() );
		if( !( inFile.is_open() ) ) { error = "Error opening input file."; }

		// If no error occurred opening the text file, read data from it.
		if( error == "" ) {

			// Read the first line, and put it into a string stream to count
			// how many probabilities are in each matrix row.
			getline( inFile, line );
			stringstream countStream( stringstream::in | stringstream::out );
			countStream << line;
			double nextProbability;
			int probabilityCount = 0;
			while( countStream >> nextProbability ) { probabilityCount++; }

			// Move the file pointer back to the beginning of the file.
			// Then, build the dot plot handler using the probabilities count.
			inFile.seekg( ios::beg );
			plotHandler = new DotPlotHandler( description, probabilityCount, false );

			// Add all possible dots to the dot plot.
			for( int i = 1; i <= probabilityCount; i++ ) {

				// Read in the text line.
				getline( inFile, line );
				stringstream lineStream( stringstream::in | stringstream::out );
				lineStream << line;

				// Read each value from the line and put it in the dots array.
				double value;
				for( int j = 1; j <= probabilityCount; j++ ) {
					lineStream >> value;
					if( ( value == 0.0 ) || ( value == -0.0 ) ) { value = numeric_limits<double>::infinity(); }
					else { value = -log10( value ); }
					plotHandler->addDotValue( i, j, value );
				}
			}

			// Close the file stream.
			inFile.close();
		}
	}

	// If the input file is a standard dot plot file, read in log10 base pair probabilities.
	else if( logPlot ) {

		// Declare a variable for line reading.
		string line;

		// Attempt to open the dot plot file.
		// If it couldn't be opened, return an error.
		ifstream in( inputFile.c_str() );
		if( !in.is_open() ) { error = "Error opening dot plot file."; }

		// Read in the first line to get the length of the dot plot sequence, and initialize the dot plot handler.
		// If the sequence length cannot be read, set an error.
		int length;
		if( error == "" ) {
			getline( in, line );
			stringstream lengthStream( line );
			if( lengthStream >> length ) { plotHandler = new DotPlotHandler( description, length ); }
			else { error = "Error reading dot plot file."; }
		}

		// Read the rest of the file if no error happened.
		if( error == "" ) {

			// Skip the second line, which has column headers.
			getline( in, line );

			// For all the other lines, read them in.
			int i, j;
			double value;
			while( !in.eof() ) {

				// Get the next line, and put it into a stream.
				getline( in, line );
				stringstream lineStream( line );
				if( !( lineStream.str() == "" ) ) {

					// Read in the first index and set an error if it didn't happen correctly.
					if( !( lineStream >> i ) ) {
						in.close();
						error = "Error reading annotation file.";
					}

					// Read in the second index and set an error if it didn't happen correctly.
					if( !( lineStream >> j ) ) {
						in.close();
						error = "Error reading annotation file.";
					}

					// Read in the probability value and set an error if it didn't happen correctly.
					if( !( lineStream >> value ) ) {
						in.close();
						error = "Error reading annotation file.";
					}

					// Add the value to the dot plot handler.
					plotHandler->addDotValue( i, j, value );
				}
			}
		}

		// Close the data file.
		in.close();
	}

	// Otherwise, read in data from a partition function save file using an RNA strand.
	else {

		// Create the RNA strand.
		RNA* strand = new RNA( inputFile.c_str(), FILE_PFS );
		ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );

		// If no error occurred creating the RNA strand, read data from it.
		if( ( error = checker->returnError() ) == "" ) {

			// Build the dot plot handler using the strand's sequence length.
			int length = strand->GetSequenceLength();
			plotHandler = new DotPlotHandler( description, length );

			// Add all possible dots to the dot plot.
			for( int i = 1; i <= length; i++ ) {
				for( int j = i+1; j <= length; j++ ) {
					double value = strand->GetPairProbability( i, j );
					if( ( value == 0.0 ) || ( value == -0.0 ) ) { value = numeric_limits<double>::infinity(); }
					else { value = -log10( value ); }
					plotHandler->addDotValue( i, j, value );
				}
			}
		}

		// Delete the error checker and RNA strand.
		delete checker;
		delete strand;
	}

	// If no errors occurred, print a message that says the dot plot data was read.
	if( error == "" ) { cout << "done." << endl; }

	/*
	 * Set the plot divider, bounds, and legend.
	 * The full legend is set last so the minimum, maximum, and divider only need to be considered once.
	 */

	if( error == "" ) {
		plotHandler->setLegendDivider( DIVIDER_PROBABILITY );
		plotHandler->setLegendMinimum( minBound );
		plotHandler->setLegendMaximum( maxBound );
		plotHandler->setLegend( entries );
	}

	/*
	 * Write the dot plot file.
	 */

	if( error == "" ) {

		// Print a message saying that the dot plot file is being written.
		if( isSVG ) { cout << "Writing SVG image..." << flush; }
		else if( writeText ) { cout << "Writing text file..." << flush; }
		else { cout << "Writing Postscript image..." << flush; }

		// Write an output file, based on the type of output the user wants.
		if( isSVG ) { plotHandler->writeSVGImage( outputFile ); }
		else if( writeText ) { plotHandler->writeTextFile( outputFile ); }
		else { plotHandler->writePostscriptImage( outputFile ); }

		// Print a message saying that the dot plot file has been written.
		cout << "done." << endl;
	}

	/*
	 * Clean up resources and show the result of the dot plot run.
	 */

	// Delete the dot plot handler, if necessary.
	if( plotHandler != 0 ) { delete plotHandler; }

	// Print confirmation of run finishing.
	if( error == "" ) { cout << calcType << " complete." << endl; }
	else { cerr << error << endl << calcType << " complete with errors." << endl; }
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {
	ProbabilityPlot* runner = new ProbabilityPlot();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
}

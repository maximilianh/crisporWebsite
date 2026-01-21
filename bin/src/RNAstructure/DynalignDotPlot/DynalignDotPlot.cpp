/*
 * A program that calculates a Dynalign energy dot plot.
 * This class can write image output to Postscript or SVG.
 * It can also write to a dot plot text file.
 *
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter.
 */

#include "DynalignDotPlot.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
DynalignDotPlot::DynalignDotPlot() {

	// Initialize the calculation type description.
	calcType = "Dynalign dot plot";

	// Initialize sequence 2 plotting, SVG image writing, and text file writing to false.
	seq2Plot = false;
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
bool DynalignDotPlot::parse( int argc, char* argv[] ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "DynalignDotPlot" );
	parser->addParameterDescription( "Dynalign save file", "A binary save file resulting from a Dynalign folding calculation." );
	parser->addParameterDescription( "output file", "The name of a file to which output will be written. Depending on the options selected, this may be one of the following file types. 1) A Postscript image file. 2) An SVG image file. 3) A plain text file." );

	// Add the entries option.
	vector<string> entriesOptions;
	entriesOptions.push_back( "-e" );
	entriesOptions.push_back( "-E" );
	entriesOptions.push_back( "--entries" );
	parser->addOptionFlagsWithParameters( entriesOptions, "Specifies the number of colors in the dot plot. Default is " + ENTRIES_DEFAULT_STRING + " colors. Minimum is " + ENTRIES_MINIMUM_STRING + " colors. Maximum is " + ENTRIES_MAXIMUM_STRING + " colors." );

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

	// Add the seccond sequence option.
	vector<string> seq2Options;
	seq2Options.push_back( "-s2" );
	seq2Options.push_back( "-S2" );
	seq2Options.push_back( "--sequence2" );
	parser->addOptionFlagsNoParameters( seq2Options, "Specifies that the dot plot should be the second sequence. If no sequence is specified, the plot is the first sequence." );

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

	// Get the plot bounds options.
	if( !parser->isError() ) {
		parser->setOptionDouble( minimumOptions, minBound );
		parser->setOptionDouble( maximumOptions, maxBound );
		if( minBound > maxBound ) {
			cerr << "Minimum plot value cannot be greater than maximum plot value." << endl;
			parser->setError();
		}
	}

	// Get the sequence number option.
	if( !parser->isError() ) { seq2Plot = parser->contains( seq2Options ); }

	// Get the SVG option.
	if( !parser->isError() ) { isSVG = parser->contains( svgOptions ); }

	// Get the text option.
	if( !parser->isError() ) { writeText = parser->contains( textOptions ); }

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void DynalignDotPlot::run() {

	/*
	 * Make sure the dot plot data can be accessed.
	 */

	// Print a message that says the dot plot data file is being checked.
	cout << "Checking dot plot input file..." << flush;

	// Initialize the dot plot handler and an error string.
	DotPlotHandler* plotHandler = 0;
	string error = "";

	// Create the Dynalign object that holds the dot plot data.
	Dynalign_object* object = new Dynalign_object( inputFile.c_str() );
	ErrorChecker<Dynalign_object>* checker = new ErrorChecker<Dynalign_object>( object );
	error = checker->returnError();
	if( error == "" ) { cout << "done." << endl; }
	else {

		cerr << error;
	}

	/*
	 * Prepare the dot plot data.
	 */

	if( error == "" ) {

		// Print a message saying that the dot plot data is being prepared.
		cout << "Preparing dot plot data..." << flush;

		// Determine the proper RNA strand and strand number to use from the Dynalign object.
		int number = ( seq2Plot ) ? 2 : 1;
		RNA* strand = ( seq2Plot ) ? object->GetRNA2() : object->GetRNA1();

		// Build the dot plot handler using the strand's sequence length.
		int length = strand->GetSequenceLength();
		plotHandler = new DotPlotHandler( inputFile, length );

		// Add all possible dots to the dot plot.
		for( int i = 1; i <= length; i++ ) {
			for( int j = i; j <= length; j++ ) {
				double energy = object->GetBestPairEnergy( number, i, j );
				if( energy > 0.0 ) { energy = numeric_limits<double>::infinity(); }
				else {
					stringstream stream( stringstream::in | stringstream::out );
					stream << fixed << setprecision( 1 ) << energy;
					stream >> energy;
				}
				plotHandler->addDotValue( i, j, energy );
			}
		}

		// Print a message saying that the dot plot data has been prepared.
		cout << "done." << endl;
	}

	/*
	 * Set the plot divider, bounds, and legend.
	 * The full legend is set last so the minimum, maximum, and divider only need to be considered once.
	 */

	if( error == "" ) {
		plotHandler->setLegendDivider( DIVIDER_ENERGY );
		plotHandler->setLegendMinimum( minBound );
		plotHandler->setLegendMaximum( maxBound );
		plotHandler->setLegend( entries );
	}

	/*
	 * Set the description based on the strand number.
	 */

	if( error == "" ) {
		RNA* strand = ( seq2Plot ) ? object->GetRNA2() : object->GetRNA1();
		string seqDescription = strand->GetCommentString();
		string seqString = ( seq2Plot ) ? "2" : "1";
		string fullDescription =
			inputFile + " -- Sequence " + seqString + ": " + seqDescription;
		plotHandler->setDescription( fullDescription );
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

	// Delete the error checker and data structure.
	// Also delete the dot plot handler, if necessary.
	delete checker;
	delete object;
	if( plotHandler != 0 ) { delete plotHandler; }

	// Print confirmation of run finishing.
	if( error == "" ) { cout << calcType << " complete." << endl; }
	else { cerr << endl << calcType << " complete with errors." << endl; }
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	DynalignDotPlot* runner = new DynalignDotPlot();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
}

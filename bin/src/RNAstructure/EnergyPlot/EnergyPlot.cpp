/*
 * A program that calculates an energy dot plot.
 * This class can write image output to Postscript or SVG.
 * It can also write to a dot plot text file.
 *
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter.
 */

#include "EnergyPlot.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
EnergyPlot::EnergyPlot() {

	// Initialize the calculation type description.
	calcType = "Energy dot plot";

	// Initialize both SVG image writing and text file writing to false.
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
bool EnergyPlot::parse( int argc, char* argv[] ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "EnergyPlot" );
	parser->addParameterDescription( "folding save file", "A binary save file resulting from a structure folding calculation." );
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

	// Get the plot bounds options.
	if( !parser->isError() ) {
		parser->setOptionDouble( minimumOptions, minBound );
		parser->setOptionDouble( maximumOptions, maxBound );
		if( minBound > maxBound ) {
			cerr << "Minimum plot value cannot be greater than maximum plot value." << endl;
			parser->setError();
		}
	}

	// Get the description option.
	if( !parser->isError() && parser->contains( descOptions ) )
		descriptionSettings.parse(parser->getOptionString(descOptions, false)); 

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
void EnergyPlot::run() {

	/*
	 * Make sure the dot plot data can be accessed.
	 */

	// Print a message that says the dot plot data file is being checked.
	cout << "Checking dot plot input file..." << flush;

	// Initialize the dot plot handler and an error string.
	DotPlotHandler* plotHandler = 0;
	string error = "";

	// Create the RNA strand that holds the dot plot data.
	RNA* strand = new RNA( inputFile.c_str(), FILE_SAV );
	ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );
	error = checker->returnError();
	if( error == "" ) { cout << "done." << endl; }

	/*
	 * Prepare the dot plot data.
	 */

	if( error == "" ) {

		// Print a message saying that the dot plot data is being prepared.
		cout << "Preparing dot plot data..." << flush;

		// Build the dot plot handler using the strand's sequence length.
		int length = strand->GetSequenceLength();
		plotHandler = new DotPlotHandler( descriptionSettings.getDescription(inputFile, true), length );

		// Add all possible dots to the dot plot.
		for( int i = 1; i <= length; i++ ) {
			for( int j = i; j <= length; j++ ) {
				double energy = strand->GetPairEnergy( i, j );
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
	delete strand;
	if( plotHandler != 0 ) { delete plotHandler; }

	// Print confirmation of run finishing.
	if( error == "" ) { cout << calcType << " complete." << endl; }
	else { cerr << endl << calcType << " complete with errors: " << error << endl; }
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	EnergyPlot* runner = new EnergyPlot();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
}

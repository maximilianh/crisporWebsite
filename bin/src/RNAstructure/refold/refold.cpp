/*
 * A program that folds a strand of nucleic acids from a previously folded save file.
 * This strand of nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#include "refold.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
//////////////////////////////////////////////////////////////////////////////
refold::refold() {

	// Initialize the calculation type description.
	calcType = "Refolding";

	// Initialize the maximum percent energy difference.
	percent = 10;

	// Initialize the maximum number of structures.
	maxStructures = 20;

	// Initialize the folding window size.
	windowSize = -1;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool refold::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "refold" );
	parser->addParameterDescription( "save file", "The name of a folding save file that contains structure data." );
	parser->addParameterDescription( "ct file", "The name of a CT file to which output will be written." );

	// Add the maximum number of structures option.
	vector<string> maxStructuresOptions;
	maxStructuresOptions.push_back( "-m" );
	maxStructuresOptions.push_back( "-M" );
	maxStructuresOptions.push_back( "--maximum" );
	parser->addOptionFlagsWithParameters( maxStructuresOptions, "Specify a maximum number of structures. Default is 20 structures." );

	// Add the percent energy difference option.
	vector<string> percentOptions;
	percentOptions.push_back( "-p" );
	percentOptions.push_back( "-P" );
	percentOptions.push_back( "--percent" );
	parser->addOptionFlagsWithParameters( percentOptions, "Specify a maximum percent energy difference. Default is 10 percent (specified as 10, not 0.1)." );

	// Add the window size option.
	vector<string> windowOptions;
	windowOptions.push_back( "-w" );
	windowOptions.push_back( "-W" );
	windowOptions.push_back( "--window" );
	parser->addOptionFlagsWithParameters( windowOptions, "Specify a window size. Default is determined by the length of the sequence." );

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		saveFile = parser->getParameter( 1 );
		ctFile = parser->getParameter( 2 );
	}

	// Get the maximum number of structures option.
	if( !parser->isError() ) {
		parser->setOptionInteger( maxStructuresOptions, maxStructures );
		if( maxStructures <= 0 ) { parser->setError( "maximum number of structures" ); }
	}

	// Get the percent energy difference option.
	if( !parser->isError() ) {
		parser->setOptionDouble( percentOptions, percent );
		if( percent < 0 ) { parser->setError( "percent energy difference" ); }
	}

	// Get the window size option.
	if( !parser->isError() ) {
		parser->setOptionInteger( windowOptions, windowSize );
		bool badWindowSize =
		  ( windowSize < 0 ) &&
		  ( windowSize != -1 );
		if( badWindowSize ) { parser->setError( "window size" ); }
	}

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void refold::run() {

	// Create a variable that handles errors.
	int error = 0;

	/*
	 * Use the constructor for RNA that specifies a filename.
	 * Specify type = 4 (folding save file).
	 * The save file handles the type of nucleic acid being folded, so it can be set as the default (RNA) in the constructor.
	 *
	 * After construction of the strand data structure, create the error checker which monitors for errors.  
	 * Throughout, the error status of the calculation is checked with a variant of the isErrorStatus method, which returns 0 if no error occurred.
	 * The calculation proceeds as long as error = 0.
	 */
	cout << "Initializing nucleic acids..." << flush;
	RNA* strand = new RNA( saveFile.c_str(), FILE_SAV );
	ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );
	error = checker->isErrorStatus();
	if( error == 0 ) { cout << "done." << endl; }

	/*
	 * Set the window size, based on the length of the sequence given as input.
	 * Only do this if window size hasn't been set on the command line.
	 * Use method GetSequenceLength to get the length of the sequence.
	 * The window sizes in relation to the length are hardcoded values.
	 */
	if( windowSize == -1 && error == 0 ) {
		int length = strand->GetSequenceLength();
		windowSize =
			( length > 1200 ) ? 20 :
			( length > 800 ) ? 15 :
			( length > 500 ) ? 11 :
			( length > 300 ) ? 7 :
			( length > 120 ) ? 5 :
			( length > 50 ) ? 3 :
			2;
	}

	/*
	 * Refold the regenerated strand using the ReFoldSingleStrand method.
	 * After the main calculation is complete, use the error checker's isErrorStatus method to check for errors.
	 */
	if( error == 0 ) {

		// Show a message saying that the main calculation has started.
		cout << "Refolding nucleic acids..." << flush;

		// Do the main calculation and check for errors.
		int mainCalcError = strand->ReFoldSingleStrand( percent, maxStructures, windowSize );
		error = checker->isErrorStatus( mainCalcError );

		// If no error occurred, print a message saying that the main calculation is done.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * Write a CT output file using the WriteCt method.
	 * After writing is complete, use the error checker's isErrorStatus method to check for errors.
	 */
	if( error == 0 ) {

		// Show a message saying that the CT file is being written.
		cout << "Writing output ct file..." << flush;

		// Write the CT file and check for errors.
		int writeError = strand->WriteCt( ctFile.c_str() );
		error = checker->isErrorStatus( writeError );

		// If no errors occurred, show a CT file writing completion message.
		if( error == 0 ) { cout << "done." << endl; }
	}

	// Delete the error checker and data structure.
	delete checker;
	delete strand;

	// Print confirmation of run finishing.
	if( error == 0 ) { cout << calcType << " complete." << endl; }
	else { cerr << calcType << " complete with errors." << endl; }
}

////////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	refold* runner = new refold();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}

/*
 * A program that calculates thermodynamic parameters from a set of oligonucleotides.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#include "oligoscreen.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
oligoscreen::oligoscreen() {

	// Initialize the calculation type description.
	calcType = "OligoScreen";

	// Initialize the nucleic acid type.
	isRNA = true;

	// Initialize the calculation temperature.
	temperature = 310.15;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool oligoscreen::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "oligoscreen" );
	parser->addParameterDescription( "list file", "The name of a file that holds the oligo list." );
	parser->addParameterDescription( "report file", "The name of a report file to which output will be written." );

	// Add the DNA option.
	vector<string> dnaOptions;
	dnaOptions.push_back( "-d" );
	dnaOptions.push_back( "-D" );
	dnaOptions.push_back( "--DNA" );
	parser->addOptionFlagsNoParameters( dnaOptions, "This flag only matters if the input file is a sequence file and has been specified as such. Specify that the sequence is DNA, and DNA parameters are to be used. Default is to use RNA parameters." );

	// Add the temperature option.
	vector<string> tempOptions;
	tempOptions.push_back( "-t" );
	tempOptions.push_back( "-T" );
	tempOptions.push_back( "--temperature" );
	parser->addOptionFlagsWithParameters( tempOptions, "Specify the temperature at which calculation takes place in Kelvin. Default is 310.15 K, which is 37 degrees C." );

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		listFile = parser->getParameter( 1 );
		reportFile = parser->getParameter( 2 );
	}

	// Get the DNA option.
	if( !parser->isError() ) { isRNA = !parser->contains( dnaOptions ); }

	// Get the temperature option.
	if( !parser->isError() ) {
		parser->setOptionDouble( tempOptions, temperature );
		if( temperature < 0 ) { parser->setError( "temperature" ); }
	}

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void oligoscreen::run() {

	// Create a variable that handles errors.
	int error = 0;

	/*
	 * Use the constructor for Oligowalk_object that specifies a nucleic acid type.
	 * isRNA identifies whether the strand is RNA (true) or DNA (false).
	 *
	 * After construction of the strand data structure, create the error checker which monitors for errors.  
	 * Throughout, the error status of the calculation is checked with a variant of the isErrorStatus method, which returns 0 if no error occurred.
	 * The calculation proceeds as long as error = 0.
	 */
	cout << "Initializing nucleic acids..." << flush;
	Oligowalk_object* strand = new Oligowalk_object( isRNA );
	ErrorChecker<Oligowalk_object>* checker = new ErrorChecker<Oligowalk_object>( strand );
	error = checker->isErrorStatus();
	if( error == 0 ) { cout << "done." << endl; }

	/*
	 * Set the temperature using the SetTemperature method.
	 * Only set the temperature if a given temperature doesn't equal the default.
	 * If the temperature does need to be set, use the error checker's isErrorStatus method to check for errors.
	 */
	if( ( error == 0 ) && ( temperature != 310.15 ) ) {

		// Show a message saying that the temperature is being set.
		cout << "Setting temperature..." << flush;

		// Set the temperature and check for errors.
		int tempError = strand->SetTemperature( temperature );
		error = checker->isErrorStatus( tempError );

		// If no error occurred, print a message saying that temperature is set.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * Run OligoScreen using the OligoScreen method.
	 * After the main calculation is complete, use the error checker's isErrorStatus method to check for errors.
	 */
	if( error == 0 ) {

		// Show a message saying that the main calculation has started.
		cout << "Running OligoScreen..." << flush;

		// Do the main calculation and check for errors.
		int mainCalcError = strand->OligoScreen( listFile.c_str(), reportFile.c_str() );
		error = checker->isErrorStatus( mainCalcError );

		// If no error occurred, print message that main calculation is done.
		if( error == 0 ) { cout << "done." << endl; }
	}

	// Delete the error checker and data structure.
	delete checker;
	delete strand;

	// Print confirmation of run finishing.
	if( error == 0 ) { cout << calcType << " complete." << endl; }
	else { cerr << calcType << " complete with errors." << endl; }
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	oligoscreen* runner = new oligoscreen();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}

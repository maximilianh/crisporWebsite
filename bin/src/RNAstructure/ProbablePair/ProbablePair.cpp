/*
 * A program that creates structures based on levels of probable pairs.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#include "ProbablePair.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
ProbablePair::ProbablePair() {

	// Initialize the calculation type description.
	calcType = "Calculation of probable structures";

	// Initialize the nucleic acid type.
	isRNA = true;

	// Initialize the sequence flag to false.
	isSequence = false;

	// Initialize the probable pairing threshold.
	threshold = 0.0;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool ProbablePair::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "ProbablePair" );
	parser->addParameterDescription( "input file", "The name of the input file. Depending on the options selected, this may be one of the following file types. 1) Partition function save file (holds probability data). 2) Sequence file (holds raw sequence: .seq or .fasta). Note that in order to use a sequence file, the \"--sequence\" flag must be specified." );
	parser->addParameterDescription( "ct file", "The name of a CT file to which output will be written." );

	// Add the DNA option.
	vector<string> dnaOptions;
	dnaOptions.push_back( "-d" );
	dnaOptions.push_back( "-D" );
	dnaOptions.push_back( "--DNA" );
	parser->addOptionFlagsNoParameters( dnaOptions, "This flag only matters if the input file is a sequence file and has been specified as such. Specify that the sequence is DNA, and DNA parameters are to be used. Default is to use RNA parameters." );

	// Add the sequence option.
	vector<string> sequenceOptions;
	sequenceOptions.push_back( "--sequence" );
	parser->addOptionFlagsNoParameters( sequenceOptions, "Identify the input file format as a sequence file." );

	// Add the threshold option.
	vector<string> thresholdOptions;
	thresholdOptions.push_back( "-t" );
	thresholdOptions.push_back( "-T" );
	thresholdOptions.push_back( "--threshold" );
	parser->addOptionFlagsWithParameters( thresholdOptions, "The threshold at which pairs should be included in a structure. This should be expressed as a number: 0.5 <= x <= 1.0. Default is 0, which signifies that structures should be generated at multiple thresholds: >= 0.99, >= 0.97, >= 0.95, >= 0.90, >= 0.80, >= 0.70, >= 0.60, and >= 0.50." );

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		input = parser->getParameter( 1 );
		ctFile = parser->getParameter( 2 );
	}

	// Get the sequence flag.
	isSequence = parser->contains( sequenceOptions );

	// Get the DNA option.
	// Only do this if the sequence flag was specified.
	if( isSequence ) { isRNA = !parser->contains( dnaOptions ); }

	// Get the threshold option.
	if( !parser->isError() ) {
		parser->setOptionFloat( thresholdOptions, threshold );
		if( threshold < 0.0 ) { parser->setError( "pairing threshold" ); }
	}

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void ProbablePair::run() {

	// Create a variable to handle errors.
	int error = 0;

	/*
	 * Use the constructor for RNA that specifies a filename.
	 * If the input file is a pfs file, specify type = 3 (pfs file).
	 * If the input file is a sequence file, specify type = 2 (sequence file).
	 *
	 * After construction of the strand, create the error checker which monitors the strand for errors.
	 * Then, check for errors with the isErrorStatus function, which returns 0 if no error occurs.
	 * Throughout, the calculation proceeds as long as error = 0.
	 */
	cout << "Initializing nucleic acids..." << flush;
	int type = ( !isSequence ) ? 3 : 2;
	RNA* strand = new RNA( input.c_str(), type, isRNA );
	ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );
	error = checker->isErrorStatus();
	if( error == 0 ) { cout << "done." << endl; }

	/*
	 * If the input file is a sequence file, calculate the partition function first before doing the main calculation.
	 */
	if( error == 0 && isSequence ) {

		// Print a message saying that the partition function has started.
		cout << "Calculating partition function..." << flush;

		// Run the partition function, then check the error status.
		int partError = strand->PartitionFunction();
		error = checker->isErrorStatus( partError );

		// Print out a message saying that partition function is done.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * Predict probable pairs using the PredictProbablePairs method.
	 * After the main calculation is complete, use the error checker's isErrorStatus method to check for errors.
	 */
	if( error == 0 ) {

		// Show a message saying that the main calculation has started.
		cout << "Calculating maximum expected accuracy structures..." << flush;

		// Do the main calculation and check for errors.
		int mainCalcError = strand->PredictProbablePairs( threshold );
		error = checker->isErrorStatus( mainCalcError );

		// If no error occurred, print message that main calculation is done.
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

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	ProbablePair* runner = new ProbablePair();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}

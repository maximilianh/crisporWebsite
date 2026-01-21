/*
 * A program that prints out an ensemble energy of nucleic acid input.
 *
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#include "EnsembleEnergy_Interface.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
EnsembleEnergy_Interface::EnsembleEnergy_Interface() {

	// Initialize the type description.
	calcType = "Ensemble energy calculation";

	// Initialize the nucleic acid type.
	isRNA = true;

	// Initialize the sequence flag to false.
	isSequence = false;

	// Initialize the silent flag to false.
	isSilent = false;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool EnsembleEnergy_Interface::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "EnsembleEnergy" );
	parser->addParameterDescription( "input file", "The name of the input file. Depending on the options selected, this may be one of the following file types. 1) Partition function save file (holds probability data). 2) Sequence file (holds raw sequence: .seq or .fasta). Note that in order to use a sequence file, the \"--sequence\" flag must be specified." );

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

	// Add the silent option.
	vector<string> silentOptions;
	silentOptions.push_back( "-s" );
	silentOptions.push_back( "-S" );
	silentOptions.push_back( "--silent" );
	parser->addOptionFlagsNoParameters( silentOptions, "Suppress all progress messages except the final ensemble energy result. Note that this does NOT suppress errors." );

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		input = parser->getParameter( 1 );
	}

	// Get the sequence flag.
	isSequence = parser->contains( sequenceOptions );

	// Get the DNA option.
	// Only do this if the sequence flag was specified.
	if( isSequence ) { isRNA = !parser->contains( dnaOptions ); }

	// Get the silent option.
	isSilent = parser->contains( silentOptions );

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void EnsembleEnergy_Interface::run() {

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
	if( !isSilent ) { cout << "Initializing nucleic acids..." << flush; }
	int type = ( !isSequence ) ? 3 : 2;
	RNA* strand = new RNA( input.c_str(), type, isRNA );
	ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );
	error = checker->isErrorStatus();
	if( ( !isSilent ) && ( error == 0 ) ) { cout << "done." << endl; }

	/*
	 * If the input file is a sequence file, calculate the partition function first before doing the main calculation.
	 */
	if( error == 0 && isSequence ) {

		// Print a message saying that the partition function has started.
		// Only do this if suppression of messages is not enabled.
		if( !isSilent ) { cout << "Calculating partition function..." << flush; }

		// Run the partition function, then check the error status.
		int partError = strand->PartitionFunction();
		error = checker->isErrorStatus( partError );

		// Print out a message saying that partition function is done.
		// Only do this if suppression of messages is not enabled.
		if( ( !isSilent ) && ( error == 0 ) ) { cout << "done." << endl; }
	}

	/*
	 * Calculate the ensemble energy using the GetEnsembleEnergy method.
	 * After the main calculation is complete, use the error checker's isErrorStatus method to check for errors.
	 */
	if( error == 0 ) {

		// Show a message that the main calculation has started.
		// Only do this if suppression of messages is not enabled.
		if( !isSilent ) { cout << "Calculating ensemble energy..." << flush; }

		// Calculate ensemble energy and check for errors.
		double energy = strand->GetEnsembleEnergy();
		error = checker->isErrorStatus();

		// Print a message that the main calculation is done.
		// Only do this if suppression of messages is not enabled.
		if( ( !isSilent ) && ( error == 0 ) ) { cout << "done." << endl << endl; }

		// If no error occurred, print out the ensemble energy.
		if( error == 0 ) { cout << "Ensemble energy for " << input << ": " << energy << " kcal/mol" << endl << endl; }
	}

	// Delete the error checker and data structure.
	delete checker;
	delete strand;

	// Print confirmation of run finishing.
	if( error == 0 ) {
		if( !isSilent ) { cout << calcType << " complete." << endl; }
	} else { cerr << calcType << " complete with errors." << endl; }
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	EnsembleEnergy_Interface* runner = new EnsembleEnergy_Interface();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}

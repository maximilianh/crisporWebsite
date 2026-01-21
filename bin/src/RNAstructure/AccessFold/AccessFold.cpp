/*
 * A program that folds two strands of nucleic acids.
 * These strands of nucleic acids can be composed of either DNA or RNA.
 *

 *	This algorithm folds each sequence independently to determine accessibility to binding.


 * (c) 2014 Mathews Lab, University of Rochester Medical Center.
 * Written by David H. Mathews, based on bifold by Jessica S. Reuter
 */

#include "AccessFold.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
accessfold::accessfold() {

	// Initialize the calculation type description.
	calcType = "AccessFold";

	// Initialize the nucleic acid type.
	alphabet = DT_RNA;

	// Initialize the maximum internal bulge loop size.
	maxLoop = 30;

	// Initialize the maximum number of structures.
	maxStructures = 20;

	// Initialize the maximum percent energy difference.
	percent = 50.0;

	// Initialize the calculation temperature.
	gamma = 0.4;

	// Initialize the calculation temperature.
	temperature = 310.15;

	// Initialize the folding window size.
	windowSize = 0;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool accessfold::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( calcType );
	parser->addParameterDescription( "seq file 1", "The name of a file containing a first input sequence." );
	parser->addParameterDescription( "seq file 2", "The name of a file containing a second input sequence." );
	parser->addParameterDescription( "ct file", "The name of a CT file to which output will be written." );

	// Add the DNA option.
	vector<string> dnaOptions;
	dnaOptions.push_back( "-d" );
	dnaOptions.push_back( "-D" );
	dnaOptions.push_back( "--DNA" );
	parser->addOptionFlagsNoParameters( dnaOptions, "Specify that the sequence is DNA, and DNA parameters are to be used. Default is to use RNA parameters." );

	
	// Add the maximum loop size option.
	vector<string> loopOptions;
	loopOptions.push_back( "-l" );
	loopOptions.push_back( "-L" );
	loopOptions.push_back( "--loop" );
	parser->addOptionFlagsWithParameters( loopOptions, "Specify a maximum internal/bulge loop size. Default is 30 unpaired numcleotides." );

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

	
	// Add the temperature option.
	vector<string> gammaOptions;
	gammaOptions.push_back( "-g" );
	gammaOptions.push_back( "-G" );
	gammaOptions.push_back( "--gamma" );
	parser->addOptionFlagsWithParameters( gammaOptions, "Specify gamma, the scaling factor for accessibility information.  The default is 0.4 ." );

	// Add the temperature option.
	vector<string> tempOptions;
	tempOptions.push_back( "-t" );
	tempOptions.push_back( "-T" );
	tempOptions.push_back( "--temperature" );
	parser->addOptionFlagsWithParameters( tempOptions, "Specify the temperature at which calculation takes place in Kelvin. Default is 310.15 K, which is 37 degrees C." );

	// Add the window size option.
	vector<string> windowOptions;
	windowOptions.push_back( "-w" );
	windowOptions.push_back( "-W" );
	windowOptions.push_back( "--window" );
	parser->addOptionFlagsWithParameters( windowOptions, "Specify a window size. Default is 0 nucleotides." );

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		seqFile1 = parser->getParameter( 1 );
		seqFile2 = parser->getParameter( 2 );
		ctFile = parser->getParameter( 3 );
	}

	// Get the DNA option.
	if( !parser->isError() && parser->contains( dnaOptions ) )
		alphabet = DT_DNA; // use DNA (unless overridden by alphabet)

	

	// Get the maximum loop size option.
	if( !parser->isError() ) {
		parser->setOptionInteger( loopOptions, maxLoop );
		if( maxLoop < 0 ) { parser->setError( "maximum loop size" ); }
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

	

	// Get the gamma option.
	if( !parser->isError() ) {
		parser->setOptionDouble( gammaOptions, gamma );
	}
	

	// Get the temperature option.
	if( !parser->isError() ) {
		parser->setOptionDouble( tempOptions, temperature );
		if( temperature < 0 ) { parser->setError( "temperature" ); }
	}

	// Get the window size option.
	if( !parser->isError() ) {
		parser->setOptionInteger( windowOptions, windowSize );
		if( windowSize < 0 ) { parser->setError( "window size" ); }
	}

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void accessfold::run() {

	// Create a variable that handles errors.
	int error = 0;

	/*
	 * Use the constructor for HybridRNA that specifies two filenames and types.
	 * For both sequences, specify type = 2 (sequence file).
	 * isRNA identifies whether the strand is RNA (true) or DNA (false).
	 *
	 * After construction of the strand data structure, create the error checker which monitors for errors.
	 * Throughout, the error status of the calculation is checked with a variant of the isErrorStatus method, which returns 0 if no error occurred.
	 * The calculation proceeds as long as error = 0.
	 */
	cout << "Initializing nucleic acids..." << flush;
	HybridRNA* strand = new HybridRNA( seqFile1.c_str(), FILE_SEQ, seqFile2.c_str(), FILE_SEQ, alphabet.c_str() );
	ErrorChecker<HybridRNA>* checker = new ErrorChecker<HybridRNA>( strand );
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
	 * Fold the hybrid strand using the FoldBimolecular method.
	 * If a save file name has been specified, the FoldBimolecular method also has the ability to write a folding save file.
	 * During calculation, monitor progress using the TProgressDialog class and the Start/StopProgress methods of the RNA class.
	 * Neither of these methods require any error checking.
	 * After the main calculation is complete, use the error checker's isErrorStatus method to check for errors.
	 */
	if( error == 0 ) {

		// Show a message saying that the main calculation has started.
		cout << "Folding two strands..." << endl;

		// Don't Create the progress monitor for now.
		//TProgressDialog* progress = new TProgressDialog();
		//strand->SetProgress( *progress );

		// Do the main calculation and check for errors.
		int mainCalcError = strand->AccessFold( gamma, percent, maxStructures, windowSize, maxLoop );
		error = checker->isErrorStatus( mainCalcError );

		// Delete the progress monitor.
		//strand->StopProgress();
		//delete progress;

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

	accessfold* runner = new accessfold();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}

/*
 * A program that finds all suboptimal structures within a predefined small increment for a strand of nucleic acids.
 * This strand of nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#include "AllSub.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
AllSub::AllSub() {

	// Initialize the calculation type description.
	calcType = "Generation of suboptimal structures";

	// Initialize the absolute energy difference.
	absolute = -1.0;



	// Initialize the nucleic acid type.
	alphabet = DT_RNA;

	// Initialize the maximum percent energy difference.
	percent = -1.0;

	// Initialize the calculation temperature.
	temperature = 310.15;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool AllSub::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "AllSub" );
	parser->addParameterDescription( "seq file", "The name of a file containing an input sequence." );
	parser->addParameterDescription( "ct file", "The name of a CT file to which output will be written." );

	// Add the absolute energy difference option.
	vector<string> absoluteOptions;
	absoluteOptions.push_back( "-a" );
	absoluteOptions.push_back( "-A" );
	absoluteOptions.push_back( "--absolute" );
	parser->addOptionFlagsWithParameters( absoluteOptions, "Specify a maximum absolute energy difference. Default is determined by the length of the sequence." );

	// Add the constraint file option.
	vector<string> constraintOptions;
	constraintOptions.push_back( "-c" );
	constraintOptions.push_back( "-C" );
	constraintOptions.push_back( "--constraint" );
	parser->addOptionFlagsWithParameters( constraintOptions, "Specify a constraints file to be applied. Default is to have no constraints applied." );

	// Add the DNA option.
	vector<string> dnaOptions;
	dnaOptions.push_back( "-d" );
	dnaOptions.push_back( "-D" );
	dnaOptions.push_back( "--DNA" );
	parser->addOptionFlagsNoParameters( dnaOptions, "Specify that the sequence is DNA, and DNA parameters are to be used. Default is to use RNA parameters." );

	// Add the Alphabet option.
	vector<string> alphabetOptions;
	alphabetOptions.push_back("--alphabet");
	parser->addOptionFlagsWithParameters(alphabetOptions, "Specify the name of a folding alphabet and associated nearest neighbor parameters. The alphabet is the prefix for the thermodynamic parameter files, e.g. \"rna\" for RNA parameters or \"dna\" for DNA parameters or a custom extended/modified alphabet. The thermodynamic parameters need to reside in the at the location indicated by environment variable DATAPATH. The default is \"rna\" (i.e. use RNA parameters). This option overrides the --DNA flag.");


	// Add the percent energy difference option.
	vector<string> percentOptions;
	percentOptions.push_back( "-p" );
	percentOptions.push_back( "-P" );
	percentOptions.push_back( "--percent" );
	parser->addOptionFlagsWithParameters( percentOptions, "Specify a maximum percent energy difference. Default is determined by the length of the sequence." );

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
		seqFile = parser->getParameter( 1 );
		ctFile = parser->getParameter( 2 );
	}

	// Get the absolute energy difference option.
	if( !parser->isError() ) {
		parser->setOptionDouble( absoluteOptions, absolute );
		bool badAbsolute =
			( absolute != -1 ) &&
			( absolute < 0 );
		if( badAbsolute ) { parser->setError( "absolute energy difference" ); }
	}

	// Get the constraint file option.
	if( !parser->isError() ) { constraintFile = parser->getOptionString( constraintOptions, true ); }

	// Get the DNA option.
	if (!parser->isError() && parser->contains(dnaOptions))
		alphabet = DT_DNA; // use DNA (unless overridden by alphabet)

	// Get the Alphabet option.
	if (!parser->isError() && parser->contains(alphabetOptions))
		alphabet = parser->getOptionString(alphabetOptions, false).c_str();

	// Get the percent energy difference option.
	if( !parser->isError() ) {
		parser->setOptionDouble( percentOptions, percent );
		bool badPercent =
			( percent != -1 ) &&
			( percent < 0 );
		if( badPercent ) { parser->setError( "percent energy difference" ); }
	}

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
void AllSub::run() {

	// Create a variable that handles errors.
	int error = 0;

	/*
	 * Use the constructor for RNA that specifies a filename.
	 * Specify type = 2 (sequence file).
	 * isRNA identifies whether the strand is RNA (true) or DNA (false).
	 *
	 * After construction of the strand data structure, create the error checker which monitors for errors.  
	 * Throughout, the error status of the calculation is checked with a variant of the isErrorStatus method, which returns 0 if no error occurred.
	 * The calculation proceeds as long as error = 0.
	 */
	cout << "Initializing nucleic acids..." << flush;
	RNA* strand = new RNA( seqFile.c_str(), FILE_SEQ, alphabet.c_str());
	ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );
	error = checker->isErrorStatus();
	if( error == 0 ) { cout << "done." << endl; }

	/*
	 * Set the percent difference and the maximum absolute energy difference.
	 * Both differences are based on the length of the sequence, and their values in relation to the sequence length are hardcoded.
	 * Only set these values if they were not set on the command line.
	 * As one or the other of them may or may not be set, each one is checked individually.
	 *
	 * Get the length of the sequence using the GetSequenceLength method.
	 * Since this method only returns a length, error checking is not necessary.
	 */
	if( error == 0 ) {

		// Get the sequence length to identify the hardcoded thresholds.
		int length = strand->GetSequenceLength();

		// Set the maximum percent difference, if applicable.
		if( percent == -1 ) {
			percent =
				( length > 1200 ) ? 5 :
				( length > 800 ) ? 8 :
				( length > 500 ) ? 10 :
				( length > 300 ) ? 15 :
				( length > 120 ) ? 20 :
				( length > 50 ) ? 25 :
				50;
		}

		// Set the absolute energy difference, if applicable.
		if( absolute == -1 ) {
			absolute =
				( length > 1200 ) ? 0.25 :
				( length > 800 ) ? 0.5 :
				( length > 500 ) ? 0.75 :
				( length > 300 ) ? 1 :
				( length > 120 ) ? 1.5 :
				( length > 50 ) ? 3 :
				10;
		}
	}

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
	 * Add constraints if a file has been given for their inclusion.
	 * Read in these constraints with method ReadConstraints.
	 * After constraints are read, use the error checker's isErrorStatus method to check for errors.
	 */
	if( error == 0 && constraintFile != "" ) {

		// Show a message saying that constraints are being applied.
		cout << "Applying constraints..." << flush;

		// Apply constraints and check for errors.
		int constraintError = strand->ReadConstraints( constraintFile.c_str() );
		error = checker->isErrorStatus( constraintError );

		// If no error occurred, print a message saying constraints were included.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * Generate structures using the GenerateAllSuboptimalStructures method.
	 * During calculation, monitor progress using the TProgressDialog class and the Start/StopProgress methods of the RNA class.
	 * Neither of these methods require any error checking.
	 * After the main calculation is complete, use the error checker's isErrorStatus method to check for errors.
	 */
	if( error == 0 ) {

		// Show a message saying that the main calculation has started.
		cout << "Generating suboptimal structures..." << flush;

		// Create the progress monitor.
		TProgressDialog* progress = new TProgressDialog();
		strand->SetProgress( *progress );

		// Do the main calculation and check for errors.
		int mainCalcError = strand->GenerateAllSuboptimalStructures( (float)percent, absolute );
		error = checker->isErrorStatus( mainCalcError );

		// Delete the progress monitor.
		strand->StopProgress();
		delete progress;

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

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	AllSub* runner = new AllSub();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}

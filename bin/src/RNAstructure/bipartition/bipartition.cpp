/*
 * A program that calculates the partition function for two strands of nucleic acids.
 * These strands of nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#include "bipartition.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
bipartition::bipartition() {

	// Initialize the calculation type description.
	calcType = "Bimolecular partition function";


	// Initialize the calculation temperature.
	temperature = 310.15;

	// Initialize the nucleic acid type.
	alphabet = DT_RNA;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool bipartition::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "bipartition" );
	parser->addParameterDescription( "seq file 1", "The name of a file containing a first input sequence." );
	parser->addParameterDescription( "seq file 2", "The name of a file containing a second input sequence." );
	parser->addParameterDescription( "pfs file", "The name of a partition function save file to which output will be written." );

	// Add the DNA option.
	vector<string> dnaOptions;
	dnaOptions.push_back( "-d" );
	dnaOptions.push_back( "-D" );
	dnaOptions.push_back( "--DNA" );
	parser->addOptionFlagsNoParameters( dnaOptions, "Specify that the sequence is DNA, and DNA parameters are to be used. Default is to use RNA parameters." );


	// Add the Alphabet option.
	vector<string> alphabetOptions;
	alphabetOptions.push_back("-a");
	alphabetOptions.push_back("--alphabet");
	parser->addOptionFlagsWithParameters(alphabetOptions, "Specify the name of a folding alphabet and associated nearest neighbor parameters. The alphabet is the prefix for the thermodynamic parameter files, e.g. \"rna\" for RNA parameters or \"dna\" for DNA parameters or a custom extended/modified alphabet. The thermodynamic parameters need to reside in the at the location indicated by environment variable DATAPATH. The default is \"rna\" (i.e. use RNA parameters). This option overrides the --DNA flag.");

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
		seqFile1 = parser->getParameter( 1 );
		seqFile2 = parser->getParameter( 2 );
		pfsFile = parser->getParameter( 3 );
	}

	// Get the DNA option.
	// Get the DNA option.
	if( !parser->isError() && parser->contains( dnaOptions ) )
		alphabet = DT_DNA; // use DNA (unless overridden by alphabet)

	// Get the Alphabet option.
	if (!parser->isError() && parser->contains(alphabetOptions))
		alphabet = parser->getOptionString(alphabetOptions, false).c_str();

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
void bipartition::run() {

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
		// The first sequence, RNA1 will propogate the thermodynamics.
		int tempError = strand->SetTemperature( temperature );
			error = checker->isErrorStatus( tempError );


		// If no error occurred, print a message saying that temperature is set.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * Calculate the bimolecular partition function using the PartitionFunctionBimolecular method.
	 * During calculation, monitor progress using the TProgressDialog class and the Start/StopProgress methods of the RNA class.
	 * Neither of these methods require any error checking.
	 * After the main calculation is complete, use the error checker's isErrorStatus method to check for errors.
	 */
	if( error == 0 ) {

		// Show a message saying that the main calculation has started.
		cout << "Calculating bimolecular partition function..." << endl;

		// Create the progress monitor.
		TProgressDialog* progress = new TProgressDialog();
		strand->SetProgress( *progress );

		// Do the main calculation and check for errors.
		int mainCalcError = strand->PartitionFunctionBimolecular( pfsFile.c_str() );
		error = checker->isErrorStatus( mainCalcError );

		// Delete the progress monitor.
		strand->StopProgress();
		delete progress;

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

	bipartition* runner = new bipartition();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}

/*
 * A program that predicts structures composed of probable base pairs and single-stranded nucleotides.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter.
 */

#include "MaxExpect.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
MaxExpect::MaxExpect() {

	// Initialize the type description.
	calcType = "Calculation of maximum expected accuracy structures";

	// Initialize the weight given to pairs.
	gamma = 1.0;

	// Initialize the nucleic acid type to "rna"
	alphabet = DT_RNA;

	// Initialize the sequence flag to false.
	isSequence = false;

	// Initialize the maximum number of structures.
	maxStructures = 1000;

	// Initialize the maximum percent energy difference.
	percent = 50;

	// Initialize the folding window size.
	windowSize = 5;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool MaxExpect::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "MaxExpect" );
	parser->addParameterDescription( "input file", "The name of the input file. Depending on the options selected, this may be one of the following file types. 1) Partition function save file (holds probability data). 2) Sequence file (holds raw sequence: .seq or .fasta). Note that in order to use a sequence file, the \"--sequence\" flag must be specified." );
	parser->addParameterDescription( "ct file", "The name of a CT file to which output will be written." );

	// Add the DNA option.
	vector<string> dnaOptions;
	dnaOptions.push_back( "-d" );
	dnaOptions.push_back( "-D" );
	dnaOptions.push_back( "--DNA" );
	parser->addOptionFlagsNoParameters( dnaOptions, "This flag only matters if the input file is a sequence file and has been specified as such. Specify that the sequence is DNA, and DNA parameters are to be used. Default is to use RNA parameters." );

	// Add the Alphabet option.
	vector<string> alphabetOptions;
	alphabetOptions.push_back( "-a" );
	alphabetOptions.push_back( "--alphabet" );
	parser->addOptionFlagsWithParameters( alphabetOptions, "Specify the name of a folding alphabet and associated nearest neighbor parameters. The alphabet is the prefix for the thermodynamic parameter files, e.g. \"rna\" for RNA parameters or \"dna\" for DNA parameters or a custom extended/modified alphabet. The thermodynamic parameters need to reside in the at the location indicated by environment variable DATAPATH. The default is \"rna\" (i.e. use RNA parameters). This option overrides the --DNA flag." );

	// Add the gamma option.
	vector<string> gammaOptions;
	gammaOptions.push_back( "-g" );
	gammaOptions.push_back( "-g" );
	gammaOptions.push_back( "--gamma" );
	parser->addOptionFlagsWithParameters( gammaOptions, "Specify the weight which should be placed on base pairs. Default is 1.0." );

	// Add the percent energy difference option.
	vector<string> percentOptions;
	percentOptions.push_back( "-p" );
	percentOptions.push_back( "-P" );
	percentOptions.push_back( "--percent" );
	parser->addOptionFlagsWithParameters( percentOptions, "Specify a maximum percent energy difference. Default is 50 (ie, 50, not 0.5)." );

	// Add the maximum number of structures option.
	vector<string> maxStructuresOptions;
	maxStructuresOptions.push_back( "-s" );
	maxStructuresOptions.push_back( "-S" );
	maxStructuresOptions.push_back( "--structures" );
	parser->addOptionFlagsWithParameters( maxStructuresOptions, "Specify a maximum number of structures. Default is 1000 structures." );

	// Add the sequence option.
	vector<string> sequenceOptions;
	sequenceOptions.push_back( "--sequence" );
	parser->addOptionFlagsNoParameters( sequenceOptions, "Identify the input file format as a sequence file." );

	// Add the window size option.
	vector<string> windowOptions;
	windowOptions.push_back( "-w" );
	windowOptions.push_back( "-W" );
	windowOptions.push_back( "--window" );
	parser->addOptionFlagsWithParameters( windowOptions, "Specify a window size. Default is 5 nucleotides." );

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
	if( !parser->isError() && parser->contains( dnaOptions ) )
		alphabet = DT_DNA; // use DNA (unless overridden by alphabet)

	// Get the Alphabet option.
	if (!parser->isError() && parser->contains(alphabetOptions))
		  alphabet = parser->getOptionString(alphabetOptions, false).c_str();

	// Get the gamma option.
	if( !parser->isError() ) { parser->setOptionDouble( gammaOptions, gamma ); }

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
void MaxExpect::run() {

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
	RNA* strand = new RNA( input.c_str(), type, alphabet.c_str() );
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
	 * Calculate accurate structures using the MaximizeExpectedAccuracy method.
	 * After the main calculation is complete, use the error checker's isErrorStatus method to check for errors.
	 */
	if( error == 0 ) {

		// Show a message saying that the main calculation has started.
		cout << "Calculating maximum expected accuracy structures..." << flush;

		// Do the main calculation and check for errors.
		int mainCalcError = strand->MaximizeExpectedAccuracy( percent, maxStructures, windowSize, gamma );
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

	MaxExpect* runner = new MaxExpect();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}

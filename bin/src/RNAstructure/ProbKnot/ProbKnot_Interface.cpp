/*
 * A program that finds pseudoknots in a strand of nucleic acids.
 * These nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#include "ProbKnot_Interface.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
ProbKnot::ProbKnot() {

	// Initialize the calculation type description.
	calcType = "Prediction of Pseudoknots";

	// Initialize the nucleic acid type.
	isRNA = true;

	// Initialize the sequence flag to false.
	isSequence = false;

	// Initialize the sequence flag to false.
	isEnsemble = false;

	// Initialize the number of iterations the calculation should undergo.
	iterations = 1;

	// Initialize the minimum helix length.
	minHelixLength = 3;

	//set the default threshold to zero, i.e. accept all pairs regardless of pairing probability
	threshold = 0.0;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool ProbKnot::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "ProbKnot" );
	parser->addParameterDescription( "input file", "The name of the input file. Depending on the options selected, this may be one of the following file types. 1) Partition function save file (holds probability data). 2) Sequence file (holds raw sequence: .seq or .fasta). Note that in order to use a sequence file, the \"--sequence\" flag must be specified. 3) Structure file (holds ensemble of structures). Note that in order to use ensemble of structures file, the \"--ensemble\" flag must be specified." );
	parser->addParameterDescription( "ct file", "The name of a CT file to which output will be written." );

	// Add the DNA option.
	vector<string> dnaOptions;
	dnaOptions.push_back( "-d" );
	dnaOptions.push_back( "-D" );
	dnaOptions.push_back( "--DNA" );
	parser->addOptionFlagsNoParameters( dnaOptions, "This flag only matters if the input file is a sequence file or an ensemble structure file and has been specified as such. Specify that the sequence is DNA, and DNA parameters are to be used. Default is to use RNA parameters." );

    // Add the ensemble option.
	vector<string> ensembleOptions;
	ensembleOptions.push_back( "--ensemble" );
	parser->addOptionFlagsNoParameters( ensembleOptions, "Identify the input file format as an ensemble ctructure file. NOTE: all structures in the file must have the same sequence." );

	// Add the iterations option.
	vector<string> iterationOptions;
	iterationOptions.push_back( "-i" );
	iterationOptions.push_back( "-I" );
	iterationOptions.push_back( "--iterations" );
	parser->addOptionFlagsWithParameters( iterationOptions, "Specify the number of iterations the calculation will undergo. Default is 1 iteration." );

	// Add the minimum helix length option.
	vector<string> helixOptions;
	helixOptions.push_back( "-m" );
	helixOptions.push_back( "-M" );
	helixOptions.push_back( "--minimum" );
	parser->addOptionFlagsWithParameters( helixOptions, "Specify the minimum length accepted for a helix. Default is 3 base pairs." );

	// Add the minimum helix length option.
	vector<string> thresholdOptions;
	thresholdOptions.push_back( "-t" );
	thresholdOptions.push_back( "-T" );
	thresholdOptions.push_back( "--threshold" );
	parser->addOptionFlagsWithParameters( thresholdOptions, "Specify a threshold probability to include a pair.  Default is 0, i.e. all pairs are allowed." );

	// Add the sequence option.
	vector<string> sequenceOptions;
	sequenceOptions.push_back( "--sequence" );
	parser->addOptionFlagsNoParameters( sequenceOptions, "Identify the input file format as a sequence file." );

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		inFile = parser->getParameter( 1 );
		ctFile = parser->getParameter( 2 );
	}

	// Get the sequence flag.
	isSequence = parser->contains( sequenceOptions );

    // Get the ensemble flag.
    isEnsemble = parser->contains( ensembleOptions );

	// Get the DNA option.
	// Only do this if the sequence flag was specified.
	if( isSequence||isEnsemble ) { isRNA = !parser->contains( dnaOptions ); }

	// Get the iterations option.
	if( !parser->isError() ) {
		parser->setOptionInteger( iterationOptions, iterations );
		if( iterations <= 0 ) { parser->setError( "number of iterations" ); }
	}

	// Get the minimum helix length option.
	if( !parser->isError() ) {
		parser->setOptionInteger( helixOptions, minHelixLength );
		if( minHelixLength <= 0 ) { parser->setError( "minimum helix length" ); }
	}

	//Get the thresholf option.
	if( !parser->isError() ) {
		parser->setOptionDouble( thresholdOptions, threshold );
		if( threshold < 0 ) { parser->setError( "threshold" ); }
	}


	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void ProbKnot::run() {

	// Create a variable to handle errors.
	int error = 0;

	/*
	 * Use the constructor for RNA that specifies a filename.
	 * If the input file is a pfs file, specify type = 3 (pfs file).
	 * If the input file is a sequence file, specify type = 2 (sequence file).
	 * If the input file is a ct file, specify type = 1 (ct file).
	 *
	 * After construction of the strand, create the error checker which monitors the strand for errors.
	 * Then, check for errors with the isErrorStatus function, which returns 0 if no error occurs.
	 * Throughout, the calculation proceeds as long as error = 0.
	 */

    // Make sure that both "--sequence" and "--ensemble" flags were not specified
    if( isSequence && isEnsemble ) {
        cerr << "\nERROR: Cannot specify both \"--sequence\" and \"--ensemble\" flags.\nPlease specify either \"--sequence\" flag if the input file contains a sequence,\nor \"--ensemble\" flag if the input file contains an ensemble of structures,\nor no flag if the input file is a partition function save file.\n";
        exit;
    }

	cout << "Initializing nucleic acids..." << flush;

    int type = 3;//partition function save file
    if ( isSequence ) type = 2;//sequence file
    else if ( isEnsemble ) type = 1;//structure file

	RNA* strand = new RNA( inFile.c_str(), type, isRNA );
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
     * Calculate pseudoknot using the ProbKnotFromSample method.
     * This is used when the input file is an ensemble of structures.
     * After the main calculation is complete, use the error checker's isErrorStatus method to check for errors.
     */

    if( error == 0 && isEnsemble ) {
        
        // Show a message saying that the main calculation has started.
		cout << "Calculating pseudoknots..." << flush;

        // Do the main calculation and check for errors.
		int mainCalcError = strand->ProbKnotFromSample( iterations, minHelixLength, threshold );
		error = checker->isErrorStatus( mainCalcError );

		// If no error occurred, print message that main calculation is done.
		if( error == 0 ) { cout << "done." << endl; }
    }

	/*
	 * Calculate pseudoknots using the ProbKnot method.
     * This is used when the input file is either a sequence file or a partition function save file.
	 * After the main calculation is complete, use the error checker's isErrorStatus method to check for errors.
	 */
	else {

		// Show a message saying that the main calculation has started.
		cout << "Calculating pseudoknots..." << flush;

		// Do the main calculation and check for errors.
		int mainCalcError = strand->ProbKnot( iterations, minHelixLength , threshold);
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

		// If the program was run in the --ensemble mode, delete all structures except for the first one
		if( isEnsemble ){
			while( strand->GetStructureNumber()>1 ) strand->GetStructure()->RemoveLastStructure();
		}

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

	ProbKnot* runner = new ProbKnot();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}

/*
 * A program that calculates the partition function for a strand of nucleic acids.
 * This strand of nucleic acids can be composed of either DNA or RNA.
 *
 */

#include "Rsample.h"
#include <vector>
#include <iostream>
#include <cstring>

///////////////////////////////////////////////////////////////////////////////
// Constructor.
//////////////////////////////////////////////////////////////////////////////
RsampleInterface::RsampleInterface() {

	// Initialize the SHAPE Cparam.
	Cparam = 0.5;

	// Initialize the nucleic acid type.
	isRNA = true;

	//Indicate the SHAPE distributions should be used by default
	isDMS = false;

	//Default maximum value for distributions
	max = 1000;

	// Initialize the maximum pairing distance between nucleotides.
	maxDistance = -1;

	// Offset the SHAPE slope.
	Offset = 1.1;

	// Initialize the calculation temperature.
	temperature = 310.15;

	// paired-end, paired-middle and unpaired filenames
	peFile = "";
	pmFile = "";
	upFile = "";

	numsamples = 10000;

	seed = 0; // random seed. 0 indicates that it should be seeded from the current time.

	// Initialize the nucleic acid type.
	alphabet = DT_RNA;

}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool RsampleInterface::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "Rsample" );
	parser->addParameterDescription( "seq file", "The name of a file containing an input sequence." );
	parser->addParameterDescription( "SHAPE file", "The name of a SHAPE reactivity file." );
	parser->addParameterDescription( "pfs file", "The name of a partition function save file to which output will be written." );	

	// Add the DNA option.
	vector<string> dnaOptions;
	dnaOptions.push_back( "-d" );
	dnaOptions.push_back( "-D" );
	dnaOptions.push_back( "--DNA" );
	parser->addOptionFlagsNoParameters( dnaOptions, "Specify that the sequence is DNA, and DNA parameters are to be used. Default is to use RNA parameters." );


	// Add the DMS option.
	vector<string> dmsOptions;
	dmsOptions.push_back("--DMS");
	parser->addOptionFlagsNoParameters(dmsOptions, "Specify that the data is DMS. Default is to use SHAPE parameters.");

	// Add the DMS option.
	vector<string> MaxOptions;
	MaxOptions.push_back("--max");
	parser->addOptionFlagsWithParameters(MaxOptions, "Specify that the data is capped at this value.  This will cap the distribution read from disk.  The default is a cap of 1000, effectively no cap.");

	// Add the Alphabet option.
	vector<string> alphabetOptions;
	alphabetOptions.push_back("-a");
	alphabetOptions.push_back("--alphabet");
	parser->addOptionFlagsWithParameters(alphabetOptions, "Specify the name of a folding alphabet and associated nearest neighbor parameters. The alphabet is the prefix for the thermodynamic parameter files, e.g. \"rna\" for RNA parameters or \"dna\" for DNA parameters or a custom extended/modified alphabet. The thermodynamic parameters need to reside in the at the location indicated by environment variable DATAPATH. The default is \"rna\" (i.e. use RNA parameters). This option overrides the --DNA flag.");


	// Add the Rsample Cparam option.
	vector<string> shapeCparamOptions;
	shapeCparamOptions.push_back( "-C" );
	shapeCparamOptions.push_back( "--cparam" );
	parser->addOptionFlagsWithParameters( shapeCparamOptions, "Specify a C parameter used in Rsample calculations. Default is 0.5 kcal/mol." );

	// Add the Rsample Offset option.
	vector<string> shapeOffsetOptions;
	shapeOffsetOptions.push_back( "-O" );
	shapeOffsetOptions.push_back( "--offset" );
	parser->addOptionFlagsWithParameters( shapeOffsetOptions, "Specify an Offset parameter used in Rsample calculations. Default is 1.1 kcal/mol." );

	// Add the number of samples option for stochastic sampling
	vector<string> numsamplesOptions;
	numsamplesOptions.push_back( "-ns" );
	numsamplesOptions.push_back( "--numsamples" );
	parser->addOptionFlagsWithParameters( numsamplesOptions, "Specify number of samples for stochastic sampling calculation used in Rsample. Default is 10,000." );

	// Specify random seed for stochastic sampling.
	vector<string> seedOption;
	seedOption.push_back( "-s" );
	seedOption.push_back( "--seed" );
	parser->addOptionFlagsWithParameters( seedOption, "Specify a random seed. Default is to set random seed from current time." );

	// Add the Reactivity Files option.
	vector<string> reactFileOptionsUP;
	reactFileOptionsUP.push_back( "-rUP" );
	reactFileOptionsUP.push_back( "--reactUnpaired" );
	parser->addOptionFlagsWithParameters( reactFileOptionsUP, "Give full path to file with unpaired nucleotide reactivities dataset. Default values are in rsample directory in $DATAPATH" );

	// Add the Reactivity Files option.
	vector<string> reactFileOptionsPM;
	reactFileOptionsPM.push_back( "-rPM" );
	reactFileOptionsPM.push_back( "--reactPairedMiddle" );
	parser->addOptionFlagsWithParameters( reactFileOptionsPM, "Give full path to file with middle-of-helix paired nucleotide reactivities dataset. Default values are in rsample directory in $DATAPATH" );

	// Add the Reactivity Files option.
	vector<string> reactFileOptionsPE;
	reactFileOptionsPE.push_back( "-rPE" );
	reactFileOptionsPE.push_back( "--reactPairedEnd" );
	parser->addOptionFlagsWithParameters( reactFileOptionsPE, "Give full path to file with end-of-helix paired nucleotide reactivities dataset. Default values are in rsample directory in $DATAPATH." );

	// Add the maximum pairing distance option.
	vector<string> distanceOptions;
	distanceOptions.push_back( "-md" );
	distanceOptions.push_back( "-MD" );
	distanceOptions.push_back( "--maxdistance" );
	parser->addOptionFlagsWithParameters( distanceOptions, "Specify a maximum pairing distance between nucleotides. Default is no restriction on distance between pairs." );

	// Add the temperature option.
	vector<string> tempOptions;
	tempOptions.push_back( "-t" );
	tempOptions.push_back( "-T" );
	tempOptions.push_back( "--temperature" );
	parser->addOptionFlagsWithParameters( tempOptions, "Specify the temperature at which calculation takes place in Kelvin. Default is 310.15 K, which is 37 degrees C." );
	
	// Add the constraint file option.
	vector<string> constraintOptions;
	constraintOptions.push_back( "-c" );
	constraintOptions.push_back( "-C" );
	constraintOptions.push_back( "--constraint" );
	parser->addOptionFlagsWithParameters( constraintOptions, "Specify a constraints file to be applied. Default is to have no constraints applied." );

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		seqFile = parser->getParameter( 1 );
		SHAPEFile = parser->getParameter( 2 );
		pfsFile = parser->getParameter( 3 );
		parser->setOptionInteger(seedOption, seed);
		isRNA = !parser->contains( dnaOptions );
		isDMS = parser->contains(dmsOptions);

		
	}

	// Get the Alphabet option.
	if (!parser->isError() && parser->contains(alphabetOptions))
		alphabet = parser->getOptionString(alphabetOptions, false).c_str();
     
	// Get the maximum distance option.
	if( !parser->isError() ) {
		parser->setOptionInteger( distanceOptions, maxDistance );
		if( maxDistance < 0 && maxDistance != -1 ) { parser->setError( "maximum pairing distance" ); }
	}

	// Get the R and Offset parameters
	if( !parser->isError() ) {
		if( !parser->isError() ) { parser->setOptionDouble( shapeCparamOptions, Cparam ); }
		if( !parser->isError() ) { parser->setOptionDouble( shapeOffsetOptions, Offset ); }
		if( !parser->isError() ) { parser->setOptionInteger( numsamplesOptions, numsamples ); }
	}

	if( !parser->isError() ) peFile = parser->getOptionString( reactFileOptionsPE );
	if( !parser->isError() ) pmFile = parser->getOptionString( reactFileOptionsPM );
	if( !parser->isError() ) upFile = parser->getOptionString( reactFileOptionsUP );

	// Get the constraint file option.
	if( !parser->isError() ) { constraintFile = parser->getOptionString( constraintOptions, true ); }

	// Get the temperature option.
	if( !parser->isError() ) {
		parser->setOptionDouble( tempOptions, temperature );
		if( temperature < 0 ) { parser->setError( "temperature" ); }
	}

	if (!parser->isError()) {
		parser->setOptionDouble(MaxOptions, max);
	}
	
	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
int RsampleInterface::run() {

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
	cout << "Initializing nucleic acids..." << endl;
	RNA strand(seqFile.c_str(), FILE_SEQ, alphabet.c_str());
	ErrorChecker<RNA> checker( &strand );
	error = checker.isErrorStatus();
	if( error != 0 ) return error;
        cout << "done." << endl;

	/*
	 * Set the temperature using the SetTemperature method.
	 * Only set the temperature if a given temperature doesn't equal the default.
	 * If the temperature does need to be set, use the error checker's isErrorStatus method to check for errors.
	 */
	if( temperature != 310.15 ) {
		// Show a message saying that the temperature is being set.
		cout << "Setting temperature..." << endl;

		// Set the temperature and check for errors.
		int tempError = strand.SetTemperature( temperature );
		error = checker.isErrorStatus( tempError );

		// If no error occurred, print a message saying that temperature is set.
		if( error != 0 ) return error;
	    cout << "done." << endl;
	}

	// Show a message saying that constraints are being applied.
	// Read folding constraints, if applicable.
	if( !constraintFile.empty() ) {
		cout << "Applying constraints..." << flush;
		error = strand.ReadConstraints( constraintFile.c_str() );
		if( checker.isErrorStatus( error ) ) return error;
		cout << "done." << endl;
	}

	/*
	 * Set maximum pairing distance using the ForceMaximumPairingDistance method.
	 */
	if( maxDistance != -1 ) {

		// Show a message saying that the maximum pairing distance is being set.
		cout << "Setting maximum distance between paired nucleotides..." << flush;

		// Set the maximum pairing distance and check for errors.
		int distError = strand.ForceMaximumPairingDistance( maxDistance );
		error = checker.isErrorStatus( distError );

		// If no error occurred, print a message saying that maximum pairing distance was set.
		if( error != 0 ) return error;
		cout << "done." << endl;
	}


	/*
	 * Perform Rsample calculation.
	 * During calculation, monitor progress using the TProgressDialog class and the Start/StopProgress methods of the RNA class.
	 * After the main calculation is complete, use the error checker's isErrorStatus method to check for errors.
	 */
	// Show a message saying that the main calculation has started.
	cout << "Performing Rsample calculation..." << endl;

	// Create the progress monitor.
	TProgressDialog progress;
	strand.SetProgress( progress );

	vector<double> shapedata(strand.GetSequenceLength());
	error = ReadRestraints(shapedata, SHAPEFile.c_str(),max);
	if (error != 0) {
		// Show a message saying that the maximum pairing distance is being set.
		cerr << RNA::GetErrorMessage(error) << "File: " << SHAPEFile << endl;
		return error;
	}


    RsampleData rsdata(isDMS, max, upFile.c_str(), peFile.c_str(), pmFile.c_str());
	if (rsdata.ErrorCode != 0) {
		cerr << RsampleData::GetErrorMessage(rsdata.ErrorCode) << endl;
		return 2; // failed to open file.
	}
	
	// Do the main calculation and check for errors.
	error = strand.Rsample(shapedata, rsdata, seed, pfsFile.c_str(), Cparam, Offset, numsamples);
	error = checker.isErrorStatus( error );
	
	// Delete the progress monitor.
	strand.StopProgress();

	// If no error occurred, print a message saying that the main calculation is done.
	if( error != 0 ) return error;

    cout << "done." << endl;

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {
	RsampleInterface runner;
	if (!runner.parse( argc, argv )) return 1;
	return runner.run() != 0 ? 1 : 0;
}



/*
 * A program that folds two strands of nucleic acids in a duplex.
 * These strands of nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#include "DuplexFold.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
DuplexFold::DuplexFold() {

	// Initialize the calculation type description.
	calcType = "Duplex folding";

	// Initialize the nucleic acid type.
	alphabet = DT_RNA;

	// Initialize the maximum internal bulge loop size.
	maxLoop = 6;

	// Initialize the maximum number of structures.
	maxStructures = 10;

	// Initialize the maximum percent energy difference.
	percent = 40;

	// Initialize the calculation temperature.
	temperature = 310.15;

	// Initialize the folding window size.
	windowSize = 0;

	//By default, the input is two sequences
	IsList = false;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool DuplexFold::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "DuplexFold" );
	parser->addParameterDescription( "seq file 1", "The name of a file containing a first input sequence." );
	parser->addParameterDescription( "seq file 2", "The name of a file containing a second input sequence." );
	parser->addParameterDescription( "ct file", "The name of a CT file to which output will be written." );

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

	// Add the maximum loop size option.
	vector<string> loopOptions;
	loopOptions.push_back( "-l" );
	loopOptions.push_back( "-L" );
	loopOptions.push_back( "--loop" );
	parser->addOptionFlagsWithParameters( loopOptions, "Specify a maximum internal/bulge loop size. Default is 6 unpaired numcleotides." );

	// Add the maximum number of structures option.
	vector<string> maxStructuresOptions;
	maxStructuresOptions.push_back( "-m" );
	maxStructuresOptions.push_back( "-M" );
	maxStructuresOptions.push_back( "--maximum" );
	parser->addOptionFlagsWithParameters( maxStructuresOptions, "Specify a maximum number of structures. Default is 10 structures." );

	// Add the percent energy difference option.
	vector<string> percentOptions;
	percentOptions.push_back( "-p" );
	percentOptions.push_back( "-P" );
	percentOptions.push_back( "--percent" );
	parser->addOptionFlagsWithParameters( percentOptions, "Specify a maximum percent energy difference. Default is 40 percent (specified as 40, not 0.4)." );

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

	// Add the list option.
	vector<string> listOptions;
	listOptions.push_back("--list");
	parser->addOptionFlagsNoParameters(listOptions, "Specify that the input is a list, rather than two sequences.");

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

	// Get the Alphabet option.
	if (!parser->isError() && parser->contains(alphabetOptions))
		alphabet = parser->getOptionString(alphabetOptions, false).c_str();

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

	if (!parser->isError()) {
		IsList = parser->contains(listOptions);
	}

	

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void DuplexFold::run() {

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
	
	
	
	
	if (!IsList) {
		//This is the traditional calculation where the user specifies two sequences.

		HybridRNA* strand = new HybridRNA(seqFile1.c_str(), FILE_SEQ, seqFile2.c_str(), FILE_SEQ, alphabet.c_str());
		ErrorChecker<HybridRNA>* checker = new ErrorChecker<HybridRNA>(strand);

		error = checker->isErrorStatus();
		if (error == 0) { cout << "done." << endl; }

		/*
		 * Set the temperature using the SetTemperature method.
		 * Only set the temperature if a given temperature doesn't equal the default.
		 * If the temperature does need to be set, use the error checker's isErrorStatus method to check for errors.
		 */
		if ((error == 0) && (temperature != 310.15)) {

			// Show a message saying that the temperature is being set.
			cout << "Setting temperature..." << flush;

			// Set the temperature and check for errors.
			int tempError = strand->SetTemperature(temperature);
			error = checker->isErrorStatus(tempError);

			// If no error occurred, print a message saying that temperature is set.
			if (error == 0) { cout << "done." << endl; }
		}

		/*
		 * Fold the hybrid strand using the FoldDuplex method.
		 * After the main calculation is complete, use the error checker's isErrorStatus method to check for errors.
		 */
		if (error == 0) {

			// Show a message saying that the main calculation has started.
			cout << "Folding duplex..." << endl;

			// Create the progress monitor.
			TProgressDialog* progress = new TProgressDialog();
			strand->SetProgress(*progress);

			// Do the main calculation and check for errors.
			int mainCalcError = strand->FoldDuplex(percent, maxStructures, windowSize, maxLoop);
			error = checker->isErrorStatus(mainCalcError);

			// Delete the progress monitor.
			strand->StopProgress();
			delete progress;

			// If no error occurred, print message that main calculation is done.
			if (error == 0) { cout << "done." << endl; }
		}

		/*
		 * Write a CT output file using the WriteCt method.
		 * After writing is complete, use the error checker's isErrorStatus method to check for errors.
		 */
		if (error == 0) {

			// Show a message saying that the CT file is being written.
			cout << "Writing output ct file..." << flush;

			// Write the CT file and check for errors.
			int writeError = strand->WriteCt(ctFile.c_str());
			error = checker->isErrorStatus(writeError);

			// If no errors occurred, show a CT file writing completion message.
			if (error == 0) { cout << "done." << endl; }
		}

		delete checker;
		delete strand;

	}
	else {
		//This is the new style calculation, where the user gives a list of sequences
		ifstream in;
		ofstream out;
		Thermodynamics *datatables;

		string sequence1;
		string sequence2;

		try {
			in.open(seqFile1);
		}
		catch (std::exception* ex) {
			cerr << "Error: " << seqFile1 << " is not readable.\n";
			error = 1;
		}

		out.open(ctFile);

		if (error = 0) {
			datatables = new Thermodynamics(isAlphabetRNA(alphabet.c_str()),alphabet.c_str());

			/*
			 * Set the temperature using the SetTemperature method.
			 * Only set the temperature if a given temperature doesn't equal the default.
			 * If the temperature does need to be set, use the error checker's isErrorStatus method to check for errors.
			 */
			if ((error == 0) && (temperature != 310.15)) {

				// Show a message saying that the temperature is being set.
				cout << "Setting temperature..." << flush;

				// Set the temperature and check for errors.
				error = datatables->SetTemperature(temperature);

				// If no error occurred, print a message saying that temperature is set.
				if (error == 0) { cout << "done." << endl; }
			}

		}

		// Show a message saying that the main calculation has started.
		if (error == 0) cout << "Folding duplexes..." << endl;

		in >> sequence1;
		in >> sequence2;

		while (!in.eof()&&error==0) {
			

			HybridRNA* strand = new HybridRNA(sequence1.c_str(), SEQUENCE_STRING, sequence2.c_str(), SEQUENCE_STRING, alphabet.c_str());
			ErrorChecker<HybridRNA>* checker = new ErrorChecker<HybridRNA>(strand);
			error = strand->GetErrorCode();
			if (error != 0) {
				cerr << "Error in sequence pair: " << sequence1 << " " << sequence2 << "\n";

			}
			

			

			/*
			 * Fold the hybrid strand using the FoldDuplex method.
			 * After the main calculation is complete, use the error checker's isErrorStatus method to check for errors.
			 */
			if (error == 0) {

				


				// Do the main calculation and check for errors.
				int mainCalcError = strand->FoldDuplex(percent, maxStructures, windowSize, maxLoop);
				error = checker->isErrorStatus(mainCalcError);


				
			}

			/*
			 * Write the structures and energies to the output file.
			 */
			if (error == 0) {

				out << strand->GetFreeEnergy(1) << "\n";
				delete checker;
				delete strand;
			}

			in >> sequence1;
			in >> sequence2;
		}
		

		in.close();
		out.close();
		
	}

	if (error == 0) { cout << "done." << endl; }

	// Delete the error checker and data structure.
	
	

	// Print confirmation of run finishing.
	if( error == 0 ) { cout << calcType << " complete." << endl; }
	else { cerr << calcType << " complete with errors." << endl; }
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	DuplexFold* runner = new DuplexFold();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}

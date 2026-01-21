/*
 * A program that folds two strands of nucleic acids.
 * These strands of nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#include "bifold.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
bifold::bifold() {

	// Initialize the calculation type description.
	calcType = "Bimolecular folding";

	// Initialize the nucleic acid type.
	alphabet = DT_RNA;

	// Initialize the allowance of intramolecular pairs.
	// False means that intramolecular pairs are forbidden.
	forbidIntramolecular = false;

	// Initialize the maximum internal bulge loop size.
	maxLoop = 30;

	// Initialize the maximum number of structures.
	maxStructures = 20;

	// Initialize the maximum percent energy difference.
	percent = 50.0;

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
bool bifold::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "bifold" );
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


	// Add the intramolecular option.
	vector<string> intramolecularOptions;
	intramolecularOptions.push_back( "-i" );
	intramolecularOptions.push_back( "-I" );
	intramolecularOptions.push_back( "--intramolecular" );
	parser->addOptionFlagsNoParameters( intramolecularOptions, "Forbid intramolecular pairs (pairs within the same strand). Default is to allow intramolecular pairs." );

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

	// Add the save file option.
	vector<string> saveOptions;
	saveOptions.push_back( "-s" );
	saveOptions.push_back( "-S" );
	saveOptions.push_back( "--save" );
	parser->addOptionFlagsWithParameters( saveOptions, "Specify the name of a save file, needed for dotplots and refolding. Default is not to generate a save file." );

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

	// Get the intramolecular option.
	if( !parser->isError() ) { forbidIntramolecular = parser->contains( intramolecularOptions ); }

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

	// Get the save file option.
	if( !parser->isError() ) { saveFile = parser->getOptionString( saveOptions, false ); }

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
void bifold::run() {

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
			// The first sequence, RNA1 will propogate the thermodynamics.
			int tempError = strand->SetTemperature(temperature);
			error = checker->isErrorStatus(tempError);

			// If no error occurred, print a message saying that temperature is set.
			if (error == 0) { cout << "done." << endl; }
		}

		/*
		 * Set allowance or denial of intramolecular pairs using the SetForbidIntramolecular method.
		 * Since this method only sets a flag, error checking is not necessary.
		 */
		if (error == 0) {
			strand->SetForbidIntramolecular(forbidIntramolecular);
		}

		/*
		 * Fold the hybrid strand using the FoldBimolecular method.
		 * If a save file name has been specified, the FoldBimolecular method also has the ability to write a folding save file.
		 * During calculation, monitor progress using the TProgressDialog class and the Start/StopProgress methods of the RNA class.
		 * Neither of these methods require any error checking.
		 * After the main calculation is complete, use the error checker's isErrorStatus method to check for errors.
		 */
		if (error == 0) {

			// Show a message saying that the main calculation has started.
			cout << "Folding two strands..." << endl;

			// Create the progress monitor.
			TProgressDialog* progress = new TProgressDialog();
			strand->SetProgress(*progress);

			// Do the main calculation and check for errors.
			int mainCalcError = strand->FoldBimolecular(percent, maxStructures, windowSize, saveFile.c_str(), maxLoop);
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

		// Delete the error checker and data structure.
		delete checker;
		delete strand;

	}
	else {
		//IsList is true, this is a list of sequences, rather than a pair of sequences.

		//This is the new style calculation, where the user gives a list of sequences
		ifstream in;
		ofstream out;
		Thermodynamics* datatables;

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
			datatables = new Thermodynamics(isAlphabetRNA(alphabet.c_str()), alphabet.c_str());

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

		while (!in.eof() && error == 0) {


			HybridRNA* strand = new HybridRNA(sequence1.c_str(), SEQUENCE_STRING, sequence2.c_str(), SEQUENCE_STRING, alphabet.c_str());
			ErrorChecker<HybridRNA>* checker = new ErrorChecker<HybridRNA>(strand);
			error = strand->GetErrorCode();
			if (error != 0) {
				cerr << "Error in sequence pair: " << sequence1 << " " << sequence2 << "\n";

			}

			/*
			* Set allowance or denial of intramolecular pairs using the SetForbidIntramolecular method.
			* Since this method only sets a flag, error checking is not necessary.
			*/
			if (error == 0) {
				strand->SetForbidIntramolecular(forbidIntramolecular);
			}




			/*
			 * Fold the hybrid strand using the FoldDuplex method.
			 * After the main calculation is complete, use the error checker's isErrorStatus method to check for errors.
			 */
			if (error == 0) {




				// Do the main calculation and check for errors.
				int mainCalcError = strand->FoldBimolecular(0, 1, 0, "", maxLoop); 
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

	// Print confirmation of run finishing.
	if( error == 0 ) { cout << calcType << " complete." << endl; }
	else { cerr << calcType << " complete with errors." << endl; }
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	bifold* runner = new bifold();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}

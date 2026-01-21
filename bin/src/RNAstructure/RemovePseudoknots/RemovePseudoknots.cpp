/*
 * A program that breaks pseudoknots to allow for easier analysis of a strand of nucleic acids.
 * This strand of nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 * Rewritten by Richard M. Watson
 */

#include "RemovePseudoknots.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
RemovePseudoknots::RemovePseudoknots() {

	// Initialize the calculation type description.
	calcType = "Pseudoknot breakage";

	// Initialize the nucleic acid type.
	isRNA = true;

	// Initialize optimization to true.
	minimum_energy = true;

	// Initialize the calculation temperature.
	temperature = 310.15;

	writeBracket = false;

	useMeaFill = false;

	quiet = false;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool RemovePseudoknots::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "RemovePseudoknots" );
	parser->addParameterDescription( "input ct file", "The name of an input CT file or dot-bracket file containing structure(s) with pseudoknotted pairs.\nIf the name is a hyphen (-) the file will be read from standard input (STDIN)." );
	parser->addParameterDescription( "output ct file", "The name of a CT file to which output will be written.\nIf the name is a hyphen (-), the converted file will be written to standard output (STDOUT) instead of a file." );

	// Add the DNA option.
	vector<string> dnaOptions;
	dnaOptions.push_back( "-d" );
	dnaOptions.push_back( "-D" );
	dnaOptions.push_back( "--DNA" );
	parser->addOptionFlagsNoParameters( dnaOptions, "Specify that the sequence is DNA, and DNA parameters are to be used. Default is to use RNA parameters." );

	// Add the dot-bracket option.
	vector<string> bracketOptions;
	bracketOptions.push_back( "-k" );
	bracketOptions.push_back( "--bracket" );
	parser->addOptionFlagsNoParameters( bracketOptions, "Write the output as a dot-bracket file instead of a CT file. (Note that the input file can be either a CT or bracket file, regardless of this option.)" );

		// Add the mea option
	vector<string> meaOptions;
	meaOptions.push_back( "--MEA" );
	parser->addOptionFlagsNoParameters( meaOptions, "Use the MaxExpect method (MEAFill) to break pseudoknots. This should produce identical results as the 'maximize' option (-m), but will be slower (for testing only)." );

	// Add the maximize pairs option.
	vector<string> maximizeOptions;
	maximizeOptions.push_back( "-m" );
	maximizeOptions.push_back( "-M" );
	maximizeOptions.push_back( "--maximize" );
	parser->addOptionFlagsNoParameters( maximizeOptions, "Specify that the NUMBER of remaining base pairs should be maximized in the pseudoknot-free structure in an energy-agnostic way -- "
		"i.e. with no regard for the folding free energy of the structure.\n"
		"This is usually much faster than the default behavior, which is to break base-pairs such that the final structure(s) have minimum free energy.\n"
		"The default method often maximizes the number of remaining base pairs, but not always." );

	// Add the temperature option.
	vector<string> tempOptions;
	tempOptions.push_back( "-t" );
	tempOptions.push_back( "-T" );
	tempOptions.push_back( "--temperature" );
	parser->addOptionFlagsWithParameters( tempOptions, "Specify the temperature at which calculation takes place in Kelvin. Default is 310.15 K, which is 37 degrees C." );

	vector<string> quietOption = parser->addFlag(false, "-q --quiet", "Suppress unnecessary output. This option is implied when the output file is '-' (STDOUT).");

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		input = parser->getParameter( 1 );
		output = parser->getParameter( 2 );
	}

	// Get the DNA option.
	if( !parser->isError() ) { isRNA = !parser->contains( dnaOptions ); }

	if( !parser->isError() ) { useMeaFill = parser->contains( meaOptions ); }

	// Get the maximize pairs option.
	if( !parser->isError() ) { minimum_energy = ! (useMeaFill || parser->contains(maximizeOptions) ); }

	// Get the temperature option.
	if( !parser->isError() ) {
		parser->setOptionDouble( tempOptions, temperature );
		if( temperature < 0 ) { parser->setError( "temperature" ); }
	}
	
	if( !parser->isError() ) {
		writeBracket = parser->contains(bracketOptions);
		quiet = parser->contains(quietOption) || isStdIoFile(output.c_str()); // suppress unnecessary output if --quiet option is present or if the CT file is output to stdout.
	}

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void RemovePseudoknots::run() {

	// Create a variable that handles errors.
	int error = 0;

	/*
	 * Use the constructor for RNA that specifies a filename.
	 * Specify type = 1 (CT file).
	 * isRNA identifies whether the strand is RNA (true) or DNA (false).
	 *
	 * After construction of the strand data structure, create the error checker which monitors for errors.  
	 * Throughout, the error status of the calculation is checked with a variant of the isErrorStatus method, which returns 0 if no error occurred.
	 * The calculation proceeds as long as error = 0.
	 */
	if (!quiet) cout << "Initializing nucleic acids..." << flush;
	RNA* strand = new RNA( input.c_str(), FILE_CT, isRNA?DT_RNA:DT_DNA, true, true );
	ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );
	error = checker->isErrorStatus();
	if( error == 0 ) { if (!quiet) cout << "done." << endl; }

	/*
	 * Set the temperature using the SetTemperature method.
	 * Only set the temperature if a given temperature doesn't equal the default.
	 * If the temperature does need to be set, use the error checker's isErrorStatus method to check for errors.
	 */
	if( ( error == 0 ) && ( temperature != 310.15 ) ) {

		// Show a message saying that the temperature is being set.
		if (!quiet) cout << "Setting temperature..." << flush;

		// Set the temperature and check for errors.
		int tempError = strand->SetTemperature( temperature );
		error = checker->isErrorStatus( tempError );

		// If no error occurred, print a message saying that temperature is set.
		if( error == 0 ) { if (!quiet) cout << "done." << endl; }
	}

	/*
	 * Break pseudoknots using the BreakPseudoknot method.
	 * Neither of these methods require any error checking.
	 * After the main calculation is complete, use the error checker's isErrorStatus method to check for errors.
	 */
	if( error == 0 ) {

		// Show a message saying that the main calculation has started.
		if (!quiet) cout << "Breaking pseudoknots..." << flush;

		// Do the main calculation and check for errors.
		error = strand->BreakPseudoknot( minimum_energy, 0, !useMeaFill);

		error = checker->isErrorStatus( error );

		// If no error occurred, print a message saying that the main calculation is done.
		if( error == 0 ) { if (!quiet) cout << "done." << endl; }
	}

	/*
	 * Write a CT output file using the WriteCt method.
	 * After writing is complete, use the error checker's isErrorStatus method to check for errors.
	 */
	if( error == 0 ) {

		// Show a message saying that the CT file is being written.
		if (!quiet) cout << "Writing output ct file..." << flush;

		// Write the CT file and check for errors.  (NO_CT_ENERGY_COMMENTS causes energy comments to be removed -- because they are invalid if the structure is changed).
		if (writeBracket)
			error = strand->WriteDotBracket(output.c_str(), -1, DBN_FMT_MULTI_TITLE, CTComments::None);
		else
			error = strand->WriteCt( output.c_str(), false, CTComments::None); // pass in "" for the energyLabel so that energies will NOT be written to the CT file.

		error = checker->isErrorStatus( error );

		// If no errors occurred, show a CT file writing completion message.
		if( error == 0 ) { if (!quiet) cout << "done." << endl; }
	}

	// Delete the error checker and data structure.
	delete checker;
	delete strand;

	// Print confirmation of run finishing.
	if( error == 0 ) { if (!quiet) cout << calcType << " complete." << endl; }
	else { cerr << calcType << " complete with errors." << endl; }
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {
	RemovePseudoknots* runner = new RemovePseudoknots();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}

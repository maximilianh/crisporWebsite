/*
 * A program that scores two structures and outputs their sensitivity and PPV.
 * These structures can be composed of either DNA or RNA.
 *
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#include "Scorer_Interface.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
Scorer_Interface::Scorer_Interface() {

	// Initialize the calculation type description.
	calcType = "Nucleic acid structure scoring";

	// Set the boolean flags to their default values.
	exact = false;
	print = false;

	// The specific structure in the predicted CT file that should be compared to the accepted CT.
	// Set this to -1 as the default, which compares all structures individually.
	number = -1;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool Scorer_Interface::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "scorer" );
	parser->addParameterDescription( "predicted ct", "The name of a file containing CT data for the predicted structure." );
	parser->addParameterDescription( "accepted ct", "The name of a file containing CT data for the accepted structure." );
	parser->addParameterDescription( "output file", "The name of a scores file to which output will be written." );

	// Add the exact option.
	vector<string> exactOptions;
	exactOptions.push_back( "-e" );
	exactOptions.push_back( "-E" );
	exactOptions.push_back( "--exact" );
	parser->addOptionFlagsNoParameters( exactOptions, "Specify exact comparison when structure comparison is scored. Default is to allow flexible pairings." );

	// Add the number option.
	vector<string> numberOptions;
	numberOptions.push_back( "-n" );
	numberOptions.push_back( "-N" );
	numberOptions.push_back( "--number" );
	parser->addOptionFlagsWithParameters( numberOptions, "Specify the index of a particular structure in the predicted CT to be compared with the accepted CT, one-indexed. Default is -1, which signifies all structures output to one file." );

	// Add the print option.
	vector<string> printOptions;
	printOptions.push_back( "-p" );
	printOptions.push_back( "-P" );
	printOptions.push_back( "--print" );
	parser->addOptionFlagsNoParameters( printOptions, "Prints the output file to standard output. This won't override the default behavior of writing to a file." );

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		predicted = parser->getParameter( 1 );
		accepted = parser->getParameter( 2 );
		output = parser->getParameter( 3 );
	}

	// Get the exact option.
	if( !parser->isError() ) { exact = parser->contains( exactOptions ); }

	// Get the number option.
	if( !parser->isError() ) {
		parser->setOptionInteger( numberOptions, number );
		bool badNumber =
			( number != -1 ) &&
			( number < 0 );
		if( badNumber ) { parser->setError( "structure number" ); }
	}

	// Get the print option.
	if( !parser->isError() ) { print = parser->contains( printOptions ); }

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run scoring calculations.
///////////////////////////////////////////////////////////////////////////////
void Scorer_Interface::run() {

	/*
	 * Create a variable to track errors.
	 * Throughout, the calculation proceeds as long as error = 0.
	 */
	int error = 0;

	/*
	 * Create two RNA strands, one for the predicted structure(s) and one for the accepted structure.
	 * For both strands, specify type = 1 (ct file).
	 * The type of nucleic acid does not matter, so it isn't specified.
	 */

	// Show message that intialization has begun.
	cout << "Initializing predicted and accepted structures..." << flush;

	// Create the first RNA strand.
	RNA* strand = new RNA( predicted.c_str(), FILE_CT, DT_RNA, true, true);
	ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );
	error = checker->isErrorStatus();

	 // Create the second RNA strand.
	RNA* strand2 = new RNA( accepted.c_str(), FILE_CT, DT_RNA, true, true);
	ErrorChecker<RNA>* checker2 = new ErrorChecker<RNA>( strand2 );
	if( error != 0 ) { error = checker2->isErrorStatus(); }

	// Show message that initialization is finished.
	if( error == 0 ) { cout << "done." << endl; }

	/*
	 * If no error occurred, check the accepted structure to make sure only one structure is present.
	 */

	if( error == 0 ) {

		// Get the number of structures in the accepted strand.
		int structures2 = strand2->GetStructureNumber();

		// If the number of structures in either strand is incorrect, show the appropriate error message.
		string numError = "Accepted ct must contain only one structure.";
		if( structures2 != 1 ) {
			cerr << numError << endl;
			error = 1;
		}
	}

	/*
	 * If no error occurred, check the length of the structures to make sure they are equal.
	 */

	if( error == 0 ) {

		// Get the sequence lengths.
		int length1 = strand->GetSequenceLength();
		int length2 = strand2->GetSequenceLength();

		// If the lengths aren't equal, show an error message.
		if( length1 != length2 ) {
			cerr << "Predicted and accepted structures are not the same length." << endl;
			error = 1;
		}
	}

	/*
	 * If no error occurred, calculate sensitivity and PPV, and write the data to a file.
	 */

	if( error == 0 ) {

		// Get the predicted and accepted structures that back the RNA strand.
		structure* predictedBack = strand->GetStructure();
		structure* acceptedBack = strand2->GetStructure();

		// Open the output stream to the written file.
		ofstream out( output.c_str() );

		// If no specific structure was specified, set the beginning and end of the structure range to all structures.
		// Otherwise, restrict the structure range to the specific structure.
		int start = ( number == -1 ) ? 1 : number;
		int end = ( number == -1 ) ? strand->GetStructureNumber() : number;

		// For each structure in the structure range, compare it with the accepted structure for sensitivity and PPV.
		for( int i = start; i <= end; i++ ) {

			// Print message saying sensitivity and PPV are being calculated.
			cout << "Calculating sensitivity and PPV for structure " << i << "..." << flush;

			// Calculate sensitivity.
			int pairs1 = 0;
			int score1 = 0;
			scorer( acceptedBack, predictedBack, &score1, &pairs1, i, exact );
			double percent1 = ( ( (double)score1 ) / ( (double)pairs1 ) ) * 100;

			// Calculate PPV.
			int pairs2 = 0;
			int score2 = 0;
			scorerppv( acceptedBack, predictedBack, &score2, &pairs2, i, exact );
			double percent2 = ( ( (double)score2 ) / ( (double)pairs2 ) ) * 100;

			// Create the sensitivity and PPV string streams.
			stringstream sensitivity( stringstream::in | stringstream::out );
			stringstream ppv( stringstream::in | stringstream::out );

			// Fill the string streams that hold sensitivity and PPV data.
			sensitivity << "Sensitivity: " << score1 << " / " << pairs1 << " = " << fixed << setprecision( 2 ) << percent1 << "%";
			ppv << "PPV:         " << score2 << " / " << pairs2 << " = " << fixed << setprecision( 2 ) << percent2 << "%";

			// Write the scoring data to the file.
			out << "Accepted Structure:  " << accepted << endl
			    << "Predicted Structure: " << predicted << endl
			    << "Score of Predicted Structure " << i << ":" << endl
			    << sensitivity.str() << endl
			    << ppv.str()
			    << endl << endl;

			// Print message saying calculation of a particular sensitivity and PPV are done.
			cout << "done." << endl;
		}

		// Close the output stream.
		out.close();
	}

	/*
	 * Output the sensitivity and PPV data to the screen, if the user requested.
	 * Only do this if no errors previously occurred.
	 */

	if( ( error == 0 ) && print ) {

		// Print out a spacing line on standard output.
		cout << endl;

		// Open the output file, pipe each line to standard output, then close the output file.
		ifstream in( output.c_str() );
		string line;
		while( !in.fail() ) {
			getline( in, line );
			if( line != "" ) { cout << line << endl; }
		}
		in.close();

		// Print out a line on standard output giving the name of the piped output file.
		cout << endl << "Saved in output file: " << output << endl << endl;
	}

	/*
	 * Clean up the strands and error checkers.
	 */

	// Clean up strand and error checker 1.
	delete strand;
	delete checker;

	// Clean up strand and error checker 2.
	delete strand2;
	delete checker2;

	/*
	 * Print out a confirmation of the run finishing.
	*/

	// Print confirmation of run finishing.
	if( error == 0 ) { cout << calcType << " complete." << endl; }
	else { cerr << calcType << " complete with errors." << endl; }
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	Scorer_Interface* runner = new Scorer_Interface();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}

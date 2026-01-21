/*
 * A program that predicts structures from multiple sequences using progressive iterations of Dynalign.
 *
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#include "Multilign_Interface.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
Multilign_Interface::Multilign_Interface() {

	// Initialize the calculation type description.
	calcType = "Multilign";

	// Initialize default number of sequences and number of processors to 1.
	seqNumber = 1;
	processors = 1;

	// Initialize the alignment window.
	alignmentWindow = 1;

	// Initialize the base pair window.
	basepairWindow = 2;

	// Initialize the gap penalty.
	gap = 0.4;

	// Initialize the index sequence.
	indexSeq = 1;

	// Initialize whether inserts are allowed to true.
	inserts = true;

	// Initialize the SHAPE intercept.
	intercept = -0.6;

	// Initialize the nucleic acid type to false.
	isRNA = true;

	// Intialize the number of iterations Multilign goes through.
	iterations = 2;

	// Initialize keeping of intermediate files to false.
	keepIntermediate = false;

	// Initialize whether local folding is done to false.
	local = false;

	// Initialize maximum DSV change.
	maxDsvChange = 1;

	// Initialize the maximum number of pairs allowed.
	// -1 specifies the average length of all sequences, calculated later.
	maxPairs = -1;

	// Initialize the default M separation parameter.
	maxSeparation = -99;

	// Initialize the maximum number of structures.
	maxStructures = 750;

	// Initialize the maximum percent energy difference.
	percent = 20;

	// Initialize the maximum percent difference in folding free energy change from single sequence folding for pairs that will be allowed in a subsequent calculation.
	percentSingle = 30;

	// Initialize whether sequences should be randomized to false.
	random = false;

	// Initialize the SHAPE slope.
	slope = 1.8;

	// Initialize the default temperature.
	temperature = 310.15;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool Multilign_Interface::parse( int argc, char** argv ) {

	// Determine the proper executable name, depending on if the executable is serial or SMP.
#ifndef COMPILE_SMP
	string type = "multilign";
#else
	string type = "multilign-smp";
#endif

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( type );
	parser->addParameterDescription( "configuration file", "The name of a file containing configuration data." );

	// Tell the parser that a special usage message is needed for this interface.
	parser->setSpecializedUsage();

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// If the specialized usage has been asked for, print it out.
	// An error is set here only for the purpose of preventing further parsing and calculation; it's not a real error.
	if( parser->isSpecializedUsage() ) {
		usage( parser );
		parser->setError();
		return false;
	}

	// Open the config file.
	// If the file isn't valid, delete the parser and return false.
	ConfigFile file( parser->getParameter( 1 ) );
	if( file.isValid() == false ) {
		delete parser;
		return false;
	}

	// Check if the config file contains the Alignment flag, which is always necessary for reading data.
	// If it exists, read it in.
	// If the flag isn't present, set an error, delete the parser, and return false.
	bool isReadable = file.contains( "Alignment" );
	if( isReadable ) { alignmentFile = file.getOption<string>( "Alignment" ); }
	else {
		parser->setErrorSpecialized( "Configuration file must contain the Alignment flag." );
		delete parser;
		return false;
	}

	// If SMP calculations are being done, make sure the processors flag has been specified.
	// If it has, set the number of processors if possible.
	// If it hasn't, show an error, delete the parser, and return false.
#ifdef COMPILE_SMP
	if( file.contains( "Processors" ) ) {
		processors = file.getOption<int>( "Processors" );
		if( processors < 0 ) {
			parser->setError( "number of processors" );
			delete parser;
			return false;
		}
	} else {
		parser->setErrorSpecialized( "Configuration file must contain the Processors flag." );
		delete parser;
		return false;
	}
#endif

	// Initialize the files vector.
	// Index 0: Sequence files
	// Index 1: CT files
	// Index 2: Optional folding constraint files
	// Index 3: Optional SHAPE constraint files
	for( int i = 1; i <= 4; i++ ) {
		vector<string> row;
		files.push_back( row );
	}

	// Check whether sequence files can be read singly.
	// If they can, attempt to read them into the 2D file vector.
	// An error occurs doing this, set an error, delete the parser, and return false.
	bool isReadableSingly =
		file.contains( "Seq1" ) &&
		file.contains( "CT1" ) &&
		file.contains( "SequenceNumber" );

	if( isReadableSingly ) {
		// Set the sequence number.
		// If the sequence number is less than or equal to 0, set an error, delete the parser, and return false.
		seqNumber = file.getOption<int>( "SequenceNumber" );
		if( seqNumber <= 0 ) {
			parser->setErrorSpecialized( "Sequence number must be greater than 0." );
			delete parser;
			return false;
		}

		// Read in all required files.
		// For each possible sequence number, read in the corresponding data file.
		for( int i = 1; i <= seqNumber; i++ ) {

			// Get the next sequence number as a string.
			string seqString = sfmt("Seq%i", i);
			// If an error occurs, set an error, delete the parser, and return false.
			if( file.contains( seqString ) ) { files[0].push_back( file.getOption<string>( seqString ) ); }
			else {
				parser->setErrorSpecialized( "The number of sequence files must be equal to SequenceNumber." );
				delete parser;
				return false;
			}

			// Read in the next CT file.
			// If an error occurs, set an error, delete the parser, and return false.
			string ctString = sfmt("CT%i", i);
			if( file.contains( ctString ) ) { files[1].push_back( file.getOption<string>( ctString ) ); }
			else {
				parser->setErrorSpecialized( "The number of CT files must be equal to SequenceNumber." );
				delete parser;
				return false;
			}
		}
	}
	// Check whether sequence files can be read in groups.
	// If they can, attempt to read them into the 2D file vector.
	// If an error occurs doing this, set an error, delete the parser, and return false.
	else if( file.contains( "InSeq" ) && file.contains( "OutCT" ) ) {
		// Get the sequences.
		string seqData = file.getOption<string>( "InSeq" );
		unsigned int seqLast = seqData.length() - 1;
		if( seqData[0] == '{' && seqData[seqLast] == '}' ) {
			seqData = seqData.erase( 0, 1 );
			seqData = seqData.erase( seqLast - 1, 1 );
			seqLast = seqData.length() - 1;
			if( seqData[seqLast] == ';' ) { seqData = seqData.erase( seqLast, 1 ); }
			stringstream seqStr( seqData );
			string seqFile;
			while( seqStr.good() ) {
				getline( seqStr, seqFile, ';' );
				if( seqFile != "" ) { files[0].push_back( seqFile ); }
			}
		} else {
			if( seqData[0] != '{' ) { parser->setErrorSpecialized( "Sequence group has no start bracket." ); }
			else { parser->setErrorSpecialized( "Sequence group has no end bracket." ); }
			delete parser;
			return false;
		}

		// Get the CT files.
		string ctData = file.getOption<string>( "OutCT" );
		unsigned int ctLast = ctData.length() - 1;
		if( ctData[0] == '{' && ctData[ctLast] == '}' ) {
			ctData = ctData.erase( 0, 1 );
			ctData = ctData.erase( ctLast - 1, 1 );
			ctLast = ctData.length() - 1;
			if( ctData[ctLast] == ';' ) { ctData = ctData.erase( ctLast, 1 ); }
			stringstream ctStr( ctData );
			string ctFile;
			while( ctStr.good() ) {
				getline( ctStr, ctFile, ';' );
				if( ctFile != "" ) { files[1].push_back( ctFile ); }
			}
		} else {
			if( ctData[0] != '{' ) { parser->setErrorSpecialized( "Sequence group has no start bracket." ); }
			else { parser->setErrorSpecialized( "Sequence group has no end bracket." ); }
			delete parser;
			return false;
		}

		// If the number of sequence files equals the number of CT files, set that number as the number of sequences.
		// If not, set an error, delete the parser, and return.
		if( files[0].size() == files[1].size() ) { seqNumber = files[0].size(); }
		else {
			parser->setErrorSpecialized( "Number of sequence files does not equal number of CT files." );
			delete parser;
			return false;
		}
	}

	// If the configuration file isn't set up properly to read sequence files singly or in groups, set an error, delete the parser, and return false.
	else { 
		// not readable singly or in groups.
		parser->setErrorSpecialized( "File groups are not specified properly. Please run \"" + type + " -h\" for help." );
		delete parser;
		return false;
	}

	// Check to see if optional constraint files or SHAPE files were specified.
	stringstream optionalStream( stringstream::in | stringstream::out );
	string constraintString, shapeString;
	hasSHAPE = false;
	for( int i = 1; i <= seqNumber; i++ ) {

		// If constraint file i is specified, get the name.
		optionalStream << "Constraint" << i;
		constraintString = optionalStream.str();
		if( file.contains( constraintString ) ) { files[2].push_back( file.getOption<string>( constraintString ) ); }
		else { files[2].push_back( "" ); }
		optionalStream.str( "" );

		// If SHAPE file i is specified, get the name.
		optionalStream << "SHAPE" << i;
		shapeString = optionalStream.str();
		if( file.contains( shapeString ) ) {
			files[3].push_back( file.getOption<string>( shapeString ) );
			hasSHAPE = true;
		} else { files[3].push_back( "" ); }
		optionalStream.str( "" );
	}

	// Get the nucleic acid type.
	if( !parser->isError() ) {
		if( file.contains( "DNA" ) ) { isRNA = !( file.getOption<bool>( "DNA" ) ); }
	}

	// Get the gap penalty.
	if( !parser->isError() ) {
		if( file.contains( "Gap" ) ) {
			gap = file.getOption<double>( "Gap" );
			if( gap <= 0.0 ) { parser->setError( "gap penalty" ); }
		}
	}

	// Get the index sequence.
	if( !parser->isError() ) {
		if( file.contains( "IndexSeq" ) ) {
			indexSeq = file.getOption<int>( "IndexSeq" );
			bool inRange = ( indexSeq >= 1 ) && ( indexSeq <= seqNumber );
			if( inRange == false ) { parser->setError( "index sequence number" ); }
		}
	}

	// Get whether base pairs should be inserted.
	if( !parser->isError() ) {
		if( file.contains( "Insert" ) ) { inserts = file.getOption<bool>( "Insert" ); }
	}

	// Get the number of iterations.
	if( !parser->isError() ) {
		if( file.contains( "Iterations" ) ) {
			iterations = file.getOption<int>( "Iterations" );
			if( iterations <= 0 ) { parser->setError( "number of iterations" ); }
		}
	}

	// Get whether intermediate files should be kept.
	if( !parser->isError() ) {
		if( file.contains( "KeepIntermediate" ) ) { keepIntermediate = file.getOption<bool>( "KeepIntermediate" ); }
	}

	// Get whether local folding should be done.
	if( !parser->isError() ) {
		if( file.contains( "Local" ) ) { local = file.getOption<bool>( "Local" ); }
	}

	// Get the maxdsvchange.
	if( !parser->isError() ) {
		if( file.contains( "MaxDsvChange" ) ) {
			maxDsvChange = file.getOption<double>( "MaxDsvChange" );
			bool inRange =
				( ( maxDsvChange > 0.0 && maxDsvChange < 99.0 ) ) ||
				( maxDsvChange == -1.0 );
			if( inRange == false ) { parser->setError( "maxdsvchange" ); }
		}
	}

	// Get the max pairs value.
	if( !parser->isError() ) {
		if( file.contains( "MaxPairs" ) ) {
			maxPairs = file.getOption<int>( "MaxPairs" );
			if( maxPairs <= 0 ) { parser->setError( "max pairs value" ); }
		}
	}

	// Get the maximum percent energy difference, if applicable.
	if( !parser->isError() ) {
		if( file.contains( "MaxPercent" ) ) {
			percent = file.getOption<int>( "MaxPercent" );
			if( percent <= 0 ) { parser->setError( "percent energy difference" ); }
		}
	}

	// Get the maximum percent energy difference in folding free energy change from single sequence folding.
	if( !parser->isError() ) {
		if( file.contains( "MaxPercentSingle" ) ) {
			percentSingle = file.getOption<int>( "MaxPercentSingle" );
			if( percentSingle <= 0 ) { parser->setError( "percent energy difference in folding free energy change" ); }
		}
	}

	// Get the maximum number of structures.
	if( !parser->isError() ) {
		if( file.contains( "MaxStructures" ) ) {
			maxStructures = file.getOption<int>( "MaxStructures" );
			if( maxStructures <= 0 ) { parser->setError( "maximum number of structures" ); }
		}
	}

	// Get the randomization flag.
	if( !parser->isError() ) {
		if( file.contains( "Random" ) ) { random = file.getOption<bool>( "Random" ); }
	}

	// Get the traditional M separation parameter.
	if( !parser->isError() ) {
		if( file.contains( "Separation" ) ) { maxSeparation = file.getOption<double>( "Separation" ); }
	}

	// Get the SHAPE intercept.
	if( !parser->isError() ) {
		if( hasSHAPE ) {
			if( file.contains( "SHAPEintercept" ) ) { intercept = file.getOption<double>( "SHAPEintercept" ); }
		}
	}

	// Get the SHAPE slope.
	if( !parser->isError() ) {
		if( hasSHAPE ) {
			if( file.contains( "SHAPEslope" ) ) { slope = file.getOption<double>( "SHAPEslope" ); }
		}
	}

	// Get the temperature.
	if( !parser->isError() ) {
		if( file.contains( "Temperature" ) ) {
			temperature = file.getOption<double>( "Temperature" );
			if( temperature < 0.0 ) { parser->setError( "temperature" ); }
		}
	}

	// Get the alignment window size.
	if( !parser->isError() ) {
		if( file.contains( "WindowAlign" ) ) {
			alignmentWindow = file.getOption<int>( "WindowAlign" );
			if( alignmentWindow < 0 ) { parser->setError( "alignment window" ); }
		}
	}

	// Get the base pair window size.
	if( !parser->isError() ) {
		if( file.contains( "WindowBP" ) ) {
			basepairWindow = file.getOption<int>( "WindowBP" );
			if( basepairWindow < 0 ) { parser->setError( "base pair window" ); }
		}
	}

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void Multilign_Interface::run() {

	/*
	 * Rearrange the files vectors so they're organized by sequence number, not by file type.
	 * This is done behind the scenes; no progress message is necessary here.
	 */

	vector< vector<string> > matrix;
	vector<string> filesSeq = files[0];
	vector<string> filesCT = files[1];
	vector<string> filesConstraint = files[2];
	vector<string> filesSHAPE = files[3];
	for( int i = 1; i <= seqNumber; i++ ) {
		vector<string> row;
		row.push_back( filesSeq[i-1] );
		row.push_back( filesCT[i-1] );
		row.push_back( filesConstraint[i-1] );
		row.push_back( filesSHAPE[i-1] );
		matrix.push_back( row );
	}

	// Create a variable that handles errors.
	int error = 0;

	/*
	 * Use the constructor for Multilign_object that specifies a matrix of file names and a nucleic acid type.
	 * This allows for many varied sequence files to be used as input.
	 *
	 * After construction of the data structure, create the error checker which monitors for errors.  
	 * Throughout, the error status of the calculation is checked with a variant of the isErrorStatus method, which returns 0 if no error occurred.
	 * The calculation proceeds as long as error = 0.
	 */
	cout << "Initializing nucleic acids..." << flush;
	Multilign_object* object = new Multilign_object( matrix, isRNA );
	ErrorChecker<Multilign_object>* checker = new ErrorChecker<Multilign_object>( object );
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
		object->SetTemperature( temperature );
		error = checker->isErrorStatus();

		// If no error occurred, print a message saying that temperature is set.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * Randomize the sequences using the Randomize method.
	 */
	if( ( error == 0 ) && ( random == true ) ) {

		// Show a message saying that randomization has started.
		cout << "Randomizing sequences..." << flush;

		// Randomize sequences.
		// No error can be thrown from this method.
		object->Randomize();

		// Print a message saying that randomization is done.
		cout << "done." << endl;
	}

	/*
	 * Set the index sequence using the SetIndexSeq method.
	 * Only set the index sequence if the given index sequence is not 1.
	 * After setting the index sequence, use the error checker's isErrorStatus method to check for errors.
	 */
	if( ( error == 0 ) && indexSeq != 1 ) {

		// Show a message saying that index sequence setting has started.
		cout << "Setting index sequence..." << flush;

		// Set the index sequence and check for errors.
		int seqError = object->SetIndexSeq( indexSeq );
		error = checker->isErrorStatus( seqError );

		// If no error occurred, print a message saying that sequence setting is done.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * Set the iterations using the SetIterations method.
	 */
	if( ( error == 0 ) && ( iterations != 2 ) ) {

		// Show a message saying that iteration setting has started.
		cout << "Setting iterations..." << flush;

		// Set iterations.
		int itError = object->SetIterations( iterations );
		error = checker->isErrorStatus( itError );

		// If no error occurred, print a message saying that iteration setting
		// is done.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * Set a maxdsvchange using the SetMaxDsv method.
	 */
	if( ( error == 0 ) && ( maxDsvChange != 1 ) ) {

		// Show a message saying that max DSV change setting has started.
		cout << "Setting max DSV change..." << flush;

		// Set max DSV change.
		int dsvError = object->SetMaxDsv( maxDsvChange );
		error = checker->isErrorStatus( dsvError );

		// If no error occurred, print a message saying that max dsv change setting is done.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * Set max pairs using the SetMaxPairs method.
	 */
	if( ( error == 0 ) && ( maxPairs != -1 ) ) {

		// Show a message saying that max pairs setting has started.
		cout << "Setting max pairs..." << flush;

		// Set max pairs.
		int pairsError = object->SetMaxPairs( maxPairs );
		error = checker->isErrorStatus( pairsError );

		// If no error occurred, print a message saying that max pairs setting is done.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/**
	 * Set SHAPE data using the SetSHAPEIntercept and SetSHAPESlope method.
	 * Only do this if SHAPE data has been read.
	 */
	if( ( error == 0 ) && ( hasSHAPE == true ) ) {

		// Show a message saying that SHAPE parameter setting has started.
		cout << "Setting SHAPE parameters..." << flush;

		// Set SHAPE parameters.
		object->SetSHAPEIntercept( intercept );
		object->SetSHAPESlope( slope );

		// If no error occurred, print a message saying that SHAPE parameter setting is done.
		cout << "done." << endl;
	}

	/*
	 * Run the Multilign algorithm using the ProgressiveMultilign method.
	 * During calculation, monitor progress using the TProgressDialog class and the Start/StopProgress methods of the RNA class.
	 * Neither of these methods require any error checking.
	 * After the main calculation is complete, use the error checker's isErrorStatus method to check for errors.
	 */
	if( error == 0 ) {

		// Show a message saying that the main calculation has started.
		cout << "Folding nucleic acids..." << endl;

		// Create the progress monitor.
		TProgressDialog* progress = new TProgressDialog();
		object->SetProgress( progress );
		// Run the Multilign algorithm and check for errors.
		int mainCalcError = object->ProgressiveMultilign(
			processors, true, true, maxStructures, basepairWindow,
			alignmentWindow, percent, maxSeparation, gap, inserts,
			percentSingle, local );
		error = checker->isErrorStatus( mainCalcError );

		// Delete the progress monitor.
		object->StopProgress();
		delete progress;

		// If no error occurred, print a message saying that the main calculation is done.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * Write the alignment file using the WriteAlignment method.
	 * After writing is complete, use the error checker's isErrorStatus method to check for errors.
	 */
	if( error == 0 ) {

		// Show a message saying that the alignment file is being written.
		cout << "Writing multiple alignment file..." << flush;

		// Write the alignment file and check for errors.
		int writeError = object->WriteAlignment( alignmentFile );
		error = checker->isErrorStatus( writeError );

		// If no errors occurred, print a message saying that alignment file writing is done.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * If requested by the user, clean up the intermediate files created by Multilign using the CleanupIntermediateFiles method.
	 * After files are cleaned up, use the error checker's isErrorStatus method to check for errors.
	 */
	if( keepIntermediate == false ) {
		int cleanupError = object->CleanupIntermediateFiles();
		checker->isErrorStatus( cleanupError );
		if (error == 0) error = cleanupError; // do not clobber existing error.
	}

	// Delete the error checker and data structure.
	delete checker;
	delete object;

	// Print confirmation of run finishing.
	if( error == 0 ) { cout << calcType << " complete." << endl; }
	else { cerr << calcType << " complete with errors." << endl; }
}

///////////////////////////////////////////////////////////////////////////////
// Print out a special usage message for this interface.
///////////////////////////////////////////////////////////////////////////////
void Multilign_Interface::usage( ParseCommandLine* parser ) {

	cout << "===================================" << endl
	     << "==== Configuration File Format ====" << endl
	     << "===================================" << endl
	     << "Note that configuration options are not case-sensitive." << endl
	     << "Any unrecognized options are ignored." << endl << endl
	     << "Configuration file line format:" << endl
	     << "<Option> = <Value>" << endl << endl;

	cout << "Required input when specifying file groups" << endl
	     << "------------------------------------------" << endl
	     << "The sequence and CT file groups must each specify the same number of files." << endl << endl;
	cout << "InSeq" << endl;
	parser->printDescription( "Flag that can be used to specify a group of sequence files, from Seq1 to Seq<SequenceNumber>. Only one sequence file group can be specified without overwriting files. Group format: {seq1File;seq2File;seq3File;}" );
	cout << "OutCT" << endl;
	parser->printDescription( "Flag that can be used to specify a group of CT files, from CT1 to CT<SequenceNumber>. Only one CT file group can be specified without overwriting files. Group format: {ct1File;ct2File;ct3File;}" );
	cout << "Alignment" << endl;
	parser->printDescription( "The name of the file used to hold a multiple sequence alignment." );

	cout << "Required input when specifying files singly" << endl
	     << "-------------------------------------------" << endl
	     << "Every specified sequence file must have a corresponding CT file specified." << endl << endl;
	cout << "SequenceNumber" << endl;
	parser->printDescription( "The number of sequences given as input." );
	cout << "Seq1 ... Seq<SequenceNumber>" << endl;
	parser->printDescription( "Names of sequence files used as input, from 1 to SequenceNumber." );
	cout << "CT1 ... CT<SequenceNumber>" << endl;
	parser->printDescription( "Names of CT files written as output, from 1 to SequenceNumber." );
	cout << "Alignment" << endl;
	parser->printDescription( "The name of the file used to hold a multiple sequence alignment." );

	cout << "General options" << endl
	     << "---------------" << endl;
	cout << "DNA" << endl;
	parser->printDescription( "Specifies whether Multilign is run with RNA or DNA thermodynamics. Use 0 to specify RNA thermodynamics, and 1 to specify DNA thermodynamics. Default is 0, which specifies RNA thermodynamics." );
	cout << "Constraint1 ... Constraint<SequenceNumber>" << endl;
	parser->printDescription( "Names of folding constraint files. The number at the end of the flag (1 to SequenceNumber) identifies which sequence the constraints will be applied to." );
#ifdef COMPILE_SMP
	cout << "Processors" << endl;
	parser->printDescription( "The number of processors on which the calculation runs. Default is 1." );
#endif
	cout << "SHAPE1 ... SHAPE<SequenceNumber>" << endl;
	parser->printDescription( "Names of SHAPE constraint files. The number at the end of the flag (1 to SequenceNumber) identifies which sequence the constraints will be applied to." );
	cout << "SHAPEintercept" << endl;
	parser->printDescription( "The SHAPE intercept. This value is only used when at least one SHAPE constraint file is specified. Default is 1.8 kcal/mol." );
	cout << "SHAPEslope" << endl;
	parser->printDescription( "The SHAPE slope. This value is only used when at least one SHAPE constraint file is specified. Default is -0.6 kcal/mol." );
	cout << "Temperature" << endl;
	parser->printDescription( "The temperature at which calculations are run, in Kelvin. Default is 310.15 K, which is 37 degrees C." );

	cout << "Multilign options" << endl
	     << "-----------------" << endl;
	cout << "IndexSeq" << endl;
	parser->printDescription( "The sequence chosen to be the index (main) sequence. Default is 1." );
	cout << "Iterations" << endl;
	parser->printDescription( "The number of iterations the calculation goes through. Default is 2 iterations." );
	cout << "KeepIntermediate" << endl;
	parser->printDescription( "Keep intermediate pairwise save files and alignments generated during calculations. Use 0 to specify that files should be deleted, and 1 to specify they should be kept. Default is 0, which specifies deletion of intermediate files." );
	cout << "MaxDsvChange" << endl;
	parser->printDescription( "The maximum interval above lowest energy for a base pair to be allowed. Default is 1 (percent)." );
	cout << "MaxPairs" << endl;
	parser->printDescription( "The maximum number of base pairs allowed in a structure. Default is -1, which specifies the average length of all sequences." );
	cout << "Random" << endl;
	parser->printDescription( "Specifies whether to randomize sequence order. Use 0 to specify no randomization, and 1 to specify randomization. Default is 0, which specifies no randomization." );

	cout << "Dynalign options" << endl
	     << "----------------" << endl;
	cout << "Gap" << endl;
	parser->printDescription( "The per nucleotide insert penalty for alignments. Default is 0.4." );
	cout << "Insert" << endl;
	parser->printDescription( "Specifies whether to allow single base pair inserts. Use 0 to specify no inserts, and 1 to specify inserts. Default is 1, which specifies inserts." );
	cout << "Local" << endl;
	parser->printDescription( "Specifies whether to do local or global folding. Use 0 to specify global folding, and 1 to specify local folding. Default is 0, which specifies global folding." );
	cout << "MaxPercent" << endl;
	parser->printDescription( "The maximum percent energy difference. Default is 20 (percent)." );
	cout << "MaxPercentSingle" << endl;
	parser->printDescription( "The maximum percent difference in folding free energy change from single sequence folding for pairs that will be allowed in a subsequent calculation. This pre-screens allowed pairs. Default is 30 (percent)." );
	cout << "Separation" << endl;
	parser->printDescription( "The traditional M parameter. Default is -99." );
	cout << "WindowAlign" << endl;
	parser->printDescription( "The alignment window size. Default is 1." );
	cout << "WindowBP" << endl;
	parser->printDescription( "The base pair window size. Default is 2." );
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	Multilign_Interface* runner = new Multilign_Interface();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}


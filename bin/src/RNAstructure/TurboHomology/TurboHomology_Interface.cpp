/*
 * A program that predicts structures using the TurboHomology algorithm.
 *
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#include "TurboHomology_Interface.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
TurboHomology_Interface::TurboHomology_Interface() {

	// Initialize the calculation type description.
	calcType = "TurboHomology";

	processors = 1;

	// Initialize the default mode.
	mode = "MEA";

	// Initialize general options, which are applied independent of mode.
		// Initialize the maximum pairing distance to 0, which indicates no limit.
		distance = 0;

		// Initialize the SHAPE intercept.
		intercept = -0.6;

		// Initialize the SHAPE slope.
		slope = 1.8;

		// Initialize the default temperature.
		temperature = 310.15;

		// Initialize the TurboHomology gamma.
        turboGamma = 0.2;

		// Initialize the TurboHomology iterations.
		turboIterations = 1;

		// Initialize the default alignment output filename.
        OutAln = "Output.aln";

        // Initialize the default alignment output format.
        AlnFormat = "Clustal";

        // Initialize the default alignment output max column number.
        ColumnNumber = 60;

		hasSHAPE = false;

        // Rsample parameters.
	    // The Cparam for SHAPE constraints.
	    Cparam = 0.5;
	    
	    // The Offset for SHAPE constraints.
	    Offset = 1.1;

	    // number of samples for stochastic sampling
	    rsample_numsamples = 10000;

        // Names of files for Rsample: paired-end, paired-middle and unpaired files
	    peFile = "";
    	pmFile = "";
    	upFile = "";

    	// pfsFile = "";

    	is_RSample_mode = false;

    	rsample_seed = 0; // random seed for Rsample. 0 indicates that it is initialized with time(0)

	// Initialize maximum expected accuracy mode parameters.
		// Initialize the maximum number of structures.
		maxStructures = 1000;

		// Initialize the maximum expected accuracy mode gamma.
		meaGamma = 1;

		// Initialize the maximum percent energy difference.
		percent = 50;

		// Initialize the window size.
		windowSize = 5;

	// Initialize ProbKnot mode parameters.
		// Initialize the number of ProbKnot iterations.
		pkIterations = 1;

		// Initialize the minimum helix length for a pseudoknot.
		minHelixLength = 3;

	// Initialize Threshold mode parameters.
		// Initialize the default probability threshold.
		// The default of 0 means that structures are generated at multiple thresholds.
		threshold = 0;

		ExistingAln = "";

		parameters.assign (2, 0);
}

TurboHomology_Interface::~TurboHomology_Interface() {
	//if (rsdata != NULL) delete rsdata;
}

#define ERR_BASE 60000
#define ERR_NO_START_BRACKET (ERR_BASE+1)
#define ERR_NO_END_BRACKET (ERR_BASE+2)
// splits a parameter group into a vector of strings, using semicolon as a delimiter.
// The parameter group must be surrounded by braces { }
// If the string ends in a semicolon, it is removed, so an empty string is NOT returned as the final value in this case.
// returns 0 on success, ERR_NO_START_BRACKET if the start bracket is missing, or ERR_NO_END_BRACKET if the end bracket is missing.
int split_params(const string &input, vector<string> &results) {
	if (input.empty()||input[0] != '{')
		return ERR_NO_START_BRACKET;
	
	size_t end = input.size()-1;
	if (end<1||input[end] != '}')
		return ERR_NO_END_BRACKET;
	
	if (end>1&&input[end-1]==';') end--; // move end back 1 if the string ends in a semicolon
	stringstream ss( input.substr(1, end-1) );
	string param;
	while( ss.good() ) {
		getline( ss, param, ';' );
		if (ss.fail()) break; // end of reading. Note: (specifically check for fail() because eof() might be true if we reached the end, but we should still add the remaining string.
		trim(param); // removes whitespace from start and end of param.
		results.push_back(param);
	}
	return 0;
}

// returns true if there is an error or false otherwise. (i.e. false==success)
bool test_error(ParseCommandLine& parser, const char*const paramName, int errorCode) {
	switch(errorCode) {
		case 0: return false; //success
		case ERR_NO_START_BRACKET: parser.setErrorSpecialized(sfmt("'%s' group has no start bracket.",paramName)); break;
		case ERR_NO_END_BRACKET: parser.setErrorSpecialized(sfmt("'%s' group has no end bracket.",paramName)); break;
		default: parser.setErrorSpecialized(sfmt("Unknown error %d parsing parameter '%s'.",errorCode, paramName)); break;
	}
	return true; // an error occurred
}

bool verify_files(ParseCommandLine& parser, const char*const fileDescription, const vector<string> &files, bool allowEmpty = false) {
	for(unsigned int i = 0; i < files.size(); i++) {
		if (files[i].empty()) {
			if (allowEmpty) continue;
 			parser.setErrorSpecialized(sfmt("'%s' file %i was blank.",fileDescription, i+1));
			return false;
		} else {
			if (fileExists(files[i].c_str())) continue;
			parser.setErrorSpecialized(sfmt("'%s' file %i was not found: \"%s\".",fileDescription, i+1, files[i].c_str()));
			return false;
		}
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool TurboHomology_Interface::parse( int argc, char** argv ) {

	// Determine the proper executable name, depending on if the executable is serial or SMP.
#ifndef COMPILE_SMP
	string type = "TurboHomology";
#else
	string type = "TurboHomology-smp";
#endif

	// Create the command line parser and build in its required parameters.
	ParseCommandLine parser( type );
	parser.addParameterDescription( "configuration file", "The name of a file containing configuration data." );

	// Tell the parser that a special usage message is needed for this interface.
	parser.setSpecializedUsage();

	// Parse the command line into pieces.
	parser.parseLine( argc, argv );

	// If the specialized usage has been asked for, print it out.
	// An error is set here only for the purpose of preventing further parsing and calculation; it's not a real error.
	if( parser.isSpecializedUsage() ) {
		usage( parser );
		parser.setError();
		return false;
	}

	if (parser.isError()) return false;
	
	// Open the config file.
	// If the file isn't valid, delete the parser and return false.
	ConfigFile file( parser.getParameter( 1 ) );
        conf_file = parser.getParameter( 1 );
	
	if( !file.isValid() ) return false;

        if(file.contains( "Parameters") )
        {
            string ParData = file.getOption<string>( "Parameters" );
            unsigned int ParLast = ParData.length() - 1;
            if( ParData[0] == '{' && ParData[ParLast] == '}' ) {
				ParData = ParData.erase( 0, 1 );
				ParData = ParData.erase( ParLast - 1, 1 );
				ParLast = ParData.length() - 1;
				if( ParData[ParLast] == ';' ) { ParData = ParData.erase( ParLast, 1 ); }
				stringstream ParStr( ParData );
				string tempParameter;
			                int parameterIndex = 0;
				while( ParStr.good()) {
			                    getline( ParStr, tempParameter, ';' );
			                    if( tempParameter != "" ) { 
			                        if(parameterIndex > parameters.size())
			                        {
			                            // parser->setErrorSpecialized( "Too many parameters are provided." );
			                            // delete parser;
			                            return false;
			                        }
			                        parameters[parameterIndex] = atof(tempParameter.c_str()); 
			                        ++parameterIndex;
			                    }
				}
			            } else {
				// if( ParData[0] != '{' ) { parser->setErrorSpecialized( "Parameter group has no start bracket." ); }
				// else { parser->setErrorSpecialized( "Parameter group has no end bracket." ); }
				// delete parser;
				return false;
            }
        }

    for(int i = 0; i < 2; ++i)
    {
        stringstream numStream( stringstream::in | stringstream::out );
        numStream << i;
        string num = numStream.str();

        string parName = "Parameter" + num;
        if(file.contains( parName )) {
                parameters[i] = file.getOption<double>(parName);
        }
    }  

	// Check if the config file contains the Mode flag, which is always necessary for reading data.
	// If it exists, read it in.
	// If the flag isn't present, set an error, delete the parser, and return false.
	bool isReadable = file.contains( "Mode" );
	if( isReadable ) { file.getOptionByRef("Mode", mode); }
	else {
		parser.setErrorSpecialized( "Configuration file must contain the Mode flag." );
		return false;
	}

	// If SMP calculations are being done, make sure the processors flag has been specified.
	// If it has, set the number of processors if possible.
	// If it hasn't, show an error, delete the parser, and return false.
#ifdef COMPILE_SMP
	if( file.contains( "Processors" ) ) {
		file.getOptionByRef("Processors", processors);
		if( processors < 0 ) {
			parser.setError( "number of processors" );
			return false;
		}
	} else {
		parser.setErrorSpecialized( "Configuration file must contain the Processors flag." );
		return false;
	}
#endif

	// Read in the list of sequences. It can be specified either as a group (InSeq) or as individual sequences (Seq1, Seq2 etc)
	if (file.contains( "InSeq" )) {
		// Input sequences specified in groups
		if (test_error(parser,"Sequence", split_params(file.getOption<string>( "InSeq" ), sequenceFiles))) return false;
	} else {
		// Input sequences specified individually
		if (!file.contains("SequenceNumber")) {
			parser.setErrorSpecialized( "SequenceNumber is a required parameter if sequences are not listed in group form (in InSeq)" );
			return false;
		}
		// Set the sequence number. If the sequence number is less than or equal to 0, return false.
		int size = file.getOption<int>("SequenceNumber");
		if( size <= 0 ) {
			parser.setErrorSpecialized( "Invalid SequenceNumber parameter -- number of sequences must be greater than 0." );
			return false;
		}
		sequenceFiles.resize(size);
		for(int i=0; i<size; i++) {
			string param(sfmt("Seq%d",i+1));  // e.g. "Seq1" etc.
			if( file.contains( param ) ) 
				sequenceFiles[i] = file.getOption<string>(param);
			else {
				parser.setErrorSpecialized(sfmt("The number of sequence files specified must be equal to SequenceNumber. Missing: %s",param.c_str()));
				return false;
			}
		}
	}
	const unsigned int sequenceCount = sequenceFiles.size(); // from here on, the size of outputCtFiles, outputPfsFiles, etc must match the size of the sequenceFiles vector.

	// Read in the list of CT files. It can be specified either as a group (OutCT) or as individual sequences (CT1, CT2 etc)
	if (file.contains( "OutCT" )) {
		// Input sequences specified in groups
		if (test_error(parser,"CT", split_params(file.getOption<string>( "OutCT" ), outputCtFiles))) return false;
	} else {
		outputCtFiles.reserve(sequenceCount); // set expected size.
		int i=1; string param;
		// even though we know the expected count, keep looping to gather all CT<N> parameters to show the use an error if they specified more than required.
		while(file.contains((param=sfmt("CT%d",i++))))
			outputCtFiles.push_back(file.getOption<string>(param));
	}
	// if (outputCtFiles.size()!=sequenceCount){
	// 	parser.setErrorSpecialized(sfmt("The number of output CT file names must be equal to the number of input Sequences. Found %d but expected %d.",outputCtFiles.size(),sequenceCount));
	// 	return false;
	// }

	if (file.contains( "RefCT" )) {
		if (test_error(parser,"RefCT", split_params(file.getOption<string>( "RefCT" ), refCTFiles))) return false;
	} else {
		refCTFiles.reserve(sequenceCount-1);
		int i=1; string param;
		while(file.contains((param=sfmt("refCT%d",i++))))
			refCTFiles.push_back(file.getOption<string>(param));
	}
	if (refCTFiles.size()!=sequenceCount-1){
		parser.setErrorSpecialized(sfmt("The number of known CT file names must be one less to the number of input Sequences. Found %d but expected %d.",refCTFiles.size(),sequenceCount-1));
		return false;
	}

	// Read in the list of PFS files. It can be specified either as a group (SaveFiles) or as individual sequences (Save1, Save2 etc)
	if (file.contains( "SaveFiles" )) {
		// Input sequences specified in groups
		if (test_error(parser,"SaveFiles", split_params(file.getOption<string>( "SaveFiles" ), outputPfsFiles))) return false;
	} else {
		outputPfsFiles.reserve(sequenceCount); // set expected size.
		unsigned int i=0; string param;
		// even though we know the expected count, keep looping to gather parameters to show the use an error if they specified more than required.
		while(file.contains((param=sfmt("Save%d",++i)))||i<=sequenceCount)
			outputPfsFiles.push_back(file.getOption<string>(param)); // note that if the use did NOT specify a Save file at this position, then "" will be added.
	}
	if (outputPfsFiles.size()<sequenceCount) outputPfsFiles.resize(sequenceCount); // append additional "" to the end to ensure the size matches that of sequenceFiles.
	if (outputPfsFiles.size()>sequenceCount){
		parser.setErrorSpecialized(sfmt("The number of output Save file names cannot be larger than the number of input Sequences. Found %d but expected %d.",outputPfsFiles.size(),sequenceCount));
		return false;
	}

	// Read in the list of SHAPE files. It can be specified either as a group (SaveFiles) or as individual sequences (Save1, Save2 etc)
	if (file.contains( "SHAPEFiles" )) {
		// Input sequences specified in groups
		if (test_error(parser,"SHAPEFiles", split_params(file.getOption<string>( "SHAPEFiles" ), shapeFiles))) return false;
	} else {
		shapeFiles.reserve(sequenceCount); // set expected size.
		unsigned int i=0; string param;
		// even though we know the expected count, keep looping to gather parameters to show the use an error if they specified more than required.
		while(file.contains((param=sfmt("SHAPE%d",++i)))||i<=sequenceCount)
			shapeFiles.push_back(file.getOption<string>(param)); // note that if the use did NOT specify a Save file at this position, then "" will be added.
	}
	if (shapeFiles.size()<sequenceCount) shapeFiles.resize(sequenceCount); // append additional "" to the end to ensure the size matches that of sequenceFiles.
	if (shapeFiles.size()>sequenceCount){
		parser.setErrorSpecialized(sfmt("The number of SHAPE file names cannot be larger than the number of input Sequences. Found %d but expected %d.",shapeFiles.size(),sequenceCount));
		return false;
	}

	if (!verify_files(parser, "input Sequence", sequenceFiles, false)) return false;
	if (!verify_files(parser, "SHAPE data", shapeFiles, true)) return false;

	for(unsigned int i=0;i<shapeFiles.size();i++)
		if (!shapeFiles[i].empty())
			hasSHAPE = true;

	// Get the alignment file.
	file.getOptionByRef("OutAln", OutAln);

	file.getOptionByRef("ExistingAln", ExistingAln);

	if(file.contains("UseRsample")&&file.getOption<bool>( "UseRsample" )) {
		is_RSample_mode = true;
		cout << "In Rsample mode. " << endl;
		file.getOptionByRef( "Seed", rsample_seed );

		file.getOptionByRef( "Cparam", Cparam );
		//if( Cparam < 0.0 ) parser.setError( "Cparam" );

		file.getOptionByRef( "Offset", Offset );
		//if( Offset < 0.0 ) { parser.setError( "Offset" ); }

		// Get the number of samples for stochastic sampling
		file.getOptionByRef( "numsamples", rsample_numsamples );
		if( rsample_numsamples < 0 ) parser.setError( "numsamples" );

		if (!hasSHAPE) parser.setErrorSpecialized( "UseRsample is enabled, but no SHAPE files have been listed." );
	}

	if( parser.isError() ) return false;

	// Get the TurboHomology gamma.
	file.getOptionByRef( "Gamma", turboGamma );
	if( turboGamma < 0.0 ) { parser.setError( "TurboHomology gamma" ); return false; }

	// Get the TurboHomology iterations.
	file.getOptionByRef( "Iterations", turboIterations);
	if( turboIterations < 0 ) { parser.setError( "TurboHomology iterations" ); return false; }

	// Get the maximum pairing distance.
	file.getOptionByRef("MaximumPairingDistance", distance);
	if( distance < 0 ) { parser.setError( "maximum pairing distance" ); }

	// Get the SHAPE intercept.
	file.getOptionByRef("SHAPEintercept", intercept);

	// Get the SHAPE slope.
	file.getOptionByRef("SHAPEslope", slope);

	// Get the temperature.
	file.getOptionByRef("Temperature", temperature);
	if( temperature < 0.0 ) { parser.setError( "temperature" ); return false; }

	// Get the alignment output format:AlnFormat.
	file.getOptionByRef("AlnFormat", AlnFormat);
	if( AlnFormat != "Clustal" && AlnFormat != "Fasta" ) { parser.setError( "Alignment output format" ); }

	// Get the alignment output max column number:ColumnNumber.
	file.getOptionByRef("ColumnNumber", ColumnNumber);
	if( ColumnNumber < 0.0 ) { parser.setError( "Alignment output max column number" ); }
	
	if (parser.isError()) return false;


	// Set the MEA mode options, if applicable.
	if( mode == "MEA" ) {
		// Get the maximum percent energy difference.
		if( !parser.isError() ) {
			if( file.contains( "MaxPercent" ) ) {
				file.getOptionByRef("MaxPercent", percent);
				if( percent < 0.0 ) { parser.setError( "maximum percent energy difference" ); }
			}
		}

		// Get the maximum number of structures.
		if( !parser.isError() ) {
			if( file.contains( "MaxStructures" ) ) {
				file.getOptionByRef("MaxStructures", maxStructures);
				if( maxStructures < 0 ) { parser.setError( "maximum number of structures" ); }
			}
		}

                // Get the MEA gamma.
                if( !parser.isError() ) {
		  if( file.contains( "MeaGamma" ) ) { file.getOptionByRef("MeaGamma", meaGamma); }
                }

		// Get the window size.
		if( !parser.isError() ) {
			if( file.contains( "Window" ) ) {
				file.getOptionByRef("Window", windowSize);
				if( windowSize < 0 ) { parser.setError( "window size" ); }
			}
		}
	}

	// Set the ProbKnot mode options, if applicable.
	else if( mode == "ProbKnot" ) {

		// Get the number of ProbKnot iterations.
		if( !parser.isError() ) {
			if( file.contains( "PkIterations" ) ) {
				file.getOptionByRef("PkIterations", pkIterations);
				if( pkIterations < 0 ) { parser.setError( "ProbKnot iterations" ); }
			}
		}

		// Get the minimum helix length, if applicable.
		if( !parser.isError() ) {
			if( file.contains( "MinHelixLength" ) ) {
				file.getOptionByRef("MinHelixLength", minHelixLength);
				if( minHelixLength < 0 ) { parser.setError( "minimum helix length" ); }
			}
		}
	}

	// Set the Threshold mode options, if applicable.
	else if( mode == "Threshold" ) {

		// Get the threshold for pairs.
		if( !parser.isError() ) {
			if( file.contains( "Threshold" ) ) {
				file.getOptionByRef("Threshold", threshold);
				if( threshold < 0.0 ) { parser.setError( "pairing threshold" ); }
			}
		}
	}

	// If an invalid mode was specified, show an error message.
	else {
		parser.setErrorSpecialized( "Invalid mode given; mode must be 'MEA', 'ProbKnot', or 'TurboHomology'." );
	}

	return !parser.isError();
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
int TurboHomology_Interface::run() {

	// Create a variable that handles errors.
	int error = 0;
	const int seqNumber = sequenceFiles.size();
	RsampleData *rsdata = NULL;

	//DEBUG: // show all interface data
	//DEBUG: cout << "--------------------" << endl;
	//DEBUG: cout << "Sequences: " << "\n\t" << join(sequenceFiles,"\n\t") << endl;
	//DEBUG: cout << "CTs: " << "\n\t" << join(outputCtFiles,"\n\t") << endl;
	//DEBUG: cout << "Pfs: " << "\n\t" << join(outputPfsFiles,"\n\t") << endl;
	//DEBUG: cout << "SHAPEfiles: " << "\n\t" << join(shapeFiles,"\n\t") << endl;
	//DEBUG: cout << "OutAln: " << OutAln << endl;
	//DEBUG: cout << "--------------------" << endl;

	/*
	 * Use the constructor for TurboHomology_object that specifies vectors of file names.
	 * This allows for many varied sequence files to be used as input.
	 *
	 * After construction of the data structure, create the error checker which monitors for errors.
	 * Throughout, the error status of the calculation is checked with a variant of the isErrorStatus method, which returns 0 if no error occurred.
	 * The calculation proceeds as long as error = 0.
	 */

	cout << "Initializing nucleic acids..." << flush;
	TurboFold* object = new TurboFold(&sequenceFiles, &outputPfsFiles, OutAln, ExistingAln, parameters);
	ErrorChecker<TurboFold>* checker = new ErrorChecker<TurboFold>( object );
	error = checker->isErrorStatus();
	if( error == 0 ) { cout << "done." << endl; }

	if (is_RSample_mode) {
		cout << "Setting up Rsample..." << flush;
		//For now, do not support DMS and use default max for SHAPE data
		rsdata = new RsampleData(false,1000,upFile.c_str(), peFile.c_str(), pmFile.c_str());
		error = rsdata->ErrorCode;

		if (error != 0) {
			cerr << RsampleData::GetErrorMessage(error) << endl;
		} else {
			error = object->setupRsample(&shapeFiles, rsdata, rsample_numsamples, rsample_seed, Cparam, Offset);
			if (error != 0)
				cerr << "Error reading SHAPE data: " << object->GetErrorDetails() << endl;
		}

		if( error == 0 ) { cout << "done." << endl; }
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
		int tempError = object->SetTemperature( temperature );
		error = checker->isErrorStatus( tempError );

		// If no error occurred, print a message saying temperature is set.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * Set the maximum pairing distance, if required.
	 */
	if( ( error == 0 ) && ( distance > 0 ) ) {

		// Show a message saying that maximum pairing distance is being set.
		cout << "Setting maximum pairing distance..." << flush;

		// Set the maximum pairing distance and check for errors.
		int distError = object->SetMaxPairingDistance( distance );
		error = checker->isErrorStatus( distError );

		// If no error occurred, print a message saying that maximum pairing
		// distance is set.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * Read SHAPE data where necessary for specific sequence files.
	 * Since SHAPE data must be read individually for each sequence, use the ReadSHAPE method.
	 */
	if( error == 0 && hasSHAPE && !is_RSample_mode) {
		cout << "Reading SHAPE data..." << flush;
		// For each sequence, read in SHAPE data, if applicable.
		for( int i = 1; i <= seqNumber; i++ ) {
			if(!shapeFiles[i-1].empty()) {
				// Show a message saying SHAPE data is being read. (but only if this is the first)
				int shapeError = object->ReadSHAPE( i, shapeFiles[i-1].c_str(), slope, intercept );
				if( ( error = checker->isErrorStatus( shapeError ) ) !=0 )
					break;
			}
		}
		// Show a message saying SHAPE data reading is done.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * Run the TurboHomology algorithm using the fold method.
	 * During calculation, monitor progress using the TProgressDialog class and the Start/StopProgress methods of the RNA class.
	 * Neither of these methods require any error checking.
	 * After the main calculation is complete, use the error checker's isErrorStatus method to check for errors.
	 */
	if( error == 0 ) {

		// Show a message saying that the main calculation has started.
        cout << "Folding nucleic acids..." << endl;

		// Create the progress monitor.
		TProgressDialog* progress = new TProgressDialog();
		object->SetProgress( *progress );
		// Get known structures of sequences from CT files.
		object->getKnownStr( &refCTFiles );

		// Run the TurboHomology algorithm and check for errors.
//		int mainCalcError = object->fold( turboGamma, turboIterations, processors );
		int mainCalcError = object->fold( turboGamma, turboIterations, processors, AlnFormat, ColumnNumber, ExistingAln);
		error = checker->isErrorStatus( mainCalcError );
		// Delete the progress monitor.
		object->StopProgress();
		delete progress;

		// If no error occurred, print message that main calculation is done.
		if( error == 0 ) { cout << "done." << endl; }
	}
	
	if (rsdata != NULL) delete rsdata;

	/*
	 * Resolve the structures generated by TurboHomology using specific methods, depending on the mode selected for TurboHomology.
	 * In MEA mode, use the MaximizeExpectedAccuracy method.
	 * In ProbKnot mode, use the ProbKnot method.
	 * In Threshold mode, use the PredictProbablePairs method.
	 * During calculation, monitor progress using the TProgressDialog class and the Start/StopProgress methods of the RNA class.
	 * Neither of these methods require any error checking.
	 * After the resolution calculation is complete, use the error checker's isErrorStatus method to check for errors.
	 */
	if( error == 0 ) {

		// Show a message saying that the resolving calculation has started.
		if( mode == "MEA" ) {
			cout << "Calculating maximum expected accuracy structures..." << flush;
		} else if( mode == "ProbKnot" ) {
			cout << "Predicting pseudoknots..." << flush;
		} else {
			cout << "Calculating probable pairs..." << flush;
		}

		// Resolve the structures and check for errors.
		for( int i = 1; i <= seqNumber; i++ ) {
			int resolveError = 0;
			if( mode == "MEA" ) {
				if(i == 1)
				resolveError = object->MaximizeExpectedAccuracy( i, percent, maxStructures, windowSize, meaGamma );
			} else if( mode == "ProbKnot" ) {
				resolveError = object->ProbKnot( i, pkIterations, minHelixLength );
			} else {
				resolveError = object->PredictProbablePairs( i, threshold );
			}
			error = checker->isErrorStatus( resolveError );

			if( error != 0 ) { i += seqNumber; }
		}

		// If no error occurred, print message that resolving is done.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * Write CT output files using the WriteCt method.
	 * After writing is complete, use the error checker's isErrorStatus method to check for errors.
	 */

	if( error == 0 ) {

		// Show a message saying that the CT files are being written.
		cout << "Writing output ct files..." << flush;

		// Write the CT files and check for errors.
		for( int i = 1; i <= seqNumber; i++ ) {
			if(i == 1){
				int writeError = object->WriteCt( i, outputCtFiles[i-1].c_str() );
				error = checker->isErrorStatus( writeError ); 
				if( error != 0 ) { i += seqNumber; }
			}
		}
		

		// If no errors occurred, show a CT files writing completion message.
		if( error == 0 ) { cout << "done." << endl; }
	}

	// Delete the error checker and data structure.
	delete checker;
	delete object;

	// Print confirmation of run finishing.
	if( error == 0 ) { cout << calcType << " complete." << endl; }
	else { cerr << calcType << " complete with errors." << endl; }

	return error; // return exit code of program. 0-->success, non-zero-->error.
}

///////////////////////////////////////////////////////////////////////////////
// Print out a special usage message for this interface.
///////////////////////////////////////////////////////////////////////////////
void TurboHomology_Interface::usage( ParseCommandLine &parser ) {

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
	parser.printDescription( "Flag that can be used to specify a group of sequence files, from Seq1 to Seq<SequenceNumber>. Only one sequence file group can be specified without overwriting files. Group format: {seq1File;seq2File;seq3File;}" );
	cout << "OutCT" << endl;
	parser.printDescription( "Flag that can be used to specify a group of CT files, from CT1 to CT<SequenceNumber>. Only one CT file group can be specified without overwriting files. Group format: {ct1File;ct2File;ct3File;}" );
	cout << "Mode" << endl;
	parser.printDescription( "The mode in which TurboHomology is run. A mode can be specified in the following ways. 1) MEA (Maximum expected accuracy). 2) ProbKnot (Pseudoknots). 3) Threshold (Probable pairs)." );

	cout << "Required input when specifying files singly" << endl
	     << "-------------------------------------------" << endl
	     << "Every specified sequence file must have a corresponding CT file specified." << endl << endl;
	cout << "SequenceNumber" << endl;
	parser.printDescription( "The number of sequences given as input." );
	cout << "Seq1 ... Seq<SequenceNumber>" << endl;
	parser.printDescription( "Names of sequence files used as input, from 1 to SequenceNumber." );
	cout << "CT1 ... CT<SequenceNumber>" << endl;
	parser.printDescription( "Names of CT files written as output, from 1 to SequenceNumber." );
	cout << "Mode" << endl;
	parser.printDescription( "The mode in which TurboHomology is run. A mode can be specified in the following ways. 1) MEA (Maximum expected accuracy). 2) ProbKnot (Pseudoknots). 3) Threshold (Probable pairs)." );

	cout << "General options" << endl
	     << "---------------" << endl;
	cout << "Gamma" << endl;
	parser.printDescription( "The TurboHomology gamma. Default is 0.3." );
	cout << "Iterations" << endl;
	parser.printDescription( "The number of iterations TurboHomology goes through. Default is 3 iterations." );
	cout << "MaximumPairingDistance" << endl;
	parser.printDescription( "The maximum distance between nucleotides that can pair. For nucleotide i to pair with j, [i - j| < MaximumPairingDistance. This applies to each sequence. Default is no limit." );
#ifdef COMPILE_SMP
	cout << "Processors" << endl;
	parser.printDescription( "The number of processors on which the calculation runs. Default is 1." );
#endif
	cout << "Save1 ... Save<SequenceNumber>" << endl;
	parser.printDescription( "Names of save files written to by TurboHomology, from 1 to SequenceNumber. The number at the end of the flag (1 to SequenceNumber) identifies which sequence the save file will be written for." );
	cout << "SHAPE1 ... SHAPE<SequenceNumber>" << endl;
	parser.printDescription( "Names of SHAPE constraint files. The number at the end of the flag (1 to SequenceNumber) identifies which sequence the constraints will be applied to." );
	cout << "SHAPEintercept" << endl;
	parser.printDescription( "The SHAPE intercept. This value is only used when at least one SHAPE constraint file is specified. Default is 1.8 kcal/mol." );
	cout << "SHAPEslope" << endl;
	parser.printDescription( "The SHAPE slope. This value is only used when at least one SHAPE constraint file is specified. Default is -0.6 kcal/mol." );
	cout << "Temperature" << endl;
	parser.printDescription( "The temperature at which calculations are run, in Kelvin. Default is 310.15 K, which is 37 degrees C." );
	cout << "AlnFormat" << endl;
	parser.printDescription( "The final output alignment format. Default is Clustal." );
	cout << "ColumnNumber" << endl;
	parser.printDescription( "The maximum column number in the final output alignment. Default is 60." );

	cout << "Maximum expected accuracy (MEA) mode options" << endl
	     << "--------------------------------------------" << endl;
	cout << "MaxPercent" << endl;
	parser.printDescription( "The maximum percent energy difference. Default is 50 percent (Specified as 50, not 0.5)." );
	cout << "MaxStructures" << endl;
	parser.printDescription( "The maximum number of structures. Default is 1000 structures." );
	cout << "MeaGamma" << endl;
	parser.printDescription( "The weight given to pairs. Default is 1.0." );
	cout << "Window" << endl;
	parser.printDescription( "The window size. Default is 5 nucleotides." );

	cout << "Pseudoknot (ProbKnot) mode options" << endl
	     << "----------------------------------" << endl;
	cout << "MinHelixLength" << endl;
	parser.printDescription( "The minimum helix length. Default is 3 nucleotides." );
	cout << "PkIterations" << endl;
	parser.printDescription( "The number of iterations. Default is 1 iteration." );

	cout << "Probable pairs (Threshold) mode options" << endl
	     << "---------------------------------------" << endl;
	cout << "Threshold" << endl;
	parser.printDescription( "The threshold at which pairs should be included in a structure. This should be expressed as a number: 0.5 <= x <= 1.0. Default is 0, which signifies that structures should be generated at multiple thresholds: >= 0.99, >= 0.97, >= 0.95, >= 0.90, >= 0.80, >= 0.70, >= 0.60, and >= 0.50." );
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {
	TurboHomology_Interface runner;
	if (!runner.parse( argc, argv )) return 1;
	return runner.run() == 0 ? 0 : 1;
}
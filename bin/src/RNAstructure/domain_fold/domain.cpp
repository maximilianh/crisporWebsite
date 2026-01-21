/*
 *
 */

#include "domain.h"
#include "../src/domain_class.h"


///////////////////////////////////////////////////////////////////////////////
// Constructor.
//////////////////////////////////////////////////////////////////////////////
domainInterface::domainInterface() {

	// Initialize the calculation type description.
	calcType = "Single strand partition function";

	// Initialize the "experimental" offset.
	experimentalOffset = 0.0;

	// Initialize the "experimental" scaling.
	experimentalScaling = 1.0;

	// Initialize the SHAPE intercept.
	intercept = -0.6;

	// Initialize the nucleic acid type.
	alphabet = DT_RNA;

	// Initialize the maximum pairing distance between nucleotides.
	maxDistance = -1;

	// Initialize the SHAPE slope.
	slope = 1.8;

	// Initialize the calculation temperature.
	temperature = 310.15;

	// Use coaxial stacking by default
	disablecoax = false;

	// Do not treat the input as a pfs file by default
	pfs_input = false;

	constraintOutputFile = "";   // The constraints file.

	w2 = 9;

	circ = false;

	min_size = 20;
}

	
	//string experimentalFile; // The experimental pair bonus file.
	//string SHAPEFile;        // The SHAPE constraints file.
	

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool domainInterface::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "partition" );
	parser->addParameterDescription( "seq file", "The name of a file containing an input sequence." );
	parser->addParameterDescription( "pfs file", "The name of a partition function save file to which output will be written." );

	// Add the constraint file option.
	vector<string> constraintOptions;
	constraintOptions.push_back( "-c" );
	constraintOptions.push_back( "-C" );
	constraintOptions.push_back( "--constraint" );
	parser->addOptionFlagsWithParameters( constraintOptions, "Specify a constraints file to be applied. Default is to have no constraints applied." );

	// Add the DNA option.
	vector<string> dnaOptions;
	dnaOptions.push_back( "-d" );
	dnaOptions.push_back( "-D" );
	dnaOptions.push_back( "--DNA" );
	parser->addOptionFlagsNoParameters( dnaOptions, "Specify that the sequence is DNA, and DNA parameters are to be used. Default is to use RNA parameters." );

	// Add the Alphabet option.
	vector<string> alphaOptions;
	alphaOptions.push_back( "-a" );
	alphaOptions.push_back( "--alphabet" );
	parser->addOptionFlagsWithParameters( alphaOptions, "Specify the name of the thermodynamic parameter alphabet (e.g. \"rna\", \"dna\", or a custom extended/modified alphabet). The default is to use RNA parameters. (This option overrides the --dna flag)." );

	// Add the double stranded offset option.
	vector<string> doubleOffsetOptions;
	doubleOffsetOptions.push_back( "-dso" );
	doubleOffsetOptions.push_back( "-DSO" );
	doubleOffsetOptions.push_back( "--doubleOffset" );
	parser->addOptionFlagsWithParameters( doubleOffsetOptions, "Specify a double-stranded offset file, which adds energy bonuses to particular double-stranded nucleotides. Default is to have no double-stranded offset specified." );

	// Add the maximum pairing distance option.
	vector<string> distanceOptions;
	distanceOptions.push_back( "-md" );
	distanceOptions.push_back( "-MD" );
	distanceOptions.push_back( "--maxdistance" );
	parser->addOptionFlagsWithParameters( distanceOptions, "Specify a maximum pairing distance between nucleotides. Default is no restriction on distance between pairs." );

	// Add the SHAPE option.
	vector<string> shapeOptions;
	shapeOptions.push_back( "-sh" );
	shapeOptions.push_back( "-SH" );
	shapeOptions.push_back( "--SHAPE" );
	parser->addOptionFlagsWithParameters( shapeOptions, "Specify a SHAPE constraints file to be applied. These constraints are pseudoenergy constraints. Default is to have no constraints applied." );

	// Add the SHAPE intercept option.
	vector<string> shapeInterceptOptions;
	shapeInterceptOptions.push_back( "-si" );
	shapeInterceptOptions.push_back( "-SI" );
	shapeInterceptOptions.push_back( "--SHAPEintercept" );
	parser->addOptionFlagsWithParameters( shapeInterceptOptions, "Specify an intercept used with SHAPE constraints. Default is -0.6 kcal/mol." );

	// Add the SHAPE slope option.
	vector<string> shapeSlopeOptions;
	shapeSlopeOptions.push_back( "-sm" );
	shapeSlopeOptions.push_back( "-SM" );
	shapeSlopeOptions.push_back( "--SHAPEslope" );
	parser->addOptionFlagsWithParameters( shapeSlopeOptions, "Specify a slope used with SHAPE constraints. Default is 1.8 kcal/mol." );

	// Add the temperature option.
	vector<string> tempOptions;
	tempOptions.push_back( "-t" );
	tempOptions.push_back( "-T" );
	tempOptions.push_back( "--temperature" );
	parser->addOptionFlagsWithParameters( tempOptions, "Specify the temperature at which calculation takes place in Kelvin. Default is 310.15 K, which is 37 degrees C." );

	// Add the experimental pair bonus option.
	vector<string> experimentalOptions;
	experimentalOptions.push_back( "-x" );
	experimentalOptions.push_back( "-X" );
	experimentalOptions.push_back( "--experimentalPairBonus" );
	parser->addOptionFlagsWithParameters( experimentalOptions, "Input text file with bonuses (in kcal) as a matrix. As with SHAPE, bonuses will be applied twice to internal base pairs, once to edge base pairs, and not at all to single stranded regions. Default is no experimental pair bonus file specified." );

	// Add the experimental pair bonus offset option.
	vector<string> experimentalOffsetOptions;
	experimentalOffsetOptions.push_back( "-xo" );
	parser->addOptionFlagsWithParameters( experimentalOffsetOptions, "Specify an intercept (overall offset) to use with the 2D experimental pair bonus file. Default is 0.0 (no change to input bonuses)." );

	// Add the experimental pair bonus scaling option.
	vector<string> experimentalScalingOptions;
	experimentalScalingOptions.push_back( "-xs" );
	parser->addOptionFlagsWithParameters( experimentalScalingOptions, "Specify a number to multiply the experimental pair bonus matrix by. Default is 1.0 (no change to input bonuses)." );

	// Add the disablecoax option.
	vector<string> disablecoaxOptions;
	disablecoaxOptions.push_back( "--disablecoax" );
	parser->addOptionFlagsNoParameters( disablecoaxOptions, "Specify that coaxial stacking recusions should not be used. This option uses a less realistic energy function in exchange for a faster calculation." );

	// Add the pfs input option.
	vector<string> pfsinputOptions;
	pfsinputOptions.push_back( "--pfs" );
	pfsinputOptions.push_back( "-p" );
	parser->addOptionFlagsNoParameters( pfsinputOptions, "Specify that the input is a partition save file" );

	// Add the pfs input option.
	vector<string> circOptions;
	circOptions.push_back( "--circ" );
	parser->addOptionFlagsNoParameters( circOptions, "Specify that the domains can be circular (span the 5' and 3' ends of the sequence" );

	// Add the constraint output option.
	vector<string> constraintOutput;
	constraintOutput.push_back( "-co" );
	parser->addOptionFlagsWithParameters( constraintOutput, "Specify a filepath for the constraint output." );

	// Add the pair probability input option.
	vector<string> ppInputFileOption;
	ppInputFileOption.push_back( "--pp" );
	parser->addOptionFlagsWithParameters( ppInputFileOption, "Specify a filepath for the pair probability input." );

	vector<string> minDomainOption;
	minDomainOption.push_back( "-m" );
	minDomainOption.push_back( "--minDomainSize" );	
	parser->addOptionFlagsWithParameters( minDomainOption, "Specify a minimum domain length." );

	// Add the pair probability input option.
	vector<string> PenaltyOption;
	PenaltyOption.push_back( "-o" );
	PenaltyOption.push_back( "--OuterPenalty" );
	parser->addOptionFlagsWithParameters( PenaltyOption, "Specify the wieght for the penalty for pairing probability outside the domain. The default value is 9 (outside pairing probability has 9x the weight as inner pairing probability" );
	
	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		seqFile = parser->getParameter( 1 );
		pfsFile = parser->getParameter( 2 );
	}

	// Get the constraint file option.
	if( !parser->isError() ) { constraintFile = parser->getOptionString( constraintOptions, true ); }

	// Get the DNA option.
	if( !parser->isError() && parser->contains( dnaOptions ) )
		alphabet = DT_DNA; // use DNA (unless overridden by alphabet)

	// Get the Alphabet option.
	if (!parser->isError() && parser->contains(alphaOptions))
		  alphabet = parser->getOptionString(alphaOptions, false).c_str();

	// Get the double stranded offset option.
	if( !parser->isError() ) { doubleOffsetFile = parser->getOptionString( doubleOffsetOptions, true ); }

	// Get the maximum distance option.
	if( !parser->isError() ) {
		parser->setOptionInteger( distanceOptions, maxDistance );
		bool badDistance =
		  ( maxDistance < 0 ) &&
		  ( maxDistance != -1 );
		if( badDistance ) { parser->setError( "maximum pairing distance" ); }
	}

	// Get the SHAPE data and options.
	if( !parser->isError() ) {
		SHAPEFile = parser->getOptionString( shapeOptions );
		if( !parser->isError() ) { parser->setOptionDouble( shapeInterceptOptions, intercept ); }
		if( !parser->isError() ) { parser->setOptionDouble( shapeSlopeOptions, slope ); }
	}

	// Get the temperature option.
	if( !parser->isError() ) {
		parser->setOptionDouble( tempOptions, temperature );
		if( temperature < 0 ) { parser->setError( "temperature" ); }
	}

	// Get the experimental bonus data and options.
	if( !parser->isError() ) {
		experimentalFile = parser->getOptionString( experimentalOptions );
		if( !parser->isError() ) { parser->setOptionDouble( experimentalOffsetOptions, experimentalOffset ); }
		if( !parser->isError() ) { parser->setOptionDouble( experimentalScalingOptions, experimentalScaling ); }
	}

	// Get the disablecoax option.
    // true by default, false if the flag is used
	if( !parser->isError() ) { disablecoax = parser->contains( disablecoaxOptions ); }

	// Get the pfsinput option.
	if( !parser->isError() && parser->contains( pfsinputOptions ) )
		pfs_input = true;

	if( !parser->isError() ) { constraintOutputFile = parser->getOptionString( constraintOutput, false ); }

	if( !parser->isError() ) { ppInputFile = parser->getOptionString( ppInputFileOption, true ); }
	
	if( !parser->isError() && parser->contains( circOptions ) )
		circ = true;

	if( !parser->isError() ) { parser->setOptionDouble( PenaltyOption, w2 ); }

	if( !parser->isError() ) { parser->setOptionInteger( minDomainOption, min_size ); }
		

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void domainInterface::run() {

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
	cout << "Initializing nucleic acids..." << flush;

	//RNA* strand = new RNA();
	RNA* strand = new RNA( seqFile.c_str(), FILE_SEQ, alphabet.c_str() );

	if (pfs_input) {
		RNA* strand = new RNA( seqFile.c_str(), 3, alphabet.c_str() );
	}
	else{
		//RNA* strand = new RNA( seqFile.c_str(), FILE_SEQ, alphabet.c_str() );
	}

	
	domain_calc Domain_Instance = domain_calc(circ, w2, strand->GetSequenceLength());

	ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );
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
		int tempError = strand->SetTemperature( temperature );
		error = checker->isErrorStatus( tempError );

		// If no error occurred, print a message saying that temperature is set.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * Set maximum pairing distance using the ForceMaximumPairingDistance method.
	 */
	if( error == 0 && maxDistance != -1 ) {

		// Show a message saying that the maximum pairing distance is being set.
		cout << "Setting maximum distance between paired nucleotides..." << flush;

		// Set the maximum pairing distance and check for errors.
		int distError = strand->ForceMaximumPairingDistance( maxDistance );
		error = checker->isErrorStatus( distError );

		// If no error occurred, print a message saying that maximum pairing distance was set.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * Add constraints if files have been given for their inclusion.
	 * For folding constraints, use the ReadConstraints method.
	 * For SHAPE constraints, use the ReadSHAPE method.
	 * For double strand offset, use the ReadDSO method.
	 * For experimental pair bonuses, use the ReadExperimentalPairBonus method.
	 * After each constraint type, check the error status of the strand as above.
	 */

	// Determine if constraints should be applied.
	bool applyConstraints =
		( constraintFile != "" ) ||
		( SHAPEFile != "" ) ||
		( doubleOffsetFile != "" ) ||
		( experimentalFile != "" );


	// If constraints should be applied, do so.
	if( error == 0 && applyConstraints ) {

		// Show a message saying that constraints are being applied.
		cout << "Applying constraints..." << flush;
		int constraintError = 0;

		// Read folding constraints, if applicable.
		if( constraintFile != "" ) {
			constraintError = strand->ReadConstraints( constraintFile.c_str() );
			error = checker->isErrorStatus( constraintError );
		}

		// Read SHAPE constraints, if applicable.
		if( error == 0 && SHAPEFile != "" ) {
			constraintError = strand->ReadSHAPE( SHAPEFile.c_str(), slope, intercept );
			error = checker->isErrorStatus( constraintError );
		}

		// Read double strand offset, if applicable.
		if( error == 0 && doubleOffsetFile != "" ) {
			constraintError = strand->ReadDSO( doubleOffsetFile.c_str() );
			error = checker->isErrorStatus( constraintError );
		}

		// Read experimental pair bonus constraints, if applicable.
		if( error == 0 && experimentalFile != "" ) {
			constraintError = strand->ReadExperimentalPairBonus( experimentalFile.c_str(), experimentalOffset, experimentalScaling );
			error = checker->isErrorStatus( constraintError );
		}

		// If no error occurred, print a message saying that constraints were applied.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * Calculate the partition function using the PartitionFunction method.
	 * During calculation, monitor progress using the TProgressDialog class and the Start/StopProgress methods of the RNA class.
	 * Neither of these methods require any error checking.
	 * After the main calculation is complete, use the error checker's isErrorStatus method to check for errors.
	 */
	if (!pfs_input && ppInputFile=="") {
		if( error == 0 ) {

			// Show a message saying that the main calculation has started.
			cout << "Calculating partition function..." << flush;

			// Create the progress monitor.
			TProgressDialog* progress = new TProgressDialog();
			strand->SetProgress( *progress );

			// Do the main calculation and check for errors.
			int mainCalcError = strand->PartitionFunction( pfsFile.c_str(), temperature, disablecoax);
			error = checker->isErrorStatus( mainCalcError );

			// Delete the progress monitor.
			strand->StopProgress();
			delete progress;

			// Copy pairing probability data into the domain calculator
			Domain_Instance.read_pfs_file(strand);

			// If no error occurred, print a message saying that the main calculation is done.
			if( error == 0 ) { cout << "done." << endl; }
		}
	}
	else{
		if (ppInputFile == ""){
			// The information is read in during the RNA strand initialization, nothing needs to be done here
			cout << "Pair probability information read in from pfs file" << endl;
			Domain_Instance.read_pfs_file(strand);
			cout << "Done" << endl;
		}
		else {
			cout << "Reading pair probability information from file." << endl;
			Domain_Instance.read_pp_file(ppInputFile);
			cout << "Done" << endl;
		}
	}

	/*
	 * Perform the domain find calculations
	 */
	cout << "Calculating domains" << endl;
	//cout << "Filling domain matrices" << endl;
	Domain_Instance.fill_domain_matrices();
	//cout << "Calculating domain scores" << endl;
	Domain_Instance.calc_domain_scores();
	//cout << "Finding optimum domain segmentation" << endl;
	Domain_Instance.complete_domain_assignment(min_size);
	Domain_Instance.copy_constraints(strand->GetStructure());

	if (constraintOutputFile != "")
		Domain_Instance.writeconstraints(constraintOutputFile);

	// Show a message saying that the main calculation has started.
	cout << "Calculating second partition function..." << flush;

	// Create the progress monitor.
	TProgressDialog* progress = new TProgressDialog();
	strand->SetProgress( *progress );

	// Do the main calculation and check for errors.
	int mainCalcError = strand->PartitionFunction( pfsFile.c_str(), temperature, disablecoax);
	error = checker->isErrorStatus( mainCalcError );

	// Delete the progress monitor.
	strand->StopProgress();
	delete progress;
//*/
	// Delete the error checker and data structure.
	delete checker;
	delete strand;
	//delete Domain_Instance;
	
	// Print confirmation of run finishing.
	if( error == 0 ) { cout << calcType << " complete." << endl; }
	else { cerr << calcType << " complete with errors." << endl; }
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	domainInterface* runner = new domainInterface();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}

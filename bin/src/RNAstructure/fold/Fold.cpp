/*
 * A program that folds a strand of nucleic acids.
 * This strand of nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#include "Fold.h"
#include <time.h>
///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
Fold_Interface::Fold_Interface() {

	// Initalize the calculation type description.
	calcType = "Single strand folding";

	// Initialize the "experimental" offset.
	experimentalOffset = 0.0;

	// Initialize the "experimental" scaling.
	experimentalScaling = 1.0;

	// Initialize the SHAPE intercept.
	intercept = -0.6;

	// Initialize the single stranded SHAPE intercept.
	interceptSingle = 0;

	// Initialize the nucleic acid type.
	alphabet = DT_RNA;

	// Initialize the maximum internal bulge loop size.
	maxLoop = 30;

	// Initialize the maximum pairing distance between nucleotides.
	maxDistance = -1;

	// Initialize the maximum number of structures.
	maxStructures = 20;

	// Initialize the maximum percent energy difference.
	percent = 10;

	// Initialize the SHAPE slope.
	slope = 1.8;

	//Initialize the differntial SHAPE slope.
	Dslope = 2.11;

	// Initialize the single stranded SHAPE slope.
	slopeSingle = 0;

	// Initialize the calculation temperature.
	temperature = 310.15;

	// The number of bootstraping iterations to be done.
	bootstrap = 0;

	// Initialize the folding window size.
	windowSize = -1;

	//  Initialize the quickfold (mfe only) variable.
	quickfold = false;

	//By default, filter out isolated pairs.
	isolated = false;

	//  Initialize the internal loop search mode flag
	simple_iloops = true;

	//  Initialize the flag to disable coaxial stacking recursions
	disablecoax = true;

	sequenceName = "";

	quiet = false;

	useBracketNotation = false;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool Fold_Interface::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "Fold" );
	parser->addParameterDescription( "sequence file", "The name of a file containing an input sequence. Acceptable formats include SEQ, FASTA and raw-sequence plain-text files.\nIf the name is a hyphen (-), the file will be read from standard input (STDIN)." );
	parser->addParameterDescription( "CT file", "The name of a CT file to which output will be written. If the --bracket flag is present, output will be written as a dot-bracket file.\nIf the file name is a hyphen (-), the structure will be written to standard output (STDOUT) instead of a file." );

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
	vector<string> alphabetOptions;
	alphabetOptions.push_back( "-a" );
	alphabetOptions.push_back( "--alphabet" );
	parser->addOptionFlagsWithParameters( alphabetOptions, "Specify the name of a folding alphabet and associated nearest neighbor parameters. The alphabet is the prefix for the thermodynamic parameter files, e.g. \"rna\" for RNA parameters or \"dna\" for DNA parameters or a custom extended/modified alphabet. The thermodynamic parameters need to reside in the at the location indicated by environment variable DATAPATH. The default is \"rna\" (i.e. use RNA parameters). This option overrides the --DNA flag." );

	// Add the quickfold option.
	vector<string> quickfoldOptions;
	quickfoldOptions.push_back( "-mfe" );
	quickfoldOptions.push_back( "-MFE" );
	quickfoldOptions.push_back( "--MFE" );
	parser->addOptionFlagsNoParameters( quickfoldOptions, "Specify that only the minimum free energy structure is needed.  No savefiles can be generated.  This saves nearly half the calculation time, but provides less information." );

	// Add the simple_iloops option.
	vector<string> iloopmodeOptions;
	iloopmodeOptions.push_back( "-y" );
	iloopmodeOptions.push_back( "-Y" );
	iloopmodeOptions.push_back( "--simple_iloops" );
	parser->addOptionFlagsNoParameters( iloopmodeOptions, "Specify that the O(N^3) internal loop search should be used. This speeds up the calculation if large internal loops are allowed using the -l option." );

	// Add the disablecoax option.
	vector<string> disablecoaxOptions;
	disablecoaxOptions.push_back( "--disablecoax" );
	parser->addOptionFlagsNoParameters( disablecoaxOptions, "Specify that coaxial stacking recusions should not be used. This option uses a less realistic energy function in exchange for a faster calculation." );

	// Add the double stranded offset option.
	vector<string> doubleOffsetOptions;
	doubleOffsetOptions.push_back( "-dso" );
	doubleOffsetOptions.push_back( "-DSO" );
	doubleOffsetOptions.push_back( "--doubleOffset" );
	parser->addOptionFlagsWithParameters( doubleOffsetOptions, "Specify a double-stranded offset file, which adds energy bonuses to particular double-stranded nucleotides. Default is to have no double-stranded offset specified." );

	// Add the allow isolated pair option.
	vector<string> isolatedOptions;
	isolatedOptions.push_back("-i");
	isolatedOptions.push_back("-I");
	isolatedOptions.push_back("--isolated");
	parser->addOptionFlagsNoParameters(isolatedOptions, "Allow isolated base pairs.  The default is to use a heuristic to forbid isolated base pairs.");

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

	// Add the maximum pairing distance option.
	vector<string> distanceOptions;
	distanceOptions.push_back( "-md" );
	distanceOptions.push_back( "-MD" );
	distanceOptions.push_back( "--maxdistance" );
	parser->addOptionFlagsWithParameters( distanceOptions, "Specify a maximum pairing distance between nucleotides. Default is no restriction on distance between pairs." );

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

	// Add the SHAPE option.
	vector<string> shapeOptions;
	shapeOptions.push_back( "-sh" );
	shapeOptions.push_back( "-SH" );
	shapeOptions.push_back( "--SHAPE" );
	parser->addOptionFlagsWithParameters( shapeOptions, "Specify a SHAPE restraints file to be applied. These constraints are pseudoenergy restraints. Default is to have no restraints applied." );

	// Add the SHAPE option.
	vector<string> DshapeOptions;
	DshapeOptions.push_back( "-dsh" );
	DshapeOptions.push_back( "-DSH" );
	DshapeOptions.push_back( "--DSHAPE" );
	parser->addOptionFlagsWithParameters( DshapeOptions, "Specify a differential SHAPE restraints file to be applied. These constraints are pseudoenergy restraints. Default is to have no restraints applied." );

	// Add the DMS option.
	vector<string> dmsOptions;
	dmsOptions.push_back( "-dms" );
	dmsOptions.push_back( "-DMS" );
	dmsOptions.push_back( "--DMS" );
	parser->addOptionFlagsWithParameters( dmsOptions, "Specify a DMS restraints file to be applied. These restraints are pseudoenergy constraints. Default is to have no restraints applied." );

    // Add the DMSnt option.
	vector<string> dmsntOptions;
	dmsntOptions.push_back( "-dmsnt" );
	dmsntOptions.push_back( "-DMSNT" );
    dmsntOptions.push_back( "-DMSnt" );
	dmsntOptions.push_back( "--DMSNT" );
	parser->addOptionFlagsWithParameters( dmsntOptions, "Specify a DMS NT restraint file to be applied. These restraints are pseudoenergy constraints applied with NT-specific potentials. Default is to have no restraints applied." );


	// Add the CMCT option.
	vector<string> cmctOptions;
	cmctOptions.push_back( "-cmct" );
	cmctOptions.push_back( "-CMC" );
	cmctOptions.push_back( "--CMCT" );
	parser->addOptionFlagsWithParameters( cmctOptions, "Specify a CMCT constraints file to be applied. These constraints are pseudoenergy constraints. Default is to have no constraints applied." );


	// Add the SHAPE intercept option.
	vector<string> shapeInterceptOptions;
	shapeInterceptOptions.push_back( "-si" );
	shapeInterceptOptions.push_back( "-SI" );
	shapeInterceptOptions.push_back( "--SHAPEintercept" );
	parser->addOptionFlagsWithParameters( shapeInterceptOptions, "Specify an intercept used with SHAPE restraints. Default is -0.6 kcal/mol." );

	// Add the SHAPE slope option.
	vector<string> shapeSlopeOptions;
	shapeSlopeOptions.push_back( "-sm" );
	shapeSlopeOptions.push_back( "-SM" );
	shapeSlopeOptions.push_back( "--SHAPEslope" );
	parser->addOptionFlagsWithParameters( shapeSlopeOptions, "Specify a slope used with SHAPE renstraints. Default is 1.8 kcal/mol." );

	// Add the Differentail SHAPE slope option.
	vector<string> DshapeSlopeOptions;
	DshapeSlopeOptions.push_back( "-dsm" );
	DshapeSlopeOptions.push_back( "-DSM" );
	DshapeSlopeOptions.push_back( "--DSHAPEslope" );
	parser->addOptionFlagsWithParameters( DshapeSlopeOptions, "Specify a slope used with differential SHAPE restraints. Default is 2.11 kcal/mol." );


	// Add the single stranded offset option.
	vector<string> singleOffsetOptions;
	singleOffsetOptions.push_back( "-sso" );
	singleOffsetOptions.push_back( "-SSO" );
	singleOffsetOptions.push_back( "--singleOffset" );
	parser->addOptionFlagsWithParameters( singleOffsetOptions, "Specify a single-stranded offset file, which adds energy bonuses to particular single-stranded nucleotides. Default is to have no single-stranded offset specified." );

	// Add the temperature option.
	vector<string> tempOptions;
	tempOptions.push_back( "-t" );
	tempOptions.push_back( "-T" );
	tempOptions.push_back( "--temperature" );
	parser->addOptionFlagsWithParameters( tempOptions, "Specify the temperature at which calculation takes place in Kelvin. Default is 310.15 K, which is 37 degrees C." );

	// Add the bootstrap option.
	vector<string> bootstrapOptions;
	bootstrapOptions.push_back( "-boot" );
	bootstrapOptions.push_back( "-B" );
	bootstrapOptions.push_back( "--bootstrap" );
	parser->addOptionFlagsWithParameters( bootstrapOptions, "Specify the number of bootstrap iterations to be done to retrieve base pair confidence. Defaults to no bootstrapping." );


	// Add the unpaired SHAPE intercept option.
	vector<string> shapeInterceptUnpairedOptions;
	shapeInterceptUnpairedOptions.push_back( "-usi" );
	shapeInterceptUnpairedOptions.push_back( "-USI" );
	shapeInterceptUnpairedOptions.push_back( "--unpairedSHAPEintercept" );
	parser->addOptionFlagsWithParameters( shapeInterceptUnpairedOptions, "Specify an intercept used with unpaired SHAPE constraints. Default is 0 kcal/mol." );

	// Add the unpaired SHAPE slope option.
	vector<string> shapeSlopeUnpairedOptions;
	shapeSlopeUnpairedOptions.push_back( "-usm" );
	shapeSlopeUnpairedOptions.push_back( "-USM" );
	shapeSlopeUnpairedOptions.push_back( "--unpairedSHAPEslope" );
	parser->addOptionFlagsWithParameters( shapeSlopeUnpairedOptions, "Specify a slope used with unpaired SHAPE constraints. Default is 0 kcal/mol." );

	// Add the window size option.
	vector<string> windowOptions;
	windowOptions.push_back( "-w" );
	windowOptions.push_back( "-W" );
	windowOptions.push_back( "--window" );
	parser->addOptionFlagsWithParameters( windowOptions, "Specify a window size. Default is determined by the length of the sequence." );

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


	// Add the option to change the output behavior of non-critical warnings.
	vector<string> warningsOption;
	warningsOption.push_back( "--warnings" );
	warningsOption.push_back( "--warn" );
	parser->addOptionFlagsWithParameters( warningsOption, "Set the behavior for non-critical warnings (e.g. related to invalid nucleotide positions or duplicate data points in SHAPE data). Valid values are: "
	"\n\t* on  -- Warnings are written to standard output. (default)"
	"\n\t* err -- Warnings are sent to STDERR. This can be used in automated scripts etc to detect problems."
	"\n\t* off -- Do not display warnings at all (not recommended)."
	);

	vector<string> nameOption = parser->addFlag(true, "--name", "Specify a name for the sequence. This will override the name in the sequence file.");
	vector<string> quietOption = parser->addFlag(false, "-q --quiet", "Suppress unnecessary output. This option is implied when the output file is '-' (STDOUT).");
	vector<string> bracketOption = parser->addFlag(false, "-k --bracket", "Write the predicted structure in dot-bracket notation (DBN) instead of CT format.");

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		seqFile = parser->getParameter( 1 );
		ctFile = parser->getParameter( 2 );
	}

	parser->setOptionString(nameOption, sequenceName); // allow user to override the sequence name
	quiet = parser->contains(quietOption) || isStdIoFile(ctFile.c_str()); // suppress unnecessary output if --quiet option is present or if the CT file is output to stdout.
	useBracketNotation = parser->contains(bracketOption); // write in dot-bracket notation.

	// Get the constraint file option.
	if( !parser->isError() ) { constraintFile = parser->getOptionString( constraintOptions, true ); }

	// Get the DNA option.
	if( !parser->isError() && parser->contains( dnaOptions ) )
		alphabet = DT_DNA; // use DNA (unless overridden by alphabet)

	// Get the Alphabet option.
	if (!parser->isError() && parser->contains(alphabetOptions))
		  alphabet = parser->getOptionString(alphabetOptions, false).c_str();

	// Get the quickfold option.
	if( !parser->isError() ) { quickfold = parser->contains( quickfoldOptions ); }

	// Get the islated pair option.
	if (!parser->isError()) { isolated = parser->contains(isolatedOptions); }

	// Get the simple_iloops option.
    // true by default, false if the flag is used
	if( !parser->isError() ) { simple_iloops = !parser->contains( iloopmodeOptions ); }

	// Get the disablecoax option.
    // true by default, false if the flag is used
	if( !parser->isError() ) { disablecoax = parser->contains( disablecoaxOptions ); }

	// Get the double stranded offset option.
	if( !parser->isError() ) { doubleOffsetFile = parser->getOptionString( doubleOffsetOptions, true ); }

	// Get the maximum loop size option.
	if( !parser->isError() ) {
		parser->setOptionInteger( loopOptions, maxLoop );
		if( maxLoop < 0 ) { parser->setError( "maximum loop size" ); }
	}

	// Get the maximum distance option.
	if( !parser->isError() ) {
		parser->setOptionInteger( distanceOptions, maxDistance );
		bool badDistance =
		  ( maxDistance < 0 ) &&
		  ( maxDistance != -1 );
		if( badDistance ) { parser->setError( "maximum pairing distance" ); }
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

	// Determine the behavior of ReadSHAPE when encountering multiple 
	// datapoints for the same nucleotide
	if(!parser->isError() && parser->contains(warningsOption)) {
		string value = parser->getOptionString(warningsOption, false/*not a file*/);
		toUpper(value);
		structure::ShowWarnings = value=="ON"?1:value=="OFF"?0:value=="ERR"?2:255;
		if (structure::ShowWarnings==255) parser->setError("warnings");
	}

	// Get the SHAPE, DMS, and CMCT data and options.
	if( !parser->isError() ) {
		SHAPEFile = parser->getOptionString( shapeOptions );
		DSHAPEFile = parser->getOptionString( DshapeOptions );
		if( !parser->isError() ) { DMSFile = parser->getOptionString( dmsOptions ); }
        if( !parser->isError() ) { DMSNTFile = parser->getOptionString( dmsntOptions ); }
		if( !parser->isError() ) { CMCTFile = parser->getOptionString( cmctOptions ); }
		if( !parser->isError() ) { parser->setOptionDouble( shapeInterceptOptions, intercept ); }
		if( !parser->isError() ) { parser->setOptionDouble( shapeSlopeOptions, slope ); }
		if( !parser->isError() ) { parser->setOptionDouble( DshapeSlopeOptions, Dslope ); }
		if( !parser->isError() ) { parser->setOptionDouble( shapeInterceptUnpairedOptions, interceptSingle ); }
		if( !parser->isError() ) { parser->setOptionDouble( shapeSlopeUnpairedOptions, slopeSingle ); }
	}

	// Get the single stranded offset option.
	if( !parser->isError() ) { singleOffsetFile = parser->getOptionString( singleOffsetOptions, true ); }

	// Get the temperature option.
	if( !parser->isError() ) {
		parser->setOptionDouble( tempOptions, temperature );
		if( temperature < 0 ) { parser->setError( "temperature" ); }
	}

	// Get the number of bootstraps.
	if( !parser->isError() ) {
		parser->setOptionDouble( bootstrapOptions, bootstrap );
		if( bootstrap < 0 ) { parser->setError( "bootstrap" ); }
	}

	// Get the window size option.
	if( !parser->isError() ) {
		parser->setOptionInteger( windowOptions, windowSize );
		bool badWindow =
		  ( windowSize < 0 ) &&
		  ( windowSize != -1 );
		if( badWindow ) { parser->setError( "window size" ); }
	}

	// Get the experimental bonus data and options.
	if( !parser->isError() ) {
		experimentalFile = parser->getOptionString( experimentalOptions );
		if( !parser->isError() ) { parser->setOptionDouble( experimentalOffsetOptions, experimentalOffset ); }
		if( !parser->isError() ) { parser->setOptionDouble( experimentalScalingOptions, experimentalScaling ); }
	}

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

string Fold_Interface::sample_file(string shapefile, int numnuc, int iter) {
	char numstr[21];
	ifstream infile;
	sprintf(numstr, "%d", iter);
	string outname = shapefile + "_boot" + numstr;
	ofstream outfile;
	outfile.open(outname.c_str());
	infile.open(shapefile.c_str());
	double *SHAPEdata = new double [numnuc];
	int index, ridx, i;
	double value;

	for( i=0;i<numnuc;i++ ) { SHAPEdata[i] = -999; }

	while( infile >> index >> value ) { SHAPEdata[index] = value; }

	srand( time(0) );
	for( i=0;i<numnuc;i++ ) {
		ridx = rand() % numnuc;
		if( SHAPEdata[ridx] != -999 ) { SHAPEdata[ridx] += SHAPEdata[ridx]; }
	}

	for( i=0;i<numnuc;i++ ) { outfile << i << " " << SHAPEdata[i] << "\n"; }
	
	outfile.close();
	infile.close();

	return outname;
}
///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
bool Fold_Interface::run() {

	// Create a variable that handles errors.
	int error = 0;
	char numstr[21];
	int b_iter;
	string consFile;

	/*
	 * Use the constructor for RNA that specifies a filename.
	 * Specify type = 2 (sequence file).
	 * isRNA identifies whether the strand is RNA (true) or DNA (false).
	 *
	 * After construction of the strand data structure, create the error checker which monitors for errors.  
	 * Throughout, the error status of the calculation is checked with a variant of the isErrorStatus method, which returns 0 if no error occurred.
	 * The calculation proceeds as long as error = 0.
	 */
	if (!quiet) cout << "Initializing nucleic acids..." << flush;
	RNA* strand = new RNA( seqFile.c_str(), FILE_SEQ, alphabet.c_str() );
	ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );
	error = checker->isErrorStatus();
	if( error == 0 ) { if (!quiet) cout << "done." << endl; }
	if (!sequenceName.empty()) strand->SetSequenceLabel(sequenceName);

	/*
	 * Set the window size, based on the length of the sequence given as input.
	 * Only do this if window size hasn't been set on the command line.
	 * Use method GetSequenceLength to get the length of the sequence.
	 * The window sizes in relation to the length are hardcoded values.
	 */
	if( windowSize == -1 && error == 0 ) {
		int length = strand->GetSequenceLength();
		windowSize =
			( length > 1200 ) ? 20 :
			( length > 800 ) ? 15 :
			( length > 500 ) ? 11 :
			( length > 300 ) ? 7 :
			( length > 120 ) ? 5 :
			( length > 50 ) ? 3 :
			2;
	}

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
	 * Set maximum pairing distance using the ForceMaximumPairingDistance method.
	 */
	if( error == 0 && maxDistance != -1 ) {

		// Show a message saying that the maximum pairing distance is being set.
		if (!quiet) cout << "Setting maximum distance between paired nucleotides..." << flush;

		// Set the maximum pairing distance and check for errors.
		int distError = strand->ForceMaximumPairingDistance( maxDistance );
		error = checker->isErrorStatus( distError );

		// If no error occurred, print a message saying that maximum pairing distance was set.
		if( error == 0 ) { if (!quiet) cout << "done." << endl; }
	}

	/*
	 * Add constraints if files have been given for their inclusion.
	 * For folding constraints, use the ReadConstraints method.
	 * For SHAPE constraints, use the ReadSHAPE method.
	 * For single strand offset, use the ReadSSO method.
	 * For double strand offset, use the ReadDSO method.
	 * For experimental pair bonuses, use the ReadExperimentalPairBonus method.
	 * After each constraint type, check the error status of the strand as above.
	 */

	// Determine if constraints should be applied.
	bool applyConstraints =
		( constraintFile != "" ) ||
		( SHAPEFile != "" ) ||
		( DSHAPEFile != "" ) ||
		( CMCTFile != "" ) ||
		( DMSFile != "" ) ||
        ( DMSNTFile != "") ||
		( singleOffsetFile != "" ) ||
		( doubleOffsetFile != "" ) ||
		( experimentalFile != "" );

	for( b_iter=0; b_iter <= bootstrap; b_iter++) {
		// If constraints should be applied, do so.
		if( error == 0 && applyConstraints ) {

			// Show a message saying that constraints are being applied.
			if (!quiet) cout << "Applying constraints..." << flush;
			int constraintError = 0;

			// Read folding constraints, if applicable.
			if( constraintFile != "" ) {
				constraintError = strand->ReadConstraints( constraintFile.c_str() );
				error = checker->isErrorStatus( constraintError );
			}

			// Read SHAPE constraints
			if( error == 0 && SHAPEFile != "" ) {
				if( b_iter > 0 ) { consFile = sample_file(SHAPEFile, strand->GetSequenceLength(), b_iter); }
				else { consFile = SHAPEFile; }
				constraintError = strand->ReadSHAPE( consFile.c_str(), slope, intercept, slopeSingle, interceptSingle, RESTRAINT_SHAPE );
				error = checker->isErrorStatus( constraintError );
			}

			// Read differential SHAPE constraints
			if( error == 0 && DSHAPEFile != "" ) {
				if( b_iter > 0 ) { consFile = sample_file(DSHAPEFile, strand->GetSequenceLength(), b_iter); }
				else { consFile = DSHAPEFile; }
				constraintError = strand->ReadSHAPE( consFile.c_str(), Dslope, 0, 0, 0, RESTRAINT_SHAPE_DIFF );
				error = checker->isErrorStatus( constraintError );
			}

			// Read DMS constraints.
			if( error == 0 && DMSFile != "" ) {
				if( b_iter > 0 ) { consFile = sample_file(DMSFile, strand->GetSequenceLength(), b_iter); }
				else { consFile = DMSFile; }
				constraintError = strand->ReadSHAPE( consFile.c_str(), slope, intercept, slopeSingle, interceptSingle, RESTRAINT_DMS );
				error = checker->isErrorStatus( constraintError );
			}
    		
            // Read DMSnt constraints.
			if( error == 0 && DMSNTFile != "" ) {
				if( b_iter > 0 ) { consFile = sample_file(DMSNTFile, strand->GetSequenceLength(), b_iter); }
				else { consFile = DMSNTFile; }
				constraintError = strand->ReadSHAPE( consFile.c_str(), slope, intercept, slopeSingle, interceptSingle, RESTRAINT_DMSNT );
				error = checker->isErrorStatus( constraintError );
			}
         

			// Read CMCT constraints.
			if( error == 0 && CMCTFile != "" ) {
				if( b_iter > 0 ) { consFile = sample_file(CMCTFile, strand->GetSequenceLength(), b_iter); }
				else { consFile = CMCTFile; }
				constraintError = strand->ReadSHAPE( consFile.c_str(), slope, intercept, slopeSingle, interceptSingle, RESTRAINT_CMCT );
				error = checker->isErrorStatus( constraintError );
			}


			// Read single strand offset, if applicable.
			if( error == 0 && singleOffsetFile != "" ) {
				constraintError = strand->ReadSSO( singleOffsetFile.c_str() );
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
			if( error == 0 ) { if (!quiet) cout << "done." << endl; }
		}

		//Make sure the user isn't using -mfe and -s, these are incompatible.

		if (quickfold&&saveFile!="") {

			error = 1;
			cerr << "Fold stopped.  The -s and -mfe commands are incompatible.\n";

		}

		/*
		 * Fold the single strand using the FoldSingleStrand method.
		 * During calculation, monitor progress using the TProgressDialog class and the Start/StopProgress methods of the RNA class.
		 * Neither of these methods require any error checking.
		 * After the main calculation is complete, use the error checker's isErrorStatus method to check for errors.
		 */
		if( error == 0 ) {
			// Show a message saying that the main calculation has started.
			if (!quiet) cout << "Folding single strand..." << flush;

			// Create the progress monitor.
			TProgressDialog* progress = new TProgressDialog();
			if (!quiet) strand->SetProgress( *progress );

			// Do the main calculation and check for errors.
			int mainCalcError = strand->FoldSingleStrand( percent, maxStructures, windowSize, saveFile.c_str(), maxLoop, quickfold, simple_iloops, disablecoax, isolated);
			error = checker->isErrorStatus( mainCalcError );

			// Delete the progress monitor.
			strand->StopProgress();
			delete progress;

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

			// Write the CT file and check for errors.
			sprintf(numstr, "%d", b_iter);
			int writeError;

			if( b_iter > 0 ){
				writeError = strand->WriteCt( (ctFile + "_boot" + numstr).c_str() ); 
			} else { 
				if (useBracketNotation) {
					const DotBracketFormat bracketFmt = DBN_FMT_MULTI_TITLE; //see DotBracketFormat.h
					// write DotBracket file
					writeError = strand->WriteDotBracket(ctFile.c_str(), -1/*write ALL*/, bracketFmt); 
				} else
					// write CT file
					writeError = strand->WriteCt( ctFile.c_str() ); 
			}

			error = checker->isErrorStatus( writeError );

			// If no errors occurred, show a CT file writing completion message.
			if( error == 0 ) { if (!quiet) cout << "done." << endl; }
		}
	}

	// Delete the error checker and data structure.
	delete checker;
	delete strand;

	// Print confirmation of run finishing.
	if( error == 0 ) { if (!quiet) cout << calcType << " complete." << endl; }
	else { cerr << calcType << " complete with errors." << endl; }

	return error == 0;
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {
	Fold_Interface runner;
	if( !runner.parse( argc, argv ) ) return 1;
	return runner.run() ? 0 : 1;
}

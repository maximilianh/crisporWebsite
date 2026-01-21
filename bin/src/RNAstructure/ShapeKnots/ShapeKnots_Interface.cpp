/*
 * ShapeKnots, a program that predicts RNA secondary structures with pseudoknots.
 *
 * (c) 2013
 * Mathews Lab, University of Rochester Medical Center
 * Weeks Lab, The University at North Carolina at Chapel Hill
 * Code contributors: Wayne Higgins, Stanislav Bellaousov, David H. Mathews
 */

#include "ShapeKnots_Interface.h"
#include "../src/ErrorChecker.h"
#include "../RNA_class/thermodynamics.h"
#define OUTPUT_TO_SCREEN
//#undef OUTPUT_TO_SCREEN

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
ShapeKnots_Interface::ShapeKnots_Interface() {

	//Initialize the differntial SHAPE slope.
	Dslope = 2.11;

	// Initialize the maximum number of possible pseudoknotted helices (used in function convertToHunderHelices())
	finallistSize=100;

	// Initialize the maximum number of structures to be generated internally.
	InMaxStructures=100;

	// Initialize the maximum percent energy difference for the generated structures.
	InPercent=20;

	// Initialize the SHAPE intercept (kcal/mol).
	intercept=-0.6;

	// Initialize the folding window size or how different internal suboptimal structures can be.
	InWindowSize=0;

	// Initialize the maximum number of structures to be outputted.
	OutMaxStructures=20;

	// Initialize the maximum percent energy difference for the outputted suboptimal structures.
	OutPercent=10;

	// Initialize the folding window size or hwo different internal suboptimal structures can be.
	// This value will be changed based on the length of the sequence.
	OutWindowSize=0;

	// Initialize pseudoknot energy model parameters (kcal/mol).
	P1=0.35;

	// Initialize pseudoknot energy model parameters (kcal/mol).
	P2=0.65;

	// Initialize the SHAPE slope (kcal/mol).
	slope=1.8;

	//////////
    // Initialize the "experimental" offset.
	experimentalOffset = 0.0;
    	// Initialize the "experimental" scaling.
	experimentalScaling = 1.0;
	//////////
}

bool ShapeKnots_Interface::Parse(int argc, char** argv){

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "ShapeKnots" );
	parser->addParameterDescription( "seq file", "The name of a sequence file containing input data. Note that lowercase nucleotides are forced single-stranded in structure prediction." );
	parser->addParameterDescription( "ct file", "The name of a CT file to which output will be written." );

    // Add the constraint file option.
	vector<string> constraintOptions;
	constraintOptions.push_back( "-c" );
	constraintOptions.push_back( "-C" );
	constraintOptions.push_back( "--constraint" );
	parser->addOptionFlagsWithParameters( constraintOptions, "Specify a constraints file to be applied. Default is to have no constraints applied." );

    // Add the DMS option.
	vector<string> dmsOptions;
	dmsOptions.push_back( "-dms" );
	dmsOptions.push_back( "-DMS" );
	dmsOptions.push_back( "--DMS" );
	parser->addOptionFlagsWithParameters( dmsOptions, "Specify a DMS restraints file to be applied. These restraints are pseudoenergy restraints. Default is to have no restraints applied." );
    
    // Add the DMS option.
	vector<string> dmsntOptions;
	dmsntOptions.push_back( "-dmsnt" );
	dmsntOptions.push_back( "-DMSNT" );
    dmsntOptions.push_back( "-DMSnt" );
	dmsntOptions.push_back( "--DMSNT" );
	parser->addOptionFlagsWithParameters( dmsntOptions, "Specify a DMS NT restraint file to be applied. These restraints are pseudoenergy restraints. Default is to have no restraints applied." );


	// Add the differential SHAPE restraint file option.
	vector<string> DshapeOptions;
	DshapeOptions.push_back( "-dsh" );
	DshapeOptions.push_back( "-DSH" );
	DshapeOptions.push_back( "--DSHAPE" );
	parser->addOptionFlagsWithParameters(DshapeOptions, "Specify a differential SHAPE restraints file to be applied. These constraints are pseudoenergy restraints. Default is to have no restraints applied." );

	// Add the differential SHAPEslope option.
	vector<string> DshapeSlopeOptions;
	DshapeSlopeOptions.push_back( "-dsm" );
	DshapeSlopeOptions.push_back( "-DSM" );
	DshapeSlopeOptions.push_back( "--DSHAPEslope" );
	parser->addOptionFlagsWithParameters( DshapeSlopeOptions, "Specify a slope used with differential SHAPE restraints. Default is 2.11 kcal/mol." );

	// Add the double offset file option.
	vector<string> doubleOffsetOptions;
	doubleOffsetOptions.push_back( "-dso" );
	doubleOffsetOptions.push_back( "-DSO" );
	doubleOffsetOptions.push_back( "--doubleOffset" );
	parser->addOptionFlagsWithParameters( doubleOffsetOptions, "Specify a double-stranded offset file, which adds energy bonuses to particular double-stranded nucleotides. Default is to have no double-stranded offset specified." );

	// Add the maximum number of structures option.
	vector<string> INmaxStructuresOptions;
	INmaxStructuresOptions.push_back( "-im" );
	INmaxStructuresOptions.push_back( "-IM" );
	INmaxStructuresOptions.push_back( "--IMaximum" );
	parser->addOptionFlagsWithParameters( INmaxStructuresOptions, "Specify a maximum number of internally generated structures for each call of the dynamic programming algorithm. Note that suboptimal structures are generated until either the maximum number of structures is reached or the maximum percent difference is reached (below).  This is not the maximum number of structures provided to the user, which is controlled by –m, -M, --maximum. Default is 100 structures." );

	// Add the option to change the output behavior of non-critical warnings.
	vector<string> warningsOption;
	warningsOption.push_back( "--warnings" );
	warningsOption.push_back( "--warn" );
	parser->addOptionFlagsWithParameters( warningsOption, "Set the behavior for non-critical warnings (e.g. related to invalid nucleotide positions or duplicate data points in SHAPE data). Valid values are: "
	"\n\t* on  -- Warnings are written to standard output. (default)"
	"\n\t* err -- Warnings are sent to STDERR. This can be used in automated scripts etc to detect problems."
	"\n\t* off -- Do not display warnings at all (not recommended)."
	);

	// Add the percent energy difference option.
	vector<string> INpercentOptions;
	INpercentOptions.push_back( "-ip" );
	INpercentOptions.push_back( "-IP" );
	INpercentOptions.push_back( "--IPercent" );
	parser->addOptionFlagsWithParameters( INpercentOptions, "Specify a maximum percent difference in folding free energy change for internally generated suboptimal structures for each call of the dynamic programming algorithm. For example, 20 would indicate 20%. This is not the maximum percent difference in energy for structures provided to the user, which is controlled by –p, -P, --percent. Default is 20%." );

	// Add the window size option.
	vector<string> INwindowOptions;
	INwindowOptions.push_back( "-iw" );
	INwindowOptions.push_back( "-IW" );
	INwindowOptions.push_back( "--IWindow" );
	parser->addOptionFlagsWithParameters( INwindowOptions, "Specify a window size for the internally generated suboptimal structures for each call of the dynamic programming algorithm.  This is not the window for structures provided to the user, which is controlled by –w, -W, --window. Default is determined by the length of the sequence." );

	// Add the maximum number of structures option.
	vector<string> OUTmaxStructuresOptions;
	OUTmaxStructuresOptions.push_back( "-m" );
	OUTmaxStructuresOptions.push_back( "-M" );
	OUTmaxStructuresOptions.push_back( "--maximum" );
	parser->addOptionFlagsWithParameters( OUTmaxStructuresOptions, "Specify a maximum number of structures to be outputted. Note that suboptimal structures are generated until either the maximum number of structures is reached or the maximum percent difference is reached (below). Default is 20 structures." );

	// Add the percent energy difference option.
	vector<string> OUTpercentOptions;
	OUTpercentOptions.push_back( "-p" );
	OUTpercentOptions.push_back( "-P" );
	OUTpercentOptions.push_back( "--percent" );
	parser->addOptionFlagsWithParameters( OUTpercentOptions, "Specify a maximum percent difference in folding free energy change for generating suboptimal structures in the output. For example, 10 would indicate 10%. Default is 10%." );

	// Add the Penalty1 option.
	vector<string> Penalty1Options;
	Penalty1Options.push_back( "-p1" );
	Penalty1Options.push_back( "-P1" );
	Penalty1Options.push_back( "--Penalty1" );
	parser->addOptionFlagsWithParameters( Penalty1Options, "Specify a pseudoknot penalty P1. Default is 0.35 kcal/mol." );

	// Add the Penalty2 option.
	vector<string> Penalty2Options;
	Penalty2Options.push_back( "-p2" );
	Penalty2Options.push_back( "-P2" );
	Penalty2Options.push_back( "--Penalty2" );
	parser->addOptionFlagsWithParameters( Penalty2Options, "Specify a pseudoknot penalty P2. Default is 0.65 kcal/mol." );

	// Add the finallistSize option.
	vector<string> finallistSizeOptions;
	finallistSizeOptions.push_back( "-ph" );
	finallistSizeOptions.push_back( "-PH" );
	finallistSizeOptions.push_back( "--PseudoknottedHelices" );
	parser->addOptionFlagsWithParameters( finallistSizeOptions, "Specify maximum number of helices to be processed. Default is 100 helices." );

	// Add the SHAPE restraint file option.
	vector<string> shapeOptions;
	shapeOptions.push_back( "-sh" );
	shapeOptions.push_back( "-SH" );
	shapeOptions.push_back( "--SHAPE" );
	parser->addOptionFlagsWithParameters(shapeOptions, "Specify a SHAPE restraints file to be applied. These restraints specifically use SHAPE pseudoenergy restraints. Default is no SHAPE restraint file specified." );

	// Add the SHAPEintercept option.
	vector<string> shapeInterceptOptions;
	shapeInterceptOptions.push_back( "-si" );
	shapeInterceptOptions.push_back( "-SI" );
	shapeInterceptOptions.push_back( "--SHAPEintercept" );
	parser->addOptionFlagsWithParameters( shapeInterceptOptions, "Specify an intercept used with SHAPE restraints. Default is -0.6 kcal/mol." );

	// Add the SHAPEslope option.
	vector<string> shapeSlopeOptions;
	shapeSlopeOptions.push_back( "-sm" );
	shapeSlopeOptions.push_back( "-SM" );
	shapeSlopeOptions.push_back( "--SHAPEslope" );
	parser->addOptionFlagsWithParameters( shapeSlopeOptions, "Specify an slope used with SHAPE restraints. Default is 1.8 kcal/mol." );

	// Add the single offset file option.
	vector<string> singleOffsetOptions;
	singleOffsetOptions.push_back( "-sso" );
	singleOffsetOptions.push_back( "-SSO" );
	singleOffsetOptions.push_back( "--singleOffset" );
	parser->addOptionFlagsWithParameters( singleOffsetOptions, "Specify a single-stranded offset file, which adds energy bonuses to particular single-stranded nucleotides. Default is to have no single-stranded offset specified." );

	// Add the window size option.
	vector<string> OUTwindowOptions;
	OUTwindowOptions.push_back( "-w" );
	OUTwindowOptions.push_back( "-W" );
	OUTwindowOptions.push_back( "--window" );
	parser->addOptionFlagsWithParameters( OUTwindowOptions, "Specify a window size for outputted suboptimal structures. Default is determined by the length of the sequence." );


	/////////
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
	/////////




	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		seqFile = parser->getParameter( 1 );
		ctFile = parser->getParameter( 2 );
	}

	// Get the constraint file option.
	if( !parser->isError() ) { constraintFile = parser->getOptionString( constraintOptions, true ); }

	// Get the differential SHAPE slope option.
	if( !parser->isError() ) {
		parser->setOptionDouble( DshapeSlopeOptions, Dslope );
	}

	// Get the double strand offset file option.
	if( !parser->isError() ) { doubleOffsetFile = parser->getOptionString( doubleOffsetOptions, true ); }

	// Get the maximum number of structures to be internally generated.
	if( !parser->isError() ) {
		parser->setOptionInteger( INmaxStructuresOptions, InMaxStructures );
		if( OutMaxStructures <= 0 ) { parser->setError( "maximum number of structures" ); }
	}

	// Get the percent energy difference option for internally generated structures.
	if( !parser->isError() ) {
		parser->setOptionInteger( INpercentOptions, InPercent );
		if( OutPercent < 0 ) { parser->setError( "percent energy difference" ); }
	}

	// Get the window size option for internally generated structures.
	if( !parser->isError() ) {
		parser->setOptionInteger( INwindowOptions, InWindowSize );
		if( OutWindowSize < 0 ) { parser->setError( "window size" ); }
	}

	// Get the maximum number of structures to be outputted option.
	if( !parser->isError() ) {
		parser->setOptionInteger( OUTmaxStructuresOptions, OutMaxStructures );
		if( OutMaxStructures <= 0 ) { parser->setError( "maximum number of structures" ); }
	}

	// Determine the behavior of ReadSHAPE when encountering multiple 
	// datapoints for the same nucleotide
	if(!parser->isError() && parser->contains(warningsOption)) {
		string value = parser->getOptionString(warningsOption, false/*not a file*/);
		toUpper(value);
		structure::ShowWarnings = value=="ON"?1:value=="OFF"?0:value=="ERR"?2:255;
		if (structure::ShowWarnings==255) parser->setError("warnings");
	}

	////////
    // Get the experimental bonus data and options.
	if( !parser->isError() ) {
		experimentalFile = parser->getOptionString( experimentalOptions );
		if( !parser->isError() ) { parser->setOptionDouble( experimentalOffsetOptions, experimentalOffset ); }
		if( !parser->isError() ) { parser->setOptionDouble( experimentalScalingOptions, experimentalScaling ); }
	}
	////////

	// Get the SHAPE, and DMS data and options.
	if( !parser->isError() ) {
		SHAPEFile = parser->getOptionString( shapeOptions );
		DSHAPEFile = parser->getOptionString( DshapeOptions );
        DMSFile = parser->getOptionString( dmsOptions );
        DMSNTFile = parser->getOptionString( dmsntOptions );
        if( !parser->isError() ) { parser->setOptionDouble( shapeInterceptOptions, intercept ); }
        if( !parser->isError() ) { parser->setOptionDouble( shapeSlopeOptions, slope ); }
		if( !parser->isError() ) { parser->setOptionDouble( DshapeSlopeOptions, Dslope ); }
	}

	// Get the percent energy difference option for outputted suboptimal structures.
    if( !parser->isError() ) {
        parser->setOptionInteger( OUTpercentOptions, OutPercent );
        if( OutPercent < 0 ) { parser->setError( "percent energy difference" ); }
	}

	// Get the pseudoknot penalty 1.
	if( !parser->isError() ) {
		parser->setOptionDouble( Penalty1Options, P1 );
	}

	// Get the pseudoknot penalty 2.
	if( !parser->isError() ) {
		parser->setOptionDouble( Penalty2Options, P2 );
	}

	// Get the finallistSize option.
	if( !parser->isError() ) {
		parser->setOptionInteger( finallistSizeOptions, finallistSize );
		if( finallistSize < 1 ) { parser->setError( "maximum number of possible pseudoknotted helices" ); }
	}

	// Get the single strand offset file option.
	if( !parser->isError() ) { singleOffsetFile = parser->getOptionString( singleOffsetOptions, true ); }

	// Get the window size option for outputted structures.
	if( !parser->isError() ) {
		parser->setOptionInteger( OUTwindowOptions, OutWindowSize );
		if( OutWindowSize < 0 ) { parser->setError( "window size" ); }
	}
	//Record if the windowOption was set
	ifWindowOptions=parser->contains(OUTwindowOptions);

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}


///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
bool ShapeKnots_Interface::run() {

#ifdef OUTPUT_TO_SCREEN
	cout << "############################################\n";
	cout << "         ShapeKnots Folding Started\n";
	cout << "############################################\n";
	cout << "Reading Files..." << flush ;
#endif

	RNA rna(seqFile.c_str(),FILE_SEQ,true);
	ErrorChecker<RNA> checker(&rna);
	if (checker.isErrorStatus()!=0) return false; // if an error has occurred, show a message and then exit.
    
    if (!constraintFile.empty()) {
        rna.ReadConstraints(constraintFile.c_str());
    };

	//If the user has specified a SHAPE restraints file, read the file.
	if(!SHAPEFile.empty()) rna.ReadSHAPE(SHAPEFile.c_str(),slope,intercept,RESTRAINT_SHAPE);

	//If the user has specified a DSHAPE restraints file, read the file.
	//Read in the shape data file and convert the data to a linear penalty
	if(!DSHAPEFile.empty()) rna.ReadSHAPE(SHAPEFile.c_str(),Dslope,0,RESTRAINT_SHAPE_DIFF /*diffSHAPE*/);

	//If the user has specified a DMS restraints file, read the file.
	//Read in the dms data file and convert the data to a linear penalty
	if(!DMSFile.empty()) rna.ReadDMS(DMSFile.c_str(), false);
    
    if(!DMSNTFile.empty()) {
        rna.ReadDMS(DMSNTFile.c_str(), true);
        //after DMS restraints are read in, DMSFile is used as a flag (empty or not empty) in pseudoknot 
        // to indicate that restraints have been applied. Set DMSFile=DMSNTFile to satisfy this flag
        DMSFile = DMSNTFile;     
    }

	// Read experimental pair bonus constraints, if applicable.
	if(!experimentalFile.empty()) rna.ReadExperimentalPairBonus( experimentalFile.c_str(), experimentalOffset, experimentalScaling );
	
	//Read the Double-strand offset data into 'ct'
	if(!doubleOffsetFile.empty()) rna.ReadDSO(doubleOffsetFile.c_str());//read the offset data

	//Read the Single-strand offset data into 'ct'
	if(!singleOffsetFile.empty()) rna.ReadSSO(singleOffsetFile.c_str());//read the offset data

	if (checker.isErrorStatus()!=0) return false; // if an error has occurred, show a message and then exit.
	
	//The rest of the code is located in the 'pseudoknot' function.
	pseudoknot(&rna, rna.GetDatatable(), InMaxStructures, InPercent, InWindowSize, ctFile, P1, P2, slope, intercept, DMSFile, SHAPEFile, Dslope, DSHAPEFile, doubleOffsetFile, OutPercent, OutWindowSize, OutMaxStructures, ifWindowOptions, finallistSize);

	return checker.isErrorStatus()==0;
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {
	ShapeKnots_Interface runner;
	if( !runner.Parse( argc, argv ) ) return 1;
	return runner.run() ? 0 : 1;
}

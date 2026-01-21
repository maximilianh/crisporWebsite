/*
 * A program that calculates the partition function for a strand of nucleic acids.
 * This strand of nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#include "StructureProb.h"
#include "../src/boltzmann.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
//////////////////////////////////////////////////////////////////////////////
StructureProbInterface::StructureProbInterface() {

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

	// Set empty string as input pfs file
	pfsFile = "";
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool StructureProbInterface::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "StructureProb" );
	parser->addParameterDescription( "ct file", "The name of a file containing an input structure." );
	//parser->addParameterDescription( "pfs file", "The name of a partition function save file to which output will be written." );

	// Add the constraint file option.
	vector<string> constraintOptions;
	constraintOptions.push_back( "-c" );
	constraintOptions.push_back( "-C" );
	constraintOptions.push_back( "--constraint" );
	parser->addOptionFlagsWithParameters( constraintOptions, "Specify a constraints file to be applied. Default is to have no constraints applied." );

	vector<string> quietOptions;
	quietOptions.push_back( "-q" );
	quietOptions.push_back( "--quiet" );
	parser->addOptionFlagsNoParameters( quietOptions, "Suppress progress display and other unnecessary output.");

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

	vector<string> inputPFSfileOptions;
	inputPFSfileOptions.push_back( "-pi" );
	parser->addOptionFlagsWithParameters( inputPFSfileOptions, "Specify an input partition save file.  If none is provided, the partition function will be calculated." );

	// Add the experimental pair bonus scaling option.
	vector<string> experimentalScalingOptions;
	experimentalScalingOptions.push_back( "-xs" );
	parser->addOptionFlagsWithParameters( experimentalScalingOptions, "Specify a number to multiply the experimental pair bonus matrix by. Default is 1.0 (no change to input bonuses)." );

	// Add the disablecoax option.
	vector<string> disablecoaxOptions;
	disablecoaxOptions.push_back( "--disablecoax" );
	parser->addOptionFlagsNoParameters( disablecoaxOptions, "Specify that coaxial stacking recusions should not be used. This option uses a less realistic energy function in exchange for a faster calculation." );

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		ctFile = parser->getParameter( 1 );
//		pfsFile = parser->getParameter( 2 );
	}

	quiet = parser->contains(quietOptions);

	// Get the constraint file option.
	if( !parser->isError() ) { constraintFile = parser->getOptionString( constraintOptions, true ); }

	// Get the DNA option.
	if( !parser->isError() && parser->contains( dnaOptions ) )
		alphabet = DT_DNA; // use DNA (unless overridden by alphabet)

	// Get the DNA option.
	// if( !parser->isError() ) { isRNA = !parser->contains( dnaOptions ); }
	

	// Get the Alphabet option.
	if (!parser->isError() && parser->contains(alphabetOptions))
		  alphabet = parser->getOptionString(alphabetOptions, false).c_str();

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

	if( !parser->isError() ) { pfsFile = parser->getOptionString( inputPFSfileOptions, true ); }

	// Get the disablecoax option.
    // true by default, false if the flag is used
	if( !parser->isError() ) { disablecoax = parser->contains( disablecoaxOptions ); }

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void StructureProbInterface::run() {

	// Create a variable that handles errors.
	int error = 0;
	ostream &info = quiet ? NullStream::Default : cout; // write output to a Null output stream if quiet==true.

	/*
	 * Use the constructor for RNA that specifies a filename.
	 * Specify type = 2 (sequence file).
	 * isRNA identifies whether the strand is RNA (true) or DNA (false).
	 *
	 * After construction of the strand data structure, create the error checker which monitors for errors.  
	 * Throughout, the error status of the calculation is checked with a variant of the isErrorStatus method, which returns 0 if no error occurred.
	 * The calculation proceeds as long as error = 0.
	 */
	info << "Initializing nucleic acids..." << flush;
	RNA* strand;
	RNA* strand2 = new RNA( ctFile.c_str(), FILE_CT, alphabet.c_str() );
//	RNA* strand = new RNA( seqFile.c_str(), FILE_SEQ, alphabet.c_str() );
	ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand2 );
	error = checker->isErrorStatus();
	if( error == 0 ) { info << "done." << endl; }

	if (pfsFile == ""){
		strand = new RNA( ctFile.c_str(), FILE_CT, alphabet.c_str() );
		/*
		* Set the temperature using the SetTemperature method.
		* Only set the temperature if a given temperature doesn't equal the default.
		* If the temperature does need to be set, use the error checker's isErrorStatus method to check for errors.
		*/
		if( ( error == 0 ) && ( temperature != 310.15 ) ) {

			// Show a message saying that the temperature is being set.
			info << "Setting temperature..." << flush;

			// Set the temperature and check for errors.
			int tempError = strand->SetTemperature( temperature );
			error = checker->isErrorStatus( tempError );

			// If no error occurred, print a message saying that temperature is set.
			if( error == 0 ) { info << "done." << endl; }
		}

		/*
		* Set maximum pairing distance using the ForceMaximumPairingDistance method.
		*/
		if( error == 0 && maxDistance != -1 ) {

			// Show a message saying that the maximum pairing distance is being set.
			info << "Setting maximum distance between paired nucleotides..." << flush;

			// Set the maximum pairing distance and check for errors.
			int distError = strand->ForceMaximumPairingDistance( maxDistance );
			error = checker->isErrorStatus( distError );

			// If no error occurred, print a message saying that maximum pairing distance was set.
			if( error == 0 ) { info << "done." << endl; }
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
			info << "Applying constraints..." << flush;
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
			if( error == 0 ) { info << "done." << endl; }
		}

		/*
		* Calculate the partition function using the PartitionFunction method.
		* During calculation, monitor progress using the TProgressDialog class and the Start/StopProgress methods of the RNA class.
		* Neither of these methods require any error checking.
		* After the main calculation is complete, use the error checker's isErrorStatus method to check for errors.
		*/
		if( error == 0 ) {

			// Show a message saying that the main calculation has started.
			info << "Calculating partition function..." << flush;

			// Create the progress monitor.
			// TProgressDialog* progress = new TProgressDialog();
			// if (!quiet) strand->SetProgress( *progress );

			// Do the main calculation and check for errors.
			// pfsFile should be an empty string at this point.  It should result in no save file being written.
			int mainCalcError = strand->PartitionFunction( pfsFile.c_str(), temperature, disablecoax);
			error = checker->isErrorStatus( mainCalcError );

			// Delete the progress monitor.
			// strand->StopProgress();
			// delete progress;

			// If no error occurred, print a message saying that the main calculation is done.
			if( error == 0 ) { info << "done." << endl; }

			// Print confirmation of run finishing.
			if( error == 0 ) { info << calcType << " complete." << endl; }
			else { cerr << calcType << " complete with errors." << endl; }
		}
	}
	else{
		strand = new RNA( pfsFile.c_str(), 3, alphabet.c_str() );
	}



	// Get the secondary structure
	//structure *ct = strand2->GetStructure();

	//Detect pseudoknots
	const int structureCount = strand2->GetStructureNumber();
	bool anyPK = false; // whether any structure has a pseudoknot
	vector<bool> hasPK(structureCount+1); // store whether each structure has pseudoknts or not, so we don't have to keep calculating this.

	for(int i=1; i<=structureCount; i++)
		if (strand2->ContainsPseudoknot(i))
			hasPK[i]=anyPK=true;

	// Calculate the span of the secondary structure
	int min_pair;
	int max_pair;
	int pair;
	double free_energy, prob;
	bool simple = true;

	for(int i=1; i<=structureCount; i++){
		min_pair = strand2->GetSequenceLength();
		max_pair = 0;
		// Iterate through the bases in the structure, looking for base pairs that expand the current span of the structure
		for(int j=1; j<= strand2->GetSequenceLength(); j++){
			pair = strand2->GetPair(j,i);
			if (pair > 0){
				if (min_pair > j) min_pair = j;
				if (max_pair < j) max_pair = j;
			}
		}
	
		// Calculate the free energy of the structure 
		// Still need to correct for exterior loop
		cout << "Span\t" << min_pair << "\t" << max_pair << endl;
		
		cout << "Free Energy\t" << strand2->CalculateFreeEnergy( i , simple ) << endl;
		cout << "External Correction\t" << strand2->ExteriorLoopCorrection(i, simple, min_pair, max_pair) << endl;

		free_energy = strand2->CalculateFreeEnergy( i , simple )+strand2->ExteriorLoopCorrection(i, simple, min_pair, max_pair);

		cout << "K\t" << boltzman(free_energy*conversionfactor, temperature) << endl;
//		cout << "Vprime/Q \t" << strand->GetVprimeQ(min_pair, max_pair) << endl;

		prob = PROD(boltzman(free_energy*conversionfactor, temperature),strand->GetVprimeQ(min_pair, max_pair));
		info << "Structure #" << i << ":  Probability: " << TO_LINEAR(prob) << endl;

		cout << "W\t" << strand->GetW(min_pair, max_pair) << endl;

		// cout << TO_LINEAR(DIV(boltzman(free_energy, temperature), strand->GetW(min_pair, max_pair))) << endl;
	}

	// Delete the error checker and data structure.
	delete checker;
	delete strand;
	delete strand2;

}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	StructureProbInterface* runner = new StructureProbInterface();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}

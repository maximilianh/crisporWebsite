/*
 * A program that calculates the free energy of a strand of nucleic acids.
 * This strand of nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2008 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */
#include "efn2.h"
#include "../src/PseudoParser.h"
#include <math.h>

/*
 * Overload the addition operator for vectors  This doesn't check if the vectors are the
 * same size.  This is used to accumulate the parameter usage counts so they are not lost 
 * when the counts are reset after each energy calculation.
 */
vector<double> operator+(const std::vector<double> &a, const std::vector<double> &x){
	// Initialize the output vector
	vector<double> sum;

	// Fill the output vector by performing pairwise addition
	for(int i = 0; i < a.size(); i++){
		sum.push_back(a[i]+x[i]);
	}
	
	// Return the output vector
  	return sum;
}

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
efn2Interface::efn2Interface() {

	// Initialize the calculation type description.
	calcType = "efn2";

	// Initialize the nucleic acid type to "rna"
	alphabet = DT_RNA;

	// Initialize the SHAPE intercept.
	intercept = DEFAULT_SHAPE_INTERCEPT;

	// Initialize the SHAPE slope.
	slope = DEFAULT_SHAPE_SLOPE;

	// Initialize flag that signifies streaming to standard output.
	stdPrint = false;

	// Initialize the calculation temperature.
	temperature = 310.15;

	// Initialize the flag that signifies writing a thermodynamic details file.
	writeTherm = false;

	// Initialize the simple energy function.
	simple = false;

	quiet = false;

	// parameters for calculating the energy of pseudoknot-containing structures
	pseudoP1=DEFAULT_PSEUDOKNOT_P1;
	pseudoP2=DEFAULT_PSEUDOKNOT_P2;
	pseudo_params[0]=0; // indicates params have not been read.
	pseudoFillMismatch=false;
	pseudoRemIsolatedPairs=false;
	omitErrors = false;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool efn2Interface::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "efn2" );
	parser->addParameterDescription( "ct file", "The name of a file containing structure CT data." );
	parser->addParameterDescription( "output file", "The energy file to which output is written. Depending on the options selected, this may be one of the following file types. 1) Simple list (Lists free energy for each structure, lowest first). 2) Thermodynamic details (Writes details of every substructure in each structure, and the corresponding free energy of each)." );

	// Add the DNA option.
	vector<string> dnaOptions;
	dnaOptions.push_back( "-d" );
	dnaOptions.push_back( "-D" );
	dnaOptions.push_back( "--DNA" );
	parser->addOptionFlagsNoParameters( dnaOptions, "Specify that the sequence is DNA, and DNA parameters are to be used. Default is to use RNA parameters." );

	// Add the DNA option.
	vector<string> alphabetOptions;
	alphabetOptions.push_back( "-a" );
	alphabetOptions.push_back( "--alphabet" );
	parser->addOptionFlagsWithParameters( alphabetOptions, "Specify the name of a folding alphabet and associated nearest neighbor parameters. The alphabet is the prefix for the thermodynamic parameter files, e.g. \"rna\" for RNA parameters or \"dna\" for DNA parameters or a custom extended/modified alphabet. The thermodynamic parameters need to reside in the at the location indicated by environment variable DATAPATH. The default is \"rna\" (i.e. use RNA parameters). This option overrides the --DNA flag." );

	// Add the simple eneergy function option.
	vector<string> simpleOptions;
	simpleOptions.push_back( "-s" );
	simpleOptions.push_back( "-S" );
	simpleOptions.push_back( "--simple" );
	parser->addOptionFlagsNoParameters( simpleOptions, "Specify the simple energy function for multibranch loops, used by the dynamic programming algorithms (Fold, partition, stochastic, AllSub, etc.), should be used. If this is not specified, an more sophisticated energy function is used, and the energies might not match those estimated for structures during structure prediction." );

	// Add the print option.
	vector<string> printOptions;
	printOptions.push_back( "-p" );
	printOptions.push_back( "-P" );
	printOptions.push_back( "--print" );
	parser->addOptionFlagsNoParameters( printOptions, "Print the simple list file to standard output. This won't override default behavior of writing to a file. Thermodynamic files (if written) are not piped. This option implies --quiet." );

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

	// Add the details option.
	vector<string> detailsOptions;
	detailsOptions.push_back( "-w" );
	detailsOptions.push_back( "-W" );
	detailsOptions.push_back( "--writedetails" );
	parser->addOptionFlagsNoParameters( detailsOptions, "Write a thermodynamic details file. The thermodynamic details file replaces the list file that is outputted by default." );

	// Add the option for a data count output file.  The file will contain a linear list of the total parameter 
	// usage counts used to perform the energy calculations for an entire ct file.
	vector<string> countOptions;
	countOptions.push_back( "-c" );
	countOptions.push_back( "-C" );
	countOptions.push_back( "--count" );
	parser->addOptionFlagsWithParameters( countOptions, "Specify a file where parameter usage counts will be exported" );

	// Add the no error option.
	vector<string> errorOptions;
	errorOptions.push_back( "--ne" );
	parser->addOptionFlagsNoParameters( errorOptions, "Do not calculate experimental uncertainties" );

	vector<string> quietOption = parser->addFlag(false, "-q --quiet", "Suppress unnecessary output. This option is implied when the output file is '-' (STDOUT) or when the --print flag is present.");

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		ctFile = parser->getParameter( 1 );
		outFile = parser->getParameter( 2 );
	}

	// Get the DNA option.
	if( !parser->isError() && parser->contains( dnaOptions ) )
		alphabet = DT_DNA; // use DNA (unless overridden by alphabet)

	// Get the Alphabet option.
	if (!parser->isError() && parser->contains(alphabetOptions))
		  alphabet = parser->getOptionString(alphabetOptions, false).c_str();

	// Get the simple energy rule option.
	if (!parser->isError() ) { simple = parser->contains( simpleOptions); }

	// Get the print option.
	if( !parser->isError() ) { stdPrint = parser->contains( printOptions ); }

	// if outFile == "-" then pritn to stdout and do not write a file.
	if (isStdIoFile(outFile.c_str())) {
		stdPrint = true;
		outFile.clear();
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

	// Get the write thermodynamic details option.
	if( !parser->isError() ) { writeTherm = parser->contains( detailsOptions ); }

	// Get the count output file option
	if( !parser->isError() ) {
		countFile = parser->getOptionString( countOptions , false);
	}

	quiet = parser->contains(quietOption) || stdPrint; // suppress unnecessary output if --quiet option is present or if printing results to STDOUT.

	if( !parser->isError() ) { omitErrors = parser->contains( errorOptions ); }

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void efn2Interface::run() {

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
	RNA* strand = new RNA( ctFile.c_str(), FILE_CT, alphabet.c_str() );
	ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );
	error = checker->isErrorStatus();
	if( error == 0 ) { if (!quiet) cout << "done." << endl; }

	// If there is no count output file specified, some of the following steps can be skipped
	bool exportCounts = (countFile != "" );
	bool calcErrors;
	std::vector<double> cumulative_Counts;

#ifdef COUNTING
	calcErrors = (alphabet == DT_RNA && !omitErrors);
#else
	calcErrors = false;
#endif
	

	

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
	 * Read SHAPE constraints, if applicable.
	 * When reading SHAPE constraints, use the ReadSHAPE method.
	 * After constraints are read, use the error checker's isErrorStatus method to check for errors.
	 */
	if( error == 0 && SHAPEFile != "" ) {

		// Show a message saying that SHAPE constraints are being read.
		if (!quiet) cout << "Applying SHAPE constraints..." << flush;

		// Initialize the single stranded SHAPE slope and intercept.
		// For now, these are hard-coded as 0.
		double slopeSingle = 0;
		double interceptSingle = 0;

		// Read SHAPE constraints and check for errors.
		int constraintError = strand->ReadSHAPE( SHAPEFile.c_str(), slope, intercept, slopeSingle, interceptSingle );
		error = checker->isErrorStatus( constraintError );

		// If no error occurred, print a message saying that SHAPE constraints are set.
		if( error == 0 ) { if (!quiet) cout << "done." << endl; }
	}

	const int structureCount = strand->GetStructureNumber();

	bool anyPK = false; // whether any structure has a pseudoknot
	vector<bool> hasPK(structureCount+1); // store whether each structure has pseudoknts or not, so we don't have to keep calculating this.

	#ifdef SMP
	#pragma omp parallel for
	#endif
	for(int i=1; i<=structureCount; i++)
		if (strand->ContainsPseudoknot(i))
			hasPK[i]=anyPK=true;
	
	if (error==0 && anyPK) // load the pseudoknot penalty parameters if any structure has a pseudoknot.
		error = checker->isErrorStatus(ReadPseudoParam(pseudo_params));

	/*
	 * Do the efn2 calculation.
	 * If the user wants a simple output file, get free energies for each structure using the CalculateFreeEnergy method.
	 * If the user wants a thermodynamic details file, write the file with the WriteThemodynamicDetails method.
	 */
	if( error == 0 ) {
		// Write a thermodynamic details file, if asked.
		if( writeTherm ) {
			// Show a message saying that the details file is being written.
			if (!quiet) cout << "Writing thermodynamic details file..." << flush;
			datatable *data = strand->GetDatatable();
			structure *ct = strand->GetStructure();
			ofstream out(outFile.c_str());
			if (out.good()) {
				for(int i=1;i<=structureCount;i++) {
					if (hasPK[i]) {
						CalculateFreeEnergywithPseudoknot(strand, i, pseudoP1, pseudoP2, pseudo_params, data, simple, pseudoFillMismatch, pseudoRemIsolatedPairs, &out);
						out << "Details follow for Broken (Pseudoknot-free) Structure." << endl;
						ct->BreakPseudoknots(i); // break so that efn2 doesn't enter an endless loop
					}
					efn2(data, ct, i, simple, &out);
				}
			} else
				error = 34; //failed to write output file.
			
			// Write the thermodynamic details file and check for errors.
			//int thermError = strand->WriteThermodynamicDetails( outFile.c_str(),simple );
			//error = checker->isErrorStatus( thermError );

			// Print a message saying that the details file has been written.
			if( error == 0 ) { if (!quiet) cout << "done." << endl; }
		}

		// Write a simple list file, if asked.
		else {

			// Show a message saying that the list file is being written.
			if (!quiet) cout << "Calculating free energies..." << flush;
			
			// Initialize the vector that will store the experimental uncertainty in the energy calculation for each structure
			std::vector<double> uncertainties;

			int i;
			      
			// For each structure, calculate its energy.
			#ifdef SMP
			#pragma omp parallel for
			#endif
			
			// Initialize the cumulative_Counts with the current parameter usage counts (they should all be zero)
#ifdef COUNTING			
			if (exportCounts)
				cumulative_Counts = strand->GetDataCounters();
#endif			

			for( i = 1; i <= structureCount; i++ ) { 
#ifdef COUNTING				
				// Reset the data counts (each structure gets a distinct uncertainty)
				strand->ResetDataCounters();
#endif
				if (hasPK[i])
					 CalculateFreeEnergywithPseudoknot(strand, i, pseudoP1, pseudoP2, pseudo_params, strand->GetDatatable(), simple, pseudoFillMismatch, pseudoRemIsolatedPairs);
				else
					strand->CalculateFreeEnergy( i , simple );

#ifdef COUNTING
				// Add the data counts for the structure to the cumulative totals
				if (exportCounts)
					cumulative_Counts = cumulative_Counts + strand->GetDataCounters();
				
				// Calculate the experimental uncertainty based on the parameter usage counts and add it to the vector
				if (calcErrors)
					uncertainties.push_back(strand->CalculateUncertainty());
#endif					
			}

			error = checker->isErrorStatus(strand->GetErrorCode()); // show a message if there was an error.
			
			// Print a message saying that the energies have been calculated.
			if( error == 0 ) { if (!quiet) cout << "done." << endl; }

			// If the output should be piped to standard output, then pipe it.
			if(error==0 && stdPrint) {
				for( i = 1; i <= structureCount; i++ ) {
					cout << "Structure: " << i << "   Energy = " << fixed << setprecision( conversionprecision ) << strand->GetFreeEnergy(i); //changed to GetFreeEnergy instead of writing the energies vector
					if (calcErrors)
						cout <<" \u00b1 "<< uncertainties[i-1]; 
					cout << endl;
				}
				if (!quiet) cout << endl;
			}

			// If all free energies were calculated correctly, write the output file.
			if( error == 0 && !outFile.empty() && outFile != ".") {
				//if (!quiet) cout << "Writing free energy list..." << flush;
				ofstream out( outFile.c_str() );
				if (out.good()) {
					for( i = 1; i <= structureCount; i++ ) {
						out << "Structure: " << i << "   Energy = " << fixed << setprecision( conversionprecision ) << strand->GetFreeEnergy(i); //changed to GetFreeEnergy instead of writing the energies vector
						if (calcErrors)
							out <<" \u00b1 "<< uncertainties[i-1]; 
						out << endl;
					}
				} else {
					cerr << "Error -- Failed to open output file " << outFile << "." << endl;
					error = 2;
				}
				out.close();
				//if (error==0) { if (!quiet) cout << "done." << endl ; }
			}
		}
	}

#ifdef COUNTING
	// Write out the cumulative parameter usage counts
	if( error == 0 && countFile != "" ){
		if (!quiet) cout << "Exporting Data Counts..." << flush;
		if (!strand->WriteDataCounters(countFile,cumulative_Counts)) {
			cerr << "\nError -- Failed to write Data Counts file " << countFile << "." << endl;
			error = 1;
		}
		if (error==0) { if (!quiet) cout << "done." << endl; }
	}
#endif	

	// Delete the error checker and data structure.
	
	delete checker;

	delete strand;

	// Print confirmation of run finishing.
	if( error == 0 ) { if (!quiet) cout << calcType << " complete." << endl; }
	else { cerr << calcType << " ended with errors." << endl; }
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	efn2Interface* runner = new efn2Interface();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}


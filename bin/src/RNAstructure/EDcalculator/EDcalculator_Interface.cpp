/*
 * EDcalculator: program that calculates Ensemble defect given a structure
 * Written by Richard M. Watson and Mohammad Kayedkhordeh in Dr. David Mathews Lab at the University of Rochester. 
 * This strand of nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2017 Mathews Lab, University of Rochester Medical Center.
 */

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include "EDcalculator_Interface.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {
	int error = 1;
	EDcalculator_Interface* runner = new EDcalculator_Interface();
	bool parseable = runner->parse( argc, argv );
	if( parseable ) { error = runner->run(); }
	delete runner;
	return error;
}

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
EDcalculator_Interface::EDcalculator_Interface() {
	// Initialize default parameters.

	// The input ct file which describes the structure to design. The sequence in the ct file will be ignored.
	ctFile = "";

	// Name of the alphabet (e.g. "rna" or "dna" etc).
	alphabet = DT_RNA;

	structurenumber = -1;
	
	//Set defaults for the local calculation (defaults are global).
	nuc_start = 0;
	nuc_end = 0;

	bool raw = false;

	constraintFile = "";
	nucfilename = "";

	//by default, do not allow isolated base pairs
	allow_isolated = false;

} 


///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool EDcalculator_Interface::parse( int argc, char** argv ) {
	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "EDCalculator" );

	parser->addParameterDescription( string("ct structure file"), string("The input ct file used to calculate the ensemble defect.") );

	vector<string> dnaOptions;
	dnaOptions.push_back( "-d" );
	dnaOptions.push_back( "--dna" );
	parser->addOptionFlagsNoParameters( dnaOptions, "Specify that the sequence is DNA, and DNA parameters are to be used. The default is to use RNA parameters." );

	vector<string> alphabetOptions;
	alphabetOptions.push_back( "-a" );
	alphabetOptions.push_back( "--alphabet" );
	parser->addOptionFlagsWithParameters( alphabetOptions, "Specify the name of a folding alphabet and associated nearest neighbor parameters. The alphabet is the prefix for the thermodynamic parameter files, e.g. \"rna\" for RNA parameters or \"dna\" for DNA parameters or a custom extended/modified alphabet. The thermodynamic parameters need to reside in the at the location indicated by environment variable DATAPATH. The default is \"rna\" (i.e. use RNA parameters). This option overrides the --DNA flag." );

	vector<string> NucFileOptions;
	NucFileOptions.push_back("--nucfile");
	parser->addOptionFlagsWithParameters(NucFileOptions, "Specify the name of a plain text file that will contain the defect per nucleotide.  The default is to specify no filename and to not output this file.  NOTE that the output is per nucleotide per structure and the data will be appended to an existing file.");


	// Add the number option.
	vector<string> numberOptions;
	numberOptions.push_back( "-n" );
	numberOptions.push_back( "--number" );
	parser->addOptionFlagsWithParameters( numberOptions, "Specify the index of a particular structure for which to calculate the defect. The default is -1, which means to calculate the defect for all structures." );

	// Add the raw option.
	vector<string> rawOptions;
	rawOptions.push_back( "-r" );
	rawOptions.push_back( "--raw" );
	parser->addOptionFlagsNoParameters( rawOptions, "Output just the *Normalized* ensemble defect as a pure number (with no additional description)." );

	// Add the output file option.
	vector<string> fileOption;
	fileOption.push_back( "-f" );
	fileOption.push_back( "--file" );
	parser->addOptionFlagsWithParameters( fileOption, "Output the results to the specified file instead of to the screen (stdout)." );

	// Add the constraint file option.
	vector<string> constraintOptions;
	constraintOptions.push_back( "-c" );
	constraintOptions.push_back( "-C" );
	constraintOptions.push_back( "--constraint" );
	parser->addOptionFlagsWithParameters( constraintOptions, "Specify a constraints file to be applied. Default is to have no constraints applied." );

	//Add the start and end parameters for local structure calculations
	vector<string> startOptions;
	startOptions.push_back("-s");
	startOptions.push_back("-S");
	startOptions.push_back("--start");
	parser->addOptionFlagsWithParameters(startOptions, "Specify a start nucleotide for a local calculation.");
	vector<string> endOptions;
	endOptions.push_back("-e");
	endOptions.push_back("-E");
	endOptions.push_back("--end");
	parser->addOptionFlagsWithParameters(endOptions, "Specify an end nucleotide for a local calculation.");

	vector<string> isolatedPairOptions;
	isolatedPairOptions.push_back("-i");
	isolatedPairOptions.push_back("--isolated");
	parser->addOptionFlagsNoParameters(isolatedPairOptions, "Specify that isolated pairs are allowed.  The default is to use a heuristic to attempt to forbid isolated pairs during structure prediction.");




	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) 
		ctFile = parser->getParameter( 1 );

	if (parser->contains( dnaOptions ))
	  alphabet = DT_DNA;

	if (!parser->isError())
		if (parser->contains(alphabetOptions))
		  alphabet = parser->getOptionString(alphabetOptions,false).c_str();

	if (!parser->isError())
		if (parser->contains(NucFileOptions))
			nucfilename = parser->getOptionString(NucFileOptions, false).c_str();

	// Get the constraint file option.
	if( !parser->isError() ) { constraintFile = parser->getOptionString( constraintOptions, true ); }

	// Get the number option.
	if( !parser->isError() ) {
		parser->setOptionInteger( numberOptions, structurenumber );
		if( structurenumber != -1  && structurenumber < 1 ) { parser->setError( "structure number" ); }
	}

	// Get the start option.
	if (!parser->isError()) {
		parser->setOptionInteger(startOptions, nuc_start);
		if (nuc_start < 0 ) { parser->setError("structure number"); }
	}

	// Get the number option.
	if (!parser->isError()) {
		parser->setOptionInteger(endOptions, nuc_end);
		
	}


	raw = parser->contains(rawOptions);

	if( !parser->isError() ) {
		outputFile = parser->getOptionString( fileOption,  false);
	}

	if (!parser->isError())
		if (parser->contains(isolatedPairOptions))
			allow_isolated = true;

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
int EDcalculator_Interface::run() {
	// Create a variable that handles errors.
	int errNum = 0;
	/*
	 * Use the constructor for RNA that specifies a filename.
	 * Specify type = 1 (CT file).
	 * isRNA identifies whether the strand is RNA (true) or DNA (false).
	 *
	 * After construction of the strand data structure, create the error checker which monitors for errors.  
	 * Throughout, the error status of the calculation is checked with a variant of the isErrorStatus method, which returns 0 if no error occurred.
	 * The calculation proceeds as long as error = 0.
	 */

	//Initialize class 'design'
	RNA rna(ctFile.c_str(), FILE_CT, alphabet.c_str());
	ErrorChecker<RNA> checker(&rna);
	errNum = checker.isErrorStatus();

	if (errNum) return errNum;

	int start = structurenumber;
	int end = structurenumber + 1;
	
	// Read folding constraints, if applicable.
	if( constraintFile != "" ) {
		errNum = checker.isErrorStatus( rna.ReadConstraints( constraintFile.c_str() ) );
		if (errNum) return errNum;
	}

	if (structurenumber == -1) { 
		start = 1; end = rna.GetStructureNumber() + 1; 
	}
	TProgressDialog* progress;
	if (!raw) {
		cout << "Calculating pair probabilities..." << endl;
		progress=new TProgressDialog();
		rna.SetProgress(*progress);
	}
	errNum = checker.isErrorStatus( rna.PartitionFunction("",-10.0,false,true,allow_isolated ));
	if (errNum) return errNum;
	if (!raw)  
		cout << endl << "done." << endl;

	if (!raw) {
		// Delete the progress monitor.
		rna.StopProgress();
		delete progress;
	}
	
	ofstream fout;
	if (!outputFile.empty()) {
		fout.open(outputFile.c_str());
		if (!fout.good())
			cerr << "Error opening output file " << outputFile << ". Results will be printed to the screen." << endl;
	}
	ostream &out = (!outputFile.empty()&&fout.good()) ? fout : cout; // write to cout by default, but if outputFile was provided and is valid, use it instead.

	//Calculate the extent of the sequence in case this is local:
	int length;
	int nstart;
	int nend;
	if (nuc_start == 0) {
		nstart = 1;	
	}
	else nstart = nuc_start;
	if (nuc_end == 0) {
		nend = rna.GetSequenceLength();
	}
	else nend = nuc_end;
	length = nend - nstart + 1;

	for(int n = start; n < end; n++) {
		double defect = rna.GetEnsembleDefect(n,nuc_start,nuc_end, nucfilename.c_str());
		if (raw)
			out << defect/((double) (length)) << endl;
		else
			out << "Structure " << n << ": Ensemble_Defect =\t" << defect << "\t\tNormalized_ED =\t" << defect/((double)(length)) << endl;
	}
	return rna.GetErrorCode();
}

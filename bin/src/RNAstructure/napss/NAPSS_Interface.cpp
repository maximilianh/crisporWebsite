/*
 * NAPSS, a program that predicts RNA secondary structures with pseudoknots with the aid of NMR constraints data.
 *
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 * Written by James Hart, David Mathews, Jon Chen, Stanislav Bellaousov
 *
 * Modified from the command-line interface to Dynalign, 
 * written by David Mathews; Copyright 2002, 2003, 2004, 2005, 2006
 *
 * Modified to parse command line interface by Jessica Reuter and Stanislav Bellaousov 2012
 *
 * Modified to read Triplet Constraints by Stanislav Bellaousov in August 2012
 *
 */

#include "NAPSS_Interface.h"
#include <iomanip>

// Flags for text output
#undef VERBOSE_MODE
//#define VERBOSE_MODE

// Flags that output dotplots
//#define DOTPLOTS_OUT
#undef DOTPLOTS_OUT

// Flags for output associated with Triplet Constraints
#undef TRIPLET_VERBOSE_MODE
//#define TRIPLET_VERBOSE_MODE

// Flags for debug features
#undef DEBUG_MODE
// #define DEBUG_MODE
#undef ALREADY_USED_CHECK
//#define ALREADY_USED_CHECK

//#define OUTPUT_MATCHES
#undef OUTPUT_MATCHES


//Define the ProbKnot mode, where the pairs are populated using ProbKnot method instead of dynamic
//#define PROBKNOT_MODE
#undef PROBKNOT_MODE

//Define to enable removing matches with broken pairs after "RemoveIsolatedPairs"
#define RemoveBrokenMatch 
//#undef RemoveBrokenMatch 

//Define to enables helix extension before the folding
#define ExtendBefore
//#undef ExtendBefore

//Defined to enable helix extension after the folding
#define ExtensionAfter
//#undef ExtensionAfter

//Define to allow matched pairs during 'dynamic'
#define AllowMatchedPairs
//#undef AllowMatchedPairs

//Define to enable debug execution of the code, which doesn't sort structures or doesn't remove
//duplicate structures. It also outputs the matches that produce lowest energy structures
//#define DEBUG_MATCHES
#undef DEBUG_MATCHES


///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
NAPSS_Interface::NAPSS_Interface() {
    
    // Initialize pseudoknot energy model parameters (kcal/mol).
	P1=DEFAULT_PSEUDOKNOT_P1; //  e.g. 0.35  Pseudoknot energy model parameters (kcal/mol).

	// Initialize pseudoknot energy model parameters (kcal/mol).
	P2=DEFAULT_PSEUDOKNOT_P2;  //e.g. 0.65  Pseudoknot energy model parameters (kcal/mol). 
    
    slope=DEFAULT_SHAPE_SLOPE;

    intercept=DEFAULT_SHAPE_INTERCEPT;//shape slope and shape intercept

    maxtracebacks=100;

    percent=5;//-d option

    windowsize=0;// If user doesn't specify windowsize it will get changed

    cutoff=0;//-p option

    pseudoknotFree=false;//default is to have pseudoknot prediction mode enabled

    warningLimit=1000000;//number of matches (matchVector->size) before the warning message is printed

}


bool NAPSS_Interface::Parse(int argc, char* argv[]){
	
	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "NAPSS" );
	parser->addParameterDescription( "seq file", "The name of a file containing an input sequence." );
	parser->addParameterDescription( "NMR file", "The name of an NMR file with constraints." );
	parser->addParameterDescription( "ct file", "The name of a CT file to which output will be written." );

	// Add the constraint file option.
	vector<string> constraintOptions;
	constraintOptions.push_back( "-c" );
	constraintOptions.push_back( "-C" );
	constraintOptions.push_back( "--constraint" );
	parser->addOptionFlagsWithParameters( constraintOptions, "Specify a constraints file to be applied. Default is to have no constraints applied." );

	// Add the maximum percent energy difference to consider in the dotplot option.
	vector<string> DotPlotPercentOptions;
	DotPlotPercentOptions.push_back( "-d" );
	DotPlotPercentOptions.push_back( "-D" );
	DotPlotPercentOptions.push_back( "--DotPercent" );
	parser->addOptionFlagsWithParameters( DotPlotPercentOptions, "Specify a maximum percent energy difference to consider in the dotplot. Default is 5 percent." );
	
	// Add the maximum number of structures option.                                                                                                                                                        
	vector<string> maxStructuresOptions;
	maxStructuresOptions.push_back( "-m" );
	maxStructuresOptions.push_back( "-M" );
	maxStructuresOptions.push_back( "--maximum" );
	parser->addOptionFlagsWithParameters( maxStructuresOptions, "Specify a maximum number of structures per matched constraint set. Default is 100 structures." );
	
	// Add the Penalty1 option.
	vector<string> Penalty1Options;
	Penalty1Options.push_back( "-p1" );
	Penalty1Options.push_back( "-P1" );
	Penalty1Options.push_back( "--Penalty1" );
	parser->addOptionFlagsWithParameters( Penalty1Options, sfmt("Specify a pseudoknot penalty P1. Default is %1.2f kcal/mol.", P1) );

	// Add the Penalty2 option.
	vector<string> Penalty2Options;
	Penalty2Options.push_back( "-p2" );
	Penalty2Options.push_back( "-P2" );
	Penalty2Options.push_back( "--Penalty2" );
	parser->addOptionFlagsWithParameters( Penalty2Options, sfmt("Specify a pseudoknot penalty P2. Default is %1.2f kcal/mol.", P2) );

	// Add the percent energy difference option.
	vector<string> percentOptions;
	percentOptions.push_back( "-p" );
	percentOptions.push_back( "-P" );
	percentOptions.push_back( "--percent" );
	parser->addOptionFlagsWithParameters( percentOptions, "Specify a maximum percent energy difference. Default is 0 which means that all structures are outputted." );

	// Add the positions paired output file option.
	vector<string> posPairsOptions;
	posPairsOptions.push_back( "-pp" );
	posPairsOptions.push_back( "-PP" );
	posPairsOptions.push_back( "--posPaired" );
	parser->addOptionFlagsWithParameters( posPairsOptions, "Specify the name of the positions paired style output file. Default is to have no file specified." );

	// Add the pseudoknot-free output option.
	vector<string> pseudoknotOptions;
	pseudoknotOptions.push_back( "-pf" );
	pseudoknotOptions.push_back( "-PF" );
	pseudoknotOptions.push_back( "--pseudoknotFree" );
	parser->addOptionFlagsNoParameters( pseudoknotOptions, "Specify pseudoknot-free prediction mode. Default is to predict pseudoknots." );

	// Add the SHAPE option.
	vector<string> shapeOptions;
	shapeOptions.push_back( "-sh" );
	shapeOptions.push_back( "-SH" );
	shapeOptions.push_back( "--SHAPE" );
	parser->addOptionFlagsWithParameters( shapeOptions, "Specify a SHAPE data file to be used to generate pseudoenergy restraints." );

	// Add the SHAPE intercept option.
	vector<string> shapeInterceptOptions;
	shapeInterceptOptions.push_back( "-si" );
	shapeInterceptOptions.push_back( "-SI" );
	shapeInterceptOptions.push_back( "--SHAPEintercept" );
	parser->addOptionFlagsWithParameters( shapeInterceptOptions, sfmt("Specify an intercept used with SHAPE constraints. Default is %1.2f kcal/mol.", intercept) );

	// Add the SHAPE slope option.
	vector<string> shapeSlopeOptions;
	shapeSlopeOptions.push_back( "-sm" );
	shapeSlopeOptions.push_back( "-SM" );
	shapeSlopeOptions.push_back( "--SHAPEslope" );
	parser->addOptionFlagsWithParameters( shapeSlopeOptions, sfmt("Specify a slope used with SHAPE constraints. Default is %1.2f kcal/mol.", slope) );

	// Add the window size option.                                                                                                                                                                         
	vector<string> windowOptions;
	windowOptions.push_back( "-w" );
	windowOptions.push_back( "-W" );
	windowOptions.push_back( "--window" );
	parser->addOptionFlagsWithParameters( windowOptions, "Specify a window size. Default is 0 nucleotides." );

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.                                                                                                                                                            
	if( !parser->isError() ) {
		inseq = parser->getParameter( 1 );
		inNMRconstraints = parser->getParameter( 2 );
		outct = parser->getParameter( 3 );
	}

	// Get the constraint file option.
	if( !parser->isError() ) { constraintFile = parser->getOptionString( constraintOptions, true ); }

	// Get the maximum number of structures option.
	if( !parser->isError() ) {
		parser->setOptionInteger( maxStructuresOptions, maxtracebacks );
		if( maxtracebacks <= 0 ) { parser->setError( "maximum number of structures" ); }
	}
	
	// Get the percent energy difference option.
	if( !parser->isError() ) {
		parser->setOptionInteger( percentOptions, cutoff );
		if( cutoff < 0 ) { parser->setError( "percent energy difference" ); }
	}
	
	// Get the pseudoknot penalty 1.
	if( !parser->isError() ) {
		parser->setOptionDouble( Penalty1Options, P1 );
	}

	// Get the pseudoknot penalty 2.
	if( !parser->isError() ) {
		parser->setOptionDouble( Penalty2Options, P2 );
	}

	// Get the window size option.
	if( !parser->isError() ) {
		parser->setOptionInteger( windowOptions, windowsize );
		if( windowsize < 0 ) { parser->setError( "window size" ); }
	}

	//Record if the windowOption was set
	ifwindowOptions=parser->contains(windowOptions);

	// Get the DotPLotPercent energy difference option.
	if( !parser->isError() ) {
		parser->setOptionInteger( DotPlotPercentOptions, percent );
		if( percent < 0 ) { parser->setError( "DotPlot Percent energy difference" ); }
	}

	// Get the pseudoknot-free output option.
	pseudoknotFree = parser->contains( pseudoknotOptions);

	// Get the SHAPE data and options.
	if( !parser->isError() ) {
		inSHAPEfile = parser->getOptionString( shapeOptions );
		if( !parser->isError() ) { parser->setOptionDouble( shapeInterceptOptions, intercept ); }
		if( !parser->isError() ) { parser->setOptionDouble( shapeSlopeOptions, slope ); }
	}

	// Get the outpairs file option.
	if( !parser->isError() ) { outpairs = parser->getOptionString( posPairsOptions, false ); }
	
	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
bool NAPSS_Interface::run() {
	int error = 0;	
	structure *results = NULL; // Will be set to a valid structure if napss is successful (i.e. returns 0).

	cout << "####################################################\n";
	cout << "                  NAPSS Folding Started\n";
	cout << "####################################################\n\n";
	cout << "Reading files..." << flush ;

	//Initialize RNA class with RNA sequence.
	RNA *rnaCT = new RNA(inseq.c_str(),FILE_SEQ);
	//Check for errors
	if ((error=rnaCT->GetErrorCode())!=0) {
		//If there was an error output the error and exit
		cerr << rnaCT->GetErrorMessage(rnaCT->GetErrorCode())<<"\n";
		goto EndRun;
	}

	//If the OutWindowSize option was not specified by the user on the command line
	if(!ifwindowOptions){
		//Scale windowsize based on the length of the sequence
		if(rnaCT->GetSequenceLength()>1200) windowsize=20;
		else if(rnaCT->GetSequenceLength()>800) windowsize=15;
		else if(rnaCT->GetSequenceLength()>500) windowsize=11;
		else if(rnaCT->GetSequenceLength()>300) windowsize=7;
		else if(rnaCT->GetSequenceLength()>120) windowsize=5;
		else if(rnaCT->GetSequenceLength()>50) windowsize=3;
		else windowsize=2;
	}	

	//Read the SHAPE data from the disk
	if(!inSHAPEfile.empty()) rnaCT->ReadSHAPE(inSHAPEfile.c_str(),slope,intercept);

	//Read the constraint data from the disk
	if(!constraintFile.empty()) rnaCT->ReadConstraints(constraintFile.c_str());

	if ((error=rnaCT->GetErrorCode())!=0) {
		//If there was an error output the error and exit
		cerr << rnaCT->GetErrorMessage(rnaCT->GetErrorCode())<<"\n";
		goto EndRun;
	}

	double pseud_param[16];//Array that holds pseudoknot penalty calculation constants.
	//Read pseudoknot parameters from the datapath
	error = ReadPseudoParam(pseud_param);
	if (error != 0) {
		cerr << "Could not read pseudoknot parameters: " << RNA::GetErrorMessage(error) << endl;
		goto EndRun;
	}

	error = napss(rnaCT, inNMRconstraints, 
		results,
		windowsize, cutoff, percent, 
		maxtracebacks, warningLimit,
		pseud_param, P1, P2, pseudoknotFree,
		inSHAPEfile,slope,intercept,
		constraintFile);

	//RMW: TODO: Use error number and add new function napss_error_message to display error message.
	if (error != 0) {
		cerr << "NAPSS Error: " << error << " - " << napss_error_message(error) << endl;
		goto EndRun;
	}

	error = results->ctout(outct.c_str());
	if (error != 0) {
		cerr << "Failed to write output CT file: " << RNA::GetErrorMessage(error) << endl;
		goto EndRun;
	}

	//  Optional: output paired-positions text file (for structure viewing in PseudoViewer3)
	//	Note: this outputs one concatenated file, individual structures must be manually cut from this file
	//  and placed in a new text file for PseudoViewer3 to read it.
	if(!outpairs.empty()) {
		error = pairout(results, outpairs.c_str());
		if (error != 0) {
			cerr << "Error writing paired-positions text file (" << outpairs << "): " << RNA::GetErrorMessage(error);
			goto EndRun;
		}
	}

	// if we get here, NAPSS was successful.
	cout << "\n####################################################\n";
	cout << "                        DONE\n";
	cout << "####################################################\n" << flush;


	// if an error occurs above, execution jumps here to perform cleanup and exit
	EndRun: 

	delete rnaCT;
	if (results != NULL) delete results;

	return error == 0;
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {
	NAPSS_Interface runner;
	return (runner.Parse( argc, argv ) && runner.run()) ? 0 : 1;
}

/*
 * ETEcalculator: program that calculates the distribution and mean end-to-end (ETE) distance of a given sequence in nanometers (nm).
 * Written by Mohammad Kayedkhordeh and Stanislav Bellaousov in Dr. David Mathews Lab at the University of Rochester. 
 * 
 * (c) 2019 Mathews Lab, University of Rochester Medical Center.
 */

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include "ETEcalculator_Interface.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {
	int error = 1;
	ETEcalculator_Interface* runner = new ETEcalculator_Interface();
	bool parseable = runner->parse( argc, argv );
	if( parseable ) { error = runner->run(); }
	delete runner;
	return error;
}

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
ETEcalculator_Interface::ETEcalculator_Interface() {
	// Initialize default parameters.

	// The input ct file which describes the structure to design. The sequence in the ct file will be ignored.
	seqFile = "";

	// Name of the alphabet (e.g. "rna" or "dna" etc).
	alphabet = DT_RNA;


	ensembleSize = 1000;

	// Initialize the random seed.
	seed = 1234;

	bool raw = false;


	//by default, expect a sequence
	load_ct = false;

	// Single strand value from Aalberts model in Angstrom
	ss_dist = 6.2;

	// Branch value from Aalberts model in Angstrom
	paired_dist = 15.0;


	string constraintFile = "";

} 


///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool ETEcalculator_Interface::parse( int argc, char** argv ) {
	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "ETECalculator" );

	parser->addParameterDescription( string("Sequence file"), string("The input .seq file used to calculate the average end-to-end (ETE) distance value.") );

	vector<string> dnaOptions;
	dnaOptions.push_back( "-d" );
	dnaOptions.push_back( "--dna" );
	parser->addOptionFlagsNoParameters( dnaOptions, "Specify that the sequence is DNA, and DNA parameters are to be used. The default is to use RNA parameters." );

	vector<string> CTOptions;
	CTOptions.push_back("--ct");
	parser->addOptionFlagsNoParameters(CTOptions, "Specify that the input is the name of a ct file with structures, rather than a sequence. The default is to expect a sequence.");


	vector<string> alphabetOptions;
	alphabetOptions.push_back( "-a" );
	alphabetOptions.push_back( "--alphabet" );
	parser->addOptionFlagsWithParameters( alphabetOptions, "Specify the name of a folding alphabet and associated nearest neighbor parameters. The alphabet is the prefix for the thermodynamic parameter files, e.g. \"rna\" for RNA parameters or \"dna\" for DNA parameters or a custom extended/modified alphabet. The thermodynamic parameters need to reside in the at the location indicated by environment variable DATAPATH. The default is \"rna\" (i.e. use RNA parameters). This option overrides the --DNA flag." );

	// Add the number option.
	vector<string> ensembleSizeOptions;
	ensembleSizeOptions.push_back( "-n" );
	ensembleSizeOptions.push_back( "--number" );
	parser->addOptionFlagsWithParameters( ensembleSizeOptions, "Specify the number of structures that are sampled using Stochastic program to claculate the average end-to-end (ETE) distance value." );

	// Add the random seed option.
	vector<string> seedOptions;
	seedOptions.push_back( "-s" );
	seedOptions.push_back( "-S" );
	seedOptions.push_back( "--seed" );
	parser->addOptionFlagsWithParameters( seedOptions, "Specify the random seed. Default is 1234." );

	// Add the raw option.
	vector<string> rawOptions;
	rawOptions.push_back( "-r" );
	rawOptions.push_back( "--raw" );
	parser->addOptionFlagsNoParameters( rawOptions, "Output only the average end-to-end (ETE) distance as a pure number (with no additional description)." );

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

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) 
		seqFile = parser->getParameter( 1 );

	if (parser->contains( dnaOptions ))
	  alphabet = DT_DNA;

	if (!parser->isError())
		if (parser->contains(CTOptions))
			load_ct = true;

	if (!parser->isError())
		if (parser->contains(alphabetOptions))
		  alphabet = parser->getOptionString(alphabetOptions).c_str();

	// Get the constraint file option.
	if( !parser->isError() ) { constraintFile = parser->getOptionString( constraintOptions, true ); }

	// Get the ensemble size option.
	if( !parser->isError() ) {
		parser->setOptionInteger( ensembleSizeOptions, ensembleSize );
		if( ensembleSize <= 0 ) { parser->setError( "ensemble size" ); }
	}

	// Get the random seed option.
	if( !parser->isError() ) {
		parser->setOptionInteger( seedOptions, seed );
		if( seed <= 0 ) { parser->setError( "random seed" ); }
	}


	raw = parser->contains(rawOptions);

	if( !parser->isError() ) {
		outputFile = parser->getOptionString( fileOption,  false);
	}

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
int ETEcalculator_Interface::run() {
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
	RNA* rna;
	TProgressDialog progress;
	if (load_ct) {
		rna = new RNA(seqFile.c_str(), FILE_CT, alphabet);

		ErrorChecker<RNA> checker(rna);
		errNum = checker.isErrorStatus();
		if (errNum) {
			delete rna;
			return errNum;
		}

	}
	else {
		rna = new RNA(seqFile.c_str(), FILE_SEQ, alphabet);
		ErrorChecker<RNA> checker(rna);
		errNum = checker.isErrorStatus();

		if (errNum) {
			delete rna;
			return errNum;
		}



		// Read folding constraints, if applicable.
		if (constraintFile != "") {
			errNum = checker.isErrorStatus(rna->ReadConstraints(constraintFile.c_str()));
			if (errNum) return errNum;
		}

		// Print a message saying that the partition function has started.
		if (!raw) {
			cout << "Calculating partition function..." << endl;
			
			rna->SetProgress(progress);
		}

		// Run the partition function, then check the error status.
		int partError = rna->PartitionFunction();
		errNum = checker.isErrorStatus(partError);
		if (errNum) return errNum;
		// Print out a message saying that partition function is done.
		//if( error == 0 & !raw) {
		if (errNum == 0 & !raw) {
			cout << "Done." << endl;
			cout << "Analyzing stochastic samples..." << endl;
		}


		// Show a message saying that the main calculation has started.

		int mainCalcError = rna->Stochastic(ensembleSize, seed);


	}
	if (!raw) {
		cout << "Done." << endl; 
		cout << "Calculating end-to-end distance..." << endl;
		//TProgressDialog progress;
		rna->SetProgress(progress);
	}

	//if (!raw)  
	//	cout << endl << "done." << endl;
	

	etecalculator(rna, outputFile);

	/*for(int n = 1; n < ensembleSize; n++) {
		double dist_value = GetETEDistance(n);
		if (raw)
			out << dist_value << endl;
		else
			out << "Structure " << n << ": ETE Distance =\t" << dist_value << endl;
	}*/
	int error = rna->GetErrorCode();
	delete rna;
	return error;
}

void ETEcalculator_Interface::etecalculator(RNA* rna, string outputFile) {

	ofstream fout;
	if (!outputFile.empty()) {
		fout.open(outputFile.c_str());
		if (!fout.good())
			cerr << "Error opening output file " << outputFile << ". Results will be printed to the screen." << endl;
	}
	ostream &out = (!outputFile.empty()&&fout.good()) ? fout : cout; // write to cout by default, but if outputFile was provided and is valid, use it instead.

	//Average end-to-end distance calculated over ensemble of structures and reported in nanometers (nm)
	double avg_distance=0;
    for(int j=1;j<=rna->GetStructureNumber();++j){
        //cerr << "\rStructure: " << j ;
        int distance1=0;//set distance counter to 0
        int distance2=0;//set distance2 counter to count 
        //walk down the sequence starting with the first nucleotide
        int i = 1;
        while (i<=rna->GetSequenceLength()){
        //for(int i=1;i<=rna->GetSequenceLength();++i){
            //if the nucleotide is unpaired
            if(rna->GetPair(i,j)==0){
                                   
                //CHECK FOR COAXIAL STACK
                bool stack=false;
                //check if it is possible by length
                /*if(rna->GetPair(i+1,j)+9<=rna->GetSequenceLength()){
                    if(rna->GetPair(rna->GetPair(i+1,j)+2,j)+1<=rna->GetSequenceLength()){
                        
                        //Check if it is a coaxial stack
                        if(rna->GetPair(i+1,j)!=0 && // check if there is a helix in i+1 position
                           rna->GetPair(rna->GetPair(i+1,j)+1,j)==0 && // check if there is an unpaired nucleotide between the two stacks
                           rna->GetPair(rna->GetPair(i+1,j)+2,j)!=0 && // check if there is a second helix after the unpaired nucleotide
                           rna->GetPair(rna->GetPair(rna->GetPair(i+1,j)+2,j)+1,j)==0){ //check if the nucleotide after the two helices is unpaired
                            
                            //COAXIAL STACK IF FOUND                            
                            stack=true;
                            distance2++;//add 2 to distance. Since we already added 1 above, the total is 3.
                            //out << " d21=" << distance2 << "@" << i;
                            i=rna->GetPair(rna->GetPair(i+1,j)+2,j)-1;//jump 'i' over the coaxial stack
                            //out << "->i=" << i;
                        }
                    }
                }*/
                
                //Otherwise it is a singlestranded link
                if(!stack){
                    distance1++;
                    ++i;
                    //out << " d11=" << distance1 << "@" << i;
                }
            }
                       
            
            //if the nucleotide is paired
            else if(rna->GetPair(i,j)!=0){

                bool stack=false;
                //Check if it is a 3' paired nucleotide. It means that there was either a
                //walk across the pair, or across the stack. 
                /*if(rna->GetPair(i,j)<i){
                    //add 1 link to singlestranded distance
                    distance1++;
                    //out << " d12=" << distance1 << "@" << i;
                    stack=true;
                }
                
                //CHECK FOR COAXIAL STACK
               
                //check if coaxial stack is possible by length
                else if(rna->GetPair(i,j)+7<=rna->GetSequenceLength()){

                    //check if it is a coaxial stack
                    if(rna->GetPair(rna->GetPair(i,j)+1,j)!=0){

                        //COAXIAL STACK IF FOUND
                        stack=true;
                        distance2++;
                        //out << " d22=" << distance2 << "@" << i;
                        i=rna->GetPair(rna->GetPair(i,j)+1,j)-1;//jump 'i' across the coaxial stack.
                        //out << "->i=" << i;
                    }
                }

                //Otherwise it is a walk across the pair
                if(!stack){

                    distance2++;
                    i=rna->GetPair(i,j)-1;
                    //out << " d23=" << distance2 << "@" << i;

                }
            }*/
                ++distance2;
                i = rna->GetPair(i,j)+1;
            }
        }
        double distance = sqrt ( ( ( pow (distance1,6.0/5) * pow (ss_dist,2) ) + ( pow (distance2,6.0/5) * pow (paired_dist,2) ) ) ) / 10.0;
        if (!raw)
			out << "Structure " << j << ": ETE Distance =\t" << distance << endl;
        avg_distance+= distance;
    }
    avg_distance/=rna->GetStructureNumber();
	if (raw)
			out << avg_distance << endl;
	else
			out << "Average End-to-end distance =\t" << avg_distance << endl;
}
/*
 * A program that ...TODO: enter description
 * This strand of nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2015 Mathews Lab, University of Rochester Medical Center.
 * Written by Richard M. Watson (2015)
 */

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "design_interface.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"


using namespace std;

const char* bstr(bool value) { return value ? "true" : "false"; }
const char* yn(bool value) { return value ? "yes" : "no"; }

//Set to true in order to output the values of all parameter variables before running design.
#define DEBUG_INTERFACE true

// DEFAULT_PRESELECTED_SEQUENCES: If true, the default behavior will be to use preselected 
//   sequence segments (instead of completely random sequences). The user will then have 
//	 to pass the "-r" option in order to use random sequences. 
//   If false, random sequences will be default, and the user will pass "-p" to use 
//   preselected sequences.
#define DEFAULT_PRESELECTED_SEQUENCES false
								

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {
	int error = 1;
	designInterface* runner = new designInterface();
	bool parseable = runner->parse( argc, argv );
	if (DEBUG_INTERFACE && parseable)
		runner->showParams();
	if( parseable ) { error = runner->run(); }
	delete runner;
	return error;
}

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
designInterface::designInterface() {

	// Initialize the calculation type description.
	calcType = "design";

	// Initialize default parameters.

	// The input ct file which describes the structure to design. The sequence in the ct file will be ignored.
	ctFile = "";

	// Determines the method for choosing sequence -- use pre-selected sequences (true) or randomly choose the sequence (false)
	preselectSequences = DEFAULT_PRESELECTED_SEQUENCES;

	biasDesign = false;

	//by default, do not allow isolated base pairs
	allow_isolated = false;
	
	// Name of the alphabet (e.g. "rna" or "dna" etc).
	alphabet = DT_RNA;

	//// if true, avoids performing a partition function on the complete sequence.  Default is false.
	//useHeuristicMode = false;

	// The maximum allowed ensemble defect per nucleotide.
	defectThreshold = 0.01;

	// The maximum extent to which the structure will be sub-divided in the binary decomposition.  The default is 5.
	maxDepth = 5;
	
	// The output ct sequence file.
	outFile = "";

	// The input bias probailities file
	biasProbFile = "";

	//// Turn pseudoknot design on or off
	//designPseudoknot = false;

	// The maximum number of redesigns per parent node. The default is 10.
	maxRedesign = 10;

	// The maximum number of times a nucleotide will be mutated during defect-weighted reoptimization. The default is 4.
	maxMutate = 4;
	
	// The maximum number of times the leaf is re-optimized at random. The default is 3.
	maxLeafRedesign = 3;

	// Print sequence to standard output.
	//stdPrint = false;   Not used. Instead, see docs for --output flag
	
	// Whether or not to time the design process.
	useTimer = false;

	// Random seed. Used to initialize the random number generator.
	randSeed = time(NULL); 

	setRandomSeed = false; // The program displays the random seed, unless it was specified by the user. This keeps track of whether or not they did.
} 



void designInterface::showParams() {
	cout << endl << "------------------------------ Selected Options -------------------------------" << endl;
	printf(	"  CtFile:   \t%s\n  OutFile:\t%s\n  EDThreshold:\t%g\tPreselectSeq:\t%s\tAlphabet:\t%s\n"
			"  MaxRedesign:\t%d\tMaxMutate:\t%d\tMaxLeafDesign:\t%d\n"
			"  MaxDepth:\t%d\tUseTimer:\t%s\tRandomSeed:\t%ld", //\n "  DesignPsknot:\t%s\tHeuristicMode:\t%s",
			
			ctFile.c_str(),  outFile.empty() ? "<standard output>" :  outFile.c_str(),
			defectThreshold, yn(preselectSequences), alphabet.c_str(),
			maxRedesign, maxMutate, maxLeafRedesign,
			maxDepth, yn(useTimer), randSeed //, "<unused>" /*yn(designPseudoknot)*/, "<unused>" /*yn(useHeuristicMode)*/
		);
	cout << endl << std::string(79, '-')  <<  endl;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool designInterface::parse( int argc, char** argv ) {
	char description[256];

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "design" );
	parser->addParameterDescription( string("ct structure file"), string("The input ct file which describes the structure to design. The sequence in the ct file will be ignored.") );
	//parser->addParameterDescription( "output sequence file", "The file to which the resulting designed sequence is written." );

	//parser->addParameterDescription( "defect threshold", "The maximum allowed ensemble defect per nucleotide." );
	vector<string> defectOptions;
	defectOptions.push_back( "-e" );
	defectOptions.push_back( "--error" );
	parser->addOptionFlagsWithParameters( defectOptions, "The maximum allowed ensemble defect per nucleotide.");

	vector<string> dnaOptions;
	dnaOptions.push_back( "-d" );
	dnaOptions.push_back( "--dna" );
	parser->addOptionFlagsNoParameters( dnaOptions, "Specify that the sequence is DNA, and DNA parameters are to be used. The default is to use RNA parameters." );

	vector<string> alphabetOptions;
	alphabetOptions.push_back( "-a" );
	alphabetOptions.push_back( "--alphabet" );
	parser->addOptionFlagsWithParameters( alphabetOptions, "Specify the name of a folding alphabet and associated nearest neighbor parameters. The alphabet is the prefix for the thermodynamic parameter files, e.g. \"rna\" for RNA parameters or \"dna\" for DNA parameters or a custom extended/modified alphabet. The thermodynamic parameters need to reside in the location indicated by environment variable DATAPATH. The default is \"rna\" (i.e. use RNA parameters). This option overrides the --DNA flag." );

	vector<string> preselectOptions;
#if DEFAULT_PRESELECTED_SEQUENCES    // This flag can be used to switch the default sequence generation algorithm from random to pre-selected (or vice-versa)
	preselectOptions.push_back( "-r" );
	preselectOptions.push_back( "--random" );
	parser->addOptionFlagsNoParameters( preselectOptions, "Specify that all nucleotides should be chosen at random. The default is to use pre-selected sequence segments."); // The default is to use random sequence generation." );
#else
	preselectOptions.push_back( "-p" );
	preselectOptions.push_back( "--preselect" );
	parser->addOptionFlagsNoParameters( preselectOptions, "Specify that use pre-selected sequence segments should be used. The default is that all nucleotides are chosen at random." );
#endif

	vector<string> isolatedPairOptions;
	isolatedPairOptions.push_back("-i");
	isolatedPairOptions.push_back("--isolated");
	parser->addOptionFlagsNoParameters(isolatedPairOptions, "Specify that isolated pairs are allowed.  The default is to use a heuristic to attempt to forbid isolated pairs during structure prediction.");


	vector<string> biasDesignOptions;
	biasDesignOptions.push_back( "-b" );
	biasDesignOptions.push_back( "--bias" );
	parser->addOptionFlagsWithParameters( biasDesignOptions, "Specify that all nucleotides should be chosen with specified probabilities. The default is to select nucleotides with equal probabilities.");

	// Heuristics mode is no longer supported
	//vector<string> heuristicsOptions;
	//heuristicsOptions.push_back( "-u" );
	//heuristicsOptions.push_back( "--heuristic" );
	//parser->addOptionFlagsNoParameters( heuristicsOptions, "Use heuristics to avoid performing a partition function on the complete sequence.  The default is to not use heuristics." );

	vector<string> maxDepthOptions;
	maxDepthOptions.push_back( "-md" );
	maxDepthOptions.push_back( "--maxdepth" );
	sprintf(description, "Max-depth: The maximum extent to which the structure will be sub-divided in the binary decomposition. The default is %d.", maxDepth);
	parser->addOptionFlagsWithParameters( maxDepthOptions, description);

	// Pseudoknot is no longer supported
	//vector<string> pseudoknotOptions;
	//pseudoknotOptions.push_back( "-k" );
	//pseudoknotOptions.push_back( "-K" );
	//pseudoknotOptions.push_back( "--pseudoknot" );
	//parser->addOptionFlagsNoParameters( pseudoknotOptions, "Turn pseudoknot design on. The default is off.");

	vector<string> maxRedesignOptions;
	maxRedesignOptions.push_back( "-mr" );
	maxRedesignOptions.push_back( "--maxredesign" );
	sprintf(description, "The maximum number of redesigns per parent node. The default is %d.", maxRedesign);
	parser->addOptionFlagsWithParameters( maxRedesignOptions, string(description));

	vector<string> maxMutateOptions;
	maxMutateOptions.push_back( "-mm" );
	maxMutateOptions.push_back( "--maxmutate" );
	sprintf(description, "The maximum number of times a nucleotide will be mutated during defect-weighted reoptimization. The default is %d.", maxMutate);
	parser->addOptionFlagsWithParameters( maxMutateOptions, string(description));

	vector<string> maxLeafRedesignOptions;
	maxLeafRedesignOptions.push_back( "-ml" );
	maxLeafRedesignOptions.push_back( "--maxleaf" );
	sprintf(description, "The maximum number of times a leaf can be re-optimized at random. The default is %d.", maxLeafRedesign);
	parser->addOptionFlagsWithParameters( maxLeafRedesignOptions, string(description));

	vector<string> outputFileOptions;
	outputFileOptions.push_back( "-o" );
	outputFileOptions.push_back( "--output" );
	parser->addOptionFlagsWithParameters( outputFileOptions, "Specify the output file. By default the resulting designed sequence is written to standard output only. This flag instructs the program to output the structure (in ct format) to the specified file." );

	vector<string> timerOptions;
	timerOptions.push_back( "-t" );
	timerOptions.push_back( "--timer" );
	parser->addOptionFlagsNoParameters( timerOptions, "Use a timer to measure the duration of the design process and print the elapsed time to standard output." );

	vector<string> seedOptions;
	seedOptions.push_back( "-s" );
	seedOptions.push_back( "--seed" );
	parser->addOptionFlagsWithParameters( seedOptions, "Specify a random seed. This is required to get exactly reproducible results. (The default is to use a seed based on the current system time)." );

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		ctFile = parser->getParameter( 1 );
		//outFile = parser->getParameter( 2 ); //now output file is set with -o flag (--output)
		//defectThreshold = atof(parser->getParameter( 3 ).c_str());
		//if( defectThreshold <= 0.0 ) { parser->setError( "defect threshold" ); }
	}
	// Get boolean parameters (true if present, false otherwise)
	if (!parser->isError() )
		preselectSequences = parser->contains( preselectOptions ) != DEFAULT_PRESELECTED_SEQUENCES; // If user passed the flag, they want the OPPOSITE of the default.
	
	//Get boolean parameters for biasDesign option and the name of bias probability file.
	if (!parser->isError()) {
		biasProbFile = parser->getOptionString(biasDesignOptions).c_str();
		if (parser->contains( biasDesignOptions ))
	  		biasDesign=true;
	}

	if (parser->contains( dnaOptions ))
	  alphabet = DT_DNA;
	//useHeuristicMode = parser->contains( heuristicsOptions); 
	//designPseudoknot = parser->contains( pseudoknotOptions); 
	useTimer = parser->contains( timerOptions );

	if (!parser->isError())
		if (parser->contains(isolatedPairOptions))
			allow_isolated = true;

	if (!parser->isError())
		if (parser->contains(alphabetOptions))
		  alphabet = parser->getOptionString(alphabetOptions,false).c_str();

	if (!parser->isError()) {
		outFile = parser->getOptionString(outputFileOptions, false);
	}

	// Get integer parameters: maxDepth, maxRedesign, maxMuate, maxLeafRedesign
	if( !parser->isError() ) {
		parser->setOptionDouble( defectOptions, defectThreshold);
		if( defectThreshold <= 0 ) { parser->setError( "Maximum Ensemble Defect" ); }
	}
	if( !parser->isError() ) {
		setRandomSeed = parser->setOptionLong( seedOptions, randSeed);
		//	cout << "Random Seed: " << randSeed << endl;
		//if( defectThreshold <= 0 ) { parser->setError( "Maximum Ensemble Defect" ); }
	}
	if( !parser->isError() ) {
		parser->setOptionInteger( maxDepthOptions, maxDepth);
		if( maxDepth <= 0 ) { parser->setError( "Max-Depth" ); }
	}
	if( !parser->isError() ) {
		parser->setOptionInteger( maxRedesignOptions, maxRedesign);
		if( maxRedesign <= 0 ) { parser->setError( "Max-Redesigns per Parent" ); }
	}
	if( !parser->isError() ) {
		parser->setOptionInteger( maxMutateOptions, maxMutate);
		if( maxMutate <= 0 ) { parser->setError( "Max-Mutate" ); }
	}
	if( !parser->isError() ) {
		parser->setOptionInteger( maxLeafRedesignOptions, maxLeafRedesign);
		if( maxLeafRedesign <= 0 ) { parser->setError( "Max-Leaf Redesign" ); }
	}

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
int designInterface::run() {
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

	//////////////////////////////////////////
	// Verify user unput
	//////////////////////////////////////////
	// heuristics was experimental and is no longer supported
	//if (useHeuristicMode) {
	//	cerr << "The use of heuristic mode is no longer supported.\n";
	//	return 1;
	//}

	//if(!preselectSequences && designPseudoknot) {  //The program cannot design pseudoknot using random design algorithm
	//	cerr << "Cannot design a sequence at random with pseudoknots.\n";
	//	return 1;
	//}

	//if(!preselectSequences && useHeuristicMode) {//The program cannot use heuristic algorithm with random design algorithm
	//	cerr << "Cannot design a sequence using random algorithm along with heuristics.\n";
	//	return 1;
	//}

	if(preselectSequences && alphabet!=DT_RNA) {  //The program does not yet have pre-selected sequences for DNA
		if (alphabet==DT_DNA)
			cerr << "The use of pre-selected sequences for DNA is not yet fully implemented.\n";
		else
			cerr << "The use of pre-selected sequences for custom nucleotides or extended alphabets is not yet fully implemented.\n";
		return 1;
	}


	time_t timerStart,timerEnd;
	if (useTimer) {
		//Start the timer
		time (&timerStart);
		//END: Start the timer
	}

	//////////////////////////////////////////////////
	/////RUN DESIGN///////////////////////////////////
	//

	//Initialize class 'design'
	design *des; 
	des = new design(ctFile.c_str(),  alphabet.c_str()); //Read the structure into the 'design' class 'des'
	errNum = des->GetErrorCode();
	if (errNum != 0) goto exitProgram; //At bottom of function. deletes des and returns errNum.
	
	// call design_sequence in design.cpp to design a sequence
	errNum = des->design_sequence(
		defectThreshold,!preselectSequences, biasDesign, biasProbFile.c_str(), maxDepth,
		/*useHeuristicMode*/ false,maxRedesign,maxMutate,
		maxLeafRedesign, randSeed, allow_isolated); 
	if (errNum != 0) goto exitProgram; //At bottom of function. deletes des and returns errNum.
	//
	/////END "Run Design"/////////////////////////////
	//////////////////////////////////////////////////

	//////////////////////////////////////////////////
	//////START Designing Pseudoknot//////////////////
	//
	//if(designPseudoknot){
	//	vector<vector<string> > Helices (11,vector<string>());//Initialize vector to hold helix sequences
	//	vector<vector<string> > Loops (11,vector<string>());//Initialize vector to hold single stranded sequences
	//	vector<int> pairsBroken;//Initialize a vector to hold broken pairs in the pseudoknot
	//	int sequenceLength;//stores the sequence length	
	//	//Read the structure into RNA class 'rna' and 'rna1'. One of them will have pseudoknot removed
	//	// and will be compared to the other one to see which pairs were removed.
	//	RNA *rna=new RNA(ctFile.c_str(),FILE_CT,isRNA);
	//	RNA *rna1=new RNA(ctFile.c_str(),FILE_CT,isRNA);
	//	DesignPseudoknot(isRNA,rna,rna1,Helices,Loops,pairsBroken,sequenceLength,des);//Design a pseudoknot
	//	delete rna;
	//	delete rna1;
	//}
	//
	/////END Designint Pseudoknot"////////////////////
	//////////////////////////////////////////////////

	// If the user did not explicitly set the random seed, display it (in case they want to repeat it).
	if (!setRandomSeed)
		cout << "Random Seed= " << randSeed << endl;

	// Print designed sequence to standard output if the print option (-p) was specified.
	if (outFile.empty()) {
		int numofbases = des->GetSequenceLength();
		char* nucs = des->GetStructure()->nucs;
		char* sequence = new char[numofbases+1];
		for (int i = 0; i < numofbases; i++) {
			sequence[i] = nucs[i+1];
		}
		sequence[numofbases] = 0;
		cout << "Result= " << sequence << endl;
		delete[] sequence;
	} else {
		//DO NOT calculate energies. User can do this separately if desired. --- errNum - des->WriteThermodynamicDetails(NULL, false); // re-calculates free energy of structure.
		errNum = des->WriteCt(outFile.c_str(), false, CTComments::None);//Write a .ct file with the correct sequence
		if (errNum != 0) goto exitProgram; //At bottom of function. deletes des and returns errNum.
	}
	cout << "NED= " << defectThreshold << endl;

	//STOP the timer
	if (useTimer) {
		time (&timerEnd);
		double dif = difftime (timerEnd,timerStart);
		cout << "RunTime (s)= " << dif << '\n' << flush;//Output the time to screen. 
	}

exitProgram:
	if (errNum != 0) //Check for errors. 
		std::cerr << "ERROR " << errNum << ": " << des->GetErrorMessage(errNum) << endl;

	delete des;
	return errNum;
}

// The DesignPseudoknot function is not currently used, but may be some time in the future.
//void DesignPseudoknot(bool IsRNA, RNA* rna, RNA* rna1, vector<vector<string> > &Helices, vector<vector<string> > &Loops, vector<int> &pairsBroken, int &sequenceLength, design* des, long randomSeed){
//	//Read the files with preselected helix and single stranded sequences 
//	string hf;
//	string lf;
//	if(IsRNA){//if GetBackbone returns true, read RNA sequence files
//		hf="design.RNA.Helices.dat";//Holds the file name with helix sequences.
//		lf="design.RNA.Loops.dat";//Holds the file name with loop sequences
//	}
//	else{//else, if GetBackbone returns false, read DNA files
//		hf="design.DNA.Helices.dat";//Holds the file name with helix sequences.
//		lf="design.DNA.Loops.dat";//Holds the file name with loop sequences
//	}
//	string dp(".");//Holds the path to data_tables.
//	if(getDataPath()!=NULL){//If Datapath is not NULL...
//		dp=string(getDataPath());//...set 'dp' to hold "DATAPATH".
//	}
//	
//	//Append the filename to the path 'dp' and store in 'fn'.
//	hf=dp+"/"+hf;
//	lf=dp+"/"+lf;
//	
//	string lineoftext;//String that temporarily holds strings of data from the read files.
//	//Open file
//	ifstream readHfile(hf.c_str());
//	ifstream readLfile(lf.c_str());
//	//          cerr << "lf=" << lf.c_str() << endl;
//	while(!readHfile.eof()){
//		readHfile >> lineoftext;
//		if(lineoftext.empty()) continue;
//		Helices.at((lineoftext.size()-3)/2).push_back(lineoftext);
//	}
//	while(!readLfile.eof()){
//		readLfile >> lineoftext;
//		//              cerr << "lineoftext=" << lineoftext << endl;
//		if(lineoftext.empty()) continue;
//		Loops.at(lineoftext.size()).push_back(lineoftext);
//	}
//	readHfile.close();
//	readLfile.close();
//	/////////END: read file with preselected helix and loop sequences
//	
//	//Read the structure into RNA class to evaluate the pseudoknot
//	//RNA *rna=new RNA(argv[1],FILE_CT,IsRNA);
//	//RNA *rna1=new RNA(argv[1],FILE_CT,IsRNA);
//	sequenceLength=rna->GetSequenceLength();
//	//Break the pseudoknot in 'rna'
//	rna->BreakPseudoknot(false);
//	//Set 'pairsBroken' to be the length of a sequence +1 to account for indesing 
//	//starting with 1 and not with 0
//	pairsBroken=vector<int>(rna1->GetSequenceLength()+1,0);
//	//When there is a mismatch in pairing between the structure with the broken
//	//pseudoknot and the structure with intact pseudoknot, store this mismatch in
//	//vector 'pairsBroken'
//	for (int i=1;i<=rna1->GetSequenceLength();++i){
//		if(rna1->GetPair(i)!=rna->GetPair(i)){
//			pairsBroken.at(i)=rna1->GetPair(i);
//			pairsBroken.at(rna1->GetPair(i))=i;
//		}
//	}
//	
//	//debug
//	//for(int i=0; i<pairsBroken.size();i++){
//	//	if(pairsBroken[i]!=0) cerr << i << " is paired with " << pairsBroken[i] << '\n' ;
//	//}
//	//debug
//
//	//initialize the random number generator
//	randomnumber dice;
//	dice.seed(randomSeed);//seed the random number generagor
//	double roll;//will hold the number generaged by the random number generagor
//	//Fill in the sequence
//	for (int j=1;j<=sequenceLength;j++) {
//		if(pairsBroken[j]!=0&&j<pairsBroken[j]){//If it is a 5' part of the helix
//			int Hcount=1;//set counter of helix length to 1
//			vector<int> HelixSplit;//initiate a vector 'HelixSplit' to hold the start positions of the helix
//			HelixSplit.push_back(j);//store the start position of the helix
//			while(pairsBroken[j]-1==(pairsBroken[j+1])){//while there is a stack without any interruptions
//				Hcount++;//add 1 to the helix length
//				j++;//go to the next nucleotide
//			}
//			if(Hcount>10){//if the helix is longer then 10 base pairs
//				while(Hcount>20){//while the helix is longer then 20 base pairs
//					Hcount-=10;//trim the helix by 10 base pairs
//					HelixSplit.push_back(HelixSplit[HelixSplit.size()-1]+10);//store the location where the starting location of the next helix will be
//				}
//						
//				int split=((int)floor((double)Hcount/2));//'split' splits the helix that is at most 20 base pairs long into two strands of largest
//				//while(split>10||Hcount-split>10)split=(floor((double)Hcount/2.0));//if any of the strands are longer then 10 base pairs,
//				// re-split the helix
//				HelixSplit.push_back(HelixSplit[HelixSplit.size()-1]+split);//store the start location of the second split helix
//						
//				Hcount-=split;
//				HelixSplit.push_back(HelixSplit[HelixSplit.size()-1]+Hcount);//store the location of the end of the second split helix +1
//						
//			}//END if(Hcount>10)
//			else HelixSplit.push_back(j+1);//if the helix length is less then or equal to 10, store the location of the end of the helix +1
//					
//			vector<string> seqlist;
//			for(int y=0;y<HelixSplit.size()-1;y++){
//				roll = dice.roll();
//				seqlist.push_back(Helices[HelixSplit[y+1]-HelixSplit[y]][(int)floor((double)Helices[HelixSplit[y+1]-HelixSplit[y]].size()*roll)]);
//			}
//					
//			//const char* seq=Helices[floor((double)Helices[Hcount].size()*roll)].c_str();//should I cast int instead of floor?
//					
//			int k=0;
//			for(int y=0;y<HelixSplit.size()-1;y++){
//				//cerr << "Helices:" << HelixSplit[y] << "-" << HelixSplit[y+1] << endl;
//				Hcount=HelixSplit[y+1]-HelixSplit[y];
//				int start5=HelixSplit[y];//store the biginning of the helix
//				int end5=HelixSplit[y+1];//store the end of the helix. end is not included in the helix
//				const char* seq=seqlist[y].c_str();
//				//cerr << "seq=" << seq << '\n';
//				for(int i=start5;i<end5;++i){
//					//cerr << des->GetStructure()->nucs[i] << " = "; 
//					des->GetStructure()->nucs[i]=seq[k];
//					if(seq[k]=='A') des->GetStructure()->numseq[i]=1;
//					else if(seq[k]=='C') des->GetStructure()->numseq[i]=2;
//					else if(seq[k]=='G') des->GetStructure()->numseq[i]=3;
//					else if(seq[k]=='U') des->GetStructure()->numseq[i]=4;
//					else if(seq[k]=='T') des->GetStructure()->numseq[i]=4;
//					des->GetStructure()->nucs[pairsBroken[i]]=seq[Hcount*2+2-k];//TAKING INTO ACCOUNT THE 3 SEPERATORS BETWEEN HELIX STRANDS
//					if(seq[Hcount*2+2-k]=='A') des->GetStructure()->numseq[pairsBroken[i]]=1;
//					else if(seq[Hcount*2+2-k]=='C') des->GetStructure()->numseq[pairsBroken[i]]=2;
//					else if(seq[Hcount*2+2-k]=='G') des->GetStructure()->numseq[pairsBroken[i]]=3;
//					else if(seq[Hcount*2+2-k]=='U') des->GetStructure()->numseq[pairsBroken[i]]=4;
//					else if(seq[Hcount*2+2-k]=='T') des->GetStructure()->numseq[pairsBroken[i]]=4;
//					k++;
//					//cerr << des->GetStructure()->nucs[i] << "-" << des->GetStructure()->nucs[pairsBroken[i]] << "\n";
//				}//END for(int i=start5;i<end5;++i){
//				k=0;
//			}// END for(int y=0;y<HelixSplit.size()-1;y++){
//		}//END  if(pairsBroken[j]!=0&&j<pairsBroken[j]){
//	}//END for (int j=1;j<=sequenceLength;j++)
//}//END DesignPseudoknot

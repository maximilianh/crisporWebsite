/*
 * OREGA -- [O]ptimizing [R]NA [E]nds using [G]enetic [A]lgorithm
 * [...Description?...]
 *
 * (c) 2015 Mathews Lab, University of Rochester Medical Center.
 * GeneticAlgorithm written by Mohammad Kayedkhordeh (2017-2019)
 * Interface written by Richard M. Watson and Mohammad Kayedkhordeh (2019)
 */

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "orega_interface.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
Orega_interface::Orega_interface() {

	// Initialize the calculation type description.
	calcType = "Orega";

	// The input ct file which describes the structure to design. The sequence in the ct file will be ignored.
	//sequenceFile  (REQUIRED)

	// Name of the alphabet (e.g. "rna" or "dna" etc).
	alphabet = DT_RNA; // set RNA as the default nucleic acid.

	// Number of individuals in the population 
	populationSize = 10;

	// Define the target segment in the sequence where mutation/optimization should be performed.
	//rangeStart, rangeLength;  (REQUIRED)

	// Number of mutation loops
	iterations = 1000;

	// probability that a nucleotide in the target segment should be mutated.
	mutationRate = 0.03;

	// Number of iterations before a recombination/crossover step should occur.
	// (i.e. recombination will occur every 6 iterations)
	recombinationFrequency = 6;

	// Probability that a nucleotide will be selected as a recombination marker
	// in the recombination/crossover steps.
	recombinationRate = 0.03;
    
    // Mutation site will be based on the mean base pair probability at specific nucleotide instead of random selection. 
	// This biases change towards nucleotides of higher base pairing probability.  The default is to choose any nucleotide at random for mutation
	MutationSwitch = 0;

    //Nucleotide G will not be an option during nucleotide mutation.  The default is to allow mutation to G.
	limitG = 0;

    //Sequences that have a premature AUG start codon will be eliminated from the population.  The default is to not test for AUG start codons.
	filterAUG = 0;  

	//Sequences that have a premature CUG start codon will be eliminated from the population.  The default is to not test for CUG start codons.
	filterCUG = 0;     

	//This options specified the maximum number of As in a row.  If this is exceeded in a sequence, that sequence will be eliminated from the population.  The default value is 0, which indicates no filtering will occur.
	filteroligoA = 0;

	//The constant number that linguistic complexity will be weighted in fitness score calculation. The default values is 1.
    complexity_constant = 1;
    //the threshold to determine if a nucleotide is considered base paired or not.
	threshold = 0.4;
	// Name of file where the optimized sequence is saved.
	// outputSequenceFile (REQUIRED)

	// Name of file where the most recently calculated state is saved.
	// (can be used to restart the program from a prior state)
	// outputStateFile; (OPTIONAL)

	// State file that can be input to restart the program from a prior saved state.
	//restartFile;  (OPTIONAL)

	// Name of the objective function to use.
	// Current options are OREGA_SIMPLE or OREGA_COMPLEX or OREGA_ALL (defined in .h)
	objectiveFunction = OREGA_COMPLEX;

	// Random seed. Used to initialize the random number generator.
	randSeed = time(NULL); 

	setRandomSeed = false; // The program displays the random seed, unless it was specified by the user. This keeps track of whether or not they did.
} 

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool Orega_interface::parse( int argc, char** argv ) {
	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( calcType );

	parser->addParameterDescription("sequence file", "The input sequence used in the optimization.");
	parser->addParameterDescription("start", "The number of the first nucleotide in the segment of the sequence that should be targeted for optimization/mutation.");
	parser->addParameterDescription("length", "The total number of nucleotides in the segment of the sequence that should be targeted for optimization/mutation.");
	parser->addParameterDescription( "output sequence file", "The file to which the resulting optimized sequence is written." );

	vector<string> dna_Flag = parser->addFlag(false, "-d --dna", "Specify that the sequence is DNA, and DNA parameters are to be used. The default is to use RNA parameters." );
	vector<string> alphabet_Flag = parser->addFlag(true, "-a --alphabet", "Specify the name of a folding alphabet and associated nearest neighbor parameters. The alphabet is the prefix for the thermodynamic parameter files, e.g. \"rna\" for RNA parameters or \"dna\" for DNA parameters or a custom extended/modified alphabet. The thermodynamic parameters need to reside in the at the location indicated by environment variable DATAPATH. The default is \"rna\" (i.e. use RNA parameters). This option overrides the --DNA flag.", alphabet); 

	vector<string> outputStateFile_Flag = parser->addFlag(true, "-sav --save", 
		"Specify the name of a file where intermediate results can be saved. "
		"This can be used later to restart the optimization where it left off, or to continue optimizing a previous result.");

	vector<string> restartFile_Flag = parser->addFlag(true, "--restart", 
		"Specify the name of a previous state file (created with the --save option) "
		"that should be loaded to restart an optimization from where it left off before, "
		"or continue optimizing a previous result.");

	vector<string> iterations_Flag = parser->addFlag(true, "-i --iter", 
		"Iterations -- The number of optimization rounds to run.  The default is 1000.", iterations);

	vector<string> populationSize_Flag = parser->addFlag(true, "-n --population", 
		"Population Size -- Number of concurrent sequences used in the genetic algorithm.", populationSize);

	vector<string> mutationRate_Flag = parser->addFlag(true, "-mr --mutate", 
		"Mutation Rate -- The probability that a nucleotide in the target segment should be mutated.", mutationRate );
	
	vector<string> complexity_constant_Flag = parser->addFlag(true, "--ComplexityConstant", 
		"ComplexityConstant -- The constant number that linguistic complexity will be weighted in fitness score calculation. The default values is 1.", complexity_constant );

	vector<string> recombinationFrequency_Flag = parser->addFlag(true, "-rf --recomb", 
		"Recombination Frequency -- Number of iterations that are run before a recombination/crossover step occurs.", recombinationFrequency);

	vector<string> MutationSwitch_Flag = parser->addFlag(false, "--MutationSwitch", 
		" Mutation site will be based on the mean base pair probability at specific nucleotide instead of random selection.  This biases change towards nucleotides of higher base pairing probability.  The default is to choose any nucleotide at random for mutation");

	vector<string> limitG_Flag = parser->addFlag(false, "--limitG", 
		" Nucleotide G will not be an option during nucleotide mutation.  The default is to allow mutation to G.");

	vector<string> filterAUG_Flag = parser->addFlag(false, "--filterAUG", 
		" Sequences that have a premature AUG start codon will be eliminated from the population.  The default is to not test for AUG start codons.");

	vector<string> filterCUG_Flag = parser->addFlag(false, "--filterCUG", 
		" Sequences that have a premature CUG start codon will be eliminated from the population.  The default is to not test for CUG start codons.");	

	vector<string> filteroligoA_Flag = parser->addFlag(true, "--filteroligoA", 
		"This options specified the maximum number of As in a row.  If this is exceeded in a sequence, that sequence will be eliminated from the population.  The default value is 0, which indicates no filtering will occur", filteroligoA);

    vector<string> threshold_Flag = parser->addFlag(true, "-t", 
		"This options specified the threshold to determine if a nucleotide is considered base paired or not.  The default value is 0.4", threshold);

	vector<string> recombinationRate_Flag = parser->addFlag(true, "-rr --crossover", 
		"Recombination Rate -- Probability that a nucleotide will be selected as a recombination marker.", recombinationRate);


	vector<string> randSeed_Flag = parser->addFlag(true, "-rs --seed", 
		"Specify a random seed. This is required to get exactly reproducible results. (The default is to use a seed based on the current system time)." );

	/* if many more are added, use a vector or map to hold the values and their names etc. */
	std::stringstream ssObjectiveFunctionDesc;
	ssObjectiveFunctionDesc << "Objective Function -- Choose which objective function to use. " <<
		 "\n\tSimple: " << OREGA_SIMPLE << 
		 "\n\tInclude Complexity: " << OREGA_COMPLEX << 
		 "\n\tInclude base pair number and Complexity: " << OREGA_ALL << "\n";
	vector<string> objectiveFunction_Flag = parser->addFlag(true, "-f --func", 
		ssObjectiveFunctionDesc.str(), objectiveFunction);

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		sequenceFile = parser->getParameter( 1, true );
		rangeStart = atof(parser->getParameter( 2 ).c_str());
		rangeLength = atof(parser->getParameter( 3 ).c_str());
		outputSequenceFile = parser->getParameter( 4 );
	}
	if( !parser->isError() ) {
		if( rangeStart <= 0 ) { parser->setError( "Start Position" ); }
		else if( rangeLength <= 0 ) { parser->setError( "Segment Length" ); }
	}

	// Choose Alphabet -- DNA or RNA	
	if (parser->contains( dna_Flag ))
	  alphabet = DT_DNA;
	if (!parser->isError()) {
		if (parser->contains(alphabet_Flag))
		  alphabet = parser->getOptionString(alphabet_Flag).c_str();
	}
	if( parser->contains(MutationSwitch_Flag)) {
			MutationSwitch = 1;
		}

	if( parser->contains(limitG_Flag)) {
			limitG = 1;
		}
	
	if( parser->contains(filterAUG_Flag)) {
			filterAUG = 1;
		}

	if( parser->contains(filterCUG_Flag)) {
			filterCUG = 1;
		}

	if( !parser->isError() ) {
		parser->setOptionInteger( filteroligoA_Flag, filteroligoA);
		if( filteroligoA < 0 ) { parser->setError( "filteroligoA" ); }
	}
	if( !parser->isError() ) {
		parser->setOptionDouble( threshold_Flag, threshold);
		if( threshold < 0.0 || threshold > 1.0) { parser->setError( "threshold" ); }
	}
	if (!parser->isError())
		outputStateFile = parser->getOptionString(outputStateFile_Flag, false);

	if (!parser->isError())
		restartFile = parser->getOptionString(restartFile_Flag, true);

	if( !parser->isError() ) {
		setRandomSeed = parser->setOptionLong( randSeed_Flag, randSeed);
	}

	if( !parser->isError() ) {
		parser->setOptionInteger( iterations_Flag, iterations);
		if( iterations <= 0 ) { parser->setError( "Iterations" ); }
	}

	if( !parser->isError() ) {
		parser->setOptionInteger( populationSize_Flag, populationSize);
		if( populationSize <= 0 ) { parser->setError( "Population" ); }
	}

	if( !parser->isError() ) {
		parser->setOptionDouble( mutationRate_Flag, mutationRate);
		if( mutationRate <= 0.0 || mutationRate > 1.0) { parser->setError( "Mutation Rate" ); }
	}

	if( !parser->isError() ) {
		parser->setOptionDouble( complexity_constant_Flag, complexity_constant);
		}
	

	if( !parser->isError() ) {
		parser->setOptionInteger( recombinationFrequency_Flag, recombinationFrequency);
		if( recombinationFrequency <= 0) { parser->setError( "Recombination Frequency" ); }
	}

	

    if( !parser->isError() ) {
		parser->setOptionDouble( recombinationRate_Flag, recombinationRate);
		if( recombinationRate <= 0 || recombinationRate > 1.0) { parser->setError( "Recombination Rate" ); }
	}

	if( !parser->isError() ) {
		parser->setOptionInteger( objectiveFunction_Flag, objectiveFunction);
		switch(objectiveFunction) {
			case OREGA_SIMPLE: case OREGA_COMPLEX: case OREGA_ALL:
				break;
			default: 
				parser->setError( "Objective Function" );
		}
	}

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
int Orega_interface::run() {
	// Create a variable that handles errors.
	int errNum = 0;

	errNum = orega(sequenceFile, outputSequenceFile, 
		outputStateFile, restartFile,
		iterations, alphabet, 
		rangeStart, rangeLength, 
		mutationRate, 
		recombinationFrequency, recombinationRate,
		populationSize, objectiveFunction, 
		randSeed,MutationSwitch,limitG,filterAUG,filterCUG,filteroligoA,complexity_constant,threshold); 
	if (errNum != 0) //Check for errors. 
		std::cerr << "ERROR: " << errNum << endl;

	return errNum;
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {
	int error = 1;
	Orega_interface* runner = new Orega_interface();
	bool parseable = runner->parse( argc, argv );
	if( parseable ) { error = runner->run(); }
	delete runner;
	return error;
}
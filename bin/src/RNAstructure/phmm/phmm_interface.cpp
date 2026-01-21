#include "../phmm/phmm_interface.h"
#include "../RNA_class/TwoRNA.h"
#include "../src/structure.h"
#include <iostream>

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
phmm_interface::phmm_interface() {

	//  Initialize the maximum likelihood variable
	bool ML = false;
	//  Initialize the log probability variable
	bool LP = false;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool phmm_interface::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "phmm" );
	parser->addParameterDescription( "seq 1", "The name of a file containing the first input sequence." );
	parser->addParameterDescription( "seq 2", "The name of a file containing the second input sequence." );
	parser->addParameterDescription( "out file", "The name of a file containing the output sequence." );

	// Add the Maximum Likelihood option.
	vector<string> alignmentOptions;
	alignmentOptions.push_back( "-M" );
	alignmentOptions.push_back( "-m" );
	alignmentOptions.push_back( "--maxlikelihood" );
	parser->addOptionFlagsNoParameters( alignmentOptions, "Specify that program should output a maximum likelihood alignment. Default is to output pairwise probabilities." );

	// Add the log probabilities option.
	vector<string> probOptions;
	probOptions.push_back( "-L" );
	probOptions.push_back( "-l" );
	probOptions.push_back( "--logprobability" );
	parser->addOptionFlagsNoParameters( probOptions, "Specify that program should output probabilities as logs (base 10). Default is to output probabilties." );



	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		seq1 = parser->getParameter( 1 );
		seq2 = parser->getParameter( 2 );
		outfile = parser->getParameter( 3 );
	}

	// Get the maximum likelihood option.
	if( !parser->isError() ) { ML = parser->contains( alignmentOptions ); }

	// Get the log probability option.
	if( !parser->isError() ) { LP = parser->contains( probOptions ); }


	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void phmm_interface::run() {
	std::cout << "Running phmm..\n";
	structure ct1, ct2;
	int l1;
	int l2;
	int error = 0;
	if(!ct1.openseq(seq1.c_str())){
		cerr << "ERROR: Could not open sequence file "<<seq1<<"\n"; 
		error = 1;
	}
	if(!ct2.openseq(seq2.c_str())){
		cerr << "ERROR: Could not open sequence file "<<seq2<<"\n"; 
		error = 1;
	}

	//put nucleotide sequences in vectors for use by phmm class
	vector<char> seq1_nucs;
	vector<char> seq2_nucs;
	for(int i = 1; i <= ct1.GetSequenceLength(); i++){
		seq1_nucs.push_back(ct1.nucs[i]);
	}
	for(int i = 1; i <= ct2.GetSequenceLength(); i++){
		seq2_nucs.push_back(ct2.nucs[i]);
	}

	if(error == 0){
		t_phmm_aln* phmm_aln = create_phmm_aln(seq1_nucs,seq2_nucs);						//t_phmm_aln class contains the two sequences and the hidden markov model
		l1 = phmm_aln->l1();										//get the length of the two sequences, used to write the output file
		l2 = phmm_aln->l2();
		if(!ML) {
			write_probability_array(phmm_aln->compute_posterior_probs(), outfile.c_str(),l1,l2,LP);
			} //pairwise probability calculation (default behavior)
		else {
			write_ML_alignment(phmm_aln->compute_ML_alignment(), outfile.c_str(),l1,l2,seq1.c_str(),seq2.c_str());
			}//maximum likelihood alignment calculation
		std::cout << "All done.\n";
	}
	else std::cerr << "Error initializing sequences"<<std::endl;
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	phmm_interface* runner = new phmm_interface();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}

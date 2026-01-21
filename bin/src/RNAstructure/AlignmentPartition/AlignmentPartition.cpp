/*
 * A program that calculates the partition function for a strand of nucleic acids.
 * This strand of nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Abhinav Mittal
 */

#include "AlignmentPartition.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
//////////////////////////////////////////////////////////////////////////////
AlignmentPartitionInterface::AlignmentPartitionInterface() {

	// Initialize the calculation type description.
	calcType = "Consensus partition function for multiple sequence alignment";

	// Initialize the "experimental" offset.
	experimentalOffset = 0.0;

	// Initialize the "experimental" scaling.
	experimentalScaling = 1.0;

	// Initialize the SHAPE intercept.
	intercept = -0.6;

	// Initialize the nucleic acid type.
	alphabet = DT_MSA;

	// Initialize the maximum pairing distance between nucleotides.
	maxDistance = -1;

	// Initialize the SHAPE slope.
	slope = 1.8;

	// Initialize the calculation temperature.
	temperature = 310.15;

	// Use coaxial stacking by default
	disablecoax = false;
	
	// base pairing cutoff
	bp_cutoff = 0.3;

	// maximum size of internal loop
	maxinter = 40;

}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool AlignmentPartitionInterface::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( PARTITION_EXE_NAME );
	parser->addParameterDescription( "alignment file", "The name of a file containing an multiple sequence alignment in FASTA format." );
	parser->addParameterDescription( "pfs file", "The name of a partition function save file (PFS) to which output will be written." );

	// Add the base pairing cutoff
	vector<string> bpcutoffOptions;
	bpcutoffOptions.push_back("-b");
	bpcutoffOptions.push_back("-B");
	bpcutoffOptions.push_back("--bpcutoff");
	parser->addOptionFlagsWithParameters(bpcutoffOptions, "Specify the base pairing cutoff. Default is 0.3");

	// Add the maximum loop size option.
	vector<string> loopOptions;
	loopOptions.push_back("-l");
	loopOptions.push_back("-L");
	loopOptions.push_back("--loop");
	parser->addOptionFlagsWithParameters(loopOptions, "Specify a maximum internal/bulge loop size. Default is 40 unpaired nucleotides.");

	vector<string> disablecoaxOptions;
	disablecoaxOptions.push_back("--disablecoax");
	parser->addOptionFlagsNoParameters(disablecoaxOptions, "Specify that coaxial stacking recusions should not be used. This option uses a less realistic energy function in exchange for a faster calculation.");

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );



	// Get required parameters from the parser.
	if( !parser->isError() ) {
		seqFile = parser->getParameter( 1 );
		pfsFile = parser->getParameter( 2 );
	}

	// Get the bp_cutoff option.
	if (!parser->isError())
	{
		parser->setOptionFloat(bpcutoffOptions, bp_cutoff);
		if (bp_cutoff < .1)
		{
			parser->setError("bp_cutoff");
		}
	}

	// Get the maxinter option.
	if (!parser->isError())
	{
		parser->setOptionInteger(loopOptions, maxinter);
		if (maxinter < 1)
		{
			parser->setError("maxinter");
		}
	}



	// Get the disablecoax option.
	// true by default, false if the flag is used
	if (!parser->isError()) { disablecoax = parser->contains(disablecoaxOptions); }

//	if (!parser->isError()) { maxinter = parser->contains(loopOptions); }

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void AlignmentPartitionInterface::run() {

	// Create a variable that handles errors.
	int error = 0;

	cout << "Initializing nucleic acids..." << flush;
	RNA* strand = new RNA( seqFile.c_str(), FILE_MSA, alphabet.c_str());
	ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );
	error = checker->isErrorStatus();
	if( error == 0 ) { cout << "done." << endl; }
	
	/*
	 * Calculate the partition function using the PartitionFunction method.
	 * During calculation, monitor progress using the TProgressDialog class and the Start/StopProgress methods of the RNA class.
	 * Neither of these methods require any error checking.
	 * After the main calculation is complete, use the error checker's isErrorStatus method to check for errors.
	 */
	if( error == 0 ) {

		// Show a message saying that the main calculation has started.
		cout << "Calculating partition function..." << flush;
		strand->GetStructure()->SetBasePairingCutoff(bp_cutoff);
		strand->GetStructure()->FillPairingMatrix();
		strand->GetStructure()->FillGapMatrix();


		strand->GetStructure()->tmp_ct_hp = new structure;
		strand->GetStructure()->tmp_ct_hp->SetThermodynamicDataTable(strand->GetStructure()->GetThermodynamicDataTable());

		strand->GetStructure()->tmp_ct_i = new structure;
		strand->GetStructure()->tmp_ct_i->SetThermodynamicDataTable(strand->GetStructure()->GetThermodynamicDataTable());

		// Do the main calculation and check for errors.
		int mainCalcError = strand->PartitionFunction( pfsFile.c_str(), temperature, disablecoax, true, false, maxinter);
		error = checker->isErrorStatus( mainCalcError );
		
		// If no error occurred, print a message saying that the main calculation is done.
		if( error == 0 ) { cout << "done." << endl; }
	}

	for (int i = 0; i < strand->GetStructure()->GetNumberofSequences(); i++)
	{
		delete strand->GetStructure()->sequences_from_alignment[i].sequence_without_gaps;
	}
	delete strand->GetStructure()->tmp_ct_hp;
	delete strand->GetStructure()->tmp_ct_i;

	// Delete the error checker and data structure.
	delete checker;
	delete strand;

	// Print confirmation of run finishing.
	if( error == 0 ) { cout << calcType << " complete." << endl; }
	else { cerr << calcType << " complete with errors." << endl; }
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	AlignmentPartitionInterface* runner = new AlignmentPartitionInterface();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}

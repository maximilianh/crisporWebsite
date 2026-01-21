/*
 * A program that determines consensus structure from multiple RNA sequence alignment.
 *
 * (c) 2020 Mathews Lab, University of Rochester Medical Center.
 * Written by Abhinav Mittal
 */


#include "AlignmentFold.h"
#include "../src/structure.h"
#include "../src/rna_library.h"


//! constructor
AlignmentFold_Interface::AlignmentFold_Interface()
{
	alphabet = DT_MSA;
	temperature = 310.15;
	bp_cutoff = 0.3;
	quiet = false;
	useBracketNotation = false;
	maxLoop = 40;
	disablecoax = true;
	quickfold = false;
}

bool AlignmentFold_Interface::run()
{	
	if (!quiet) cout << "Initializing nucleic acids..." << flush;
	RNA *strand = new RNA(seqFile.c_str(), FILE_MSA, alphabet.c_str());
	strand->GetStructure()->SetBasePairingCutoff(bp_cutoff);
	
	strand->GetStructure()->tmp_ct_hp = new structure;
	strand->GetStructure()->tmp_ct_hp->SetThermodynamicDataTable(strand->GetStructure()->GetThermodynamicDataTable());

	strand->GetStructure()->tmp_ct_i = new structure;
	strand->GetStructure()->tmp_ct_i->SetThermodynamicDataTable(strand->GetStructure()->GetThermodynamicDataTable());

	strand->GetStructure()->FillPairingMatrix();
	strand->GetStructure()->FillGapMatrix();
	strand->FoldSingleStrand(20, 20, 5, "", maxLoop, quickfold, true, disablecoax);
	if (useBracketNotation) strand->WriteDotBracket(ctFile.c_str());
	else strand->WriteCt(ctFile.c_str());

	for (int i = 0; i < strand->GetStructure()->GetNumberofSequences(); i++)
	{
		delete strand->GetStructure()->sequences_from_alignment[i].sequence_without_gaps;
	}
	delete strand->GetStructure()->tmp_ct_hp;
	delete strand->GetStructure()->tmp_ct_i;
	delete strand;
	return true;
}

bool AlignmentFold_Interface::parse(int argc, char** argv) 
{	
	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine("AlignmentFold");
	parser->addParameterDescription("sequence file", "The name of a file containing multiple sequence alignment of homologous RNA sequences in FASTA format.");
	parser->addParameterDescription("CT file", "The name of a CT file to which output will be written. If the --bracket flag is present, output will be written as a dot-bracket file.\nIf the file name is a hyphen (-), the structure will be written to standard output (STDOUT) instead of a file.");

	// Add the Alphabet option.
//	vector<string> alphabetOptions;
//	alphabetOptions.push_back("-a");
//	alphabetOptions.push_back("--alphabet");
//	parser->addOptionFlagsWithParameters(alphabetOptions, "Specify the name of a folding alphabet and associated nearest neighbor parameters. The alphabet is the prefix for the thermodynamic parameter files, e.g. \"rna\" for RNA parameters or \"dna\" for DNA parameters or a custom extended/modified alphabet. The thermodynamic parameters need to reside in the at the location indicated by environment variable DATAPATH. The default is \"rna\" (i.e. use RNA parameters). This option overrides the --DNA flag.");
	
	// Add the base pairing cutoff
	vector<string> bpcutoffOptions;
	bpcutoffOptions.push_back("-b");
	bpcutoffOptions.push_back("-B");
	bpcutoffOptions.push_back("--bpcutoff");
	parser->addOptionFlagsWithParameters(bpcutoffOptions, "Specify the base pairing cutoff. Default is 0.3");

	// Add the quickfold option.
	vector<string> quickfoldOptions;
	quickfoldOptions.push_back("-mfe");
	quickfoldOptions.push_back("-MFE");
	quickfoldOptions.push_back("--MFE");
	parser->addOptionFlagsNoParameters(quickfoldOptions, "Specify that only the minimum free energy structure is needed.  No savefiles can be generated.  This saves nearly half the calculation time, but provides less information.");

	// Add the maximum loop size option.
	vector<string> loopOptions;
	loopOptions.push_back("-l");
	loopOptions.push_back("-L");
	loopOptions.push_back("--loop");
	parser->addOptionFlagsWithParameters(loopOptions, "Specify a maximum internal/bulge loop size. Default is 40 unpaired nucleotides.");

	// Add the disablecoax option.
	vector<string> disablecoaxOptions;
	disablecoaxOptions.push_back("--disablecoax");
	parser->addOptionFlagsNoParameters(disablecoaxOptions, "Specify that coaxial stacking recusions should not be used. This option uses a less realistic energy function in exchange for a faster calculation.");

	vector<string> quietOption = parser->addFlag(false, "-q --quiet", "Suppress unnecessary output. This option is implied when the output file is '-' (STDOUT).");
	vector<string> bracketOption = parser->addFlag(false, "-k --bracket", "Write the predicted structure in dot-bracket notation (DBN) instead of CT format.");

	// Parse the command line into pieces.
	parser->parseLine(argc, argv);

	// Get required parameters from the parser.
	if (!parser->isError()) 
	{
		seqFile = parser->getParameter(1);
		ctFile = parser->getParameter(2);
	}

	quiet = parser->contains(quietOption) || isStdIoFile(ctFile.c_str()); // suppress unnecessary output if --quiet option is present or if the CT file is output to stdout.
	useBracketNotation = parser->contains(bracketOption); // write in dot-bracket notation.

	// Get the Alphabet option.
//	if (!parser->isError() && parser->contains(alphabetOptions))
//		alphabet = parser->getOptionString(alphabetOptions, false).c_str();

	
	// Get the bp_cutoff option.
	if (!parser->isError())
	{
		parser->setOptionFloat(bpcutoffOptions, bp_cutoff);
		if (bp_cutoff < .1)
		{
			parser->setError("bp_cutoff");
		}
	}

	// Get the maximum loop size option.
	if (!parser->isError()) {
		parser->setOptionInteger(loopOptions, maxLoop);
		if (maxLoop < 0) { parser->setError("maximum loop size"); }
	}

	// Get the quickfold option.
	if (!parser->isError()) { quickfold = parser->contains(quickfoldOptions); }

	// Get the disablecoax option.
	// true by default, false if the flag is used
	if (!parser->isError()) { disablecoax = parser->contains(disablecoaxOptions); }


	// Delete the parser and return whether the parser encountered an error.
	bool noError = (parser->isError() == false);
	delete parser;
	return noError;
}

int main(int argc, char* argv[]) 
{
	AlignmentFold_Interface runner;
	if (!runner.parse(argc, argv)) return 1;
	return runner.run() ? 0 : 1;
}

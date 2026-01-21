#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "../RNA_class/RNA.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

class AlignmentFold_Interface {
public:
	AlignmentFold_Interface();
	bool parse(int argc, char** argv);
	bool run();
	float bp_cutoff;		 // minimum fraction base pairing column

private:
	string seqFile;			 // input sequence file
	string ctFile;			 // output ct file
	string alphabet;		 // nucleic acid type (rna, dna or custom)
	int maxLoop;			 // the maximum internal bulge loop size.
	double temperature;		 // temperature at which calculation occurs.
//	float bp_cutoff;		 // minimum fraction base pairing column
	bool quiet;				 // suppress unnecessary output
	bool useBracketNotation; // write structure in dot-bracket notation.
	bool disablecoax;		 // Initialize the flag to disable coaxial stacking recursions
	bool quickfold;			 // Initialize the quickfold (mfe only) variable.
};


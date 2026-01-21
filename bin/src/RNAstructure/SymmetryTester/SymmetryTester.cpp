/*
* A program that checks the required symmetries of RNAstructure data tables.
*
* (c)2021 Mathews Lab, University of Rochester Medical Center.
* Written by David H. Mathews
*/


#include "ST.h"
///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
ST_Interface::ST_Interface() {
	// Initialize the nucleic acid type.
	alphabet = DT_RNA;
	
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool ST_Interface::parse(int argc, char** argv) {

	ParseCommandLine* parser = new ParseCommandLine("SymmetryTester");
	

	// Add the Alphabet option.
	vector<string> alphabetOptions;
	alphabetOptions.push_back("-a");
	alphabetOptions.push_back("--alphabet");
	parser->addOptionFlagsWithParameters(alphabetOptions, "Specify the name of a folding alphabet and associated nearest neighbor parameters. The alphabet is the prefix for the thermodynamic parameter files, e.g. \"rna\" for RNA parameters or \"dna\" for DNA parameters or a custom extended/modified alphabet. The thermodynamic parameters need to reside in the at the location indicated by environment variable DATAPATH. The default is \"rna\" (i.e. use RNA parameters). This option overrides the --DNA flag.");

	// Parse the command line into pieces.
	parser->parseLine(argc, argv);

	

	// Get the Alphabet option.
	if (!parser->isError() && parser->contains(alphabetOptions))
		alphabet = parser->getOptionString(alphabetOptions, false).c_str();

	// Delete the parser and return whether the parser encountered an error.
	bool noError = (parser->isError() == false);
	delete parser;
	return noError;
}


///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
bool ST_Interface::run() {

	// Create a variable that handles errors.
	int error = 0;
	char numstr[21];
	

	/*
	 * Use the constructor for RNA that specifies a filename.
	 * Specify type = 2 (sequence file).
	 * isRNA identifies whether the strand is RNA (true) or DNA (false).
	 *
	 * After construction of the strand data structure, create the error checker which monitors for errors.
	 * Throughout, the error status of the calculation is checked with a variant of the isErrorStatus method, which returns 0 if no error occurred.
	 * The calculation proceeds as long as error = 0.
	 */
	
	RNA* strand = new RNA("AAAAAAAA", SEQUENCE_STRING, alphabet.c_str());
	ErrorChecker<RNA>* checker = new ErrorChecker<RNA>(strand);
	error = checker->isErrorStatus();
	
	

	if (error == 0) {
		//Now check the symmetries in the tables
		datatable *tables = strand->GetDatatable();

		cout << "Checking Helix Stack Table\n\n";

		//check the base pair stack table:
		for (int i = 0; i < tables->alphabet.size(); ++i) {
			for (int j = 0; j < tables->alphabet.size(); ++j) {
				for (int ip = 0; ip < tables->alphabet.size(); ++ip) {
					for (int jp = 0; jp < tables->alphabet.size(); ++jp) {

						if (tables->stack[i][j][ip][jp] != tables->stack[jp][ip][j][i]) {

							cerr << "Problem in .stack.dg " << tables->alphabet[i][0] << " " << tables->alphabet[j][0] << " " << tables->alphabet[ip][0] << " " << tables->alphabet[jp][0] << "\n";

						}

					}
				}
			}
		}

		cout << "Checking 1x1 Internal Loop Table\n\n";

		//now check the 1x1 internal loops
		//check the base pair stack table:
		for (int i = 0; i < tables->alphabet.size(); ++i) {
			for (int j = 0; j < tables->alphabet.size(); ++j) {
				for (int ip = 0; ip < tables->alphabet.size(); ++ip) {
					for (int jp = 0; jp < tables->alphabet.size(); ++jp) {
						for (int a = 0; a < tables->alphabet.size(); ++a) {
							for (int b = 0; b < tables->alphabet.size(); ++b) {

								if (tables->iloop11[i][a][ip][j][b][jp] != tables->iloop11[jp][b][j][ip][a][i]) {

									cerr << "Problem in .int11.dg " << tables->alphabet[i][0] << " " << tables->alphabet[j][0] << " " << tables->alphabet[ip][0] << " " << tables->alphabet[jp][0] << " " << tables->alphabet[a][0] << " " << tables->alphabet[b][0] << "\n";
									cerr << tables->alphabet[i][0] << " " << tables->alphabet[a][0] << " " << tables->alphabet[ip][0] << " "<<tables->iloop11[i][a][ip][j][b][jp] << "\n";
									cerr << tables->alphabet[j][0] << " " << tables->alphabet[b][0] << " " << tables->alphabet[jp][0] << " "<< tables->iloop11[jp][b][j][ip][a][i] << "\n";
								}

							}
						}

					}
				}
			}
		}

		cout << "Checking 2x2 Internal Loop Table\n\n";

		//now check the 2x2 internal loops
		//check the base pair stack table:
		for (int i = 0; i < tables->alphabet.size(); ++i) {
			for (int j = 0; j < tables->alphabet.size(); ++j) {
				for (int ip = 0; ip < tables->alphabet.size(); ++ip) {
					for (int jp = 0; jp < tables->alphabet.size(); ++jp) {
						for (int a = 0; a < tables->alphabet.size(); ++a) {
							for (int b = 0; b < tables->alphabet.size(); ++b) {
								for (int c = 0; c < tables->alphabet.size(); ++c) {
									for (int d = 0; d < tables->alphabet.size(); ++d) {

										if (tables->iloop22[i][ip][j][jp][a][b][c][d] != tables->iloop22[jp][j][ip][i][d][c][b][a]) {

											cerr << "Problem in .int22.dg " << tables->alphabet[i][0] << " " << tables->alphabet[j][0] << " " << tables->alphabet[ip][0] << " " << tables->alphabet[jp][0] << " " << tables->alphabet[a][0] << " " << tables->alphabet[b][0] << " " << tables->alphabet[c][0] << " " << tables->alphabet[d][0] << "\n";
											cerr << tables->alphabet[i][0] << " " << tables->alphabet[a][0] << " " << tables->alphabet[b][0] << " "<< tables->alphabet[ip][0] << " "<< tables->iloop22[i][ip][j][jp][a][b][c][d]<<"\n";
											cerr << tables->alphabet[j][0] << " " << tables->alphabet[c][0] << " " << tables->alphabet[d][0] << " "<< tables->alphabet[jp][0] << " "<< tables->iloop22[jp][j][ip][i][d][c][b][a]<<"\n";
										}

									}
								}
							}
						}

					}
				}
			}
		}



	}

	
		
	delete strand;

	// Print confirmation of run finishing.
	if (error == 0) { cout <<  " complete." << endl; }
	else { cerr <<  " complete with errors." << endl; }

	return error == 0;
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {
	ST_Interface runner;
	if (!runner.parse(argc, argv)) return 1;
	return runner.run() ? 0 : 1;
}
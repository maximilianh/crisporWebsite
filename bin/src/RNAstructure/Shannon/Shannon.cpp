/*
 * A program that calculates the Shannon entropy of each nucleotide in a sequence.
 
 *
 * (c) 2021 Mathews Lab, University of Rochester Medical Center.
 * Written by David H. Mathews
 */

#include "Shannon.h"

 ///////////////////////////////////////////////////////////////////////////////
 // Constructor.
 ///////////////////////////////////////////////////////////////////////////////
Shannon::Shannon() {


}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool Shannon::parse(int argc, char* argv[]) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine("Shannon");
	parser->addParameterDescription("input file", "The name of the input file that holds base pairing probabilities. This file is a Partition function save file (binary file created by partition).");
	parser->addParameterDescription("output file", "The name of a text file to which output will be written.");


	// Parse the command line into pieces.
	parser->parseLine(argc, argv);

	// Get required parameters from the parser.
	if (!parser->isError()) {
		inputFile = parser->getParameter(1);
		outputFile = parser->getParameter(2);
	}





	
	// Delete the parser and return whether the parser encountered an error.
	bool noError = (parser->isError() == false);
	delete parser;
	return noError;
}


///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void Shannon::run() {

	/*
	 * Read in dot plot data.
	 */
	string error = "";
	 // Print a message that says the dot plot data is being read.
	cout << "Reading partition function data..." << flush;

	// Create the RNA strand.
	RNA* strand = new RNA(inputFile.c_str(), FILE_PFS);
	ErrorChecker<RNA>* checker = new ErrorChecker<RNA>(strand);

	// If no error occurred creating the RNA strand, read data from it.
	if ((error = checker->returnError()) == "") {

		
		int length = strand->GetSequenceLength();
		vector<double> entropies(length+1);
		// Calculate entropies for each nucleotide.
		for (int i = 1; i <= length; i++) {
			entropies[i] = 0;
			for (int j = 1; j < i; ++j) {
				double value = strand->GetPairProbability(j, i);
				if (value > 0.0) entropies[i] -= value * log10(value);
			}
			for (int j = i + 1; j <= length; j++) {
				double value = strand->GetPairProbability(i, j);
				if (value > 0.0) entropies[i] -= value * log10(value);
				
			}
		}

		// If no errors occurred, print a message that says the dot plot data was read.
		if (error == "") { cout << "done." << endl; }



		/*
		 * Write the dot plot file.
		 */

		if (error == "") {

			// Print a message saying that the dot plot file is being written.
			cout << "Writing Output File..." << flush;
			ofstream outfile;
			outfile.open(outputFile.c_str());
			outfile << "i\tShannon Entropy(i)\n";
			for (int i = 1; i <= length; i++) {
				outfile << i << "\t" << entropies[i] << "\n";

			}


			outfile.close();

		}
	}

	// Delete the error checker and RNA strand.
	delete checker;
	delete strand;
	

	

	/*
	 * Clean up resources and show the result of the dot plot run.
	 */


	// Print confirmation of run finishing.
	if (error == "") { cout << "complete." << endl; }
	else { cerr << error << endl << " complete with errors." << endl; }
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {
	Shannon* runner = new Shannon();
	bool parseable = runner->parse(argc, argv);
	if (parseable == true) { runner->run(); }
	delete runner;
}

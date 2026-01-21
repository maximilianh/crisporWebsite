/*
 * A program that calculates the distribution and mean end-to-end (ETE) distance, 
 * of a given sequence in nanometers (nm) based on a polymer model from Aalberts et al.
 *
 * (c) 2019 Mathews Lab, University of Rochester Medical Center.
 * Written by Mohammad Kayedkhordeh, Stanislav Bellaousov
 */

#ifndef ETECALCULATOR_INTERFACE_H
#define ETECALCULATOR_INTERFACE_H

#include <vector>
#include <string>
#include "../RNA_class/RNA.h"

class ETEcalculator_Interface {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	ETEcalculator_Interface();

	/*
	 * Name:        parse
	 * Description: Parses command line arguments to determine what options are required for a particular calculation.
	 * Arguments:
	 *     1.   The number of command line arguments.
	 *     2.   The command line arguments themselves.
	 * Returns:
	 *     True if parsing completed without errors, false if not.
	 */
	bool parse( int argc, char** argv );

	/*
	 * Name:        run
	 * Description: Run calculations.
	 */
	int run();

	/*
	 * Name:        etecalculator    
	 * Description: Calculates the distribution and average end-to-end distance values for a sequence using Aalberts model.
	 * Arguments:
	 *     1.   Input RNA used in the calculations
	 *     2.   The output filestream that will be used (if provided in the options) to store the results. 
	 * Returns:
	 *     Nothing (Void return).

	*/
	void etecalculator(RNA* rna, string outputFileName);

 private:
	// Private variables.

	// The input ct file which describes the structure to design. The sequence in the ct file will be ignored.
	string seqFile;           

	// Name of the alphabet (e.g. "rna" or "dna" etc).
	const char* alphabet;	

	int ensembleSize;

	int seed;

	bool raw;

	//indicate if the program should read a ct file instead of a sequence
	bool load_ct;

	//Coefficient value for single stranded (unpaired) nucleotide from Aalberts model used in calculating the end-to-end distance;
	double ss_dist;

	//Coefficient value for a branch in the exterior loop from Aalberts model used in calculating end-to-end distance
	double paired_dist;

	string constraintFile;

	string outputFile;
};

#endif /* ETECALCULATOR_INTERFACE_H */

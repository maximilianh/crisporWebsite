/*
 * NAPSS, a program that predicts RNA secondary structures with pseudoknots with the aid of NMR constraints data.
 *
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 * Written by Stanislav Bellaousov
 */

#ifndef NAPSS_Interface_H
#define NAPSS_Interface_H

#include "napss.h"

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cstring>
#include <cstdlib>

#include "../RNA_class/RNA.h"
#include "../src/PseudoParser.h"
#include "../src/ParseCommandLine.h"

//#include "stdafx.h"
//#include "../RNA_class/RNA.h"
//#include "../src/configfile.h"
//#include "../src/defines.h"
//#include "../src/platform.h"
//#include "../src/rna_library.h"
//#include "../src/structure.h"
//#include "../src/algorithm.h"
//#include "../src/ParseCommandLine.h"

//using namespace std;

class NAPSS_Interface {
 public:
	// Public constructor and methods.
	
	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	NAPSS_Interface();

	/*
	 * Name:        parse
	 * Description: Parses command line arguments to determine what options are required for a particular calculation.
	 * Arguments:
	 *     1.   The number of command line arguments.
	 *     2.   The command line arguments themselves.
	 * Returns:
	 *     True if parsing completed without errors, false if not.
	 */
	bool Parse( int argc, char** argv );

	/*
	 * Name:        run
	 * Description: Run calculations.
	 */
	bool run();

 private:
	//Define private variables.

    double P1, P2;

    string inseq;//address for the input sequence file

    string inNMRconstraints;//address for the NMR constraints file

    string outct;//address for the output ct file

    string inSHAPEfile;//address for the input SHAPE file

    string outpairs;//address for the output file with pairs

    string constraintFile;//address for the structural constraints file

    double slope,intercept;//shape slope and shape intercept

    int maxtracebacks,percent,windowsize,cutoff;

    bool pseudoknotFree;//default is to have pseudoknot prediction mode enabled

    int warningLimit;//number of matches (matchVector->size) before the warning message is printed

    bool ifwarningMessage;//keeps track if the warning message was printed. false - warning message was not printed yet, so should be printed when warning limit is reached; true - warning message has already been printed, so no need to print it again.

	bool ifwindowOptions;//Tracks if the windowsize flag was specified by the user.

};

#endif /* NAPSS_Interface_H */

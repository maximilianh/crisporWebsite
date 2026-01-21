/*
 * NAPSS, a program that predicts RNA secondary structures with pseudoknots with the aid of NMR constraints data.
 *
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 * Written by Stanislav Bellaousov
 */

#ifndef NAPSS_H
#define NAPSS_H


#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cstring>
#include <cstdlib>

#include "../RNA_class/RNA.h"
#include "../src/PseudoParser.h"

using namespace std;

typedef vector<short> conMatch;

#define NAPSS_BAD_CONSTRAINT_FILE 30001
#define NAPSS_ERR_CONSTRAINT_LEN 30002
#define NAPSS_BAD_RNA_COPY 30003
#define NAPSS_ERR_NO_MATCHES 30004


void firstDotMatch(short**,bool*,bool*,short,short,short**,short**,vector<conMatch>*,short*,short*,short**,short**,int,short, structure *ct,char**** tripletArray, int warningLimit, bool* ifwarningMessage);
// Loops through x,y of entire dotplot searching for matches to the first basepair of each constraint

void recursiveMatch(short**,bool*,short,short,short**,short**,vector<conMatch>*,short*,short*,short**,short**,int,short,bool*,structure *ct,char**** tripletArray, int warningLimit, bool* ifwarningMessage);
// Takes over for firstDotMatch to find the matches to the remaining basepairs in a constraint

void pseudodptrim(structure*,int*, int*); // Searches dotplot for complicated pseudoknot folds to exclude

int pseudoenergy(datatable *data, structure *ct, int count, bool *nopseudos, structure *ctbroken, int *init, int *brokenpairs);
// Searches ct for pseudoknots, returns energy for breaking apart the pseudoknot (including Beta penalties)

int filterByEnergy(structure *ct, double energyCutoff);

int pairout(structure *ct, const char* pairsoutfile); // Outputs PseudoViewer3 structural data

bool TripletMatch(structure *ct, char**** tripletArray, short** xCoords, short** yCoords, int currConNum, int currConPos);//Match triplet constraints

bool NucCompare(char ConNuc, int SeqNuc, char ConSign=0);//Given a constraint and a nucleotide the program returns 'true' if there is a match and 'false' if not

void HelicalExtension(structure *ct, short** convertedDotplot, int* helixExtended, short** dgArray);//Check for helical extension possibilities

//! Main napss function
int napss(RNA *rnaCT, const string inNMRconstraints, 
		/* pass in a NULL or uninitialized structure* it will be set to the address of the output structure on success. The caller must then delete the returned structure. */
		structure* &outputStructureResults, 
		const int windowsize, int cutoff, const int percent, 
		const int maxtracebacks, const int warningLimit,
		double pseud_param[16], const double P1, const double P2, const bool pseudoknotFree,
		const string inSHAPEfile, const double slope, const double intercept,
		const string constraintFile,
		const bool useSimpleMBRules=false);

/**
  * Get an error message from an integer value returned from napss()
*/
string napss_error_message(int errorcode);

#endif
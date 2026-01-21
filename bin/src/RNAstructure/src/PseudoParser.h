/*
 * Pseudoparser functions that are used by the pkEnergyCalc cunction to
 * calculate the pseudoenergy of the pseudoknot. Used in ShapeKnots.
 *
 * (c) 2013 Stanislav Bellaousov
 * Mathews Lab, University of Rochester Medical Center
 */

#ifndef PSEUDOPARSER_H
#define PSEUDOPARSER_H

#include "Pseudoknot.h"
#include "../RNA_class/RNA.h"

#define DEFAULT_PSEUDOKNOT_P1 0.35 //pseudoknot energy model parameters (kcal/mol).
#define DEFAULT_PSEUDOKNOT_P2 0.65 //pseudoknot energy model parameters (kcal/mol).

//Function 'FillMiemsatch' fills in single or tandem mismatches in rna structure 'rna'
void FillMismatch(RNA *rna, int strnum=1); 

//Function 'RemoveIsolatedPairs' removes isolated pairs in rna structure 'rna'
void RemoveIsolatedPairs(RNA *rna, int strnum=1);

//Function 'RNAcopy' creates a vector 'rnacopy' that stores a copy of all pairs from RNA structure 'rna'
vector<int> RNAcopy(RNA *rna, int strnum=1);

//Function 'BreakPseudoknot' breakes pseudoknot in a given rna from class 'RNA'.
void BreakPseudoknot(RNA *rna, int strnum=1);

//Function 'paursBroken' makes a vector that holds only the broken pairs.
vector<int> PairsBroken(RNA *rna, vector<int> &rnacopy, int strnum=1);

//Function 'PseudoParser' parses through a structure from 'RNA' class and creates a vector of 
//class 'Pseudoknot' that stores pairing values for all pseudoknotted pairs in the structure
vector<Pseudoknot> PseudoParser(RNA *rna, vector<int> &pairsBroken, int strnum=1);

//Function 'DistanceCounter' calculates the number of unpaired nucleotides, number of internal branches, 
//and length and type of spanning helixes for each pseudoknot from class 'Pseudoknot'
void DistanceCounter(RNA *rna, vector<Pseudoknot> &pk, vector<int> &rnacopy, vector<int> &UNcount, vector<int> &IBcount, vector<vector<int> > &SHLcount, int strnum=1);

//! Function 'ReadPseudoParam' reads the pseudoknot parameters
//! \return 0 on success or an error code compatible with RNA::GetErrorMessage if the file could not be read.
int ReadPseudoParam(double *pseudo_param);

//Function 'RestorePseudoknot' restores pairs that were broken when the pseudoknot was broken. 
//'pairsBroken' contains the pairs that were broken when using function BreakPseudoknot
void RestorePseudoknot(RNA *rna, vector<int> &pairsBroken, int strnum=1);

//Function 'RestoreStructure' restores any changes made to the structure. 
//'pairsChanged' contains the copy of the structure
void RestoreStructure(RNA *rna, vector<int> &pairsChanged, int strnum=1);

//Function 'RetToScreen' returns results to the screen
void RetToScreen(vector<Pseudoknot> &pk, vector<int> &UNcount, vector<int> &IBcount, vector<vector<int> > &SHLcount);

//Function that checks if two nucleotides are canonical
//Returns 'true' if canonical, 'false' if non-canonical
//THIS FUNCTION IS LOCATED IN RNAstructure/src/MaxExpect.cpp
bool isCanonical(char i, char j);

//Function 'pkEnergyCalc' calculated the energy penalty for having a 
//a pseudoknot.
double pkEnergyCalc(RNA *rna, double *a, double P1, double P2, int strnum=1);

//Function calculates the energy of the structure broken during the removal
//of the pseudoknot. The function adds the energy of stacks and  bulges/internal loops.
double CalculatePKhelixEnergy(RNA* ct, int strnum, datatable* data);
double CalculatePKhelixEnergy(structure* ct, int strnum, datatable* data);

//Function calculates the energy of the pseudoknotted structure using the ShapeKnots energy model
//Input:
//1. RNA rna - rna class that holds the structure
//2. strnum - structure number 
//3 and 4. P1 and P2 - Pseudoknot energy model parameters (kcal/mol)
//5. pseud_param - Experimental pseudoknot distances
//6. data - datatables
//7. UseSimpleMBLoopRules - input for the RNA->CalculateFreeEnergy to use simplified rules
//8. fillMismatch - indicates if filling of single and double mismatches should be performed when calculating the free energy
//9. remIsolatedPairs - indicates if isolated pairs should be removed
double CalculateFreeEnergywithPseudoknot(RNA* rna, int strnum, double P1, double P2, double* pseud_param, datatable* data, bool UseSimpleMBLoopRules, bool fillMismatch, bool remIsolatedPairs, ofstream* writeDetails = NULL);

#endif
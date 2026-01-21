/*
 * Pseudoknot Class
 * Pseudoknot class stores and accesses pseudoknotted pairs.
 * Files with Pseudoknot class are: basepair.cpp, basepair.h, Pseudoknot.cpp, Pseudoknot.h
 *
 * (c) 2013 
 * Mathews Lab, University of Rochester Medical Center
 */

#ifndef PSEUDOKNOT_H
#define PSEUDOKNOT_H

#include "basepair.h"
#include <vector>
#include <algorithm>

using namespace std;

//Class to keep track of pairs in a pseudoknot
//Pairs are in lists called broken and intact
//where broken are those pairs removed by a breakpseudoknot routine
//and intact are the crossing pairs to those from broken
class Pseudoknot {

 public:
	//Empty constructor
	Pseudoknot();
		
	//Return the number of broken pairs
	int GetNumberBroken();

	//return the number of intact pairs
	int GetNumberIntact();

	//Get the broken basepair at index i
	basepair GetBroken(const int i);

	//Get the intact basepair at index i
	basepair GetIntact(const int i);

	//Set a new broken pair, with i paired to j and i<j
	//return true when it works, return false if i > j
	bool SetBroken(const int i, const int j);

	//Set a new intact pair, with i paired to j and i<j
	//return true when it works, return false if i > j
	bool SetIntact(const int i, const int j);

	//Sort nucleotides
	void Sort();

 private:
	//store the broken pairs
	vector<basepair> broken;

	//store the intact pairs
	vector<basepair> intact;




};

#endif
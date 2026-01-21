/*
 * Pseudoknot Class
 * Pseudoknot class stores and accesses pseudoknotted pairs.
 * Files with Pseudoknot class are: basepair.cpp, basepair.h, Pseudoknot.cpp, Pseudoknot.h
 *
 * (c) 2013 
 * Mathews Lab, University of Rochester Medical Center
 */

#ifndef BASEPAIR_H
#define BASEPAIR_H

//Class to define a basepair.  It keeps track of 5' and 3' nucs in pair.
class basepair {


 public:

	//an empty constructor
	basepair();

	//a constructor that takes pairs
	//i must be 5' to j
	basepair(const int i, const int j);

	//Get the index of the 5' nuc
	int Get5pNucleotide();

	//Get the index of the 3' nuc
	int Get3pNucleotide();

	//Set the index to the 5' nuc
	void Set5pNucleotide(const int nuc);

	//Set the index to the 3' nuc
	void Set3pNucleotide(const int nuc);


 private:
	//5' and 3' nucs
	int p5,p3;

};


//set up comparison of two basepairs
bool operator==(basepair a, basepair b);

//define <
//a is less than b if the 5' nuc of a < 5' nuc of b
bool operator<(basepair a, basepair b);


//set up == on int to basepair
//i is equal to a basepair if either the 5' or 3' nuc is equal to i
bool operator==(int a, basepair b);
bool operator==(basepair b, int a);

#endif

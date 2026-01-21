/*
 * Pseudoknot Class
 * Pseudoknot class stores and accesses pseudoknotted pairs.
 * Files with Pseudoknot class are: basepair.cpp, basepair.h, Pseudoknot.cpp, Pseudoknot.h
 *
 * (c) 2013 
 * Mathews Lab, University of Rochester Medical Center
 */

#include "basepair.h"


basepair::basepair() {


}

basepair::basepair(const int i, const int j) {

	p5 = i;
	p3 = j;

}
		
int basepair::Get5pNucleotide() {

	return p5;

}
	
int basepair::Get3pNucleotide() {

	return p3;

}

void basepair::Set5pNucleotide(const int nuc) {

	p5 = nuc;

}
		
void basepair::Set3pNucleotide(const int nuc) {

	p3 = nuc;

}

//set up comparison of two basepairs
bool operator==(basepair a, basepair b) {

	return (a.Get5pNucleotide() == b.Get5pNucleotide());

}

//define <
//a is less than b if the 5' nuc of a < 5' nuc of b
bool operator<(basepair a, basepair b) {

	return (a.Get5pNucleotide() < b.Get5pNucleotide());

}


//set up == on int to basepair
//i is equal to a basepair if either the 5' or 3' nuc is equal to i
bool operator==(int a, basepair b) {

	return ( (a == b.Get5pNucleotide())	 || (a == b.Get3pNucleotide()));

}


//set up == on int to basepair
//i is equal to a basepair if either the 5' or 3' nuc is equal to i
bool operator==(basepair b, int a) {

	return ( (a == b.Get3pNucleotide()) || (a == b.Get5pNucleotide()));

}

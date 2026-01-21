/*
 * Pseudoknot Class
 * Pseudoknot class stores and accesses pseudoknotted pairs.
 * Files with Pseudoknot class are: basepair.cpp, basepair.h, Pseudoknot.cpp, Pseudoknot.h
 *
 * (c) 2013 
 * Mathews Lab, University of Rochester Medical Center
 */

#include "Pseudoknot.h"

//Empty constructor
Pseudoknot::Pseudoknot() {


}
		
//Return the number of broken pairs
int Pseudoknot::GetNumberBroken() {

	return (int) broken.size();

}

//return the number of intact pairs
int Pseudoknot::GetNumberIntact() {
	
	return (int) intact.size();

}

//Get the broken basepair at index i
basepair Pseudoknot::GetBroken(const int i) {

	return broken[i];

}

//Get the intact basepair at index i
basepair Pseudoknot::GetIntact(const int i) {

	return intact[i];

}

//Set a new broken pair, with i paired to j and i<j
//return true when it works, return false if i > j
bool Pseudoknot::SetBroken(const int i, const int j) {

	if (i>j) return false;
	broken.push_back(basepair(i,j));
	return true;

}

//Set a new intact pair, with i paired to j and i<j
//return true when it works, return false if i > j
bool Pseudoknot::SetIntact(const int i, const int j) {

	if (i>j) return false;
	intact.push_back(basepair(i,j));
	return true;

}

//Sort nucleotides 
void Pseudoknot::Sort(){
	sort(broken.begin(),broken.end());
	sort(intact.begin(),intact.end());
}




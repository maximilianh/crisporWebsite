/*class to encapsulate functionality of a helix which could be involved in a pseudoknot
  The helices can not handle bulges at this point
  This class maintains a list that is a helix in the form 1,2,3,4,5,6,7,8
  where 1 is base-paired to 2, 3 is base paired to 4, etc
  and the 5' half of the helix is 1, 3, 5, 7 (increasing numbers) and the 3' half of the helix is 2, 4, 6, 8 (decreasing numbers)
  This class also maintains the length of the helix and the energy of the helix.  The energy of the helix is calculated as the sum 
  of the au-end penalty plus the SHAPE enhanced stacking energies between base pairs
*/

#if !defined(PKHELIX_H)
#define PKHELIX_H


#include "../RNA_class/RNA.h"
#include "../src/platform.h"
#include "../src/rna_library.h"
#include "../src/structure.h"
#include "../src/algorithm.h"


#include <iomanip>
#include <iostream>

class pkHelix{

 public:

	pkHelix();
	~pkHelix();	//null constructor for lists
	pkHelix(datatable *data, structure *st, int size, short *helixList);
	//4 argument constructor
	//data and st are the datatable and structure class defined by DHM
	//helixList is the list of nucleotides that defines a helix
	//	in the arrangement 1,2,3,4,5,6,7 etc
	//	where 1 is basepaired with 2, 3 is basepaired with 4, etc
	//	5'half of helix is 1,3,5,7
	//	3'half of helix is 2,4,6,8
	//size is the total number of nucleotides in the helix, therefore the length of helix is size/2
	//This function is used to build all lists of helices, including the list of possible pseudoknots
	//	in the function findhelix, as well as to hold the list of helices in the pseudoknot free minimum
	//	free energy structure built in the function convertStructureToHelixList 
	//	and used to get rid of false positives in the function comparepkHelixListToMFE
	void setEqual(pkHelix & ExistingPKhelix);
	//Function that is used to assign a newly constructed pkhelix object an already established pkhelix object
	//		essentially an assignment constructor without overloading =.
	//Function is used in building lists of helices for the list of possible pseudoknots in find helices and
	//		in the function that resizes the list of possible pseudoknots 
	void destroy();
	//Function to explicitly deallocate dynamic memory allocated to a pkHelix object, may be unnecessary
	void removeHelixFromStructure(structure* st);
	//Function that removes a helix from a structure object.  Assigns all nucleotides in the helix to be single stranded.
    
    void reverseRemoveHelix(structure* st);
    //Function that reverses single-stranded constraints imposed by removeHelixFromStructure

	void addHelixToStructure(structure* st, int structureNumber);
	//Function that adds a helix that was forced to be single-stranded back to a structure object after folding the rest 
	//	of the RNA.	 Also increments the energy of the structure by the SHAPE corrected energy of the helix.
	int getSize();
	//Returns the total number of nucleotides in a helix.  Length of the helix is therefore size/2.
	int getEnergy();
	//Returns the SHAPE corrected energy of the helix.
	void setEnergy(int i);
	//Sets the energy of a helix to i.	Used to increment the energy of the helix in the pseudoknotFold function
	//		by the pseudoknot penalty parameter.
	int calculateEnergy(structure *st, datatable *data);
	//Function that calculates the energy of a helix using SHAPE corrected terms for base-pair stacking parameters and
	//		the au-end penalty function.  Function does not handle bulges at this time.
	short getelement(int i);
	//Function returns the nucleotide at the ith position.	
	int get5primeFirst();
	//Returns the first nucleotide in the 5' half of the helix.
	int get5primeLast();
	//Returns the last nucleotide in the 5' half of the helix.
	int get3primeFirst();
	//Returns the first nucleotide in the 3' half of the helix.
	int get3primeLast();
	//Returns the last nucleotide in the 3' half of the helix.
	bool isEqual(pkHelix & PKhelix);
   
	//Copy Constructur
	pkHelix(const pkHelix& f){ copy(f); }

	const pkHelix& operator=(const pkHelix &f){
		copy(f);
		return(*this);
	}
   

 private:

	void copy(const pkHelix& f){
		if(&f == this)return;
		length=f.length;
		energy=f.energy;
		forwardFirst=f.forwardFirst;
		forwardLast=f.forwardLast;
		reverseFirst=f.reverseFirst;
		reverseLast=f.reverseLast;
		if(f.helix!=0){
			helix=new short[f.length];
			memcpy(helix,f.helix,f.length*sizeof(short));
		}
		else helix=0;
	}


	int length;			//holds the total number of nucleotides in the helix.  Length of helix is length/2.
	short *helix;		//pointer to the list of nucleotides
	int energy;			//holds energy of the helix
	int forwardFirst;	//First nucleotide in the 5' half of the helix
	int forwardLast;	//last nucleotide in the 5' half of the helix
	int reverseFirst;	//first nucleotide in the 3' half of the helix
	int reverseLast;	//last nucleotide in the 3' half of the helix	
};

#endif

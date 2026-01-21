
#if !defined(PKHELIX_CPP)
#define PKHELIX_CPP


#include "pkHelix.h"

#include <iostream>
using namespace std;


pkHelix::pkHelix(){
	helix=NULL;
}
void pkHelix::setEqual(pkHelix	&ExistingPKhelix)
{
	length = ExistingPKhelix.getSize();
	if (helix!=NULL) delete[] helix;
	helix = new short[length];
	//0 is base-paired to 1, 2 is basepaired to 3, etc.
	//0,2,4,6 is five prime half of helix
	//1,3,5,7 is 3' half of helix
	for(int i=0; i<length; i++){
		helix[i]=ExistingPKhelix.getelement(i);
	}
	energy = ExistingPKhelix.getEnergy();
	//store the first and last nucleotides of both the 5' and 3' half of the helices
	forwardFirst = helix[0];   
	reverseLast = helix[1];
	forwardLast = helix[length-2];
	reverseFirst = helix[length-1];
}

bool pkHelix::isEqual(pkHelix & PKhelix)
//check if two objects are equal
{
	if(length !=PKhelix.getSize() || energy !=PKhelix.getEnergy())
		return false;

	for (int i=0; i<length; i++){
		if(helix[i]!=PKhelix.getelement(i))
			return false;
	}
	return true;
}
	
pkHelix::pkHelix(datatable *data, structure *st, int size, short *helixList){
	//4 argument constructor that builds a new object from a pre-existing array of shorts
	helix = new short[size];
	length = size;
	//0 is base-paired to 1, 2 is basepaired to 3, etc.
	//0,2,4,6 is five prime half of helix
	//1,3,5,7 is 3' half of helix		
	for(int i=0; i<size; i++){
		helix[i]=helixList[i];
	}
	energy = calculateEnergy(st, data);	
	//store the first and last nucleotides of both the 5' and 3' half of the helices
	forwardFirst = helix[0];
	reverseLast = helix[1];
	forwardLast = helix[size-2];
	reverseFirst = helix[size-1];
		

}

pkHelix::~pkHelix(){
	//destructor
	if (helix!=NULL) delete[] helix;
}


int pkHelix::calculateEnergy(structure *st, datatable *data){
	//this function calculates the energy of a helix using the 1)penalty function defined by DHM to assign
	//		penalty to an au-end of a helix and 2) the erg1 function defined by DHM to assign the energy of
	//		stacking a base pair onto another

	int tempenergy = 0;	 //holds the energy of the helix as it is calculated
	//get two elements of the array at a time, corresponds to 1 base pair
	for(int i=0;i<length;i+=2){
		if (i==0) {
			//for the first base-pair, only calculate the au-end penalty
			tempenergy+=penalty(helix[i],helix[i+1], st, data);
		}
		else if(i==length-2){
			//for the last base-pair, calculate the au-end penalty first
			tempenergy+=penalty(helix[i],helix[i+1], st, data);
			//then calculate the energy of stacking this base pair onto another
			tempenergy+=erg1(helix[i-2],helix[i-1],helix[i], helix[i+1], st, data);			
		}
		//check that there is not a bulge
		else if( ((helix[i-2]+1)==helix[i]) && ((helix[i-1]-1)==helix[i+1])){
			//for interior base-pairs, only calculate stacking energies
			tempenergy+=erg1(helix[i-2],helix[i-1],helix[i], helix[i+1], st, data);
		}			
	}
	return tempenergy;
}
		

	
void pkHelix::removeHelixFromStructure(structure* st){
	//Set all nucleotides in the helx to be single stranded using structure object parameters defined by DHM.
	//Used to remove a helix from an RNA before the rest of the structure is folded using the dynamic function defined by DHM.
	//st->nnopair=length;		//structure object parameter defining number of nucleotides forced to be single-stranded
	for (int i=0; i<length; i++){
		st->AddSingle(helix[i]);
		
	}
}


void pkHelix::reverseRemoveHelix(structure* st){
    //Remove single-stranded constraints from the structure set by removeHelixFromStructure
    //added by TONY      
    st->RemoveSingleStrandConstraints(length);
}
  
void pkHelix::addHelixToStructure(structure * st, int structureNumber){
	//Resets all base-pairs in a helix that was previously forced to be single stranded
	st->SetEnergy(structureNumber,st->GetEnergy(structureNumber)+energy);  //structure object parameter defining overall energy of the folded RNA
    
    //RemoveConstraints has been replaced by call to reverseRemoveHelix in pseudoknotFold
    //st->RemoveConstraints();				  //structure object parameter defining number of nucleotides forced to be single-stranded
	
    //restore two nucleotides, or 1 basepair, at a time
	for (int i=0; i<length; i+=2){
		st->SetPair(helix[i+1],helix[i],structureNumber);	//basepr is list of base pairs in structure object. Each element in this list holds its basepairing partner, ie base pair 10-20 is basepair[10]=20.
															//structureNumber idicates which structure, defined in maxtracebacks.
	}
}

		
int pkHelix::get5primeFirst(){
	return forwardFirst;
}
int pkHelix::get5primeLast(){
	return forwardLast;
}
void pkHelix::setEnergy(int i){
	energy = i;
}
int pkHelix::get3primeFirst(){
	return reverseFirst;
}
int pkHelix::get3primeLast(){
	return reverseLast;
}
	
int pkHelix::getSize(){
	return length;
}
short pkHelix::getelement(int i){
	return helix[i];
}

int pkHelix::getEnergy(){
	return energy;
}

#endif

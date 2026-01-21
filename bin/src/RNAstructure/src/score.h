//Calculate scores when comparing a predicted structure (called test) to a known structure (called correct).
#include "structure.h"

//Calculate PPV, the fraction of predicted pairs that are correct.
	//Requires pointers to two structure classes, the test and correct structures
	//score is a pointer to integer that tabulates the number of correct predictions
	//basepairs is the total number of predicted pairs
	//structure number is the number of the structure to score in test (correct can only have one structure)
	//exact is a bool that indicates whether exact matches are required, normally this is set false to allow slippage

//This function returns an int that prevides an error code
	//0 = no error
	//1 = more than one structure in correct
	//2 = two structures are different lengths
	
int scorerppv(structure* correct, structure* test, int *score, int *basepairs, const int structurenumber, const bool exact);

//Calculate sensitivity, the fraction of known pairs correctly predicted.
	//Requires pointers to two structure classes, the test and correct structures
	//score is a pointer to integer that tabulates the number of correct predictions
	//basepairs is the total number of known pairs
	//structure number is the number of the structure to score in test (correct can only have one structure)
	//exact is a bool that indicates whether exact matches are required, normally this is set false to allow slippage

//This function returns an int that prevides an error code
	//0 = no error
	//1 = more than one structure in correct
	//2 = two structures are different lengths

int scorer(structure* correct, structure* test, int *score, int *basepairs, const int structurenumber, const bool exact);


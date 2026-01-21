// OligoWalk_class.cpp : Defines the entry point for the console application.
//

#include "Oligowalk_object.h"
#include <iostream>


int main(int argc, char* argv[])
{
	//Test OligoScreen
	Oligowalk_object *oligo;
	int error;

	//Instantiate a new instance of Oligowalk_object.
	oligo = new Oligowalk_object(false);


	if (oligo->GetErrorCode()!=0) {
		//check for errors.

		std::cerr << oligo->GetErrorMessage(oligo->GetErrorCode()) << "\n";
		delete oligo;

		return 1;

	}

	error = oligo->OligoScreen("oligoscreen_example.lis","test.out");
	if (error!=0) {
		//check for errors.

		std::cerr << oligo->GetErrorMessage(error) << "\n";
		delete oligo;

		return 1;

	}

	delete oligo;


	return 0;
}


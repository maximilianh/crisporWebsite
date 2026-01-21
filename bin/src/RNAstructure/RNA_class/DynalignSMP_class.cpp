// Dynalign_class.cpp : A simple demonstrator for the Dyalign_object class.
//



#include <iostream>
#include <cstdlib>
#include "Dynalign_object.h"

using namespace std;

int main(int argc, char* argv[])
{

	Dynalign_object *dyn; 
	int error;

	//check for the correct # of parameters
	if (argc!=7) {
		std::cout << "Usage: Dynalign_class #_processors seq1.seq seq2.seq ct1.ct ct2.ct alignment.ali\n";
		return 0;
	}
	dyn = new Dynalign_object(argv[2],2,argv[3],2);

	error = dyn->GetErrorCode();
	if (error!=0) {
		cout << dyn->GetErrorMessage(error);
		delete dyn;
		return 0;
	}


	//Force an alignment as a test
	//error = dyn->ForceAlignment(1,5);
	//if (error!=0) cout << dyn->GetErrorMessage(error);

	error = dyn->Dynalign(100,0,0,20,-99,0.4,true,NULL,false,30,false,std::atoi(argv[1]));
	if (error!=0) cout << dyn->GetErrorMessage(error);

	
	else {

		dyn->GetRNA1()->WriteCt(argv[4]);
		dyn->GetRNA2()->WriteCt(argv[5]);
		dyn->WriteAlignment(argv[6]);
	}
	delete dyn;

	return 0;
}

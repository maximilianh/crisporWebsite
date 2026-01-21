// Example of RNA Class Usage.
//

//Provide example of free energy calculation.
//Program takes one parameters.  The first is the name of a ct file to read.
//The program then output is the folding free energy change to std out. 

#include "RNA.h"
#include "../src/structure.h" //get the definition of ct just to read the ct from disk
#include <iostream>
#include <cstdlib>

// example usage:
// ./RNA_class example.ct 1
int main(int argc, char* argv[])
{

	RNA *rna;
	structure ct;
	int s;
	int i,j;
	char *string;
	int error;
	double free_energy;


	if (argc!=3) {
		std::cout << "Usage: RNA_class ct_file_input structure_#_for_energy_calculation\n";
		return 0;
	}

	//convert the input for structure number to integer
	s = atoi(argv[2]);

	// construct an RNA object from the input
	rna = new RNA(argv[1], FILE_CT, true); //construct instance of RNA


	//now return the calulated energy and display it:
	free_energy = rna->CalculateFreeEnergy(s);
	error = rna->GetErrorCode();
	if (error==0) {
		//Note that when calculate energy is called the first time, RNA reads parameter files from
		//disk at the location specified by environment variable DATAPATH.

		//These are the .dat files found in the data_tables directory of RNAstructure
		std::cout << "Free energy change is: "<<rna->CalculateFreeEnergy(s) << "\n";
	}
	else {
		std::cerr << rna->GetErrorMessage(error);
	}
		


	delete rna;
	return 0;
}


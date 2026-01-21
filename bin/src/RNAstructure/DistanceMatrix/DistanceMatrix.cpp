// DistanceMatrix.cpp : This file contains the 'main' function. Program execution begins and ends there.
// Open a ct file, and calculate a distance matrix bwteen all pairs of structures.
//

#include <iostream>
#include "../RNA_class/RNA.h"


int main(int argc, char **argv)
{
	if (argc != 3) {
		std::cout << "Usage: DistanceMatrix inputctfile outputmatrixfile\n";
		exit(0);
	}

	RNA* rna;
	//Open the ct file
	rna = new RNA(argv[1], 1);

	//make an array to save the distances
	int** distances;
	distances = new int *[rna->GetStructureNumber() + 1];
	for (int i = 1; i <= rna->GetStructureNumber(); ++i) {
		distances[i]=new int[rna->GetStructureNumber() + 1];
	}

	//loop over all structures:
	for (int i = 1; i <= rna->GetStructureNumber(); ++i) {
		distances[i][i] = 0;
		for (int j = i + 1; j <= rna->GetStructureNumber(); ++j) {
			distances[i][j] = 0;
			//now loop over nucleotides
			for (int k = 1; k <= rna->GetSequenceLength(); ++k) {
				if (rna->GetPair(k, i) != rna->GetPair(k, j)) {
					++distances[i][j];
					if (rna->GetPair(k, i) != 0 && rna->GetPair(k, j)!=0) ++distances[i][j];
				}

			}
			//fill in the other triangular part of the matrix
			distances[j][i] = distances[i][j];


		}
		
	}
	std::ofstream outfile;
	outfile.open(argv[2]);
	for (int i = 1; i <= rna->GetStructureNumber(); ++i) {
		for (int j = 1; j <= rna->GetStructureNumber(); ++j) {
			if (j== rna->GetStructureNumber()) outfile << distances[i][j] << "\n";
			else outfile << distances[i][j] << "\t";

		}
	}



	
	for (int i = 1; i <= rna->GetStructureNumber(); ++i) {
		delete[] distances[i];
	}
	delete[] distances;
	delete rna;


}



#include <vector>
#include <complex>
#include <valarray>


#include "./RnaContainerClass.h"

#include "../RNA_class/RNA.h"
#include "../RNA_class/thermodynamics.h"
#include "../RNA_class/HybridRNA.h"
#include "../src/defines.h"
#include "../src/ErrorChecker.h"

static const double t1 = 29+273.15;//Set the first temperature to be 29C
static const double t2 = 33+273.15;//Set the second temperature to be 33C
const double PI = 3.141592653589793238460;
typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

// Enum used to determine which objective function to use.
enum OregaObjectiveFunction {
	OREGA_SIMPLE=1,
	OREGA_COMPLEX=2,
        OREGA_ALL=3
};

void fft(CArray& x);
//Declare the function the calculates the complexity of a sequences according to the model in
//Gabrielian, Andrei, and Alexander Bolshoy. "Sequence complexity and DNA curvature." Computers & chemistry 23.3-4 (1999): 263-274.
double complexity_calc(const string& sequence,int filterAUG, int filterCUG, int filteroligoA);

int orega(std::string InputFile,
	std::string OutputFile,
	std::string savefile,
	std::string restartFile,
	int NumIterations,
	const char* const alphabet,
    int MutationStartNuc,//position at which to start nucleotide mutations
    int NumberOfNucs,
    double MutationRate,
    int RecombinationFrequency,
    double RecombinationRate,
    int NumberOfSequences,
	int objectiveFunction,
    int randomseed,
    int MutationSwitch,
    int limitG,
    int filterAUG,
    int filterCUG,
    int filteroligoA,
    double complexity_constant,double threshold);

//Declare the function that calculates and returns the RMSD
double CalculateRMSD(RNA* rnaWT, RNA* rna);

double CalculateFFT(RNA* rna, int start_mutation_site, int number_of_nucs_to_mutate, int objectiveFunction, int filterAUG, int filterCUG, int filteroligoA, double complexity_constant,int total_length, double thershold);

//Declare the function that cutates a nucleotide
//returns 1 if the
//int MutateNuc(RNA* rna);
int MutateNuc(RNA* rna, int MutationStartNuc, int NumberOfNucs, double MutationRate, int MutationSwitch,int limitG, randomnumber *randomnum,int filterAUG,int filterCUG,int filteroligoA);

//Declare the function that converts nucleotide to number. A=a=1, C=c=2, G=g=3, T=U=t=u=4
int nuc2num(char nuc);

//Declare the function that cycles through codons keeping the Amino Acid identity
//Returns "err" in case of error, or "NNN" in case no change to codon is possible.
string CycleCodon(std::string Codon);

//int Recombine(RNA* rna1, RNA* rna2, RNA* rnaStore1, RNA* rnaStore2);
int Recombine(RNA* rna1, RNA* rna2, RNA* rnaStore1, RNA* rnaStore2, int NumberOfNucs, int MutationStartNuc, double RecombinationRate, randomnumber *randomnum);

//Function that checks if the recombination position has been repeated
bool CheckForRepeat(std::vector<std::vector<int> > RecombinationPair, int position);


bool AllSwapped(std::vector<bool> &KeepSequences,int NumberOfSequences);


void CopySequence(RNA* rnaCopyTo, RNA* rnaCopyFrom);

void ObjectiveFunction(std::vector<vector<double> > &CalculatedRMSDs, int i);

//Function that calculates the total base pair probabilities of the "length" number of nucleotides from "start" nucleotide
double totalPairwiseProbability(RNA* rna, int start, int length);

double Basepair_fraction(RNA* rna, int start, int length,double threshold);

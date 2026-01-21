#define MAX_INTER 30
#include "DynProgArray.h"
#include "TProgressDialog.h"

#include <vector>
#include "structure.h"
#include "../src/phmm/utils/xmath/log/xlog_math.h"

//static variables
static double DOUBLE_INFINITY = 1e300;

// core functionality methods

// this function runs MaxExpect using data written to disk
void bpMatch(structure *ct, char* pfsfile, double gamma, double maxPercent, int maxStructures, int Window, ProgressHandler *progress=NULL);

//This function sets up the fill routine and runs the traceback routine
void MaxExpectFill(structure *ct, DynProgArray<PFPRECISION> *v, PFPRECISION *w5, pfdatatable *pfdata, bool *lfce, bool *mod, forceclass *fce, double maxPercent, int maxStructures, int Window, double gamma=1.0, ProgressHandler *progress=NULL);

//This is actual fill routine
void MEAFill(structure *ct, double **bpProbArray, double *bpSSProbArray, double **vwArray, double **vwPArray, double *w5Array, double *w3Array, std::vector< std::vector<bool> >* pfdata, double gamma, double maxPercent, ProgressHandler *progress, bool OnlyCanonical=true);

// execute the recursion function based on probabilities from the partition function
void bpProbRecursion(double **bpProbArray, double **vwArray, structure *ct, char* pfsfile);

// execute the traceback utilizing the v and w Arrays - internal fragments ( nucs i to j, inclusive)
void traceBack(structure *ct, double **vwArray, double **bpProbArray, double gamma, int ip, int jp);

// execute the traceback utilizing the v and w Arrays - external fragments ( nucs 1 to i and j to N)
void traceBackExternal(structure *ct, double **vwArray, double **vwPArray, double **bpProbArray, double gamma, int ip, int jp);

//Coordinates traceback of suboptimal structures
void trace(structure *ct, double **vwArray, double **vwPArray, double **bpProbArray, double gamma, double maxPercent, int maxStructures, int Window);

// compares 2 double values for equality
bool doubleEqual(double double1, double double2);

// check to see if the pair is a canonical
bool isCanonical(char i, char j);


// get the location of a passed branch value looking for a structure value
bool getStructure(int i, int j, double branchValue,
                         double **vwArray, int *branchPt);

// gets the max from a array of numbers of a defined size
void getMax(double *max, double *valueArray, int size);

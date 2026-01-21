#ifndef _TURBOFOLD_PROBABILISTICMODEL_
#define _TURBOFOLD_PROBABILISTICMODEL_

#include <list>
#include <cmath>
#include <cstdio>
#include "SafeVector.h"
#include "SparseMatrix.h"
#include "MultiSequence.h"

#ifdef TURBOHOMOLOGY
	#define ChooseBestOfThree ChooseBestOfThreeturbohomology
	#define ProbabilisticModel ProbabilisticModelturbohomology
#endif


using namespace std;

//! ProbabilisticModel Class
/*!
    The ProbabilisticModel Class stores the parameters of a probabilistic model.

*/

#define NumMatrixTypes 3   // One match state and two insert states.

//! Store the largest of three values x1, x2, and x3 in *x.  
//! If xi is the largest value, then store bi in *b.
void ChooseBestOfThree(float x1, float x2, float x3, char b1, char b2, char b3, float *x, char *b);

class ProbabilisticModel{

public:


    //! Computes an alignment based on given posterior matrix.
    //! This is done by finding the maximum summing path (or maximum weight trace) through the posterior matrix.
    //! The final alignment is returned as a pair consisting of:
    //! (1) a string (e.g., XXXBBXXXBBBBBBYYYYBBB) where X's and denote insertions in one of the two sequences and B's denote that both sequences are present (i.e. matches).
    //! (2) a float indicating the sum achieved.
    pair<SafeVector<char> *, float> ComputeAlignment(int seq1Length, int seq2Length, const SafeVector<float> &posterior) const ;

    //! Builds a posterior probability matrix needed to align a pair of alignments.
#ifdef TURBOHOMOLOGY
    SafeVector<float> *BuildPosterior(MultiSequence *align1, MultiSequence *align2,
                       const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices, const vector<int> &mappingSeqIndex = vector<int>(), float cutoff = 0.0f) const ;


    SafeVector<float> *BuildAlnScore(MultiSequence *align1, MultiSequence *align2,
                       const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices, 
                       const vector<int> &mappingSeqIndex = vector<int>(), 
                       const vector<vector<vector<double> > > basePairScore = vector<vector<vector<double> > >(),
                       const vector<double> similarity_list = vector<double>(),
                       float cutoff = 0.0f) const ;
#else
    SafeVector<float> *BuildPosterior(MultiSequence *align1, MultiSequence *align2,
                     const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices, float cutoff = 0.0f) const ;

#endif
};

#endif
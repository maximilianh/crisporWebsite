#include "SparseMatrix.h"
#include <iostream>

using namespace std;

//! SparseMatrix
/*! 
    Class for sparse matrix computations.
*/

//! Private constructor.
SparseMatrix::SparseMatrix(){}

//! Constructor. 
//! Builds a sparse matrix from a posterior matrix.
//! Note that the expected format for the posterior matrix is as a (seq1Length+1) x (seq2Length+1) matrix 
//! where the 0th row and 0th column are ignored (they should contain all zeroes).
SparseMatrix::SparseMatrix(int seq1Length, int seq2Length, const SafeVector<float> &posterior) :
    seq1Length(seq1Length), seq2Length(seq2Length) {

    int numCells = 0;

    assert (seq1Length > 0);
    assert (seq2Length > 0);

    // calculate memory required; count the number of cells in the
    // posterior matrix above the threshold
    SafeVector<float>::const_iterator postPtr = posterior.begin();
    for (int i = 0; i <= seq1Length; i++){
        for (int j = 0; j <= seq2Length; j++){
            if (*(postPtr++) >= POSTERIOR_CUTOFF){
                assert (i != 0 && j != 0);
                numCells++;
            }
        }
    }
    
    // allocate memory
    data.resize(numCells);
    rowSize.resize (seq1Length + 1); rowSize[0] = -1;
    rowPtrs.resize (seq1Length + 1); rowPtrs[0] = data.end();

    // build sparse matrix
    postPtr = posterior.begin() + seq2Length + 1;           // Note that we're skipping the first row here
    SafeVector<PIF>::iterator dataPtr = data.begin();
    for (int i = 1; i <= seq1Length; i++){
        postPtr++;                                          // Skipping the first column of each row
        rowPtrs[i] = dataPtr;
        for (int j = 1; j <= seq2Length; j++){
            if (*postPtr >= POSTERIOR_CUTOFF){
                dataPtr->first = j;
                dataPtr->second = *postPtr;
                dataPtr++;
            }
            postPtr++;
        }
        rowSize[i] = dataPtr - rowPtrs[i];
    }
}

//! Returns the pointer to a particular row in the sparse matrix.
SafeVector<PIF>::iterator SparseMatrix::GetRowPtr(int row) const {
    assert (row >= 1 && row <= seq1Length);
    return rowPtrs[row];
}

//! Returns the number of entries in a particular row.
int SparseMatrix::GetRowSize (int row) const {
    assert (row >= 1 && row <= seq1Length);
    return rowSize[row];
}

//! Returns the first dimension of the matrix.
int SparseMatrix::GetSeq1Length () const {
    return seq1Length;
}

//! Returns the second dimension of the matrix.
int SparseMatrix::GetSeq2Length () const {
    return seq2Length;
}


//! Returns a new sparse matrix containing the transpose of the current matrix.
SparseMatrix * SparseMatrix::ComputeTranspose() const {

    // create a new sparse matrix
    SparseMatrix *ret = new SparseMatrix();
    int numCells = data.size();

    ret->seq1Length = seq2Length;
    ret->seq2Length = seq1Length;

    // allocate memory
    ret->data.resize (numCells);
    ret->rowSize.resize (seq2Length + 1); ret->rowSize[0] = -1;
    ret->rowPtrs.resize (seq2Length + 1); ret->rowPtrs[0] = ret->data.end();

    // compute row sizes
    for (int i = 1; i <= seq2Length; i++) ret->rowSize[i] = 0;
    for (int i = 0; i < numCells; i++)
        ret->rowSize[data[i].first]++;

    // compute row ptrs
    for (int i = 1; i <= seq2Length; i++){
        ret->rowPtrs[i] = (i == 1) ? ret->data.begin() : ret->rowPtrs[i-1] + ret->rowSize[i-1];
    }

    // now fill in data
    SafeVector<SafeVector<PIF>::iterator> currPtrs = ret->rowPtrs;

    for (int i = 1; i <= seq1Length; i++){
        SafeVector<PIF>::iterator row = rowPtrs[i];
        for (int j = 0; j < rowSize[i]; j++){
            currPtrs[row[j].first]->first = i;
            currPtrs[row[j].first]->second = row[j].second;
            currPtrs[row[j].first]++;
        }
    }

    return ret;
}


//! Return the posterior representation of the sparse matrix.
SafeVector<float> * SparseMatrix::GetPosterior() const {

    // create a new posterior matrix
    SafeVector<float> *posteriorPtr = new SafeVector<float>((seq1Length+1) * (seq2Length+1)); 
    assert (posteriorPtr);
    SafeVector<float> &posterior = *posteriorPtr;

    // build the posterior matrix
    for (int i = 0; i < (seq1Length+1) * (seq2Length+1); i++) posterior[i] = 0;
    for (int i = 1; i <= seq1Length; i++){
        SafeVector<float>::iterator postPtr = posterior.begin() + i * (seq2Length+1);
        for (int j = 0; j < rowSize[i]; j++){
            postPtr[rowPtrs[i][j].first] = rowPtrs[i][j].second;
        }
    }
    return posteriorPtr;
}


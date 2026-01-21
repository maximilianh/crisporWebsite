#ifndef _SPARSEMATRIX_
#define _SPARSEMATRIX_

#include <iostream>
#include "SafeVector.h"
#include "SparseMatrix.h"

using namespace std;

#define POSTERIOR_CUTOFF 0.01         // minimum posterior probability
                                             // value that is maintained in the
                                             // sparse matrix representation

typedef pair<int,float> PIF;                 // Sparse matrix entry type
                                             //   first --> column
                                             //   second --> value

//! SparseMatrix
/*! 
    Class for sparse matrix computations.
*/

class SparseMatrix {

    int seq1Length, seq2Length;                     // dimensions of matrix
    SafeVector<int> rowSize;                        // rowSize[i] = # of cells in row i
    SafeVector<PIF> data;                           // data values
    SafeVector<SafeVector<PIF>::iterator> rowPtrs;  // pointers to the beginning of each row

    //! Private constructor.
    SparseMatrix();

public:
    //! Constructor. 
    //! Builds a sparse matrix from a posterior matrix.
    //! Note that the expected format for the posterior matrix is as a (seq1Length+1) x (seq2Length+1) matrix 
    //! where the 0th row and 0th column are ignored (they should contain all zeroes).
    SparseMatrix(int seq1Length, int seq2Length, const SafeVector<float> &posterior);

    //! Returns the pointer to a particular row in the sparse matrix.
    SafeVector<PIF>::iterator GetRowPtr(int row) const ;

    //! Returns the number of entries in a particular row.
    int GetRowSize (int row) const ;
    
    //! Returns the first dimension of the matrix.
    int GetSeq1Length () const ;

    //! Returns the second dimension of the matrix.
    int GetSeq2Length () const ;

    //! Returns a new sparse matrix containing the transpose of the current matrix.
    SparseMatrix *ComputeTranspose() const ;


    //! Return the posterior representation of the sparse matrix.
    SafeVector<float> *GetPosterior() const ;
};

#endif
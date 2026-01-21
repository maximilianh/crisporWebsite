#include "ProbabilisticModel.h"
#include <list>
#include <cmath>
#include <cstdio>


using namespace std;

//! ProbabilisticModel Class
/*!
    The ProbabilisticModel Class stores the parameters of a probabilistic model.

*/


//! Store the largest of three values x1, x2, and x3 in *x.  
//! If xi is the largest value, then store bi in *b.
void ChooseBestOfThree(float x1, float x2, float x3, char b1, char b2, char b3, float *x, char *b){
    if (x1 >= x2){
        if (x1 >= x3){
            *x = x1;
            *b = b1;
            return;
        }
        *x = x3;
        *b = b3;
        return;
    }
    if (x2 >= x3){
        *x = x2;
        *b = b2;
        return;
    }
    *x = x3;
    *b = b3;
}

//! Computes an alignment based on given posterior matrix.
//! This is done by finding the maximum summing path (or maximum weight trace) through the posterior matrix.
//! The final alignment is returned as a pair consisting of:
//! (1) a string (e.g., XXXBBXXXBBBBBBYYYYBBB) where X's and denote insertions in one of the two sequences and B's denote that both sequences are present (i.e. matches).
//! (2) a float indicating the sum achieved.
pair<SafeVector<char> *, float> ProbabilisticModel::ComputeAlignment(int seq1Length, int seq2Length, const SafeVector<float> &posterior) const {

    float *twoRows = new float[(seq2Length+1)*2]; 
    assert (twoRows);
    float *oldRow = twoRows;
    float *newRow = twoRows + seq2Length + 1;

    char *tracebackMatrix = new char[(seq1Length+1)*(seq2Length+1)]; 
    assert (tracebackMatrix);
    char *tracebackPtr = tracebackMatrix;

    SafeVector<float>::const_iterator posteriorPtr = posterior.begin() + seq2Length + 1;

    // Initialization. Top row are all follow left.
    for(int i = 0; i <= seq2Length; i++){
        oldRow[i] = 0;
        *(tracebackPtr++) = 'L';
    }

    // Fills the matrix.
    for(int i = 1; i <= seq1Length; i++){
        // Initializes the left column.
        newRow[0] = 0;
        posteriorPtr++;
        *(tracebackPtr++) = 'U';

        // Fills the matrix.
        for(int j = 1; j <= seq2Length; j++){
            ChooseBestOfThree(*(posteriorPtr++) + oldRow[j-1], newRow[j-1], oldRow[j], 'D', 'L', 'U', &newRow[j], tracebackPtr++);
        }

        // Swaps the old and new rows.
        float *temp = oldRow;
        oldRow = newRow;
        newRow = temp;
    }
    // store best score
    float total = oldRow[seq2Length];
    delete [] twoRows;

    // compute traceback
    SafeVector<char> *alignment = new SafeVector<char>; 
    assert (alignment);
    
    int r = seq1Length, c = seq2Length;
    while (r != 0 || c != 0){
        char ch = tracebackMatrix[r*(seq2Length+1) + c];
        switch (ch){
        case 'L': c--; alignment->push_back ('Y'); break;
        case 'U': r--; alignment->push_back ('X'); break;
        case 'D': c--; r--; alignment->push_back ('B'); break;
        default: assert (false);
        }
    }

    delete [] tracebackMatrix;
    reverse(alignment->begin(), alignment->end());

    return make_pair(alignment, total);
}

//! Builds a posterior probability matrix needed to align a pair of alignments.

#ifdef TURBOHOMOLOGY
SafeVector<float> * ProbabilisticModel::BuildPosterior(MultiSequence *align1, MultiSequence *align2,
                   const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices, const vector<int> &mappingSeqIndex, float cutoff) const {
#else
SafeVector<float> * ProbabilisticModel::BuildPosterior(MultiSequence *align1, MultiSequence *align2,
                   const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices, float cutoff) const {
#endif
    const int seq1Length = align1->GetSequence(0)->GetLength();
    const int seq2Length = align2->GetSequence(0)->GetLength();

    SafeVector<float> *posteriorPtr = new SafeVector<float>((seq1Length+1) * (seq2Length+1), 0); assert (posteriorPtr);
    SafeVector<float> &posterior = *posteriorPtr;
    SafeVector<float>::iterator postPtr = posterior.begin();

    // Loops through align1
    for (int i = 0; i < align1->GetNumSequences(); i++){
        int first = align1->GetSequence(i)->GetLabel();
        SafeVector<int> *mapping1 = align1->GetSequence(i)->GetMapping();

        // Loops through align2
        for (int j = 0; j < align2->GetNumSequences(); j++){
            int second = align2->GetSequence(j)->GetLabel();
            SafeVector<int> *mapping2 = align2->GetSequence(j)->GetMapping();

#ifdef TURBOHOMOLOGY
            bool isMap = false;
            // cout << "size: " << mappingSeqIndex[] << endl;
            for (int index = 0; index < mappingSeqIndex.size(); index++){
                // cout << "mappingSeqIndex[index]: " << mappingSeqIndex[index] << endl;
                if (mappingSeqIndex[index] == j)
                    isMap = true;
            }
            if (!isMap){
                // cout << "skip index: " << j << endl;
                continue;
            }
            else {
                // cout << "index: " << j << endl;
#endif

                if(first < second){
                    // get the associated sparse matrix
                    SparseMatrix *matrix = sparseMatrices[first][second];
                    
                    for (int ii = 1; ii <= matrix->GetSeq1Length(); ii++){
                        SafeVector<PIF>::iterator row = matrix->GetRowPtr(ii);
                        int base = (*mapping1)[ii] * (seq2Length+1);
                        int rowSize = matrix->GetRowSize(ii);
                      
                        // add in all relevant values
                        for (int jj = 0; jj < rowSize; jj++)
                            posterior[base + (*mapping2)[row[jj].first]] += row[jj].second;
                      
                        // subtract cutoff 
                        for (int jj = 0; jj < matrix->GetSeq2Length(); jj++)
                            posterior[base + (*mapping2)[jj]] -= cutoff;
                    }
                }
                else {
                    // get the associated sparse matrix
                    SparseMatrix *matrix = sparseMatrices[second][first];
                    
                    for (int jj = 1; jj <= matrix->GetSeq1Length(); jj++){
                        SafeVector<PIF>::iterator row = matrix->GetRowPtr(jj);
                        int base = (*mapping2)[jj];
                        int rowSize = matrix->GetRowSize(jj);
                      
                        // add in all relevant values
                        for (int ii = 0; ii < rowSize; ii++)
                            posterior[base + (*mapping1)[row[ii].first] * (seq2Length + 1)] += row[ii].second;
                      
                        // subtract cutoff 
                        for (int ii = 0; ii < matrix->GetSeq2Length(); ii++)
                            posterior[base + (*mapping1)[ii] * (seq2Length + 1)] -= cutoff;
                    }
                }
#ifdef TURBOHOMOLOGY
            }
#endif
            delete mapping2;
        }

        delete mapping1;
    }

    return posteriorPtr;
}

#ifdef TURBOHOMOLOGY
SafeVector<float> * ProbabilisticModel::BuildAlnScore(MultiSequence *align1, MultiSequence *align2,
                   const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices, 
                   const vector<int> &mappingSeqIndex, 
                   const vector<vector<vector<double> > > basePairScore, 
                   const vector<double> similarity_list,
                   float cutoff) const {

    const int seq1Length = align1->GetSequence(0)->GetLength();
    const int seq2Length = align2->GetSequence(0)->GetLength();

    SafeVector<float> *posteriorPtr = new SafeVector<float>((seq1Length+1) * (seq2Length+1), 0); assert (posteriorPtr);
    SafeVector<float> &posterior = *posteriorPtr;
    SafeVector<float>::iterator postPtr = posterior.begin();
    // Loops through align1
    for (int i = 0; i < align1->GetNumSequences(); i++){
        int first = align1->GetSequence(i)->GetLabel();
        SafeVector<int> *mapping1 = align1->GetSequence(i)->GetMapping();

        // Loops through align2
        for (int j = 0; j < align2->GetNumSequences(); j++){
            int second = align2->GetSequence(j)->GetLabel();
            SafeVector<int> *mapping2 = align2->GetSequence(j)->GetMapping();

            bool isMap = false;
            for (int index = 0; index < mappingSeqIndex.size(); index++){
                if (mappingSeqIndex[index] == j+1)
                    isMap = true;
            }
            // isMap = true;
            if (!isMap){
                // continue;
            }
            else {
                cout << "index: " << j << endl;
                if(first < second){
                    // get the associated sparse matrix
                    SparseMatrix *matrix = sparseMatrices[first][second];
                    for (int n = 0; n < matrix->GetSeq2Length(); n++){
                    }   
                    for (int ii = 1; ii <= matrix->GetSeq1Length(); ii++){
                        SafeVector<PIF>::iterator row = matrix->GetRowPtr(ii);
                        int base = (*mapping1)[ii] * (seq2Length+1);
                        int rowSize = matrix->GetRowSize(ii);

                        for (int jj = 0; jj < rowSize; jj++){
                            posterior[base + (*mapping2)[row[jj].first]] += row[jj].second * 0.5;

                            posterior[base + (*mapping2)[row[jj].first]] += basePairScore[j][ii][(*mapping2)[row[jj].first]] * 0.5;
                            
                        }
                    }
                }
                else {
                    // get the associated sparse matrix
                    SparseMatrix *matrix = sparseMatrices[second][first];
                    for (int jj = 1; jj <= matrix->GetSeq1Length(); jj++){
                        SafeVector<PIF>::iterator row = matrix->GetRowPtr(jj);
                        int base = (*mapping2)[jj];
                        int rowSize = matrix->GetRowSize(jj);
                        // add in all relevant values
                        for (int ii = 0; ii < rowSize; ii++){
                            posterior[base + (*mapping1)[row[ii].first] * (seq2Length + 1)] += row[ii].second;
                        }
                    }
                }
            }
            
            delete mapping2;
            
        }

        delete mapping1;
    }

    return posteriorPtr;
}
#endif
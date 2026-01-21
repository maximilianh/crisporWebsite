#include "Alignment.h"

#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <set>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <iomanip>

//! This function takes two multiple sequence alignments as input.
//! Returns the alignment of the two MultiSequence objects.
//! \param align1 and \param align2 are two multiple sequence alignments.
//! \param sparseMatrices is matrices storing posterior probability for each pair of sequences.
//! \param model is storing the parameters of a probabilistic model. (initDistrib, gapOpen, gapExtend, emitPairs, emitSingle) 
MultiSequence *AlignAlignments (MultiSequence *align1, MultiSequence *align2,
                                const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                                const ProbabilisticModel &model){

    // Print some info about the alignment
    SafeVector<float> *posterior = model.BuildPosterior (align1, align2, sparseMatrices);
    pair<SafeVector<char> *, float> alignment;

    // Choose the alignment routine depending on the "cosmetic" gap penalties used
    alignment = model.ComputeAlignment (align1->GetSequence(0)->GetLength(), align2->GetSequence(0)->GetLength(), *posterior);

    delete posterior;

    // Build final alignment
    MultiSequence *result = new MultiSequence();
    for (int i = 0; i < align1->GetNumSequences(); i++)
        result->AddSequence (align1->GetSequence(i)->AddGaps(alignment.first, 'X'));
    for (int i = 0; i < align2->GetNumSequences(); i++)
        result->AddSequence (align2->GetSequence(i)->AddGaps(alignment.first, 'Y'));
    result->SortByLabel();

    // Free temporary alignment
    delete alignment.first;

    return result;
}

#ifdef TURBOHOMOLOGY
MultiSequence *AlignAlignmentsMappingSeq (MultiSequence *align1, MultiSequence *align2,
                                const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                                const ProbabilisticModel &model, vector<int> mappingSeqIndex,
                                const vector<vector<vector<double> > > basePairScore,
                                const vector<double> similarity_list){

    // Print some info about the alignment
    // SafeVector<float> *posterior = model.BuildPosterior (align1, align2, sparseMatrices, mappingSeqIndex);
    SafeVector<float> *posterior = model.BuildAlnScore (align1, align2, sparseMatrices, mappingSeqIndex, basePairScore, similarity_list);
    pair<SafeVector<char> *, float> alignment;
    // Choose the alignment routine depending on the "cosmetic" gap penalties used
    alignment = model.ComputeAlignment (align1->GetSequence(0)->GetLength(), align2->GetSequence(0)->GetLength(), *posterior);
    delete posterior;

    // Build final alignment
    MultiSequence *result = new MultiSequence();
    for (int i = 0; i < align1->GetNumSequences(); i++)
        result->AddSequence (align1->GetSequence(i)->AddGaps(alignment.first, 'X'));
    for (int i = 0; i < align2->GetNumSequences(); i++)
        result->AddSequence (align2->GetSequence(i)->AddGaps(alignment.first, 'Y'));
    result->SortByLabel();
    // Free temporary alignment
    delete alignment.first;

    return result;
}

//! AlignProfile()
//! This function takes a MSA of H sequences and a new sequence as input.
//! Returns the aligned sequences of H+1 sequences.
MultiSequence *AlignProfile (Sequence *newSeq, MultiSequence *align, 
                              const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                              const ProbabilisticModel &model, vector<int> mappingSeqIndex,
                              const vector<vector<vector<double> > > basePairScore, 
                              const vector<double> similarity_list){
    // Build align from known MSA

    // BuildPosterior

    // ComputeAlignment

    // Build final alignment

    // Or Build alignment from known MSA
    MultiSequence *new_seq_to_aln = new MultiSequence();
    new_seq_to_aln->AddSequence(newSeq->Clone());

    // AlignAlignments
    // MultiSequence *alignment = AlignAlignments(new_seq_to_aln, align, sparseMatrices, model);
    MultiSequence *alignment = AlignAlignmentsMappingSeq(new_seq_to_aln, align, sparseMatrices, model, mappingSeqIndex, basePairScore, similarity_list);

    delete new_seq_to_aln;
    return alignment;
}
#endif

//! This function computes the consistency transformation for sequence Z by taking two posterior probabilities matrices of alignments between X-Z and Z-Y.
//! For the case that sequence Z's index is larger than that of X.
//! The transformed matrix is added to \param posterior.
void ConsistencyTransform (SparseMatrix *matXZ, SparseMatrix *matZY, SafeVector<float> &posterior){

    assert (matXZ);
    assert (matZY);

    int lengthX = matXZ->GetSeq1Length();
    int lengthY = matZY->GetSeq2Length();
    assert (matXZ->GetSeq2Length() == matZY->GetSeq1Length());

    // For every x[i]
    for (int i = 1; i <= lengthX; i++){
        SafeVector<PIF>::iterator XZptr = matXZ->GetRowPtr(i);
        SafeVector<PIF>::iterator XZend = XZptr + matXZ->GetRowSize(i);

        SafeVector<float>::iterator base = posterior.begin() + i * (lengthY + 1);

        // Iterate through all x[i]-z[k]
        while (XZptr != XZend){
            SafeVector<PIF>::iterator ZYptr = matZY->GetRowPtr(XZptr->first);
            SafeVector<PIF>::iterator ZYend = ZYptr + matZY->GetRowSize(XZptr->first);
            const float XZval = XZptr->second;

            // Iterate through all z[k]-y[j]
            while (ZYptr != ZYend){
                base[ZYptr->first] += XZval * ZYptr->second;
                ZYptr++;
            }
            XZptr++;
        }
    }
}

//! Same as the above, but for the case that Z's index is smaller than that of X.
void ConsistencyTransform1 (SparseMatrix *matZX, SparseMatrix *matZY, SafeVector<float> &posterior){

    assert (matZX);
    assert (matZY);

    int lengthZ = matZX->GetSeq1Length();
    int lengthY = matZY->GetSeq2Length();

    // For every z[k]
    for (int k = 1; k <= lengthZ; k++){
        SafeVector<PIF>::iterator ZXptr = matZX->GetRowPtr(k);
        SafeVector<PIF>::iterator ZXend = ZXptr + matZX->GetRowSize(k);

        // Iterate through all z[k]-x[i]
        while (ZXptr != ZXend){
            SafeVector<PIF>::iterator ZYptr = matZY->GetRowPtr(k);
            SafeVector<PIF>::iterator ZYend = ZYptr + matZY->GetRowSize(k);
            const float ZXval = ZXptr->second;
            SafeVector<float>::iterator base = posterior.begin() + ZXptr->first * (lengthY + 1);

            // Iterate through all z[k]-y[j]
            while (ZYptr != ZYend){
                base[ZYptr->first] += ZXval * ZYptr->second;
                ZYptr++;
            }
            ZXptr++;
        }
    }
}

//! This function takes multiple sequences and posterior probability matices to perform three-way probabilistic consistency transformation.
//! Returns new re-estimated alignment score matrices. 
//! The formula is: P'(x[i]-y[j])=(1/|S|)*sum_z_in_S{ sum_k{ P(x[i]-z[k]) * P(z[k]-y[j]) } }
//! \param sequences is multiple sequences.
//! \param sparseMatrices is posterior probaility from pairwise HMM alignment.
SafeVector<SafeVector<SparseMatrix *> > MultiConsistencyTransform(MultiSequence *sequences, 
                                                      SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices){
    
    const int numSeqs = sequences->GetNumSequences();
    SafeVector<SafeVector<SparseMatrix *> > newSparseMatrices (numSeqs, SafeVector<SparseMatrix *>(numSeqs, NULL));

    // For every pair of sequences
    for (int i = 0; i < numSeqs; i++){
        for (int j = i+1; j < numSeqs; j++){
            Sequence *seq1 = sequences->GetSequence (i);
            Sequence *seq2 = sequences->GetSequence (j);

            // Get the original posterior matrix
            SafeVector<float> *posteriorPtr = sparseMatrices[i][j]->GetPosterior(); assert (posteriorPtr);
            SafeVector<float> &posterior = *posteriorPtr;

            const int seq1Length = seq1->GetLength();
            const int seq2Length = seq2->GetLength();

            // Contribution from the summation where z = x and z = y
            for (int k = 0; k < (seq1Length+1) * (seq2Length+1); k++) posterior[k] += posterior[k];

            // Contribution from all other sequences
            for (int k = 0; k < numSeqs; k++) if (k != i && k != j){
                if (k < i)
                    ConsistencyTransform1 (sparseMatrices[k][i], sparseMatrices[k][j], posterior);
                else if (k > i && k < j)
                    ConsistencyTransform (sparseMatrices[i][k], sparseMatrices[k][j], posterior);
                else {
                    SparseMatrix *temp = sparseMatrices[j][k]->ComputeTranspose();
                    ConsistencyTransform (sparseMatrices[i][k], temp, posterior);
                    delete temp;
                }
            }

            // Renormalization
            for (int k = 0; k < (seq1Length+1) * (seq2Length+1); k++) posterior[k] /= numSeqs;

            // Mask out positions not originally in the posterior matrix
            SparseMatrix *matXY = sparseMatrices[i][j];
            for (int y = 0; y <= seq2Length; y++) posterior[y] = 0;
            for (int x = 1; x <= seq1Length; x++){
                SafeVector<PIF>::iterator XYptr = matXY->GetRowPtr(x);
                SafeVector<PIF>::iterator XYend = XYptr + matXY->GetRowSize(x);
                SafeVector<float>::iterator base = posterior.begin() + x * (seq2Length + 1);
                int curr = 0;
                while (XYptr != XYend){
                    // Zero out all cells until the first filled column
                    while (curr < XYptr->first){
                    base[curr] = 0;
                    curr++;
                    }

                // Skip over this column
                curr++;
                ++XYptr;
                }
        
                // Zero out cells after last column
                while (curr <= seq2Length){
                    base[curr] = 0;
                    curr++;
                }
            }

            // Save the new posterior matrix
            newSparseMatrices[i][j] = new SparseMatrix (seq1->GetLength(), seq2->GetLength(), posterior);
            newSparseMatrices[j][i] = NULL;

            delete posteriorPtr;
        }
    }
    return newSparseMatrices;
}

//! This function performs randomized partitioning iterative refinement. 
//! Taking posterior probability matrices, parameters of probabilistic model, and multiple sequence alignments.
//! Returns a new multiple sequence alignment.
void DoIterativeRefinement (const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                            const ProbabilisticModel &model, MultiSequence* &alignment, int i){
    set<int> groupOne, groupTwo;
    randomnumber rn;
    rn.seed(1234+i);
    // Create two separate groups
    for (int i = 0; i < alignment->GetNumSequences(); i++){
        int x = rn.roll_int(1,10);
        //cout << "rand: " << x << '\t';
        if (x % 2) {
          groupOne.insert (i);
        }
        else
          groupTwo.insert (i);
    }

    if (groupOne.empty() || groupTwo.empty()) return;

    // Project into the two groups
    MultiSequence *groupOneSeqs = alignment->Project (groupOne); assert (groupOneSeqs);
    MultiSequence *groupTwoSeqs = alignment->Project (groupTwo); assert (groupTwoSeqs);
    delete alignment;

    // Realign
    alignment = AlignAlignments (groupOneSeqs, groupTwoSeqs, sparseMatrices, model);

    delete groupOneSeqs;
    delete groupTwoSeqs;
}

//! This function takes guide tree (computed by distance) as input.
//! Returns the aligned sequences corresponding to a node or leaf of a guide tree.
MultiSequence *ProcessTree (const TreeNode *tree, MultiSequence *sequences,
                            const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                            const ProbabilisticModel &model){
  MultiSequence *result;

  // Check if this is a node of the alignment tree
  if (tree->GetSequenceLabel() == -1){
    MultiSequence *alignLeft = ProcessTree (tree->GetLeftChild(), sequences, sparseMatrices, model);
    MultiSequence *alignRight = ProcessTree (tree->GetRightChild(), sequences, sparseMatrices, model);

    assert (alignLeft);
    assert (alignRight);

    result = AlignAlignments (alignLeft, alignRight, sparseMatrices, model);
    assert (result);

    delete alignLeft;
    delete alignRight;
  }

  // Otherwise, this is a leaf of the alignment tree
  else {
    result = new MultiSequence(); assert (result);
    result->AddSequence (sequences->GetSequence(tree->GetSequenceLabel())->Clone());
  }

  return result;
}

//! This function computes the final alignment by calling ProcessTree() and performing iterative refinement.
MultiSequence *ComputeFinalAlignment (const TreeNode *tree, MultiSequence *sequences,
                                      const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                                      const ProbabilisticModel &model){

    MultiSequence *alignment = ProcessTree (tree, sequences, sparseMatrices, model);

    // Iterative refinement
    for (int i = 0; i < numIterativeRefinementReps; i++)
        DoIterativeRefinement (sparseMatrices, model, alignment, i);

    // Return final alignment
    return alignment;
}
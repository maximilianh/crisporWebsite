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


#include "SafeVector.h"
#include "MultiSequence.h"
#include "ProbabilisticModel.h"
#include "GuideTree.h"
#include "SparseMatrix.h"
#include "../src/random.h"

#define numConsistencyReps 2     // Number of probabilistic consistency transformation.
#define numIterativeRefinementReps 100   // Number of iterative refinement of the multiple alignment.

#ifdef TURBOHOMOLOGY
	#define AlignAlignments AlignAlignmentsturbohomology
	#define AlignAlignmentsMappingSeq AlignAlignmentsMappingSeqturbohomology
	#define AlignProfile AlignProfileturbohomology
	#define ConsistencyTransform ConsistencyTransformturbohomology
	#define ConsistencyTransform1  ConsistencyTransform1turbohomology
	#define MultiConsistencyTransfor MultiConsistencyTransforturbohomology
	#define DoIterativeRefinement DoIterativeRefinementturbohomology
	#define ProcessTree ProcessTreeturbohomology
	#define ComputeFinalAlignment ComputeFinalAlignmentturbohomology

#endif

//! This function takes two multiple sequence alignments as input.
//! Returns the alignment of the two MultiSequence objects.
//! \param align1 and \param align2 are two multiple sequence alignments.
//! \param sparseMatrices is matrices storing posterior probability for each pair of sequences.
//! \param model is storing the parameters of a probabilistic model. (initDistrib, gapOpen, gapExtend, emitPairs, emitSingle) 
MultiSequence *AlignAlignments (MultiSequence *align1, MultiSequence *align2,
                                const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                                const ProbabilisticModel &model);

#ifdef TURBOHOMOLOGY
MultiSequence *AlignAlignmentsMappingSeq (MultiSequence *align1, MultiSequence *align2,
                                const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                                const ProbabilisticModel &model, vector<int> mappingSeqIndex,
                                const vector<vector<vector<double> > > basePairScore,
                                const vector<double> similarity_list);

MultiSequence *AlignProfile (Sequence *newSeq, MultiSequence *align, 
                              const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                              const ProbabilisticModel &model, vector<int> mappingSeqIndex,
                              const vector<vector<vector<double> > > basePairScore, 
                              const vector<double> similarity_list);
#endif

//! This function computes the consistency transformation for sequence Z by taking two posterior probabilities matrices of alignments between X-Z and Z-Y.
//! For the case that sequence Z's index is larger than that of X.
//! The transformed matrix is added to \param posterior.
void ConsistencyTransform (SparseMatrix *matXZ, SparseMatrix *matZY, SafeVector<float> &posterior);

//! Same as the above, but for the case that Z's index is smaller than that of X.
void ConsistencyTransform1 (SparseMatrix *matZX, SparseMatrix *matZY, SafeVector<float> &posterior);

//! This function takes multiple sequences and posterior probability matices to perform three-way probabilistic consistency transformation.
//! Returns new re-estimated alignment score matrices. 
//! The formula is: P'(x[i]-y[j])=(1/|S|)*sum_z_in_S{ sum_k{ P(x[i]-z[k]) * P(z[k]-y[j]) } }
//! \param sequences is multiple sequences.
//! \param sparseMatrices is posterior probaility from pairwise HMM alignment.
SafeVector<SafeVector<SparseMatrix *> > MultiConsistencyTransform(MultiSequence *sequences, 
                                                      SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices);

//! This function performs randomized partitioning iterative refinement. 
//! Taking posterior probability matrices, parameters of probabilistic model, and multiple sequence alignments.
//! Returns a new multiple sequence alignment.
void DoIterativeRefinement (const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                            const ProbabilisticModel &model, MultiSequence* &alignment, int i);

//! This function takes guide tree (computed by distance) as input.
//! Returns the aligned sequences corresponding to a node or leaf of a guide tree.
MultiSequence *ProcessTree (const TreeNode *tree, MultiSequence *sequences,
                            const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                            const ProbabilisticModel &model);

//! This function computes the final alignment by calling ProcessTree() and performing iterative refinement.
MultiSequence *ComputeFinalAlignment (const TreeNode *tree, MultiSequence *sequences,
                                      const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                                      const ProbabilisticModel &model);
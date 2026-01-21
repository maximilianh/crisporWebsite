
#include "GuideTree.h"
#include <string>
#include <list>
#include <stdio.h>
#include "SafeVector.h"
#include "MultiSequence.h"
#include "Sequence.h"

using namespace std;

//! TreeNode Class
/*!
    The TreeNode Class provides unit for representing an alignment tree. 
*/

//! Constructor for a tree node. 
//! Note that sequenceLabel = -1 implies that the current node is not a leaf in the tree.
TreeNode::TreeNode(int sequenceLabel) : sequenceLabel(sequenceLabel), left(NULL), right(NULL), parent(NULL) {
    assert(sequenceLabel >= -1);
}

//! Destructor for a tree node.  Recursively deletes all children.
TreeNode::~TreeNode() {
    if(left){
        delete left; 
        left = NULL;
    }
    if(right){
        delete right;
        right = NULL;
    }
    parent = NULL;
}

// Getters
int TreeNode::GetSequenceLabel() const { return sequenceLabel; }
TreeNode * TreeNode::GetLeftChild() const { return left; }
TreeNode * TreeNode::GetRightChild() const { return right; }
TreeNode * TreeNode::GetParent() const { return parent; }

// Setters
void TreeNode::SetSequenceLabel(int sequenceLabel){ 
    this->sequenceLabel = sequenceLabel; 
    assert (sequenceLabel >= -1); 
}
void TreeNode::SetLeftChild(TreeNode *left){ this->left = left; }
void TreeNode::SetRightChild(TreeNode *right){ this->right = right; }
void TreeNode::SetParent(TreeNode *parent){ this->parent = parent; }

//! Computes a guide tree based on a given distance matrix. 
TreeNode * TreeNode::ComputeTree(const SafeVector<SafeVector<float> > &distMatrix){

    // Number of sequences
    int numSeqs = distMatrix.size();
    // A copy of distance matrix.
    SafeVector<SafeVector<float> > distances(numSeqs, SafeVector<float> (numSeqs));
    SafeVector<TreeNode *> nodes (numSeqs, NULL);
    // valid[i] tells whether or not the ith nodes in the distances and nodes array are valid.
    SafeVector<int> valid (numSeqs, 1);

    // initialization: make a copy of the distance matrix
    for (int i = 0; i < numSeqs; i++)
        for (int j = 0; j < numSeqs; j++)
            distances[i][j] = distMatrix[i][j];

    // initialization: create all the leaf nodes
    for (int i = 0; i < numSeqs; i++){
        nodes[i] = new TreeNode (i);
        assert (nodes[i]);
    }

    // repeat until only a single node left
    for (int numNodesLeft = numSeqs; numNodesLeft > 1; numNodesLeft--){
        float bestProb = -1;
        pair<int,int> bestPair;

        // find the closest pair
        for (int i = 0; i < numSeqs; i++) if (valid[i]){
            for (int j = i+1; j < numSeqs; j++) if (valid[j]){
                if (distances[i][j] > bestProb){
                    bestProb = distances[i][j];
                    bestPair = make_pair(i, j);
                }
            }
        }

        // merge the closest pair
        TreeNode *newParent = new TreeNode (-1);
        newParent->SetLeftChild (nodes[bestPair.first]);
        newParent->SetRightChild (nodes[bestPair.second]);
        nodes[bestPair.first]->SetParent (newParent);
        nodes[bestPair.second]->SetParent (newParent);
        nodes[bestPair.first] = newParent;
        nodes[bestPair.second] = NULL;

        // now update the distance matrix
        for (int i = 0; i < numSeqs; i++) if (valid[i]){
            distances[bestPair.first][i] = distances[i][bestPair.first]
                = (distances[i][bestPair.first] + distances[i][bestPair.second]) * bestProb / 2;
        }

        // mark the second node entry as no longer valid
        valid[bestPair.second] = 0;
    }

    assert (nodes[0]);
    return nodes[0];
}


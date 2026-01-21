#ifndef _GUIDETREE_
#define _GUIDETREE_

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

class TreeNode {

    int sequenceLabel;                  // sequence label
    TreeNode *left, *right, *parent;    // pointers to left, right children

public:

    //! Constructor for a tree node. 
    //! Note that sequenceLabel = -1 implies that the current node is not a leaf in the tree.
    TreeNode(int sequenceLabel);
    ~TreeNode();


    // Getters
    int GetSequenceLabel() const ;
    TreeNode *GetLeftChild() const ;
    TreeNode *GetRightChild() const ;
    TreeNode *GetParent() const ;

    // Setters
    void SetSequenceLabel(int sequenceLabel);
    void SetLeftChild(TreeNode *left);
    void SetRightChild(TreeNode *right);
    void SetParent(TreeNode *parent);

    //! Computes a guide tree based on a given distance matrix. 
    static TreeNode *ComputeTree(const SafeVector<SafeVector<float> > &distMatrix);
};

#endif
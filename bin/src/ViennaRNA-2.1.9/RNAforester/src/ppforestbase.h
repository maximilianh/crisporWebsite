/*
  Copyright by Matthias Hoechsmann (C) 2002-2004
  =====================================
  You may use, copy and distribute this file freely as long as you
  - do not change the file,
  - leave this copyright notice in the file,
  - do not make any profit with the distribution of this file
  - give credit where credit is due
  You are not allowed to copy or distribute this file otherwise
  The commercial usage and distribution of this file is prohibited
  Please report bugs and suggestions to <mhoechsm@TechFak.Uni-Bielefeld.DE>
*/

#ifndef _PPFOREST_BASE_H
#define _PPFOREST_BASE_H

#include <algorithm>
#include <cassert>
#include <cstring>

#include "misc.h"
#include "types.h"

using namespace std;

/** PPForestBase is the base class of the template class PPForest<L>.
 *  To reduce the size of compiled programs functions and variables that are
 *  independent of the labelling are implemented in this class.
 */
class PPForestBase
{
  //template <class R, class LX, class AL>
  //friend class Alignment;

 public:
  typedef Uint size_type;

 private:
  bool m_isSumUpCSFValid;
  bool m_isRMBValid;

  size_type getNumRightBrothers(size_type i) const;
  size_type countLeaves(size_type i,size_type k) const;

protected:
  size_type m_size;         /**< size is the number of nodes of a tree */
  size_type *m_rb;          /**< m_rb[i] stores the preorder index of the rightbrother node of the ith node, or 0 if there is none. */
  size_type *m_noc;         /**< m_noc[i] stores the  number of children of the ith node */
  size_type *m_sumUpCSF;    /**< m_sumUpCSF[i] stores the sum of non empty closed subforests of nodes k<i. This is used to map index pairs to single indexes. */
  size_type *m_rmb;         /**< m_rmb[i] stores the index of the rightmost brother of the ith node. */

  /** Calculate m_sumUpCSF from m_rb and m_noc. */
  void calcSumUpCSF();
  /** Calculate m_rmb from m_rb and m_noc */
  void calcRMB();

  /** @name Alignment construction functions for backtrack routine of class Alignment */
  //@{

  /** Allocate memory and initialize variables */
  void initialize(size_type size);

  inline void setSize(size_type size) {m_size=size;};
  inline void setNumChildren(size_type node,size_type num) {m_noc[node]=num;};
  inline void setRightBrotherIndex(size_type node,size_type index) {m_rb[node]=index;};

  //@}

 public:
  /** Default constructor creates an empty forest */
  PPForestBase() : m_isSumUpCSFValid(false), m_isRMBValid(false), m_size(0),m_rb(NULL),m_noc(NULL),m_sumUpCSF(NULL),m_rmb(NULL){};
  /** Allocates the memory for a forest of "size" nodes */
  PPForestBase(size_type size);
  /** Copy constructor */
  PPForestBase(const PPForestBase &ppfBase);
  /** Frees allocated memory */
  ~PPForestBase();

  inline size_type size() const {return m_size;};               /**< Returns m_size */
  inline size_type noc(size_type i) const {return m_noc[i];};   /**< Returns m_noc[i] */
  inline size_type rb(size_type i) const {return m_rb[i];};     /**< Returns m_rb[i] */
  size_type rb(size_type i, size_type k) const;                 /**< Returns the kth brother of i */


  inline bool isLeave(size_type i) const
    {
      if(m_noc[i])
	return false;
      else
	return true;
    };

  inline bool isInternalNode(size_type i) const {return !isLeave(i);};

  /** Calculates the maximal length of a closed subforest with node index i in constant time.
   *  len(i)=sumUpCSF[i+1]-sumUpCSF[i], if i<size(F) and 1 otherwise.
   *  Requires that function calcSumUpCSF() was already called.
   */
  inline size_type getMaxLength(size_type i) const
    {
      assert(m_isSumUpCSFValid);

      if(m_size==1)
	return 1;
      else
	return m_sumUpCSF[i+1]-m_sumUpCSF[i]; /* len(i)=sumUpCSF[i+1]-sumUpCSF[i], if i<nodes(F)*/
    }

  /** Returns the number of all non empty closed subforests in a forest
   *  Requires that function calcSumUpCSF() was already called.
   */
  inline size_type getNumCSFs() const
    {
      assert(m_isSumUpCSFValid);
      return m_sumUpCSF[m_size];  // max node is m_size-1, hence m_sumUpCSF[m_size] stores the number of all csfs
    };

  /** Returns the index if the rightmost brother node of i in constant time.
   *  Requires that function calcRMB() was already called.
   */
  inline size_type getRightmostBrotherIndex(size_type i) const
    {
      assert(m_isRMBValid);
      assert(i<m_size);
      return m_rmb[i];
    };

  /** Returns the index if the rightmost brother node of i.
   *  Does not require that function calcRMB() was already called and runs linear in the number of brother nodes.
   */
  inline size_type getRightmostBrotherIndex2(size_type i) const
    {
      size_type h=i;
      while(m_rb[h])
	h=m_rb[h];

      return h;
    };

  size_type maxDegree() const;
  size_type numLeaves() const;
  size_type maxDepth() const;

  /** Returns the number of leaf nodes until the leftmost leave node of i.
   *  This is useful to calculate the nucleotide positions of local alignments.
   */
  size_type countLeaves(size_type i) const;

  /** @name Index transition functions */
  //@{

  /** Calculates one dimensional index for a closed subforest (i,j) */
  inline size_type indexpos(size_type i,size_type j) const
    {
      assert(m_isSumUpCSFValid);

      if(j==0)
	return 0;
      else
	return m_sumUpCSF[i]+j-1;
    }

  /** Returns indexpos(i+1,noc(i)) for a csf (i,j) */
  inline size_type down(size_type i) const {return indexpos(i+1,m_noc[i]);};
  /** Returns indexpos(rb(i),j-1) for a csf (i,j) */
  inline size_type over(size_type i,size_type j) const {return indexpos(m_rb[i],j-1);}
  /** Returns indexpos(rb(i),getMaxLength(i)-1) for a node index i */
  inline size_type over(size_type i) const {return indexpos(m_rb[i],getMaxLength(i)-1);}

  /** Returns the mapped index of csf (i2,..,in-1) where i1,...,in are the children of i.
    * This transition is important for the extended alignment model for RNA structures where
    * a P-node and the pairing bases can be replaced by a single edit operation
    */
  inline size_type mdown(size_type i) const
    {
      if(m_noc[i]<=2)
	return 0;
      else
	return indexpos(i+1+1,m_noc[i]-2);
    };

  //@}


  /** Tests if (i,j) and (i2,j2) are disjoint.
   *  Intuitively, (i,j) and (i2,j2) are disjoint if they do not include
   *  each other and dont intersect.
   */
  inline bool isDisjoint(size_type i, size_type j, size_type i2, size_type j2) const
    {
      size_type max_node;

      // the empty forest is included in any forest
      if(j==0 || j2==0)
	return false;

      // (i2,j2) included in (i,j)

      max_node=getRightmostBrotherIndex(i);
      if(noc(max_node))
	max_node=getRightmostBrotherIndex(max_node+1);

      if(i2>=i && i2<=max_node)
	return false;

      if(rb(i2,j2-1)>=i && rb(i2,j2-1)<=max_node)
	return false;

      // (i,j) included in (i2,j2)
      max_node=getRightmostBrotherIndex(i2);
	if(noc(max_node))
	  max_node=getRightmostBrotherIndex(max_node+1);


	if(i>=i2 && i<=max_node)
	  return false;

	if(rb(i,j-1)>=i2 && rb(i,j-1)<=max_node)
	  return false;

	return true;
    }
};

#endif

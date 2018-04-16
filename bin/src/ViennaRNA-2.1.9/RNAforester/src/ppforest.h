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

#ifndef _PPFOREST_H_
#define _PPFOREST_H_

#include <string>
#include <iostream>

#include "algebra.h"
//#include "alignment.h"
#include "misc.h"
#include "ppforestbase.h"

using namespace std;

//template<class R,class L,class AL>
//class Alignment;

/** The template class PPForest<L> is the basic class that is handled by the alignmnent
 *  algorithms of the RNAforester template library. The template parameter L is the
 *  datatype of the forest labels.
 */
template <class L>
class PPForest : public PPForestBase
{
  template <class R, class LX, class AL>
  friend class Alignment;

  public:
  typedef typename PPForestBase::size_type size_type;

private:
	/** internal function used by '<<' operator to print bracket notation of the forest*/
	void print(ostream &s,size_type node,bool brackets) const;

	/** Calculate maximal score of csf (i,j) for a certain Algebra.
	*  It is assumed that a perfect match of the forest obtains the best
	*  possible score.
	*/
	template <class R>
	R maxScore(const Algebra<R,L> &alg, size_type i, size_type j) const;
	/** Calculate maximal score of csf (i,j) for a certain RNAAlgebra.
	*  It is assumed that a perfect match of the forest obtains the best
	*  possible score.
	*/
	template <class R>
	R maxScore(const RNA_Algebra<R,L> &alg, size_type i, size_type j) const;

	/** Function showLabel is used by print routines of PPForest */
//	virtual void showLabel(ostream &s,char c) const {s << c;};
	virtual void showLabel(ostream &s,L l) const {s << 'X';};

	/** Function makeLabel is used by function buildForest */
//	virtual void makeLabel(char &a,char c) {a=c;};

	void makeLabel(L &a,char c) {};


protected:
	L *m_lb;	

	void initialize(size_type size);

public:
	typedef L label_type;

	/** Default Constructor */
	PPForest() : PPForestBase(), m_lb(NULL) {};
	/** Construct a PPForest with 'size' nodes */
	PPForest(size_type size);
	/** Copy constructor */
	PPForest(const PPForest<L> &ppf);
	/** Read PPForest from binary file */
	PPForest(istream &s);

	virtual ~PPForest();	

	/** returns label of node i */
	inline L label(size_type i) const {return m_lb[i];};

	/** Calculate maximal score of a forest alignment against itself for a certain Algebra.
	*  It is assumed that a perfect match of the forest obtains the best
	*  possible score.
	*/
	template <class R>
	inline R maxScore(const Algebra<R,L> &alg) const
	{
		return maxScore(alg,0,getMaxLength(0));
	}  

	/** Calculate maximal score of a forest alignment against itself for a certain RNA_Algebra.
	*  It is assumed that a perfect match of the forest obtains the best
	*  possible score.
	*/
	template <class R>
	inline R maxScore(const RNA_Algebra<R,L> &alg) const
	{
		return maxScore(alg,0,getMaxLength(0));
	}  

	/** Print forest in GraphViz format */
	void printDot(ostream &s) const;

	/** Save forest to binary file */
	void save(ostream &s) const;
};

/** Stream output operator for PPForest<L> */
template <class L>
ostream& operator<< (ostream &s, PPForest<L> &ppf);

#endif


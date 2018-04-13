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

#ifndef _PPFORESTALI_H_
#define _PPFORESTALI_H_

//#include "alignment.h"
#include "ppforest.h"

//template<class R,class L,class AL>
//class Alignment;

/** PPForestAli is the Base class of forest Alignemnts.
 *  The pure virtual functions must be implemented when 
 *  inherited to allow the construction of alignments.
 */

template <class L,class AL>
class PPForestAli : public PPForest<AL>
{
public:
  typedef typename PPForest<AL>::size_type size_type;

  //private:
public:
	/** @name Interface functions for alignment construction */
	//@{ 

	virtual void makeRepLabel(size_type node,L a, L b) = 0;
	virtual void makeDelLabel(size_type node,L a) = 0;
	virtual void makeInsLabel(size_type node,L b) = 0;

	//@}

public:
	PPForestAli() : PPForest<AL>(){};
	PPForestAli(size_type size) : PPForest<AL>(size) {};

	//  virtual ~PPForestAli(){};
};

#endif

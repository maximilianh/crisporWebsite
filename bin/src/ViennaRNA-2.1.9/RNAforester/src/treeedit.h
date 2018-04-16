/*
  Copyright by Matthias Hoechsmann (C) 2002
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

#ifndef _TREE_EDIT_H_
#define _TREE_EDIT_H_

#include <cassert>
#include <list>

#include "types.h"
#include "algebra.h"
#include "ppforestsz.h"
//#include "ppforestali.h"

template<class R,class L>
class Mapping
{
    R *m_mtrxTD;
	R *m_mtrxFD;
    Ulong *m_rowStartTD;
	Ulong *m_rowStartFD;
    Ulong m_mtrxSizeTD;
	Ulong m_mtrxSizeFD;
	R m_optimum;

    PPForestSZ<L> *m_ppfx;
    PPForestSZ<L> *m_ppfy;
    const SZAlgebra<R,L> *m_alg;

    inline R getMtrxTDVal(Ulong i,Ulong j) const          
      {
		assert(m_rowStartTD[i]+j<m_mtrxSizeTD);
		return m_mtrxTD[m_rowStartTD[i]+j];
      };

    inline void setMtrxTDVal(Ulong i,Ulong j,R val)
      {
		assert(m_rowStartTD[i]+j<m_mtrxSizeTD);
		m_mtrxTD[m_rowStartTD[i]+j]=val;
      };

   inline R getMtrxFDVal(Ulong i,Ulong j) const          
      {
		assert(m_rowStartFD[i]+j<m_mtrxSizeFD);
		return m_mtrxFD[m_rowStartFD[i]+j];
      };

    inline void setMtrxFDVal(Ulong i,Ulong j,R val)
      {
		assert(m_rowStartFD[i]+j<m_mtrxSizeFD);
		m_mtrxFD[m_rowStartFD[i]+j]=val;
      };

//    Uint backtrack(PPForestAli<L,AL> &ppf,Uint i, Uint j, Uint k, Uint l, Uint &node);
public:    
	Mapping(const PPForestSZ<L> *ppfx, const PPForestSZ<L> *ppfy,const SZAlgebra<R,L> &alg);
    ~Mapping();
    
	R getGlobalOptimum() {return m_optimum;};
};

#endif

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

#ifndef _ALIGNMENT_H_
#define _ALIGNMENT_H_

#include <cassert>
#include <list>

#include "types.h"
#include "algebra.h"
#include "ppforest.h"
#include "ppforestali.h"

template<class R,class L,class AL>
class Alignment
{
private:
	struct CSFPair
	{
		Uint i;
		Uint j;
		Uint k;
		Uint l;

		CSFPair(Uint i, Uint j, Uint k, Uint l) : i(i),j(j),k(k),l(l) {};
	};

	R *m_mtrx;
	Ulong *m_rowStart;
    Ulong m_mtrxSize;

    PPForest<L> *m_ppfx;
    PPForest<L> *m_ppfy;

    const Algebra<R,L> *m_alg;
    const RNA_Algebra<R,L> *m_rnaAlg;

    R m_localOptimum;
    R m_localSubOptimum;
    std::list<CSFPair> m_localAlis;	    // holds which alignments are already produced by getOptLocalAlignment
    double m_suboptimalsPercent;            // percentage is stores as a value between 0 and 1

    inline R getMtrxVal(Ulong i,Ulong j) const          
      {
	assert(m_rowStart[i]+j<m_mtrxSize);
	return m_mtrx[m_rowStart[i]+j];
      };

    inline void setMtrxVal(Ulong i,Ulong j,R val)
      {
	assert(m_rowStart[i]+j<m_mtrxSize);

	m_mtrx[m_rowStart[i]+j]=val;
	
	m_localOptimum=m_alg->choice(m_localOptimum,val);   // here we calculate local similarity on the fly
      };

    void calculateLocal(const PPForest<L> *ppfx, const PPForest<L> *ppfy,const Algebra<R,L> &alg, bool noSpeedup=false);
    void calculateGlobal(const PPForest<L> *ppfx, const PPForest<L> *ppfy,const Algebra<R,L> &alg, bool noSpeedup=false);
    Uint backtrack(PPForestAli<L,AL> &ppf,Uint i, Uint j, Uint k, Uint l, Uint &node);
public:    
    Alignment(const PPForest<L> *ppfx, const PPForest<L> *ppfy,const Algebra<R,L> &alg, bool local, bool noSpeedup=false);
    Alignment(const PPForest<L> *ppfx, const PPForest<L> *ppfy,const RNA_Algebra<R,L> &rnaAlg);
    ~Alignment();
    
    R getGlobalOptimum();

    /** This function should only be used in conjunction with
	 *  similarity based algebras
	 */
    double getGlobalOptimumRelative();

    inline R getLocalOptimum(){return m_localSubOptimum;};
    R getSILOptimum();
    void getOptGlobalAlignment(PPForestAli<L,AL> &ppfali);

    void resetOptLocalAlignment(int suboptimalsPercent=100);
    bool nextLocalSuboptimum();
    void getOptLocalAlignment(PPForestAli<L,AL> &ppfali,Uint &xbasepos, Uint &ybasepos);
    void getOptSILAlignment(PPForestAli<L,AL> &ppfali,Uint &ybasepos);
};

#endif

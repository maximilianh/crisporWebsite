/*
  Copyright by Matthias Hoechsmann (C) 2002-2003
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

#ifndef _TREE_EDIT_T_CPP_
#define _TREE_EDIT_T_CPP_

#include <algorithm>
#include <cassert>

#include "alignment.h"
#include "debug.h"

#include "misc.t.cpp"
#include "ppforestsz.t.cpp"

template<class R,class L>
Mapping<R,L>::Mapping(const PPForestSZ<L> *ppfx, const PPForestSZ<L> *ppfy, const SZAlgebra<R,L> &alg)
{
  assert(ppfx != NULL);
  assert(ppfy != NULL);

  Ulong m,n,h,cols;
  Ulong i,j,k,l,r,s;
  R score,h_score;


  // alloc space for the score matrices
  m=ppfx->size();
  n=ppfy->size();

  m_mtrxSizeTD=(m)*(n);
  m_mtrxSizeFD=(m+1)*(n+1);	// +1 because index 0 means the empty forest
  m_mtrxTD=new R[m_mtrxSizeTD];
  m_mtrxFD=new R[m_mtrxSizeFD];
  m_rowStartTD=new Ulong[m];
  m_rowStartFD=new Ulong[m+1];

  m_ppfx = new PPForestSZ<L>(*ppfx);
  m_ppfy = new PPForestSZ<L>(*ppfy);
  m_alg=&alg;
  
  cols=n;
  m_rowStartTD[0]=0;
  for(h=1;h<m;h++)
    m_rowStartTD[h]=m_rowStartTD[h-1]+cols;

  cols=n+1;
  m_rowStartFD[0]=0;
  for(h=1;h<m+1;h++)
    m_rowStartFD[h]=m_rowStartFD[h-1]+cols;

  // calculate edit distance
  for(i=0;i<m;i++)						// for all pairs of subtrees
  {
    if(!m_ppfx->keyroot(i))
		continue;

	for(j=0;j<n;j++)
	{
		if(!m_ppfy->keyroot(j))
			continue;

		// calculate forest distance
		
		// edit to empty forest
		setMtrxFDVal(0,0,0);
		
		for(k=m_ppfx->lml(i),r=1;k<=i;k++,r++)		// r is the indexpos of k
			setMtrxFDVal(r,0,m_alg->del(m_ppfx->label(k),getMtrxFDVal(r-1,0)));
	
		for(l=m_ppfy->lml(j),s=1;l<=j;l++,s++)
			setMtrxFDVal(0,s,m_alg->insert(getMtrxFDVal(0,s-1),m_ppfy->label(l)));	// s is the indexpos of j

		for(k=m_ppfx->lml(i),r=1;k<=i;k++,r++)
		{
			for(l=m_ppfy->lml(j),s=1;l<=j;l++,s++)
			{
				// fdist(k,i,l,j)
				// lml(k)==lml(i) && lml(l)==lml(j)
				if(m_ppfx->lml(k)==m_ppfx->lml(i) && m_ppfy->lml(l)==m_ppfy->lml(j))
				{
					score=m_alg->replace(m_ppfx->label(k),getMtrxFDVal(r-1,s-1),m_ppfy->label(l));

					h_score=m_alg->del(m_ppfx->label(k),getMtrxFDVal(r-1,s));
					score=alg.choice(score,h_score);

					h_score=m_alg->insert(getMtrxFDVal(r,s-1),m_ppfy->label(l));
					score=alg.choice(score,h_score);

					setMtrxFDVal(r,s,score);
					setMtrxTDVal(k,l,score);
				}
				else
				{
					long p,q;
					p=m_ppfx->lml(k) - m_ppfx->lml(i);
					q=m_ppfy->lml(l) - m_ppfy->lml(j);

					score=getMtrxFDVal(p,q) + getMtrxTDVal(k,l);

					h_score=m_alg->del(m_ppfx->label(k),getMtrxFDVal(r-1,s));
					score=alg.choice(score,h_score);

					h_score=m_alg->insert(getMtrxFDVal(r,s-1),m_ppfy->label(l));
					score=alg.choice(score,h_score);

					setMtrxFDVal(r,s,score);
				}
			}
		}
	}
  }

  	// until here the original tree edit distance was calculated
	// to allow edit distances between forests, the distance for a virtual root node is calculated
	// leftmost leaf of the root node is the first node in postorder traversal
	// The distances are calculated until the last tree node in postorder traversal which then holds the result
	// for a virtual root node that is matched with zero cost

	// edit to empty forest
	setMtrxFDVal(0,0,0);
	
	for(k=0,r=1;k<m;k++,r++)		// r is the indexpos of k
		setMtrxFDVal(r,0,m_alg->del(m_ppfx->label(k),getMtrxFDVal(r-1,0)));

	for(l=0,s=1;l<n;l++,s++)
		setMtrxFDVal(0,s,m_alg->insert(getMtrxFDVal(0,s-1),m_ppfy->label(l)));	// s is the indexpos of j

	for(k=0,r=1;k<m;k++,r++)	
	{
		for(l=0,s=1;l<n;l++,s++)
		{

		      long p,q;
		      p=m_ppfx->lml(k);
		      q=m_ppfy->lml(l);
		      
		      score=getMtrxFDVal(p,q) + getMtrxTDVal(k,l);
		      
		      h_score=m_alg->del(m_ppfx->label(k),getMtrxFDVal(r-1,s));
		      score=alg.choice(score,h_score);
		      
		      h_score=m_alg->insert(getMtrxFDVal(r,s-1),m_ppfy->label(l));
		      score=alg.choice(score,h_score);
		      
		      setMtrxFDVal(r,s,score);						    
		}
	}

	m_optimum=getMtrxFDVal(m,n);
}

template<class R,class L>
Mapping<R,L>::~Mapping()
{
	DELETE(m_ppfx);
	DELETE(m_ppfy);

	delete[] m_mtrxTD;
	delete[] m_mtrxFD;
	delete[] m_rowStartTD;
	delete[] m_rowStartFD;
}
#endif

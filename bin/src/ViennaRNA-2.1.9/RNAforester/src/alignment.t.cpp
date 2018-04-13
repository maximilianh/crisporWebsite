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

#ifndef _ALIGNMENT_T_CPP_
#define _ALIGNMENT_T_CPP_

#include <algorithm>
#include <cassert>
#include <new>
#include <exception>
#include <limits.h>

#include "alignment.h"
#include "debug.h"

#include "misc.t.cpp"
#include "ppforest.t.cpp"

/* ****************************************** */
/*    Constructor and Destruktor functions    */
/* ****************************************** */

template<class R,class L,class AL>
Alignment<R,L,AL>::Alignment(const PPForest<L> *ppfx, const PPForest<L> *ppfy, const Algebra<R,L> &alg, bool local, bool noSpeedup)
: m_suboptimalsPercent(100)
{
  if(local)
    calculateLocal(ppfx,ppfy,alg,noSpeedup);
  else
    calculateGlobal(ppfx,ppfy,alg,noSpeedup);
}

template<class R,class L,class AL>
Alignment<R,L,AL>::Alignment(const PPForest<L> *ppfx, const PPForest<L> *ppfy, const RNA_Algebra<R,L> &rnaAlg)
: m_suboptimalsPercent(100)
{
	assert(ppfx != NULL);
	assert(ppfy != NULL);

	Ulong m,n,h,cols;
	long i,k;
	Uint j,l,r;

	R score,h_score;

	// alloc space for the score matrix, backtrack structure  and , if wanted, for the calculation-order-matrix
	// check for an overflow
	if (ppfx->getNumCSFs() > ULONG_MAX / ppfy->getNumCSFs()) { 
	  cerr << "Error: Overflow in calculation matrix multiplication. Calculation terminated." << endl;
	  exit(EXIT_FAILURE);
	}

	m_mtrxSize=ppfx->getNumCSFs()*ppfy->getNumCSFs();

	// maximum array size is 2GB
	if(m_mtrxSize > 2000000000 || ppfx->getNumCSFs() > 2000000000) {
	  cerr << "Error: Maximum array size of 2GB exceeded due to large input data. Calculation terminated." << endl;
	  exit(EXIT_FAILURE);
	}

	m_mtrx = new R[m_mtrxSize];
		
	m_rowStart=new Ulong[ppfx->getNumCSFs()];
	m_ppfx = new PPForest<L>(*ppfx);			// copy the ppforests
	m_ppfy = new PPForest<L>(*ppfy);
	
	m_rnaAlg=&rnaAlg;
	m_alg=(const Algebra<R,L>*)&rnaAlg;
	m_localOptimum=rnaAlg.worst_score();

	// initialize variables
	m=ppfx->size();
	n=ppfy->size();
	cols=ppfy->getNumCSFs();
	 
	m_rowStart[0]=0;
	for(h=1;h<ppfx->getNumCSFs();h++){
	  m_rowStart[h]=m_rowStart[h-1]+cols;
	}
		
	// align forests fx and fy

	// the easiest case .. 
	setMtrxVal(0,0,rnaAlg.empty());

	// align fx to the empty forest (fill first row of array)
	for(i=m-1;i>=0;i--)  // for all nodes in fx
	{
	  
		for(j=1;j<=ppfx->getMaxLength(i);j++)  // for all non empty csfs induced by i
		{
		  score = rnaAlg.del(ppfx->label(i),
				               getMtrxVal(ppfx->down(i),0),
				               getMtrxVal(ppfx->over(i,j),0));	  

			
			setMtrxVal(ppfx->indexpos(i,j),0,score);
			
		}
	}

	

	// align fy to the empty forest (fill first column of array)
	for(k=n-1;k>=0;k--)  // for all nodes in fx
	{
		for(l=1;l<=ppfy->getMaxLength(k);l++)  // for all non empty csfs induced by k
		{
			score = rnaAlg.insert(getMtrxVal(0,ppfy->down(k)),
					              ppfy->label(k),
					              getMtrxVal(0,ppfy->over(k,l)));

			
			setMtrxVal(0,ppfy->indexpos(k,l),score);
		}
	}

	// align the rest
	for(i=m-1;i>=0;i--)  // for all nodes in fx  
		for(k=n-1;k>=0;k--)  // for all nodes in fx
			for(j=1;j<=ppfx->getMaxLength(i);j++)    // for all non empty csfs induced by i
				for(l=1;l<=ppfy->getMaxLength(k);l++)  // for all non empty csfs induced by k
				{
					// basepair replacement
					if(ppfx->down(i) && ppfy->down(k))
					{
						// must be two P nodes !!
						h_score = rnaAlg.replacepair(ppfx->label(i+1),
							                         ppfy->label(k+1),
							                         getMtrxVal(ppfx->mdown(i),ppfy->mdown(k)),
							                         ppfx->label(ppfx->getRightmostBrotherIndex(i+1)),
							                         ppfy->label(ppfy->getRightmostBrotherIndex(k+1)),
							                         getMtrxVal(ppfx->over(i,j),ppfy->over(k,l)));
					}
					else
					{
						h_score = rnaAlg.replace(ppfx->label(i),
							                     getMtrxVal(ppfx->down(i),ppfy->down(k)),
							                     ppfy->label(k),
							                     getMtrxVal(ppfx->over(i,j),ppfy->over(k,l)));
					}

					score=h_score;

					// delete
					h=k;               // h is the node where the suffix of the split begins
					for(r=0;r<=l;r++)  // for all splits of fy
					{
						h_score = rnaAlg.del(ppfx->label(i),
							                 getMtrxVal(ppfx->down(i),ppfy->indexpos(k,r)),
							                 getMtrxVal(ppfx->over(i,j),ppfy->indexpos(h,l-r)));

						score=rnaAlg.choice(score,h_score);
						h=ppfy->rb(h);
					}

					// insert
					h=i;
					for(r=0;r<=j;r++) // for all splits of fx
					{
						h_score = rnaAlg.insert(getMtrxVal(ppfx->indexpos(i,r),ppfy->down(k)),
							                    ppfy->label(k),
							                    getMtrxVal(ppfx->indexpos(h,j-r),ppfy->over(k,l)));

						score=rnaAlg.choice(score,h_score);
						h=ppfx->rb(h);
					}

					setMtrxVal(ppfx->indexpos(i,j),ppfy->indexpos(k,l),score);
				}

	//      showArray(m_mtrx,ppfx->getNumCSFs(),ppfy->getNumCSFs());
	resetOptLocalAlignment(100);
}



template<class R,class L,class AL>
Alignment<R,L,AL>::~Alignment()
{
    delete[] m_mtrx;
    delete[] m_rowStart;
    delete m_ppfx;
    delete m_ppfy;
}



/* ****************************************** */
/*            Private functions               */
/* ****************************************** */

template<class R,class L,class AL>
void Alignment<R,L,AL>::calculateLocal(const PPForest<L> *ppfx, const PPForest<L> *ppfy, const Algebra<R,L> &alg, bool noSpeedup)
{
	assert(ppfx != NULL);
	assert(ppfy != NULL);

	Ulong m,n,h,cols;
	long i,k;
	Uint j,l,r;

	R score,h_score;

	// alloc space for the score matrix, backtrack structure  and , if wanted, for the calculation-order-matrix
	if (ppfx->getNumCSFs() > ULONG_MAX / ppfy->getNumCSFs()) { 
	  cerr << "Error: Overflow in calculation matrix multiplication. Calculation terminated." << endl;
	  exit(EXIT_FAILURE);
	}

	m_mtrxSize=ppfx->getNumCSFs()*ppfy->getNumCSFs();

	if(m_mtrxSize > 2000000000 || ppfx->getNumCSFs() > 2000000000) {
	  cerr << "Error: Maximum array size of 2GB exceeded due to large input data. Calculation terminated." << endl;
	  exit(EXIT_FAILURE);
	}

	m_mtrx=new R[m_mtrxSize];
	m_rowStart=new Ulong[ppfx->getNumCSFs()];
	m_ppfx = new PPForest<L>(*ppfx);				// copy the ppforests
	m_ppfy = new PPForest<L>(*ppfy);
	m_alg=&alg;
	m_rnaAlg=NULL;
	m_localOptimum=alg.worst_score();


	// initialize variables
	m=ppfx->size();
	n=ppfy->size();
	cols=ppfy->getNumCSFs();

	m_rowStart[0]=0;
	for(h=1;h<ppfx->getNumCSFs();h++)
		m_rowStart[h]=m_rowStart[h-1]+cols;

	// align forests fx and fy

	// the easiest case .. 
	setMtrxVal(0,0,alg.empty());

	// align fx to the empty forest (fill first row of array)
	for(i=m-1;i>=0;i--)  // for all nodes in fx
	{
		for(j=1;j<=ppfx->getMaxLength(i);j++)  // for all non empty csfs induced by i
		{
			score = alg.del(ppfx->label(i),
				            getMtrxVal(ppfx->down(i),0),
				            getMtrxVal(ppfx->over(i,j),0));	  

			setMtrxVal(ppfx->indexpos(i,j),0,score);
		}
	}

	// align fy to the empty forest (fill first column of array)
	for(k=n-1;k>=0;k--)  // for all nodes in fx
	{
		for(l=1;l<=ppfy->getMaxLength(k);l++)  // for all non empty csfs induced by k
		{
			score = alg.insert(getMtrxVal(0,ppfy->down(k)),
					            ppfy->label(k),
					            getMtrxVal(0,ppfy->over(k,l)));

			setMtrxVal(0,ppfy->indexpos(k,l),score);
		}
	}

	// align the rest
	for(i=m-1;i>=0;i--)  // for all nodes in fx  
		for(k=n-1;k>=0;k--)  // for all nodes in fx
			for(j=1;j<=ppfx->getMaxLength(i);j++)    // for all non empty csfs induced by i
				for(l=1;l<=ppfy->getMaxLength(k);l++)  // for all non empty csfs induced by k
				{
					// replace
					score = alg.replace(ppfx->label(i),
							            getMtrxVal(ppfx->down(i),ppfy->down(k)),
							            ppfy->label(k),
							            getMtrxVal(ppfx->over(i,j),ppfy->over(k,l)));
					
					// delete				       
					if(ppfx->noc(i)==0 && !noSpeedup)  // no child
					  {
					    h_score = alg.del(ppfx->label(i),0,
							              getMtrxVal(ppfx->over(i,j),ppfy->indexpos(k,l)));

						score=alg.choice(score,h_score);
					  }
					else					  
					  {	
					    if(ppfx->rb(i)==0 && !noSpeedup) // no right brother
					      {
						h_score = alg.del(ppfx->label(i),
								  getMtrxVal(ppfx->down(i),ppfy->indexpos(k,l)),
								  0);
						
						score=alg.choice(score,h_score);
					      }
					    else
					      {
						h=k;               // h is the node where the suffix of the split begins		   					
						for(r=0;r<=l;r++)  // for all splits of fy
						  {
						    h_score = alg.del(ppfx->label(i),
							              getMtrxVal(ppfx->down(i),ppfy->indexpos(k,r)),
							              getMtrxVal(ppfx->over(i,j),ppfy->indexpos(h,l-r)));
						    
						    score=alg.choice(score,h_score);
						    h=ppfy->rb(h);
						  }
					      }
					  }

					// insert
					if(ppfy->noc(k)==0 && !noSpeedup) // no child
					  {
					    h_score = alg.insert(0,
							                 ppfy->label(k),
							                 getMtrxVal(ppfx->indexpos(i,j),ppfy->over(k,l)));

						score=alg.choice(score,h_score);
					  }
					else
					  {
					    if(ppfy->rb(k)==0 && !noSpeedup) // no right brother
					      {
						h_score = alg.insert(getMtrxVal(ppfx->indexpos(i,j),ppfy->down(k)),
								     ppfy->label(k),
								     0);
						
						score=alg.choice(score,h_score);
					      }					 
					    else
					      {
						h=i;
						for(r=0;r<=j;r++) // for all splits of fx
						  {
						    h_score = alg.insert(getMtrxVal(ppfx->indexpos(i,r),ppfy->down(k)),
							                 ppfy->label(k),
									 getMtrxVal(ppfx->indexpos(h,j-r),ppfy->over(k,l)));
						    
						    score=alg.choice(score,h_score);
						    h=ppfx->rb(h);
						  }
					      }
					  }

					// set value
					setMtrxVal(ppfx->indexpos(i,j),ppfy->indexpos(k,l),score);
				}

	/*
					// delete
					h=k;               // h is the node where the suffix of the split begins
					for(r=0;r<=l;r++)  // for all splits of fy
					{
						h_score = alg.del(ppfx->label(i),
							              getMtrxVal(ppfx->down(i),ppfy->indexpos(k,r)),
							              getMtrxVal(ppfx->over(i,j),ppfy->indexpos(h,l-r)));

						score=alg.choice(score,h_score);
						h=ppfy->rb(h);
					}

					// insert
					h=i;
					for(r=0;r<=j;r++) // for all splits of fx
					{
						h_score = alg.insert(getMtrxVal(ppfx->indexpos(i,r),ppfy->down(k)),
							                 ppfy->label(k),
							                 getMtrxVal(ppfx->indexpos(h,j-r),ppfy->over(k,l)));

						score=alg.choice(score,h_score);
						h=ppfx->rb(h);
					}

					// set value
					setMtrxVal(ppfx->indexpos(i,j),ppfy->indexpos(k,l),score);
	*/


	//      showArray(m_mtrx,ppfx->getNumCSFs(),ppfy->getNumCSFs());
	resetOptLocalAlignment(100);
}

template<class R,class L,class AL>
void Alignment<R,L,AL>::calculateGlobal(const PPForest<L> *ppfx, const PPForest<L> *ppfy, const Algebra<R,L> &alg, bool noSpeedup)
{
	assert(ppfx != NULL);
	assert(ppfy != NULL);

	Ulong m,n,h,cols;
	long i,k;
	Uint j,l,r;

	R score,h_score;

	// alloc space for the score matrix, backtrack structure  and , if wanted, for the calculation-order-matrix
	if (ppfx->getNumCSFs() > ULONG_MAX / ppfy->getNumCSFs()) { 
	  cerr << "Error: Overflow in calculation matrix multiplication. Calculation terminated." << endl;
	  exit(EXIT_FAILURE);
	}

	m_mtrxSize=ppfx->getNumCSFs()*ppfy->getNumCSFs();

	if(m_mtrxSize > 2000000000 || ppfx->getNumCSFs() > 2000000000) {
	  cerr << "Error: Maximum array size of 2GB exceeded due to large input data. Calculation terminated." << endl;
	  exit(EXIT_FAILURE);
	}

	m_mtrx=new R[m_mtrxSize];
	m_rowStart=new Ulong[ppfx->getNumCSFs()];
	  m_ppfx = new PPForest<L>(*ppfx);				// copy the ppforests
	  m_ppfy = new PPForest<L>(*ppfy); 
	m_alg=&alg;
	m_rnaAlg=NULL;
	m_localOptimum=alg.worst_score();


	// initialize variables
	m=ppfx->size();
	n=ppfy->size();
	cols=ppfy->getNumCSFs();

	m_rowStart[0]=0;
	for(h=1;h<ppfx->getNumCSFs();h++)
		m_rowStart[h]=m_rowStart[h-1]+cols;

	// align forests fx and fy

	// the easiest case .. 
	setMtrxVal(0,0,alg.empty());

	// align fx to the empty forest (fill first row of array)
	for(i=m-1;i>=0;i--)  // for all nodes in fx
	{
		for(j=1;j<=ppfx->getMaxLength(i);j++)  // for all non empty csfs induced by i
		{
			score = alg.del(ppfx->label(i),
				            getMtrxVal(ppfx->down(i),0),
				            getMtrxVal(ppfx->over(i,j),0));	  

			setMtrxVal(ppfx->indexpos(i,j),0,score);
		}
	}

	// align fy to the empty forest (fill first column of array)
	for(k=n-1;k>=0;k--)  // for all nodes in fx
	{
		for(l=1;l<=ppfy->getMaxLength(k);l++)  // for all non empty csfs induced by k
		{
			score = alg.insert(getMtrxVal(0,ppfy->down(k)),
					            ppfy->label(k),
					            getMtrxVal(0,ppfy->over(k,l)));

			setMtrxVal(0,ppfy->indexpos(k,l),score);
		}
	}

	// align the rest
	for(i=m-1;i>=0;i--)  // for all nodes in fx  
		for(k=n-1;k>=0;k--)  // for all nodes in fx
		        {
			        j=ppfx->getMaxLength(i);
				for(l=1;l<=ppfy->getMaxLength(k);l++)  // for all non empty csfs induced by k
				{
					// replace
		  
					score = alg.replace(ppfx->label(i),
							            getMtrxVal(ppfx->down(i),ppfy->down(k)),
							            ppfy->label(k),
							            getMtrxVal(ppfx->over(i,j),ppfy->over(k,l)));
					
					// delete				       
					if(ppfx->noc(i)==0 && !noSpeedup)  // no child
					  {
					    h_score = alg.del(ppfx->label(i),0,
							              getMtrxVal(ppfx->over(i,j),ppfy->indexpos(k,l)));

						score=alg.choice(score,h_score);
					  }
					else					  
					  {	
					    if(ppfx->rb(i)==0 && !noSpeedup) // no right brother
					      {
						h_score = alg.del(ppfx->label(i),
								  getMtrxVal(ppfx->down(i),ppfy->indexpos(k,l)),
								  0);
						
						score=alg.choice(score,h_score);
					      }
					    else
					      {
						h=k;               // h is the node where the suffix of the split begins		   					
						for(r=0;r<=l;r++)  // for all splits of fy
						  {
						    h_score = alg.del(ppfx->label(i),
							              getMtrxVal(ppfx->down(i),ppfy->indexpos(k,r)),
							              getMtrxVal(ppfx->over(i,j),ppfy->indexpos(h,l-r)));
						    
						    score=alg.choice(score,h_score);
						    h=ppfy->rb(h);
						  }
					      }
					  }

					// insert
					if(ppfy->noc(k)==0 && !noSpeedup) // no child
					  {
					    h_score = alg.insert(0,
							                 ppfy->label(k),
							                 getMtrxVal(ppfx->indexpos(i,j),ppfy->over(k,l)));

						score=alg.choice(score,h_score);
					  }
					else
					  {
					    if(ppfy->rb(k)==0 && !noSpeedup) // no right brother
					      {
						h_score = alg.insert(getMtrxVal(ppfx->indexpos(i,j),ppfy->down(k)),
								     ppfy->label(k),
								     0);
						
						score=alg.choice(score,h_score);
					      }					 
					    else
					      {
						h=i;
						for(r=0;r<=j;r++) // for all splits of fx
						  {
						    h_score = alg.insert(getMtrxVal(ppfx->indexpos(i,r),ppfy->down(k)),
							                 ppfy->label(k),
									 getMtrxVal(ppfx->indexpos(h,j-r),ppfy->over(k,l)));
						    
						    score=alg.choice(score,h_score);
						    h=ppfx->rb(h);
						  }
					      }
					  }

					// set value
					setMtrxVal(ppfx->indexpos(i,j),ppfy->indexpos(k,l),score);
				}

				l=ppfy->getMaxLength(k);
				for(j=1;j<=ppfx->getMaxLength(i);j++)    // for all non empty csfs induced by i
				{
					// replace
					score = alg.replace(ppfx->label(i),
							            getMtrxVal(ppfx->down(i),ppfy->down(k)),
							            ppfy->label(k),
							            getMtrxVal(ppfx->over(i,j),ppfy->over(k,l)));

					// delete				       
					if(ppfx->noc(i)==0 && !noSpeedup)  // no child
					  {
					    h_score = alg.del(ppfx->label(i),0,
							              getMtrxVal(ppfx->over(i,j),ppfy->indexpos(k,l)));

						score=alg.choice(score,h_score);
					  }
					else					  
					  {	
					    if(ppfx->rb(i)==0 && !noSpeedup) // no right brother
					      {
						h_score = alg.del(ppfx->label(i),
								  getMtrxVal(ppfx->down(i),ppfy->indexpos(k,l)),
								  0);
						
						score=alg.choice(score,h_score);
					      }
					    else
					      {
						h=k;               // h is the node where the suffix of the split begins		   					
						for(r=0;r<=l;r++)  // for all splits of fy
						  {
						    h_score = alg.del(ppfx->label(i),
							              getMtrxVal(ppfx->down(i),ppfy->indexpos(k,r)),
							              getMtrxVal(ppfx->over(i,j),ppfy->indexpos(h,l-r)));
						    
						    score=alg.choice(score,h_score);
						    h=ppfy->rb(h);
						  }
					      }
					  }

					// insert
					if(ppfy->noc(k)==0 && !noSpeedup) // no child
					  {
					    h_score = alg.insert(0,
							                 ppfy->label(k),
							                 getMtrxVal(ppfx->indexpos(i,j),ppfy->over(k,l)));

						score=alg.choice(score,h_score);
					  }
					else
					  {
					    if(ppfy->rb(k)==0 && !noSpeedup) // no right brother
					      {
						h_score = alg.insert(getMtrxVal(ppfx->indexpos(i,j),ppfy->down(k)),
								     ppfy->label(k),
								     0);
						
						score=alg.choice(score,h_score);
					      }					 
					    else
					      {
						h=i;
						for(r=0;r<=j;r++) // for all splits of fx
						  {
						    h_score = alg.insert(getMtrxVal(ppfx->indexpos(i,r),ppfy->down(k)),
							                 ppfy->label(k),
									 getMtrxVal(ppfx->indexpos(h,j-r),ppfy->over(k,l)));
						    
						    score=alg.choice(score,h_score);
						    h=ppfx->rb(h);
						  }
					      }
					  }

					// set value
					setMtrxVal(ppfx->indexpos(i,j),ppfy->indexpos(k,l),score);
				}
			}

	//      showArray(m_mtrx,ppfx->getNumCSFs(),ppfy->getNumCSFs());
	resetOptLocalAlignment(100);
}

template<class R,class L,class AL>
Uint Alignment<R,L,AL>::backtrack(PPForestAli<L,AL> &ppf,Uint i, Uint j, Uint k, Uint l, Uint &node)
{
	R score, h_score;
	Uint p_node,rb_node,noc,rbs,r,h;

	// empty alignment
	if(j==0 && l==0)
	{
		return 0;
	}

	score = getMtrxVal(m_ppfx->indexpos(i,j),m_ppfy->indexpos(k,l));
	p_node=node;
	node++;

	// could it be a replacement
	if(j>0 && l>0)
	{
		// check for basepair replacement only if Algebra is of type RNA_Algebra
		if(m_rnaAlg && m_ppfx->down(i) && m_ppfy->down(k))
		{
			h_score = m_rnaAlg->replacepair(m_ppfx->label(i+1),
				                            m_ppfy->label(k+1),
				                            getMtrxVal(m_ppfx->mdown(i),m_ppfy->mdown(k)),
				                            m_ppfx->label(m_ppfx->getRightmostBrotherIndex2(i+1)),
				                            m_ppfy->label(m_ppfy->getRightmostBrotherIndex2(k+1)),
				                            getMtrxVal(m_ppfx->over(i,j),m_ppfy->over(k,l)));      

			if(score == h_score)
			{
				// it is a basepair replacement
				ppf.makeRepLabel(p_node,m_ppfx->label(i),m_ppfy->label(k));   // P

				// set labels of leftmost child
				ppf.makeRepLabel(p_node+1,m_ppfx->label(i+1),m_ppfy->label(k+1));
				ppf.setRightBrotherIndex(p_node+1,p_node+1+1);
				ppf.setNumChildren(p_node+1,0);  // base node has no children
				node++;

				// down alignment
				assert(m_ppfx->noc(i)>=2);
				assert(m_ppfy->noc(k)>=2);
				noc=backtrack(ppf,i+1+1,m_ppfx->noc(i)-2,k+1+1,m_ppfy->noc(k)-2,node);  // !! mdown !!
				ppf.setNumChildren(p_node,noc+2);

				if(noc==0)
				{
					ppf.setRightBrotherIndex(p_node+1,p_node+1+1);
					ppf.setRightBrotherIndex(p_node+1+1,0);
				}
				else
				{
					ppf.setRightBrotherIndex(ppf.getRightmostBrotherIndex2(p_node+1+1),node);
				}


				// set labels of leftmost child
				ppf.makeRepLabel(node,m_ppfx->label(m_ppfx->getRightmostBrotherIndex2(i+1+1)),m_ppfy->label(m_ppfy->getRightmostBrotherIndex2(k+1+1)));
				ppf.setRightBrotherIndex(node,0);
				ppf.setNumChildren(node,0);  // base node has no children	  
				node++;
				rb_node=node; // !!	  

				// right alignment
				rbs=backtrack(ppf,m_ppfx->rb(i),j-1,m_ppfy->rb(k),l-1,node);
				if(rbs)
					ppf.setRightBrotherIndex(p_node,rb_node);
				else
					ppf.setRightBrotherIndex(p_node,0);

				return rbs+1;
			}
		}
		else
		{
			h_score = m_alg->replace(m_ppfx->label(i),
				                     getMtrxVal(m_ppfx->down(i),m_ppfy->down(k)),
				                     m_ppfy->label(k),
				                     getMtrxVal(m_ppfx->over(i,j),m_ppfy->over(k,l)));

			if(score == h_score)
			{
				// it is a replacement
				ppf.makeRepLabel(p_node,m_ppfx->label(i),m_ppfy->label(k));

				// down alignment
				noc=backtrack(ppf,i+1,m_ppfx->noc(i),k+1,m_ppfy->noc(k),node);
				ppf.setNumChildren(p_node,noc);
				rb_node=node;

				// right alignment
				rbs=backtrack(ppf,m_ppfx->rb(i),j-1,m_ppfy->rb(k),l-1,node);
				if(rbs)
					ppf.setRightBrotherIndex(p_node,rb_node);
				else
					ppf.setRightBrotherIndex(p_node,0);

				return rbs+1;
			}
		}
	}

	// could it be a deletion 
	if(j>0)
	{
		h=k;               // h is the node where the suffix of the split begins
		for(r=0;r<=l;r++)  // for all splits of fy
		{
			h_score = m_alg->del(m_ppfx->label(i),
				                 getMtrxVal(m_ppfx->down(i),m_ppfy->indexpos(k,r)),
				                 getMtrxVal(m_ppfx->over(i,j),m_ppfy->indexpos(h,l-r)));

			if(score == h_score)
			{
				// it is a deletion
				ppf.makeDelLabel(p_node,m_ppfx->label(i));

				// down alignment
				noc=backtrack(ppf,i+1,m_ppfx->noc(i),k,r,node);
				ppf.setNumChildren(p_node,noc);
				rb_node=node;

				// right alignment
				rbs=backtrack(ppf,m_ppfx->rb(i),j-1,h,l-r,node);
				if(rbs)
					ppf.setRightBrotherIndex(p_node,rb_node);
				else
					ppf.setRightBrotherIndex(p_node,0);

				return rbs+1;
			}

			//	  if(r<l)       // do not calculate rightbrother of h=0=no-rightbrother
			h=m_ppfy->rb(h);
		}
	}

	// could it be an insertion
	if(l>0)
	{
		h=i;
		for(r=0;r<=j;r++) // for all splits of fx
		{
			h_score = m_alg->insert(getMtrxVal(m_ppfx->indexpos(i,r),m_ppfy->down(k)),
				                    m_ppfy->label(k),
				                    getMtrxVal(m_ppfx->indexpos(h,j-r),m_ppfy->over(k,l)));

			if(score == h_score)
			{
				// it is an insertion
				ppf.makeInsLabel(p_node,m_ppfy->label(k));

				// down alignment
				noc=backtrack(ppf,i,r,k+1,m_ppfy->noc(k),node);
				ppf.setNumChildren(p_node,noc);
				rb_node=node;

				// right alignment
				rbs=backtrack(ppf,h,j-r,m_ppfy->rb(k),l-1,node);
				if(rbs)
					ppf.setRightBrotherIndex(p_node,rb_node);
				else
					ppf.setRightBrotherIndex(p_node,0);

				return rbs+1;      
			}

			//	  if(r<j)
			h=m_ppfx->rb(h);
		}
	}

	// you should never be here
	cerr << "Strange things happening in backtrack" << endl;
	exit(EXIT_FAILURE);
} 

/* ****************************************** */
/*             Public functions               */
/* ****************************************** */

template<class R,class L,class AL>
R Alignment<R,L,AL>::getGlobalOptimum()
{
  return getMtrxVal(m_ppfx->indexpos(0,m_ppfx->getMaxLength(0)),m_ppfy->indexpos(0,m_ppfy->getMaxLength(0)));
};

template<class R,class L,class AL>
double Alignment<R,L,AL>::getGlobalOptimumRelative()
{
  double opt;
  double max_x,max_y;

  opt=(double)getGlobalOptimum();

  if(m_rnaAlg)
    {
      max_x=(double)m_ppfx->maxScore(*m_rnaAlg);
      max_y=(double)m_ppfy->maxScore(*m_rnaAlg);
    }
  else
    {
      max_x=(double)m_ppfx->maxScore(*m_alg);
      max_y=(double)m_ppfy->maxScore(*m_alg);
    }

  assert(max_x+max_y>0);	
  opt=2*opt/(max_x+max_y);  

  return opt;
};

template<class R,class L,class AL>
R Alignment<R,L,AL>::getSILOptimum()
{
  R silOptimum;
  int k=0;
  Uint j,l=0;
  Ulong n;
  
  n=m_ppfy->size();

  if(m_rnaAlg)
    silOptimum=m_rnaAlg->worst_score();
  else
    silOptimum=m_alg->worst_score();

  j=m_ppfx->getMaxLength(0);

  // find the best match
  for(k=n-1;k>=0;k--)  // for all nodes in fx
    for(l=1;l<=m_ppfy->getMaxLength(k);l++)  // for all non empty csfs induced by k
      {
	silOptimum=m_alg->choice(silOptimum,getMtrxVal(m_ppfx->indexpos(0,j),m_ppfy->indexpos(k,l)));
      }

  return silOptimum;
}

template<class R,class L,class AL>
void Alignment<R,L,AL>::getOptGlobalAlignment(PPForestAli<L,AL> &ppfali)
{
  Uint node=0;

  // allocate a forest of the maximal size that a forest alignment can have
  ppfali.initialize(m_ppfx->size()+m_ppfy->size());
  backtrack(ppfali,0,m_ppfx->getMaxLength(0),0,m_ppfy->getMaxLength(0),node);
  ppfali.setSize(node);
  ppfali.calcSumUpCSF();
  ppfali.calcRMB();
}

template<class R,class L,class AL>
void Alignment<R,L,AL>::resetOptLocalAlignment(int suboptimalsPercent)
{
	m_localAlis.clear();
	m_suboptimalsPercent=suboptimalsPercent/100.0;
	m_localSubOptimum=m_localOptimum;
};

// calculate the score of the next best local alignment that is not "included" 
// in a local alignment returned by getOptLocalAlignment before 
template<class R,class L,class AL>
bool Alignment<R,L,AL>::nextLocalSuboptimum()
{
	Ulong m,n;
	int i=0,k=0;
	Uint j=0,l=0;
	
  m_localSubOptimum=m_alg->worst_score();
  m=m_ppfx->size();
  n=m_ppfy->size();

  // find a matrix element that is optimal
  for(i=m-1;i>=0;i--)  // for all nodes in fx  
    for(k=n-1;k>=0;k--)  // for all nodes in fx
      for(j=1;j<=m_ppfx->getMaxLength(i);j++)    // for all non empty csfs induced by i
		for(l=1;l<=m_ppfy->getMaxLength(k);l++)  // for all non empty csfs induced by k
		{
			bool disjoint=true;

			// check if i,j,k,l is included 
			typename list<CSFPair>::const_iterator it;
			for(it=m_localAlis.begin();it!=m_localAlis.end();it++)
			{
				if(!m_ppfx->isDisjoint(it->i,it->j,i,j) || !m_ppfy->isDisjoint(k,l,it->k,it->l))
				{
					disjoint=false;
					break;
				}
				
			}

			if(disjoint)
			{
				if(getMtrxVal(m_ppfx->indexpos(i,j),m_ppfy->indexpos(k,l))>=m_suboptimalsPercent*m_localOptimum)
					m_localSubOptimum=m_alg->choice(m_localSubOptimum,getMtrxVal(m_ppfx->indexpos(i,j),m_ppfy->indexpos(k,l)));
			}
		}

	return m_localSubOptimum!=m_alg->worst_score();	
}

template<class R,class L,class AL>
void Alignment<R,L,AL>::getOptLocalAlignment(PPForestAli<L,AL> &ppfali,Uint &xbasepos, Uint &ybasepos)
{
  Ulong m,n;
  int i=0,k=0;
  Uint j=0,l=0,node=0;

  m=m_ppfx->size();
  n=m_ppfy->size();
  // allocate a forest of the maximal size that a forest alignment can have
  ppfali.initialize(m_ppfx->size()+m_ppfy->size());

  // find a matrix element that is optimal
  for(i=m-1;i>=0;i--)  // for all nodes in fx  
    for(k=n-1;k>=0;k--)  // for all nodes in fx
      for(j=1;j<=m_ppfx->getMaxLength(i);j++)    // for all non empty csfs induced by i
		for(l=1;l<=m_ppfy->getMaxLength(k);l++)  // for all non empty csfs induced by k
		{
			bool disjoint=true;

			// check if i,j,k,l is included 
			typename list<CSFPair>::const_iterator it;
			for(it=m_localAlis.begin();it!=m_localAlis.end();it++)
			{
				if(!m_ppfx->isDisjoint(it->i,it->j,i,j) || !m_ppfy->isDisjoint(k,l,it->k,it->l))
				{
					disjoint=false;
					break;
				}
				
			}

			if(disjoint && getMtrxVal(m_ppfx->indexpos(i,j),m_ppfy->indexpos(k,l)) == m_localSubOptimum)
				goto found;
		}

 found: 
  backtrack(ppfali,i,j,k,l,node);
  ppfali.setSize(node);
  ppfali.calcSumUpCSF();
  ppfali.calcRMB();  

  m_localAlis.push_back(CSFPair(i,j,k,l));
  xbasepos=m_ppfx->countLeaves(i);
  ybasepos=m_ppfy->countLeaves(k);
}

template<class R,class L,class AL>
void Alignment<R,L,AL>::getOptSILAlignment(PPForestAli<L,AL> &ppfali,Uint &ybasepos)
{
  R silOptimum;
  Ulong m,n;
  int k=0;
  Uint j=0,l=0,node=0;

  m=m_ppfx->size();
  n=m_ppfy->size();
  // allocate a forest of the maximal size that a forest alignment can have
  ppfali.initialize(m_ppfx->size()+m_ppfy->size());

  silOptimum=getSILOptimum();
  j=m_ppfx->getMaxLength(0);
  
  // find a matrix element that is optimal
  for(k=n-1;k>=0;k--)  // for all nodes in fx
    for(l=1;l<=m_ppfy->getMaxLength(k);l++)  // for all non empty csfs induced by k
      {
	if(silOptimum == getMtrxVal(m_ppfx->indexpos(0,j),m_ppfy->indexpos(k,l)))
	  goto found;
      }

 found: 
  backtrack(ppfali,0,j,k,l,node);
  ppfali.setSize(node);
  ppfali.calcSumUpCSF();
  ppfali.calcRMB();  

  ybasepos=m_ppfy->countLeaves(k);
}


#endif























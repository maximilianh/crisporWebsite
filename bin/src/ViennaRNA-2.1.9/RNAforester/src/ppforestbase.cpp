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

#include "misc.h"
#include "ppforestbase.h"
#include "debug.h"

#include <algorithm>
#include <stdio.h>
#include <stddef.h>
#include <string.h>

#include "misc.t.cpp"

/* ****************************************** */
/*    Constructor and Destructor functions    */
/* ****************************************** */

PPForestBase::PPForestBase(size_type size)
{
  initialize(size);
};

PPForestBase::PPForestBase(const PPForestBase &ppfBase)
{
  initialize(ppfBase.size());

  memcpy(m_rb,ppfBase.m_rb,sizeof(size_type)*m_size);
  memcpy(m_noc,ppfBase.m_noc,sizeof(size_type)*m_size);
  memcpy(m_sumUpCSF,ppfBase.m_sumUpCSF,sizeof(size_type)*(m_size+1));
  memcpy(m_rmb,ppfBase.m_rmb,sizeof(size_type)*m_size);
  
  m_isSumUpCSFValid=ppfBase.m_isSumUpCSFValid; 
  m_isRMBValid=ppfBase.m_isRMBValid;
}

PPForestBase::~PPForestBase()
{
  DELETE_ARRAY(m_rb);
  DELETE_ARRAY(m_noc);
  DELETE_ARRAY(m_sumUpCSF);
  DELETE_ARRAY(m_rmb);
};

/* ****************************************** */
/*            Private functions               */
/* ****************************************** */

PPForestBase::size_type PPForestBase::countLeaves(size_type i) const
{
  size_type numLeaves=0;
  
  for(size_type k=0;k<i;k++)
    {
      if(isLeave(k))
	numLeaves++;
    }
  
  return numLeaves;
};


/* ****************************************** */
/*            Protected function              */
/* ****************************************** */

void PPForestBase::calcSumUpCSF()
{	
  // sum CSFs
  m_sumUpCSF[0] = 1;
  for(size_type i=1; i<m_size; i++)
  {
    m_sumUpCSF[i] = m_sumUpCSF[i-1]+getNumRightBrothers(i-1)+1;
  }
  
  m_sumUpCSF[m_size] = m_sumUpCSF[m_size-1]+1;  // last node is a leaf, so the number of siblings is one	
  
  m_isSumUpCSFValid=true;
}

void PPForestBase::calcRMB()
{
  for(int i=m_size-1;i>=0;i--)
    {
      if(m_rb[i])
	m_rmb[i]=m_rmb[m_rb[i]];
      else
	m_rmb[i]=i;
    }
  
  m_isRMBValid=true;
}

void PPForestBase::initialize(size_type size)
{ 
  m_size=size;
  m_rb=new size_type[size];
  m_noc=new size_type[size];
  m_sumUpCSF=new size_type[size+1];
  m_rmb=new size_type[size];

  m_isSumUpCSFValid=false;
  m_isRMBValid=false;
}

/* ****************************************** */
/*             Public functions               */
/* ****************************************** */

PPForestBase::size_type PPForestBase::rb(size_type i, size_type k) const
{
  if(k==0)
    return i;
  else
    return rb(rb(i),k-1);
}

PPForestBase::size_type PPForestBase::maxDegree() const
{
  size_type degree=0;
  
  for(size_type i=0;i<m_size;i++)
      degree=max(degree,getMaxLength(i));
  
  return degree;
}
  
PPForestBase::size_type PPForestBase::numLeaves() const
{
  size_type count=0;
  
  for(size_type i=0;i<m_size;i++)
    if(isLeave(i))
      count++;
  
  return count;
}

PPForestBase::size_type PPForestBase::getNumRightBrothers(size_type node) const
{
  size_type numRbs=0;
  
  while(m_rb[node])
    {
      numRbs++;  
      node=m_rb[node];
    }
  
  return numRbs;
}

PPForestBase::size_type PPForestBase::maxDepth() const
  {
    size_type *maxDepth;
    size_type m,h,r,i;
    
    maxDepth=new size_type[m_size];
    
    for(r=0, i=m_size-1;r<m_size;r++,i--)
      {
	if(isLeave(i))
	  maxDepth[i]=1;
	else
	  {
	    m=0;
	    h=i+1;		// h = first child
	    do
	      {
		m=max(m,maxDepth[h]);
		h=rb(h);
	      }
	    while(h);
	    
	    maxDepth[i]=1+m;
	  }
      }
    
    // calculate maxDepth for the forest
    m=0;
    h=0;		// h = first root
    do
      {
	m=max(m,maxDepth[h]);
	h=rb(h);
      }
    while(h);
    
    delete maxDepth;
    return m;
  }




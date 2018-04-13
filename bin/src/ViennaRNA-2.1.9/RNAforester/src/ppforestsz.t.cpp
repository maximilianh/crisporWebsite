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

#ifndef _PPFORESTSZ_T_CPP_
#define _PPFORESTSZ_T_CPP_

#include <map>

#include "misc.h"
#include "ppforestsz.h"


// PPForest<T>

template <class L> 
PPForestSZ<L>::PPForestSZ(Uint nrOfNodes)
: m_size(nrOfNodes)
{
  m_lml=new Uint[nrOfNodes];	
  m_lb=new L[nrOfNodes];
  m_keyroot=new bool[nrOfNodes];
}

template <class L>
PPForestSZ<L>::PPForestSZ(const PPForestSZ<L> &ppf)
{
  m_size=ppf.size();	

  m_lml=new Uint[ppf.size()];	
  m_lb=new L[ppf.size()];
  m_keyroot=new bool[ppf.size()];

  memcpy(m_lml,ppf.m_lml,sizeof(Uint)*m_size);
  memcpy(m_keyroot,ppf.m_keyroot,sizeof(bool)*m_size);

  for(Uint i=0;i<m_size;i++)
	m_lb[i]=ppf.m_lb[i];
} 

template <class L> 
PPForestSZ<L>::~PPForestSZ()
{
  DELETE_ARRAY(m_lml);
  DELETE_ARRAY(m_lb);
  DELETE_ARRAY(m_keyroot);
}

template <class L> 
void PPForestSZ<L>::calcKeyroots()
{
  std::map<Uint,Uint> keyrootMap;

  for(Uint i=0;i<m_size;i++)
  {
    m_keyroot[i]=false;
    keyrootMap[lml(i)]=i;
  }

  std::map<Uint,Uint>::const_iterator it;
  for(it=keyrootMap.begin();it!=keyrootMap.end();it++)
  {
	m_keyroot[it->second]=true;
  }

}

#endif

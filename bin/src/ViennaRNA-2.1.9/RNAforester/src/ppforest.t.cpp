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

#ifndef _PPFOREST_T_CPP_
#define _PPFOREST_T_CPP_

#include "misc.h"
#include "ppforest.h"

/* ****************************************** */
/*    Constructor and Destruktor functions    */
/* ****************************************** */

template <class L> 
PPForest<L>::PPForest(size_type size): PPForestBase(size)
{
	m_lb=new L[size];
}

template <class L>
PPForest<L>::PPForest(const PPForest<L> &ppf) : PPForestBase(ppf)
{
	m_lb=new L[ppf.size()];

	// copy labels
	for(size_type i=0;i<size();i++)
		m_lb[i]=ppf.m_lb[i];
} 

template <class L>
PPForest<L>::PPForest(istream &s)
{
  /*
  // read ppforest from file

  // first of all save the size
  s.read(static_cast<char*>(m_size),sizeof(PPForestBase::size_type));

  // initialize structures
  ::PPForestBase(m_size);  
  m_lb=new L[ppf.size()];

  // save the arrays
  s.read(static_cast<char*>(m_lb),sizeof(label_type)*m_size);
  s.read(static_cast<char*>(m_rb),sizeof(PPForestBase::size_type)*m_size);
  s.read(static_cast<char*>(m_noc),sizeof(PPForestBase::size_type)*m_size);
  s.read(static_cast<char*>(m_sumUpCSF),sizeof(PPForestBase::size_type)*m_size);
  s.read(static_cast<char*>(m_rmb),sizeof(PPForestBase::size_type)*m_size);    
  */
}

template <class L> 
PPForest<L>::~PPForest()
{
   DELETE_ARRAY(m_lb);
}

/* ****************************************** */
/*            Private functions               */
/* ****************************************** */

template<class L>
void PPForest<L>::print(ostream &s,size_type node,bool brackets) const
{
	if(brackets)
		s << "(";

	s <<  "{" << node << "}";
	showLabel(s,m_lb[node]);
	if(m_noc[node])
		print(s,node+1,true);
	if(m_rb[node])
		print(s,m_rb[node],false);

	if(brackets)
		s << ")";	
}

/** Calculate maximal score of csf (i,j) for a certain Algebra.
 *  It is assumed that a perfect match of the forest obtains the best
 *  possible score.
 */
template <class L>
template <class R>
R PPForest<L>::maxScore(const Algebra<R,L> &alg, size_type i, size_type j) const
{
	R down, over;

	if(j==0)
		return 0;

	if(isLeave(i))
	{
		over=maxScore(alg,rb(i),j-1);
		return alg.replace(label(i),0,label(i),over);
	}
	else
	{
		down=maxScore(alg,i+1,noc(i));
		over=maxScore(alg,rb(i),j-1);
		return alg.replace(label(i),down,label(i),over);
	}
}  

/** Calculate maximal score of csf (i,j) for a certain RNA_Algebra.
 *  It is assumed that a perfect match of the forest obtains the best
 *  possible score.
 */
template <class L>
template <class R>
R PPForest<L>::maxScore(const RNA_Algebra<R,L> &alg, Uint i, Uint j) const
{
	R down, over;

	if(j==0)
		return 0;

	if(isLeave(i))
	{
		over=maxScore(alg,rb(i),j-1);
		return alg.replace(label(i),0,label(i),over);
	}
	else
	{
		down=maxScore(alg,i+1+1,noc(i)-2);
		over=maxScore(alg,rb(i),j-1);
		return alg.replacepair(label(i+1),label(i+1),down,label(getRightmostBrotherIndex(i+1)),label(getRightmostBrotherIndex(i+1)),over);
	}
}  

/* ****************************************** */
/*            Protected function              */
/* ****************************************** */

template <class L>
void PPForest<L>::initialize(size_type size)
{
	PPForestBase::initialize(size);
	m_lb=new L[size];
}

/* ****************************************** */
/*             Public functions               */
/* ****************************************** */


template<class L>
void PPForest<L>::printDot(ostream &s) const
{
	size_type i,h;

	s << "digraph forest" << endl << "{" << endl;

	// edges
	for(i=0;i<m_size;i++)
	{
		s << i;
		if(m_noc[i])
		{
			s << " -> {";
			h=i+1;
			for(Uint r=0;r<getMaxLength(i+1);r++)
			{
				s << h << " ";
				h=m_rb[h];
			}
			s << "}";	  	     	    
		}
		s << endl;
	}

	// labels
	s << endl << endl;
	for(i=0;i<m_size;i++)
	{
		s << i << "[label=\"";
		showLabel(s,m_lb[i]);
		s << "\"]" << endl;
	}

	s << "}";
}

template<class L>
void PPForest<L>::save(ostream &s) const
{
  /*	
  // save the pforest to stream in binary format
  
  // first of all save the size
  s.write(reinterpret_cast<char*>(&m_size),sizeof(PPForestBase::size_type));
  // save the arrays
  s.write(reinterpret_cast<char*>(m_lb),sizeof(label_type)*m_size);
  s.write(reinterpret_cast<char*>(m_rb),sizeof(PPForestBase::size_type)*m_size);
  s.write(reinterpret_cast<char*>(m_noc),sizeof(PPForestBase::size_type)*m_size);
  s.write(reinterpret_cast<char*>(m_sumUpCSF),sizeof(PPForestBase::size_type)*m_size);
  s.write(reinterpret_cast<char*>(m_rmb),sizeof(PPForestBase::size_type)*m_size);
  */
}
 
template <class L>
ostream& operator<< (ostream &s, PPForest<L> &ppf)
{
  ppf.print(s,0,true);	
  return s;
}

#endif









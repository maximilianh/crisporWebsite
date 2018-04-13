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

#include "rnaforestsz.h"
#include "rnafuncs.h"

#include "ppforestsz.t.cpp"

RNAForestSZ::RNAForestSZ(const string &baseStr, const string &viennaStr, const string &name)
 : PPForestSZ<RNA_Alphabet>(RNAFuncs::treeSize(viennaStr)),
   m_baseStr(baseStr), m_viennaStr(viennaStr)
{  	
  Uint pos,node;
  Ulong basePairCount,maxDepth;	

  assert(RNAFuncs::isRNAString(baseStr));
  RNAFuncs::isViennaString(viennaStr,basePairCount,maxDepth);

  // check if base string and vienna string have the same length
//  if(baseStr.length() > 0)
//    if(baseStr.length() != m_viennaStr.length())
//      throw RNAForestExceptionInput(RNAForestExceptionInput::Error_BaseStringAndViennaStringIncompatible);

  if(name.empty())
    m_name="unknown";
  else
    m_name=name;

  pos=0;
  node=0;
  buildForest(pos,node);
  calcKeyroots();

  // debug 
  /*
  int i;
  cout << "m_lb:" << endl;
  for(i=0;i<m_size;i++)
	  cout << i << ": " << RNA_Alpha2alpha(m_lb[i]) << endl;

  cout << "m_lml:" << endl;
  for(i=0;i<m_size;i++)
	  cout << i << ": " << m_lml[i] << endl; */
}

void RNAForestSZ::buildForest(Uint &pos, Uint &node)
{				
  Uint node2;	

  switch(m_viennaStr[pos])
  {
  case '.':
	  // set label 
	  if(m_baseStr.length())
		m_lb[node]=alpha2RNA_Alpha(m_baseStr[pos]);
	  else
	    m_lb[node]=ALPHA_BASE;

	  m_lml[node]=node;

	  break;
  case '(':   
	  // left paring base
	  if(m_baseStr.length())
		m_lb[node]=alpha2RNA_Alpha(m_baseStr[pos]);
	  else
	    m_lb[node]=ALPHA_BASE;

	  m_lml[node]=node;
	  
	  // build subforest right to (
	  node2=node;
	  node++;
	  pos++;
	  buildForest(pos,node);
	  

	  // right bracket
	  assert(m_viennaStr[pos]==')');

	  // right pairing base
	  if(m_baseStr.length())
		m_lb[node]=alpha2RNA_Alpha(m_baseStr[pos]);
	  else
	    m_lb[node]=ALPHA_BASE;

	  m_lml[node]=node;
	  	

	  node++;
	  m_lb[node]=ALPHA_BASEPAIR;
	  m_lml[node]=node2;
	  break;	  	
  
  case ')':
	  return;
	  break;
  }

  pos++;
  node++;	
  if(pos<m_viennaStr.length())
	buildForest(pos,node);			  
}

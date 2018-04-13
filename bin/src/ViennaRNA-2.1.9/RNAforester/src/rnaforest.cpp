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

#ifndef WIN32
#include "config.h"
#endif

#ifndef HAVE_LIBRNA
#undef HAVE_LIBG2
#endif

#include <string.h>

#include "rnaforest.h"
#include "rnafuncs.h"

#include "ppforest.t.cpp"

/* ****************************************** */
/*    Constructor and Destruktor functions    */
/* ****************************************** */

RNAForest::RNAForest(const string &baseStr, const string &viennaStr, const string &name)
: PPForest<RNA_Alphabet>(RNAFuncs::treeSize(viennaStr)),
  m_baseStr(baseStr),
  m_viennaStr(viennaStr)
{	
  if(name.empty())
    m_name="unknown";
  else
    m_name=name;

  buildForest(baseStr,viennaStr);
};

/* ****************************************** */
/*            Private functions               */
/* ****************************************** */

inline void RNAForest::makeLabel(RNA_Alphabet &a,char c)
{
  a=c;
}

inline void RNAForest::showLabel(ostream &s,RNA_Alphabet a) const
{
  s << a;
}

void RNAForest::buildForest(const string &baseStr, const string &viennaStr)
{
	Ulong basePairCount=0,maxDepth=0,stackPtr;
	Uint baseStrLen,viennaStrLen,node,*nodeStack, *numChildrenStack;

	assert(RNAFuncs::isRNAString(baseStr));

	RNAFuncs::isViennaString(viennaStr,basePairCount,maxDepth);

	baseStrLen=baseStr.length();
	viennaStrLen=viennaStr.length();

	// check if base string and vienna string have the same length
	//  if(baseStr.length() > 0)
	//    if(baseStr.length() != viennaStrLen)
	//      throw RNAForestExceptionInput(RNAForestExceptionInput::Error_BaseStringAndViennaStringIncompatible);

	nodeStack=new Uint[maxDepth+1];
	numChildrenStack=new Uint[maxDepth+1];	 
	memset(nodeStack,0,sizeof(Uint)*maxDepth+1);
	memset(numChildrenStack,0,sizeof(Uint)*maxDepth+1);

	// fill PPForest structure
	stackPtr=0;
	node=0;
	for(Uint i=0;i<viennaStrLen;i++)
	{		
		switch(viennaStr[i])
		{
		case '.':
			// set label
			if(baseStrLen)
				makeLabel(m_lb[node],baseStr[i]);
			else
				makeLabel(m_lb[node],'B');				

			// set right brother
			if(node==size()-1)
				setRightBrotherIndex(node,0);
			else
				setRightBrotherIndex(node,node+1);

			// set num children
			setNumChildren(node,0);

			// increase stack values
			numChildrenStack[stackPtr]++;

			node++;
			break;

		case '(':
			// set label
			makeLabel(m_lb[node],'P');

			// increase stack values
			numChildrenStack[stackPtr]++;

			// push 
			stackPtr++;
			nodeStack[stackPtr]=node;
			numChildrenStack[stackPtr]=1;

			node++;

			// set label
			if(baseStrLen)
				makeLabel(m_lb[node],baseStr[i]);
			else
				makeLabel(m_lb[node],'B');

			// set right brother
			setRightBrotherIndex(node,node+1);

			// set num children
			setNumChildren(node,0);

			node++;				
			break;
		case ')':
			// set label
			if(baseStrLen)
				makeLabel(m_lb[node],baseStr[i]);
			else
				makeLabel(m_lb[node],'B');  

			// set right brother
			setRightBrotherIndex(node,0);    

			// set num children
			setNumChildren(node,0);  

			// pop
			if(node==size()-1)
				setRightBrotherIndex(nodeStack[stackPtr],0);
			else
				setRightBrotherIndex(nodeStack[stackPtr],node+1);  

			setNumChildren(nodeStack[stackPtr],numChildrenStack[stackPtr]+1);
			stackPtr--;

			node++;
			break;							
		}		
	}	

	delete[] nodeStack;
	delete[] numChildrenStack;
	calcSumUpCSF();
	calcRMB();

	assert (m_noc[m_size-1]==0);
}

/* ****************************************** */
/*             Public functions               */
/* ****************************************** */

#ifdef HAVE_LIBG2  // This features require the g2 library
void RNAForest::plot2d(const string &filename_prefix, const list<pair<Uint,Uint> > &regions, const RNAFuncs::SquigglePlotOptions &sqOptions) const
{
  RNAFuncs::drawRNAStructure(m_baseStr,m_viennaStr,filename_prefix,m_name,regions,sqOptions);
}
#endif

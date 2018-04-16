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

#include <algorithm>
#include <string>

#ifndef WIN32
#include "config.h"
#endif

#include "misc.h"
#include "rna_alignment.h"
#include "rnafuncs.h"
#include "utils.h"

/* ****************************************** */
/*            Private functions               */
/* ****************************************** */

void RNA_Alignment::makePairTable(map<Uint,Uint> &pairs, bool first) const
{
        pair<int,int> *baseIndex;   // stores for each node in the alignment the index position of the left and right base (or -1)
	RNA_Alphabet c;

	baseIndex=new pair<int,int>[m_size];

	// initialize pairs, all gaps
	for(int i=size()-1;i>=0;i--)
	  {
	    baseIndex[i].first=-1;
	    baseIndex[i].second=-1;
	  }

	for(int i=size()-1;i>=0;i--)
	{
		// first or second component
		if(first)
			c=m_lb[i].a;
		else
			c=m_lb[i].b;

		if(isLeave(i))
		  {
			if(c!=ALPHA_GAP)
			{
			    baseIndex[i].first=i;
			    baseIndex[i].second=i;
			}
		  }
		else
		{
			// internal node
			// leftmost and rightmost base
			bool lmBaseFound=false;
			for(size_type r=0,h=i+1;r<getMaxLength(i+1);r++,h=rb(h))
			{					  
				// leftmost base
				if(!lmBaseFound && baseIndex[h].first != -1)
				{
					baseIndex[i].first=baseIndex[h].first;
					lmBaseFound=true;
				}

				// rightmost base
				if(baseIndex[h].second != -1)
				{
					baseIndex[i].second=baseIndex[h].second;
				}
			}

			// report pairing bases if P node
			if(c==ALPHA_BASEPAIR)
			{
			        assert(baseIndex[i].first != -1);
			        assert(baseIndex[i].second != -1);
				assert(baseIndex[i].first < baseIndex[i].second);

				pairs[baseIndex[i].first]=baseIndex[i].second;
				pairs[baseIndex[i].second]=baseIndex[i].first;
			}
		}
	}

	delete[] baseIndex;
}

/* ****************************************** */
/*             Public functions               */
/* ****************************************** */

void RNA_Alignment::getSequenceAlignments(string &s1, string &s2) const
{
  s1="";
  s2="";

  // generate base strings
  for(size_type i=0;i<size();i++)
    {
      if(m_lb[i].a != ALPHA_BASEPAIR && m_lb[i].b != ALPHA_BASEPAIR)
	  {
		s1 += m_lb[i].a;
		s2 += m_lb[i].b;
	  }
    }
}

void RNA_Alignment::getStructureAlignment(string &s, bool first) const
{
  s="";

  map<Uint,Uint> pairs;
  makePairTable(pairs, first);

  // iterate through leaves nodes and use information of pairs
  for(size_type i=0;i<size();i++)
    {
	  RNA_Alphabet c;
	  
	  if(first)
		c=m_lb[i].a;
	  else
		c=m_lb[i].b;

      if(isLeave(i))
	  {
		if(c != ALPHA_GAP)
		{
			if(pairs.find(i) != pairs.end())	// is base paired ?
			{
				size_type j=pairs[i];
				if(i<j)
					s+='(';
				else
					s+=')';
			}
			else
				s+='.';
		}
		else
			s+='-';
	  }
    }
}

#ifdef HAVE_LIBG2  // This features require the g2 library
void RNA_Alignment::squigglePlot(const string &filename_suffix, const RNAFuncs::SquigglePlotOptions &options) const
{
  string str1,str2,seq1,seq2;
  string filename;

  getStructureAlignment(str1,true);
  getStructureAlignment(str2,false);
  getSequenceAlignments(seq1,seq2);

  filename = "x_" + filename_suffix;
  RNAFuncs::drawRNAAlignment(str1,str2,seq1,seq2,m_strname1,m_strname2,filename,true,options);
  filename = "y_" + filename_suffix;
  RNAFuncs::drawRNAAlignment(str2,str1,seq1,seq2,m_strname1,m_strname2,filename,false,options);
}
#endif

void RNA_Alignment::setStructureNames(const string &s1,const string &s2)
{
  m_strname1=s1;
  m_strname2=s2;
}

#ifdef HAVE_LIBRNA
void RNA_Alignment::generateXML(ostream &s) const
{
  string str1,str2,seq1,seq2;

  getStructureAlignment(str1,true);
  getStructureAlignment(str2,false);
  getSequenceAlignments(seq1,seq2);

  RNAFuncs::generateRNAAlignmentXML(str1,str2,seq1,seq2,m_strname1,m_strname2,s);
}
#endif



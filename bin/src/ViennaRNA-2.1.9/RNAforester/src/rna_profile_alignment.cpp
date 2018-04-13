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

#include <cmath>

#ifdef HAVE_LIBG2
#include <g2.h>
#include <g2_PS.h>
#endif

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string.h>

#ifdef HAVE_LIBXMLPLUSPLUS
#include <libxml++/libxml++.h>
#endif

#include "matrix.h"
#include "rna_profile_alignment.h"
#include "rnafuncs.h"
#include "utils.h"

#ifdef HAVE_LIBRNA  // This features require the ViennaRNA library
#include "part_func.h"
#include "fold_vars.h"

extern "C"{
#include "fold.h"
}

#endif

#include "ppforest.t.cpp"

ostream& operator <<(ostream &s,RNA_Alphabet_Profile p)
{
    
  s.precision(2);
  s  << p.p[ALPHA_PRO_BASE_A] << "\\n";
  s.precision(2);
  s << p.p[ALPHA_PRO_BASE_C] << "\\n";
  s.precision(2);
  s << p.p[ALPHA_PRO_BASE_G] << "\\n";
  s.precision(2);
  s << p.p[ALPHA_PRO_BASE_U] << "\\n";
  s.precision(2);
  s << p.p[ALPHA_PRO_BASEPAIR] << "\\n";
  s.precision(2);
  s << p.p[ALPHA_PRO_GAP];
    s.precision(2);
    s << p.p[ALPHA_PRO_BASE];

  //  for(Uint i=0;i<p.columnStr.length();i++)
  //    s << p.columnStr[i] << "\\n";

  return s;
}

/* ****************************************** */
/*    Constructor and Destructor functions    */
/* ****************************************** */

RNAProfileAlignment::RNAProfileAlignment(const string &baseStr, const string &viennaStr, const string &name)
: PPForestAli<RNA_Alphabet_Profile,RNA_Alphabet_Profile>(RNAFuncs::treeSize(viennaStr)),
  m_name(name),
  m_numStructures(1)
{
	buildForest(baseStr,viennaStr);
	addStrName(name);
};

#ifdef HAVE_LIBRNA  // This features require the ViennaRNA library
RNAProfileAlignment::RNAProfileAlignment(const string &baseStr, const string &name, const string &constraint, double t)
  : PPForestAli<RNA_Alphabet_Profile,RNA_Alphabet_Profile>(2*baseStr.length()),
    m_name(name),
    m_numStructures(1)
{
  char *viennaStr=NULL;
  
  // calculate partition function for the sequence
  do_backtrack=1;
  init_pf_fold(baseStr.length());

  //if(constraint.length()>0)
  //pf_fold((char*)baseStr.c_str(),(char*)constraint.c_str());    // expicit conversion to non-const value, but pf_fold does not alter baseStr
  //else
  pf_fold((char*)baseStr.c_str(),NULL);    // expicit conversion to non-const value, but pf_fold does not alter baseStr

  viennaStr=new char[baseStr.length()+1];
  dangles=2;
  fold((char*)baseStr.c_str(),viennaStr);

  setSize(RNAFuncs::treeSize(viennaStr));
  buildForest(baseStr,viennaStr,true);
  
  free_pf_arrays();
  delete[] viennaStr;
  
  //  hasSequence=true;
  addStrName(name);
}
#endif

RNAProfileAlignment::RNAProfileAlignment(const string &filename)
{
  // read ppforest from file
  size_type size;
  ifstream s;
  char *colStr;

  s.open(filename.c_str());

  // first of all save the size
  s.read((char*)&size,sizeof(size_type));
  s.read((char*)&m_numStructures,sizeof(Uint));
  
  colStr=new char[m_numStructures+1];
  colStr[m_numStructures]=0;   // terminate string
  
  // initialize structures
  //  ::PPForestBase(m_size);  
  //  m_lb=new L[ppf.size()];
  initialize(size);

  // save the arrays
    for(Uint i=0;i<size;i++)
      {
	label_type l;
	
	for(int r=0;r<RNA_ALPHABET_SIZE;r++)
	  {
	    double d;
	    s.read(reinterpret_cast<char*>(&d),sizeof(double));
	    l.p[r]=d;
	  }

	s.read(colStr,sizeof(char)*m_numStructures);

	l.columnStr=colStr;
	m_lb[i]=l;
    }

  s.read(reinterpret_cast<char*>(m_rb),sizeof(size_type)*size);
  s.read(reinterpret_cast<char*>(m_noc),sizeof(size_type)*size);
  //  s.read((char*)m_sumUpCSF,sizeof(size_type)*size);
  //  s.read((char*)m_rmb,sizeof(size_type)*size);    

  for(Uint r=0;r<m_numStructures;r++)
    {
      char str[21];
      str[20]=0;
      
      s.read(str,sizeof(char)*20);
      m_strNames.push_back(string(str));
    }


  calcSumUpCSF();
  calcRMB();  
}

/*
RNAProfileForest::RNAProfileForest(const Profile_RNA_Alignment_Forest &ppf)
  : RNAForestBase<RNA_Alphabet_Profile>(ppf)
{
  m_numStructures=ppf.m_numStructuresX+ppf.m_numStructuresY;
}*/

/* ****************************************** */
/*            Private functions               */
/* ****************************************** */

void RNAProfileAlignment::makeRepLabel(size_type node, RNA_Alphabet_Profile a,  RNA_Alphabet_Profile b)
    {
	  double p,q;
	  double m;

	  p=m_numStructuresX;	// weight number of structures in profile
	  q=m_numStructuresY;
	  m=p+q;		
	  p/=m;
	  q/=m;

	  // profile
      for(int i=0;i<RNA_ALPHABET_SIZE;i++)
		m_lb[node].p[i]=p*a.p[i]+q*b.p[i];

	  // alignment column
	  m_lb[node].columnStr=a.columnStr + b.columnStr;
    };
  
void RNAProfileAlignment::makeDelLabel(size_type node, RNA_Alphabet_Profile a)
    {
	  double p,q;
	  double m;

	  p=m_numStructuresX;	// weight number of structures in profile
	  q=1;
	  m=p+q;		
	  p/=m;
	  q/=m;

	  // profile
	  m_lb[node]=a;
	  //m_lb[node].p[ALPHA_GAP]=(1+m_lb[node].p[ALPHA_GAP])/2.0;	

	  for(int i=0;i<RNA_ALPHABET_SIZE;i++)
	    {
	      m_lb[node].p[i]*=p;
	    }
	  m_lb[node].p[ALPHA_PRO_GAP]+=q;

	  // alignment column
	  m_lb[node].columnStr=a.columnStr + string(m_numStructuresY,ALPHA_GAP);
    };
  
void RNAProfileAlignment::makeInsLabel(size_type node, RNA_Alphabet_Profile b)
    {
	  double p,q;
	  double m;

	  p=1;						// weight number of structures in profile
	  q=m_numStructuresY;
	  m=p+q;		
	  p/=m;
	  q/=m;

	  // profile
	  m_lb[node]=b;
	  //m_lb[node].p[ALPHA_GAP]=(1+m_lb[node].p[ALPHA_GAP])/2.0;	

	  for(int i=0;i<RNA_ALPHABET_SIZE;i++)
	    {
	      m_lb[node].p[i]*=q;
	    }
	  m_lb[node].p[ALPHA_PRO_GAP]+=p;
	  
	  //	  m_lb[node].p[ALPHA_PRO_GAP]=(p+q*m_lb[node].p[ALPHA_PRO_GAP])/2.0;

	  // alignment column
	  m_lb[node].columnStr=string(m_numStructuresX,ALPHA_GAP) + b.columnStr;
    };
  
void RNAProfileAlignment::showLabel(ostream &s,RNA_Alphabet_Profile p) const
  {
    s << p;
  };


void RNAProfileAlignment::makeLabel(RNA_Alphabet_Profile &p,char c)
{ 
  int i;

  // initialize profile entries to zero
  for(i=0;i<RNA_ALPHABET_SIZE;i++)
    p.p[i]=0.0;

  i=alpha2RNA_Alpha(c);
  p.p[i]=1.0;

  // if it is acgu it is also a base 
  if(i<=ALPHA_PRO_BASE_U)
    p.p[ALPHA_PRO_BASE]=1.0;

  // set column string
  p.columnStr=c;
}

void RNAProfileAlignment::buildForest(const string &baseStr, const string &viennaStr, bool use_bp_prob)
{
	Ulong basePairCount=0,maxDepth=0,stackPtr;
	Uint baseStrLen,viennaStrLen,node,*nodeStack, *numChildrenStack, *baseposStack;

	assert(RNAFuncs::isRNAString(baseStr));

	RNAFuncs::isViennaString(viennaStr,basePairCount,maxDepth);

	baseStrLen=baseStr.length();
	viennaStrLen=viennaStr.length();

	if(baseStrLen)
		hasSequence=true;
	else
		hasSequence=false;	

	// check if base string and vienna string have the same length
	//  if(baseStr.length() > 0)
	//    if(baseStr.length() != viennaStrLen)
	//      throw RNAForestExceptionInput(RNAForestExceptionInput::Error_BaseStringAndViennaStringIncompatible);

	nodeStack=new Uint[maxDepth+1];
	numChildrenStack=new Uint[maxDepth+1];
	baseposStack=new Uint[maxDepth+1];
	memset(nodeStack,0,sizeof(Uint)*maxDepth+1);
	memset(numChildrenStack,0,sizeof(Uint)*maxDepth+1);
	memset(baseposStack,0,sizeof(Uint)*maxDepth+1);

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
			baseposStack[stackPtr]=i+1;

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

#ifdef HAVE_LIBRNA			
			// set basepair probability
			if(use_bp_prob)
			  {			    
			    m_lb[nodeStack[stackPtr]].p[ALPHA_PRO_BASEPAIR]=pr[iindx[baseposStack[stackPtr]]-(i+1)];
			    m_lb[nodeStack[stackPtr]].p[ALPHA_PRO_GAP]=1-m_lb[nodeStack[stackPtr]].p[ALPHA_PRO_BASEPAIR];
			  }			
#endif

			
			node++;
			break;							
		}		
	}	
	
	delete[] nodeStack;
	delete[] numChildrenStack;
	delete[] baseposStack;
	calcSumUpCSF();

	calcRMB();
	
	assert (m_noc[m_size-1]==0);
}

void RNAProfileAlignment::getStructureAlignmentFromCSF(string &s, deque<double> &pairprob, double t,size_type i, size_type j) const
{
  size_type h;
  double bestPairScore;
  //  list<Uint> leftPairList;
  //list<Uint> rightPairList;
  //Uint bestLeftIndex,bestRightIndex;
  //Uint lastLeftIndex,lastRightIndex;
  
  //  QWATCH(i);
  //  QWATCH(j);	

  if(j==0)
    return;
  
  if(isPair(i))
    {
      // TRACE(DBG_GET_PROFILE_STRUCTURE,"Profile_RNA_Alignment_Forest::getStructureAlignmentFromCSF","basepair");
      // WATCH(DBG_GET_PROFILE_STRUCTURE,"Profile_RNA_Alignment_Forest::getStructureAlignmentFromCSF",m_lb[i].p[ALPHA_BASEPAIR]);

      bestPairScore=bestPairs(i);

      // backtrack best pairs
      if(bestPairScore == bestPairs(i+1) + bestPairs(getRightmostBrotherIndex(i+1)) || m_lb[i].p[ALPHA_PRO_BASEPAIR] < t)
      {
	// i pairs not
	//	cout << "unpaired" << endl;
	getStructureAlignmentFromCSF(s,pairprob,t,i+1,noc(i));
	//	cout << "back to:" << endl;
	//	QWATCH(i);
	//	QWATCH(j);
      }
      else
      {
	//	cout << "paired" << endl;
	// i pairs
	s += '(';
	pairprob.push_back(m_lb[i].p[ALPHA_PRO_BASEPAIR]);
	
	// left path - righthand best pairs
	h=i+1;
	while(h < size() && isPair(h))
	{
	  //	  cout << "left" << endl;
	  //	  QWATCH(h);

	  assert((int)noc(h)-1>=0);
	  getStructureAlignmentFromCSF(s,pairprob,t,rb(h+1),noc(h)-1);
	  h=h+1;
	}


	assert((int)noc(i)-2>=0);
	getStructureAlignmentFromCSF(s,pairprob,t,rb(i+1),noc(i)-2);
	//	cout << "back to:" << endl;
	//	QWATCH(i);
	//	QWATCH(j);

	// right path - lefthand best pairs
	h=getRightmostBrotherIndex(i+1);
	while(h < size() && isPair(h))
	{
	  //	  cout << "right" << endl;
	  //	  QWATCH(h);

	  assert((int)noc(h)-1>=0);
	  getStructureAlignmentFromCSF(s,pairprob,t,h+1,noc(h)-1);
	  //h=h+1;
	  h=getRightmostBrotherIndex(h+1);
	}      
 

	s += ')';
	pairprob.push_back(m_lb[i].p[ALPHA_PRO_BASEPAIR]);
      }
    }
  else
    {
      s+= '.';
      pairprob.push_back(m_lb[i].p[ALPHA_PRO_BASE]);
    }
  
  // right forest
  getStructureAlignmentFromCSF(s,pairprob,t,rb(i),j-1);
}

double RNAProfileAlignment::bestPairs(size_type node) const
{
  size_type i=node;
  double d1,d2;

  WATCH(DBG_GET_PROFILE_STRUCTURE,"RNAProfileForest::getStructureAlignment",node);

  if(isBase(node))
    return 0;
  else
    {
      // node pairs not
      d1=bestPairs(node+1) + bestPairs(getRightmostBrotherIndex(node+1));

      // node pairs
      d2=label(node).p[ALPHA_PRO_BASEPAIR];

      // left path - righthand best pairs
      i=node+1;
      while(i < size() && isPair(i))
	{
	  d2+=bestPairs(getRightmostBrotherIndex(i+1));
	  i=i+1;
	}

      // right path - lefthand best pairs
      i=getRightmostBrotherIndex(node+1);
      while(isPair(i) && i < size())
	{
	  d2+=bestPairs(i+1);
	  i=getRightmostBrotherIndex(i+1);
	}

      return max(d1,d2);
    }
} 

#ifdef HAVE_LIBG2  // This features require the g2 library
void RNAProfileAlignment::drawBaseCircles(int device_id,const BaseProbs &bp,double center_x,double center_y) const
{
  const double box_size=6.0;
  const double max_radius=3.0;

  double xpos,ypos;
  int color;
  //  double dashes=0.5;

  // draw the base probabilities as circles on the edges of the square and the gap probability
  // as the center square

  // upper left corner = a
  xpos=center_x-box_size/2.0;
  ypos=center_y+box_size/2.0;
  
  color=g2_ink(device_id,1,0,0);
  g2_pen(device_id,color);
  g2_filled_circle(device_id,xpos,ypos,max_radius*bp.a);

  // upper right corner = c
  xpos=center_x+box_size/2.0;
  ypos=center_y+box_size/2.0;
  
  color=g2_ink(device_id,0,1,0);
  g2_pen(device_id,color);
  g2_filled_circle(device_id,xpos,ypos,max_radius*bp.c);

  // lower right corner = g
  xpos=center_x+box_size/2.0;
  ypos=center_y-box_size/2.0;
  
  color=g2_ink(device_id,0,0,1);
  g2_pen(device_id,color);
  g2_filled_circle(device_id,xpos,ypos,max_radius*bp.g);

  // lower right corner = u
  xpos=center_x-box_size/2.0;
  ypos=center_y-box_size/2.0;
  
  color=g2_ink(device_id,1,0,1);
  g2_pen(device_id,color);
  g2_filled_circle(device_id,xpos,ypos,max_radius*bp.u);

  // gap in the center  
  color=g2_ink(device_id,0,0,0);
  g2_pen(device_id,color);
  g2_filled_circle(device_id,center_x,center_y,max_radius*bp.gap);  

  // draw rectangle for orientation
  //g2_set_dash(device_id,1,&dashes);
  color=g2_ink(device_id,0.75,0.75,0.75);
  g2_pen(device_id,color);
  g2_rectangle(device_id,center_x-box_size/2,center_y-box_size/2,center_x+box_size/2,center_y+box_size/2);
}
#endif

inline double RNAProfileAlignment::getMlBaseFreq(const BaseProbs &bp) const
{
  return max(max(max(bp.a,bp.c),bp.g),bp.u);
}

char RNAProfileAlignment::getMlBase(const BaseProbs &bp) const
{
  double p = getMlBaseFreq(bp);
  if(bp.a==p)
    return ALPHA_BASE_A;
  if(bp.c==p)
    return ALPHA_BASE_C;
  if(bp.g==p)
    return ALPHA_BASE_G;
  if(bp.u==p)
    return ALPHA_BASE_U;

  assert(false); // you should never get so far ...
  return 'x';  
}

void RNAProfileAlignment::getSeqAli(string &seq,Uint row, Uint i, Uint j) const
{
	if(i==0 && j==getMaxLength(0))
		seq="";

	if(j==0)
		return;

	// basepair=internal node
	if(isPair(i))
	{
		getSeqAli(seq,row,i+1,noc(i));
		getSeqAli(seq,row,rb(i),j-1);		
	}
	else
	{
		// base=leaf
		seq+=m_lb[i].columnStr[row];
		getSeqAli(seq,row,rb(i),j-1);
	}
}

void RNAProfileAlignment::getStructAli(string &s,Uint row) const
{
  s="";

  map<Uint,Uint> pairs;
  makePairTable(pairs, row);

  // iterate through leaves nodes and use information of pairs
  for(Uint i=0;i<m_size;i++)
    {
	  RNA_Alphabet c;
	  
	  c=m_lb[i].columnStr[row];

      if(isBase(i))
	  {
		if(c != '-')
		{
			if(pairs.find(i) != pairs.end())	// is base paired ?
			{
				Uint j=pairs[i];
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

void RNAProfileAlignment::makePairTable(map<Uint,Uint> &pairs, Uint row) const
{
	pair<int,int> *baseIndex;		// left and right base
	RNA_Alphabet c;

	assert(row<m_numStructures);
	baseIndex=new pair<int,int>[m_size];

	// initialize pairs, all gaps
	for(int i=size()-1;i>=0;i--)
	  {
	    baseIndex[i].first=-1;
	    baseIndex[i].second=-1;
	  }

	for(int i=size()-1;i>=0;i--)
	{
	        c=m_lb[i].columnStr[row];

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

void RNAProfileAlignment::filterConsensus(string &structure, deque<double> &pairprob, deque<BaseProbs> &baseprobs, double minFreq) const
{
	deque<double>::iterator itPP;
	//	deque<BaseProbs>::iterator itBP;
	string::iterator it;

	for(itPP=pairprob.begin(),it=structure.begin();itPP!=pairprob.end();itPP++,it++)
	{
		if(*itPP<minFreq)
		{
			*itPP=0.0;
			*it='.';
		}
	}

	/*
	for(it=structure.begin(),itPP=pairprob.begin(),itBP=baseprobs.begin();it!=structure.end();)
	{
		if(itBP->base<minFreq)
		{				
			structure.erase(it);
			baseprobs.erase(itBP);
			pairprob.erase(itPP);

			it=structure.begin();
			itPP=pairprob.begin();
			itBP=baseprobs.begin();			
		}
		else
		{
			it++;
			itPP++;
			itBP++;
		}
	}
	*/
}

/* ****************************************** */
/*             Public functions               */
/* ****************************************** */

void RNAProfileAlignment::getSequenceAlignment(deque<BaseProbs> &baseprob) const
{
	BaseProbs bp;

	// generate base strings
	for(size_type i=0;i<size();i++)
	{
		if(isBase(i))
		{	
			bp.a=label(i).p[ALPHA_PRO_BASE_A];
			bp.c=label(i).p[ALPHA_PRO_BASE_C];
			bp.g=label(i).p[ALPHA_PRO_BASE_G];
			bp.u=label(i).p[ALPHA_PRO_BASE_U];
			bp.gap=label(i).p[ALPHA_PRO_GAP];
			bp.base=bp.a+bp.c+bp.g+bp.u;

			baseprob.push_back(bp);
		}
	} 
}

void RNAProfileAlignment::getStructureAlignment(double t,string &s, deque<double> &pairprob) const
{
	WATCH(DBG_GET_PROFILE_STRUCTURE,"Profile_RNA_Alignment_Forest::getStructureAlignment",size());

	s="";
	getStructureAlignmentFromCSF(s,pairprob,t,0,getMaxLength(0));
}

#ifdef HAVE_LIBG2  // This features require the g2 library
void RNAProfileAlignment::squigglePlot(const string &filename, SquigglePlotOptions &options) const
{
	const double base_fontsize=8;
	const Uint num_grey_colors=100;
	const double min_grey_color=1.0;

	string seq,structure;
	string base,structname;
	float *X,*Y,min_X=0,max_X=0,min_Y=0,max_Y=0;
	Uint i;
	short *pair_table;
	int id_PS,id;
	int ps_grey_colors[num_grey_colors];
	int ps_color_red;
	int ps_color_black;
	double xpos,ypos;

	deque<double> pairprob;
	deque<BaseProbs> baseprobs;

	getStructureAlignment(options.minPairProb,structure,pairprob);
	getSequenceAlignment(baseprobs);

	//  filterConsensus(structure,pairprob,baseprobs,0.5);

	//assert(baseprobs.size() == structure.size());
	if(baseprobs.size() != structure.size())
		cerr <<  "Error in resolving consensus structure!" << endl;

	X = new float[structure.size()];
	Y = new float[structure.size()];

	pair_table = make_pair_table(structure.c_str());
	i = naview_xy_coordinates(pair_table, X, Y);
	if(i!=structure.size())
		cerr << "strange things happening in squigglePlot ..." << endl;

	// calculate image dimesions
	for(i=0;i<structure.size();i++)
	{
		min_X=min(min_X,X[i]);
		max_X=max(max_X,X[i]);
		min_Y=min(min_Y,Y[i]);
		max_Y=max(max_Y,Y[i]);
	}

	//  id_PS  = g2_open_PS("ali.ps", g2_A4, g2_PS_port);
	id_PS  = g2_open_EPSF((char*)filename.c_str());
	id     = g2_open_vd();
	g2_attach(id,id_PS);

	//  cout << "min_X: " << min_X <<",max_X: " << max_X << ",min_Y: " << min_Y << "max_Y: " << max_Y << endl; 
	g2_set_coordinate_system(id_PS,595/2.0,842/2.0,0.5,0.5);
	g2_set_line_width(id,0.2);


	// set colors
	double intv=min_grey_color/(double)num_grey_colors;
	for(i=0;i<num_grey_colors;i++)
	{
		double grey_color=min_grey_color-i*intv;
		ps_grey_colors[i]=g2_ink(id_PS,grey_color,grey_color,grey_color);
	}

	ps_color_black=g2_ink(id_PS,0,0,0);
	if(options.greyColors)
		ps_color_red=g2_ink(id_PS,0,0,0);
	else
		ps_color_red=g2_ink(id_PS,1,0,0);

	// draw sequence
	g2_set_font_size(id,base_fontsize);
	for(i=0;i<structure.size();i++)
	{

		if(options.mostLikelySequence)
		{
			double p=getMlBaseFreq(baseprobs[i]);

			//base color
			if(p==1)
				g2_pen(id,ps_color_red);
			else
				g2_pen(id,ps_grey_colors[(int)floor(p*num_grey_colors-1)]);
			  
			base=getMlBase(baseprobs[i]);

			xpos=X[i]-base.length()*base_fontsize/2.0;
			ypos=Y[i]-4;
			g2_string(id,xpos,ypos,(char*)base.c_str());     
		}
		else
		{
			drawBaseCircles(id_PS,baseprobs[i],X[i],Y[i]);
		}

		// connection to next base
		if(i<structure.size()-1)
		{
			if((1-baseprobs[i].gap)*(1-baseprobs[i+1].gap)==1)
				g2_pen(id,ps_color_red);
			else
				g2_pen(id,ps_grey_colors[(int)floor((1-baseprobs[i].gap)*(1-baseprobs[i+1].gap)*num_grey_colors-1)]);

			g2_line(id,X[i],Y[i],X[i+1],Y[i+1]);
		}
	}

	// draw pairings
	// !!! pair_table indexing begins at 1 !!!
	for(i=0;i<structure.size();i++)
	{
		if((unsigned short)pair_table[i+1]>i+1)
		{	    	    
			// pairs in both structures
			if(pairprob[i]==1)
				g2_pen(id,ps_color_red);
			else
				g2_pen(id,ps_grey_colors[(int)floor(pairprob[i]*num_grey_colors-1)]);	   	    

			g2_line(id,X[i],Y[i],X[pair_table[i+1]-1],Y[pair_table[i+1]-1]);
		}
	}

	g2_flush(id);
	g2_close(id);

	free(pair_table);
	DELETE(X);
	DELETE(Y);
}
#endif

void RNAProfileAlignment::printSeqAli() const
{
	Uint i,l;
	deque<string> seqs;
	string seq,info;

	// get alignment rows
	for(i=0;i<m_numStructures;i++)
	{
	  
		getSeqAli(seq,i,0,getMaxLength(0));
		seqs.push_back(seq);
	}

	l=seq.length();

	// sequence 
	//	if(hasSequence)
	//	{
		// calculate info line
		for(i=0;i<l;i++)
		{
				bool equal=true;
			
				for(Uint r=1;r<m_numStructures;r++)
			{
				if((seqs[r])[i]!=(seqs[r-1])[i])
				{
				equal=false;
				break;
				}
			}

			info += equal ? "*" : " ";		  
		}  

		// print it
		// sequences
		for(i=0;i<l;i+=55)
		{
			for(Uint r=0;r<m_numStructures;r++)
				cout << setw(20) << setfill(' ') << left << m_strNames[r].substr(0,20) << setw(5) << " " << RNAFuncs::UpperCase(seqs[r].substr(i,55)) << endl;

			cout << setw(25) << " " << info.substr(i,55) << endl;
			cout << endl;
		}
		//	}
}

deque<pair<string,string> > RNAProfileAlignment::getSeqAli() const
{
	Uint i;
	deque<pair<string,string> > seqAli;
	pair<string,string> p;
	string seq;

	 // get alignment rows^M
         for(i=0;i<m_numStructures;i++)
         {
           getSeqAli(seq,i,0,getMaxLength(0));
	   p.first = m_strNames[i];
	   p.second = seq;
           seqAli.push_back(p);
         }

	 return seqAli; 
}

void RNAProfileAlignment::printStrAli() const
{
        Uint i,l;
	deque<string> strs;
	string str,info;

	// get alignment rows
	for(i=0;i<m_numStructures;i++)
	{
		getStructAli(str,i);

		strs.push_back(str);
	}

	l=str.length();


	for(i=0;i<l;i++)
	  {
	    bool equal=true;
	    
	    for(Uint r=1;r<m_numStructures;r++)
	      {
		if((strs[r])[i]!=(strs[r-1])[i])
		  {
		    equal=false;
		    break;
		  }
              }
	    
	    info += equal ? "*" : " ";		  
	  }  
	
	// structures
	for(i=0;i<l;i+=55)
	{
		for(Uint r=0;r<m_numStructures;r++)
			cout << setw(20) << setfill(' ') << left << m_strNames[r].substr(0,20) << setw(5) << " " << strs[r].substr(i,55) << endl;

		cout << setw(25) << " " <<  info.substr(i,55) << endl;
		cout << endl;
	} 
  
}

deque<string> RNAProfileAlignment::getStrAli() const
{
	Uint i;
	deque<string> strs;
        string str;

        // get alignment rows^M
        for(i=0;i<m_numStructures;i++)
        {
          getStructAli(str,i);
          strs.push_back(str);
        }

	return strs;
}

string RNAProfileAlignment::getConsSeq() const
{
	string seq;
	deque<BaseProbs> baseprobs;
	getSequenceAlignment(baseprobs);

	 // build sequence
	 deque<BaseProbs>::const_iterator it;
	 for(it=baseprobs.begin();it!=baseprobs.end();it++)
	 {
	   seq += getMlBase(*it);
	 }

	 return seq;
}

string RNAProfileAlignment::getConsStr(double minPairProb) const
{
	string str;
	deque<double> pairprob;
	getStructureAlignment(minPairProb,str,pairprob);

	return str;
}

deque<double> RNAProfileAlignment::getBaseProb() const
{
	deque<BaseProbs> baseprobs;
	deque<double> bestBaseprob;

	getSequenceAlignment(baseprobs);
        // build sequence
        deque<BaseProbs>::const_iterator it;
        for(it=baseprobs.begin();it!=baseprobs.end();it++)
        {
           bestBaseprob.push_back(getMlBaseFreq(*it));
        }

	return bestBaseprob;
}

//deque<pair<pair<int,int>,double> > RNAProfileAlignment::getPairProb(double minPairProb)
void RNAProfileAlignment::getPairProb(double &minPairProb, deque<pair<pair<int,int>,double> > &pairprobs)
{
  //deque<pair<pair<int,int>,double> > structureProbs;
	deque<double> pairprob;
	string str;
	
	getStructureAlignment(minPairProb,str,pairprob);

	// store the position of open brackets in an integer deque 
	deque<int> openBrackets;
	// iterator for structureProbs
	pair<int,int> tmpPair;
	pair<pair<int,int>,double> tmpPair2;
	
	// iterate through structure str	
	for(int i=0;i<=str.length();i++)
	{
	  if(str[i]=='(') {
	    openBrackets.push_back(i);
	  } else if(str[i]==')') {
	    // insert the int value of the corresponding opening bracket into structureProbs
	    // insert position of the actual closing bracket
	    // insert probability according to the actual bases
	    tmpPair = pair<int,int> (openBrackets.back(),i);
	    tmpPair2 = pair<pair<int,int>,double> (tmpPair,pairprob[i]);
	    pairprobs.push_back(tmpPair2);
	    // delete position of the last opening bracket
	    openBrackets.pop_back();
	  }
	}
}

void RNAProfileAlignment::printFastaAli(bool noStructure) const
{
	string str,seq;

	// get alignment rows
	for(Uint i=0;i<m_numStructures;i++)
	{
		getStructAli(str,i);	  
		getSeqAli(seq,i,0,getMaxLength(0));

		cout << ">" << m_strNames[i] << endl;
		cout << seq << endl;
		if(!noStructure)
		  cout << str << endl;
	}
}

void RNAProfileAlignment::printConsensus(double minPairProb) const
{
	const int mountain_high=10;
	int j;

	string seq,structure;

	deque<BaseProbs> baseprobs;
	deque<double> pairprob;

	getSequenceAlignment(baseprobs);
	getStructureAlignment(minPairProb,structure,pairprob);

	// build sequence
	deque<BaseProbs>::const_iterator it;
	for(it=baseprobs.begin();it!=baseprobs.end();it++)
	{
		seq += getMlBase(*it);
	}

//	assert(seq.size() == structure.size());
	if(seq.size() != structure.size())
		cerr <<  "Error in resolving consensus structure!" << endl;

	for(Uint i=0;i<seq.length();i+=55)
	{
		// show structure frequency mountain
		for(j=0;j<mountain_high;j++)
		{
			cout << setw(20) << setfill(' ') << " ";
			cout << setw(3) << right << 100 - (100 / mountain_high) * j << "% ";

			for(int k=0;k<55;k++)
			{
				if(i+k<seq.length())
				{
					if(getMlBaseFreq(baseprobs[i+k])>=1-(1/(double)mountain_high)*(double)j)
						cout << "*";
					else
						cout << " ";
				}
			}
			cout << endl;
		}

		cout << setw(25) << " " << RNAFuncs::UpperCase(seq.substr(i,55)) << endl;
		cout << setw(25) << " " << structure.substr(i,55) << endl;

		// show structure frequency mountain
		for(j=0;j<mountain_high;j++)
		{
			cout << setw(20) << setfill(' ') << " ";
			cout << setw(3) << right << (100 / mountain_high) * (j+1) << "% ";

			for(int k=0;k<55;k++)
			{		
				if(i+k<seq.length())
				{
					if(pairprob[i+k]>=(1/(double)mountain_high)*(double)j)
						cout << "*";
					else
						cout << " ";
				}
			}
			cout << endl;
		}

		cout << endl;
	}
}

void RNAProfileAlignment::addStrNames(const deque<string>& strNames)
{
	deque<string>::const_iterator it;
	for(it=strNames.begin();it!=strNames.end();it++)
		m_strNames.push_back(*it);
}

void RNAProfileAlignment::save(const string &filename)
{
  ofstream s(filename.c_str());

  // save the pforest to stream in binary format
  
  // first of all save the size
  s.write(reinterpret_cast<const char*>(&m_size),sizeof(size_type));
  s.write(reinterpret_cast<const char*>(&m_numStructures),sizeof(Uint));
  
  // save the arrays
  for(Uint i=0;i<m_size;i++)
    {
    for(int r=0;r<RNA_ALPHABET_SIZE;r++)
      s.write(reinterpret_cast<const char*>(&m_lb[i].p[r]),sizeof(double));

    s.write(reinterpret_cast<const char*>(m_lb[i].columnStr.c_str()),sizeof(char)*m_numStructures);
    }
  
  
  //  for(Uint i=0;i<m_size;i++)
  //  {
  //    s.write(reinterpret_cast<char*>(&m_lb[i]),sizeof(label_type)*m_size);
  //  }

  //  s.write(reinterpret_cast<char*>(m_lb),sizeof(label_type)*m_size);
  s.write(reinterpret_cast<const char*>(m_rb),sizeof(size_type)*m_size);
  s.write(reinterpret_cast<const char*>(m_noc),sizeof(size_type)*m_size);
  
  for(Uint r=0;r<m_numStructures;r++)
    {
      ostringstream ss;

      ss << setw(20) << setfill(' ') << left <<  m_strNames[r];
      s.write(reinterpret_cast<const char*>(ss.str().c_str()),sizeof(char)*20);
    }
  

  //  s.write(reinterpret_cast<char*>(m_sumUpCSF),sizeof(size_type)*m_size);
  //  s.write(reinterpret_cast<char*>(m_rmb),sizeof(size_type)*m_size);
  
  //  s.write(m_name;
  //Uint m_numStructures;
  //Uint m_numStructuresX;
  //Uint m_numStructuresY;
  //deque<string> m_strNames;
  //bool hasSequence;
}

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

#ifndef _RNA_PROFILE_ALIGNMENT_H
#define _RNA_PROFILE_ALIGNMENT_H

#include <deque>
#include <iostream>
#include <map>

#include "alignment.h"
#include "matrix.h"
#include "ppforestali.h"
#include "rna_algebra.h"
#include "types.h"

//#include "alignment.t.cpp"

using namespace std;

class RNAProfileAlignment
	: public PPForestAli<RNA_Alphabet_Profile,RNA_Alphabet_Profile>
{
 public:
  //  typedef typename PPForestAli<RNA_Alphabet_Profile,RNA_Alphabet_Profile>::size_type size_type;

  struct SquigglePlotOptions
  {
    bool hideBaseNumbers;
    Uint baseNumInterval;
    bool greyColors;
    bool mostLikelySequence;
    double minPairProb;	
  };

  struct BaseProbs
  {
    double a;
    double c;
    double g;
    double u;
    double gap;
    double base;
  };

 private:
  string m_name;
  Uint m_numStructures;
  Uint m_numStructuresX;
  Uint m_numStructuresY;
  deque<string> m_strNames;
  bool hasSequence;

  void makeRepLabel(size_type node, RNA_Alphabet_Profile a,  RNA_Alphabet_Profile b); 
  void makeDelLabel(size_type node, RNA_Alphabet_Profile a); 
  void makeInsLabel(size_type node, RNA_Alphabet_Profile b);
  
  void showLabel(ostream &s,RNA_Alphabet_Profile p) const;
  void makeLabel(RNA_Alphabet_Profile &a,char c);

  /** Build forest from structure, sequence pair. */
  void buildForest(const string &baseStr, const string &viennaStr, bool use_bp_prob=false);

  void getStructureAlignmentFromCSF(string &s, deque<double> &pairprop, double t,size_type i, size_type j) const; 
  double bestPairs(size_type node) const;
  void drawBaseCircles(int device_id,const BaseProbs &bp,double center_x,double center_y) const;
  double getMlBaseFreq(const BaseProbs &baseprobs) const;
  char getMlBase(const BaseProbs &bp) const;
  void getSeqAli(string &seq,Uint row,Uint i,Uint j) const;
  void getStructAli(string &str,Uint row) const;
  void makePairTable(map<Uint,Uint> &pairs, Uint row) const;
  void filterConsensus(string &structure, deque<double> &pairprob, deque<BaseProbs> &baseprobs, double minFreq) const; 
  void addStrName(const string &strName) {m_strNames.push_back(strName);};

  inline bool isBase(size_type node) const
  {
	  if(label(node).p[ALPHA_PRO_BASE] > 0)
		  return true;
	  else
		  return false;
  };

  inline bool isPair(size_type node) const {return !isBase(node);};

public:
	RNAProfileAlignment(const string &baseStr, const string &viennaStr, const string &name);
	RNAProfileAlignment(const string &baseStr, const string &constraint, const string &name, double t);
	RNAProfileAlignment(const string &filename);

	void printSeqAli() const;
	deque<pair<string,string> > getSeqAli() const;
	deque<string> getStrAli() const;
	string getConsSeq() const;
	string getConsStr(double minPairProb) const;
	deque<double> getBaseProb() const;
	//deque<pair<pair<int,int>,double> > getPairProb(double minPairProb);
        void getPairProb(double &minPairProb, deque<pair<pair<int,int>,double> > &pairprobs);
	void printStrAli() const;
	void printConsensus(double minPairProb) const;
	void printFastaAli(bool noStructure=false) const;

	void squigglePlot(const string &filename, SquigglePlotOptions &options) const;
	void getSequenceAlignment(deque<BaseProbs> &baseprob) const;
	void getStructureAlignment(double t, string &s, deque<double> &pairprob) const;
	string getName() const {return m_name;};
	void setName(const string &name) {m_name=name;};
	Uint getNumStructures() const {return m_numStructures;};
	const deque<string>& getStrNames() const {return m_strNames;};
	void addStrNames(const deque<string>& strNames);       
	void save(const string &filename);
	
	RNAProfileAlignment(Uint numStructuresX,Uint numStructuresY) : 
	   PPForestAli<RNA_Alphabet_Profile,RNA_Alphabet_Profile>(),
	   m_numStructures(numStructuresX + numStructuresY),
	   m_numStructuresX(numStructuresX),
	   m_numStructuresY(numStructuresY){};
};

/*
class Profile_RNA_Alignment: public Alignment<double,RNA_Alphabet_Profile,RNA_Alphabet_Profile>
{
 private:
  string m_name;

 public:
  Profile_RNA_Alignment(const RNAProfileForest *ppfx, const RNAProfileForest *ppfy, const DoubleScoreProfileAlgebraType &alg)
      : Alignment<double,RNA_Alphabet_Profile,RNA_Alphabet_Profile>((const PPForest<RNA_Alphabet_Profile>*)ppfx,(const PPForest<RNA_Alphabet_Profile>*)ppfy,alg)
    {
      m_name=ppfx->getName();
      m_name+=".";
      m_name+=ppfy->getName();
    };

  string getName() const {return m_name;};
  };*/

#endif

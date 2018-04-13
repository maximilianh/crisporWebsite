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

#ifndef _RNA_ALIGNMENT_H
#define _RNA_ALIGNMENT_H

#include <iostream>
#include <map>
#include <string>

#include "alignment.h"
#include "ppforestali.h"
#include "rna_algebra.h"
#include "rnaforest.h"
#include "rnafuncs.h"
#include "types.h"

#include "alignment.t.cpp"

using namespace std;

/** Alignment of RNAForests */
class RNA_Alignment : public PPForestAli<RNA_Alphabet,RNA_AlphaPair>
{
private:
	string m_strname1;
	string m_strname2;

	/** @name Implementations of virtual functions */	
	//@{ 

	void makeRepLabel(size_type node, RNA_Alphabet a,  RNA_Alphabet b)
	{
		m_lb[node].a=a;
		m_lb[node].b=b;
	};

	void makeDelLabel(size_type node, RNA_Alphabet a)
	{
		m_lb[node].a=a;
		m_lb[node].b=ALPHA_GAP;
	};

	void makeInsLabel(size_type node, RNA_Alphabet b)
	{
		m_lb[node].a=ALPHA_GAP;
		m_lb[node].b=b;
	};

	void showLabel(ostream &s,RNA_AlphaPair p) const
	{
		s << "[" << p.a << p.b << "]";
	};

	//@}

	void makePairTable(map<Uint,Uint> &pairs, bool first) const;

public:
	void setStructureNames(const string &s1,const string &s2);
	const string& getStructureNameX() const {return m_strname1;};
	const string& getStructureNameY() const {return m_strname2;};
	void getSequenceAlignments(string &s1, string &s2) const;
	void getStructureAlignment(string &s, bool first) const;
	void generateXML(ostream &s) const;
	void squigglePlot(const string &filename_suffix, const RNAFuncs::SquigglePlotOptions &options) const;
};

/*
class PW_RNA_Alignment: public Alignment<int,RNA_Alphabet,RNA_AlphaPair>
{
 public:
  PW_RNA_Alignment(const RNAForest *ppfx, const RNAForest *ppfy, const IntScoreRNA_AlgebraType &alg)
      : Alignment<int,RNA_Alphabet,RNA_AlphaPair>((const PPForest<RNA_Alphabet>*)ppfx,(const PPForest<RNA_Alphabet>*)ppfy,alg) {};
};*/

#endif

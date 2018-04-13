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

#ifndef _RNA_FOREST_H_
#define _RNA_FOREST_H_

#include <iostream>

#include "algebra.h"
#include "ppforest.h"
#include "rna_alphabet.h"
#include "rnafuncs.h"


using namespace std;


/** RNAForest encapsulates the forest representation of an RNA */
class RNAForest
	: public PPForest<RNA_Alphabet>
{
private:
	typedef PPForest<RNA_Alphabet> PPForestRNAType;

	string m_name;
	string m_baseStr;
	string m_viennaStr;

	void showLabel(ostream &s,RNA_Alphabet a) const;	// virtual function	of PPForest
	void makeLabel(RNA_Alphabet &a,char c);

	/** Build	forest from	structure, sequence	pair. */
	void buildForest(const string	&baseStr, const	string &viennaStr);

public:
	/**	Build forest from structure	sequence pair.
	*	This allows	for	direct use of Vienna RNA strings.
	*	e.g.: str="(..(...))", seq="accguuucu"	 
	*/
	RNAForest(const string &baseStr, const string	&viennaStr, const string &name);

	const	string&	getName() const	{return	m_name;};
	const	string&	getBaseStr() const {return m_baseStr;};
	const	string&	getViennaStr() const {return m_viennaStr;};

	void plot2d(const string &filename_prefix, const list<pair<Uint,Uint> > &regions, const RNAFuncs::SquigglePlotOptions &sqOptions) const;	
};

#endif



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

#ifndef _RNA_FORESTSZ_H_
#define _RNA_FORESTSZ_H_

#include <iostream>

#include "ppforestsz.h"
#include "rna_alphabet.h"

using namespace std;


// RNAForestSZ

class RNAForestSZ : public PPForestSZ<RNA_Alphabet>
{
 private:
  string m_name;
  string m_baseStr;
  string m_viennaStr;

//  void showLabel(ostream &s,RNA_Alphabet a);	// virtual function of PPForest
//  void makeLabel(RNA_Alphabet &a,char c);
    void buildForest(Uint &pos, Uint &node);
 public:
  RNAForestSZ(const string &baseStr, const string &viennaStr, const string &name);

  const string& getName() const {return m_name;};
  const string& getBaseStr() const {return m_baseStr;};
  const string& getViennaStr() const {return m_viennaStr;};

//  void plot2d(const string &filename_prefix, const list<pair<Uint,Uint> > &regions, const RNAFuncs::SquigglePlotOptions &sqOptions);
};

#endif



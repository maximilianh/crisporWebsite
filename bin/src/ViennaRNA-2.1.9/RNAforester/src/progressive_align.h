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

#include "graphtypes.h"
#include "matrix.h"
#include "rnaforester_options.h"
#include "rna_profile_alignment.h"
#include "wmatch.h"

/* ****************************************** */
/*          Definitions and typedefs          */
/* ****************************************** */

typedef map<long,RNAProfileAlignment*> RNAProfileAliMapType;
typedef pair<long,RNAProfileAlignment*> RNAProfileAliKeyPairType;

/* ****************************************** */
/*            Function prototypes             */
/* ****************************************** */

void progressiveAlign(deque<RNAProfileAlignment*> &inputList, deque<pair<double,RNAProfileAlignment*> > &resultList, const DoubleScoreProfileAlgebraType *alg,const RNAforesterOptions &options);
Graph makePairsGraph(const RNAProfileAliMapType &inputListProfile, const DoubleScoreProfileAlgebraType *alg, const Matrix<double> *score_mtrx, double threshold);

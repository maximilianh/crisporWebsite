#ifndef _RNA_ALGEBRA_H_
#define _RNA_ALGEBRA_H_

#include <assert.h>
#include <algorithm>
#include <climits>
#include <cstring>

#include "algebra.h"
#include "debug.h"
#include "misc.h"
#include "rna_alphabet.h"
#include "rnaforester_options.h"

/* ****************************************** */
/*          Definitions and typedefs          */
/* ****************************************** */

const double DBL_NEG = -100000000.0;		// the values from limits.h caused problems..
const double DBL_POS = 100000000.0;

typedef Algebra<double,RNA_Alphabet_Profile> DoubleScoreProfileAlgebraType;
typedef Algebra<int,RNA_Alphabet> IntScore_AlgebraType;
typedef RNA_Algebra<int,RNA_Alphabet> IntScoreRNA_AlgebraType;
typedef SZAlgebra<int,RNA_Alphabet> IntScoreSZAlgebraType;

/* ****************************************** */
/*                 Class Score                */
/* ****************************************** */

/** Class Score reads scoring parameters from RNAforester command line. */
class Score
{
 private:
  bool m_isDistance;
  bool m_isLocal;
  bool m_isRIBOSUM;

 public:
  int m_bp_rep_score;
  int m_bp_del_score;
  int m_b_match_score;
  int m_b_rep_score;
  int m_b_del_score;

  Score(RNAforesterOptions &options);
  Score(const Score &s);
  void print();
};

/* ****************************************** */
/*            RNA Algebra Classes             */
/* ****************************************** */

/** Similarity algebra for RNA forests */
class IntSimiRNA_Algebra : public IntScoreRNA_AlgebraType
{
 private:
  Score m_s;
  
 public:
  int empty() const {return 0;};
  int replacepair(RNA_Alphabet la, RNA_Alphabet lb, int down, RNA_Alphabet ra, RNA_Alphabet rb, int over) const
    {
      return m_s.m_bp_rep_score+down+over;
    };

  int replace(RNA_Alphabet a,int down, RNA_Alphabet b, int over) const
  {
    if(a==ALPHA_BASEPAIR && b == ALPHA_BASEPAIR)
      return m_s.m_bp_rep_score+down+over;
    else
      {
	if(a==ALPHA_BASEPAIR || b==ALPHA_BASEPAIR)
	  return INT_MIN/4;
	else
	  {
	    if(a==b)
	      return m_s.m_b_match_score+down+over;
	    else
	      return m_s.m_b_rep_score+down+over;
	  }
      }	     
  };

  int del(RNA_Alphabet a,int down, int over) const
  {
    if(a==ALPHA_BASEPAIR)
      return m_s.m_bp_del_score+down+over;
    else
      return m_s.m_b_del_score+down+over;
  };

  int insert(int down,RNA_Alphabet b,int over) const
  {
    if(b==ALPHA_BASEPAIR)
     return m_s.m_bp_del_score+down+over;
    else
      return m_s.m_b_del_score+down+over;
  };

  int choice(int a, int  b) const
  {
	  return max(a,b);
  };

  int worst_score() const
  {
    return INT_MIN;
  };

  IntSimiRNA_Algebra(const Score &s)
    : m_s(s) {};
};

/** Distance algebra for RNA forests */
class IntDistRNA_Algebra : public IntScoreRNA_AlgebraType
{
 private:
  Score m_s;

 public:
  int empty() const {return 0;};
  int replacepair(RNA_Alphabet la, RNA_Alphabet lb, int down, RNA_Alphabet ra, RNA_Alphabet rb, int over) const
    {
      return m_s.m_bp_rep_score+down+over;
    };

  int replace(RNA_Alphabet a,int down, RNA_Alphabet b, int over) const
  {
    if(a==ALPHA_BASEPAIR && b == ALPHA_BASEPAIR)
      return m_s.m_bp_rep_score+down+over;
    else
      {
	if(a==ALPHA_BASEPAIR || b==ALPHA_BASEPAIR)
	  return INT_MAX/4;
	else
	  {
	    if(a==b)
	      return m_s.m_b_match_score+down+over;
	    else
	      return m_s.m_b_rep_score+down+over;
	  }
      }	     
  };

  int del(RNA_Alphabet a,int down,int over) const
  {
    if(a==ALPHA_BASEPAIR)
      return m_s.m_bp_del_score+down+over;
    else
      return m_s.m_b_del_score+down+over;
  };

  int insert(int down,RNA_Alphabet b,int over) const
  {
    if(b==ALPHA_BASEPAIR)
      return m_s.m_bp_del_score+down+over;
    else
      return m_s.m_b_del_score+down+over;
  };

  int choice(int a, int  b) const
  {
	  return min(a,b);
  };

  int worst_score() const
  {
    return INT_MAX;
  };

  IntDistRNA_Algebra(const Score &s)
    : m_s(s) {};
};

/** RIBOSUM85-60 matrix published in RSEARCH: Finding homologs of single structured RNA sequences
 *  R. Klein and S. Eddy, BMC Bioinformatics 2003 Vol.4
 */
class RIBOSUM8560 : public IntScoreRNA_AlgebraType
{
 private:
  Score m_s;
  int m_baseSubstMtrx[4][4];
  int m_basepairSubstMtrx[4][4][4][4];

 public:
  int empty() const {return 0;};
  int replacepair(RNA_Alphabet la, RNA_Alphabet lb, int down, RNA_Alphabet ra, RNA_Alphabet rb, int over) const
    {
      int i,j,k,l;
      i=alpha2RNA_Alpha(la);
      j=alpha2RNA_Alpha(ra);
      k=alpha2RNA_Alpha(lb);
      l=alpha2RNA_Alpha(rb);

      return m_basepairSubstMtrx[i][j][k][l]+down+over;
    };

  int replace(RNA_Alphabet a,int down, RNA_Alphabet b, int over) const
  {
    assert(!(a==ALPHA_BASEPAIR && b==ALPHA_BASEPAIR));
    
    if(a==ALPHA_BASEPAIR || b==ALPHA_BASEPAIR)
      return INT_MIN/4;
    else
      {
	int i,j;
	i=alpha2RNA_Alpha(a);
	j=alpha2RNA_Alpha(b);
	
	return m_baseSubstMtrx[i][j]+down+over;
      }
  };

  int del(RNA_Alphabet a,int down, int over) const
  {
    if(a==ALPHA_BASEPAIR)
      return m_s.m_bp_del_score+down+over;
    else
      return m_s.m_b_del_score+down+over;
  };

  int insert(int down,RNA_Alphabet b,int over) const
  {
    if(b==ALPHA_BASEPAIR)
     return m_s.m_bp_del_score+down+over;
    else
      return m_s.m_b_del_score+down+over;
  };

  int choice(int a, int  b) const
  {
    return max(a,b);
  };

  int worst_score() const
  {
    return INT_MIN;
  };

  RIBOSUM8560(const Score &s)
    : m_s(s)
    {
      int i,j,k,l;

      // set substitution matrices

      // base replacement
      m_baseSubstMtrx[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_A]=222;
      m_baseSubstMtrx[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_C]=-186;
      m_baseSubstMtrx[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_G]=-146;
      m_baseSubstMtrx[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U]=-139;

      m_baseSubstMtrx[ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_C]=116;
      m_baseSubstMtrx[ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G]=-248;
      m_baseSubstMtrx[ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_U]=-105;

      m_baseSubstMtrx[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_G]=103;
      m_baseSubstMtrx[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U]=-174;

      m_baseSubstMtrx[ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_U]=165;

      // copy triangle
      for(i=0;i<=ALPHA_PRO_BASE_U;i++)
	for(j=0;j<i;j++)
	  m_baseSubstMtrx[i][j]=m_baseSubstMtrx[j][i];
	  
      // basepair replacement

      // set default score. This score should never be used since the scores for canonical basepairs are defined later
      for(i=0;i<=ALPHA_PRO_BASE_U;i++)
	for(j=0;j<=ALPHA_PRO_BASE_U;j++)
	  for(k=i;k<=ALPHA_PRO_BASE_U;k++)
	    for(l=j;l<=ALPHA_PRO_BASE_U;l++)
	      m_basepairSubstMtrx[i][j][k][l]=-1000;

      m_basepairSubstMtrx[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U]=449;
      m_basepairSubstMtrx[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G]=167;
      m_basepairSubstMtrx[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_C]=270;
      m_basepairSubstMtrx[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U]=59;
      m_basepairSubstMtrx[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_A]=161;
      m_basepairSubstMtrx[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G]=-51;      

      m_basepairSubstMtrx[ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G]=536;
      m_basepairSubstMtrx[ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_C]=211;
      m_basepairSubstMtrx[ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U]=-27;
      m_basepairSubstMtrx[ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_A]=275;
      m_basepairSubstMtrx[ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G]=132;

      m_basepairSubstMtrx[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_C]=562;
      m_basepairSubstMtrx[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U]=121;
      m_basepairSubstMtrx[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_A]=167;
      m_basepairSubstMtrx[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G]=-8;
      
      m_basepairSubstMtrx[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U]=347;
      m_basepairSubstMtrx[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_A]=-57;
      m_basepairSubstMtrx[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G]=-209;

      m_basepairSubstMtrx[ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_A]=497;
      m_basepairSubstMtrx[ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G]=114;

      m_basepairSubstMtrx[ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G]=336;

      // copy triangle
      for(i=0;i<=ALPHA_PRO_BASE_U;i++)
	for(j=0;j<=ALPHA_PRO_BASE_U;j++)
	  for(k=0;k<=ALPHA_PRO_BASE_U;k++)
	    for(l=0;l<=ALPHA_PRO_BASE_U;l++)
	      if(k<i || (k==i && l<j))
		m_basepairSubstMtrx[i][j][k][l]=m_basepairSubstMtrx[k][l][i][j];
    };
};

/* ****************************************** */
/*        RNA Profile Algebra Classes         */
/* ****************************************** */

/** Similarity algebra for RNA profile forests */
class DoubleSimiProfileAlgebra : public DoubleScoreProfileAlgebraType
{
 private:
  Score m_s;

 public:
  double empty() const {return 0.0;};
  double replace(RNA_Alphabet_Profile a,double down, RNA_Alphabet_Profile b, double over) const
  {
    if(a.p[ALPHA_PRO_BASEPAIR]>0 && b.p[ALPHA_PRO_BASEPAIR]>0)
      {
	// pair replacement
	return a.p[ALPHA_PRO_BASEPAIR]*b.p[ALPHA_PRO_BASEPAIR]*m_s.m_bp_rep_score +
	       down+over;
      }
    else
      {
	if(a.p[ALPHA_PRO_BASE]>0 && b.p[ALPHA_PRO_BASE]>0)
	  {
		double s=0;

	    // base replacement  
		for(int i=ALPHA_PRO_BASE_A;i<=ALPHA_PRO_BASE_U;i++)
			for(int j=ALPHA_PRO_BASE_A;j<=ALPHA_PRO_BASE_U;j++)
				s+= i==j ? a.p[i]*b.p[j]*m_s.m_b_match_score : a.p[i]*b.p[j]*m_s.m_b_rep_score;

		if(s==0) // no sequence information
			s=a.p[ALPHA_PRO_BASE]*b.p[ALPHA_PRO_BASE]*m_s.m_b_rep_score;

		return s+down+over;
	  }
	else
	  {
	    // undefined operation (replace base by basepair ??)
	    return DBL_NEG/4;
	  }	  
      }	     
  };

  double del(RNA_Alphabet_Profile a,double down, double over) const
  {
    if(a.p[ALPHA_PRO_BASEPAIR]>0)
      return a.p[ALPHA_PRO_BASEPAIR]*m_s.m_bp_del_score+down+over;
    else
      return a.p[ALPHA_PRO_BASE]*m_s.m_b_del_score+down+over;
  };

  double insert(double down,RNA_Alphabet_Profile b,double over) const
  {
    if(b.p[ALPHA_PRO_BASEPAIR]>0)
      return b.p[ALPHA_PRO_BASEPAIR]*m_s.m_bp_del_score+down+over;
    else
      return b.p[ALPHA_PRO_BASE]*m_s.m_b_del_score+down+over;
  };

  double choice(double a, double  b) const
  {
    return max(a,b);
  };

  double worst_score() const
  {
    return DBL_NEG;
  };

  DoubleSimiProfileAlgebra(const Score &s)
    : m_s(s) {};
};

/** Distance algebra for RNA profile forests */
class DoubleDistProfileAlgebra : public DoubleScoreProfileAlgebraType
{
 private:
  Score m_s;

 public:
  double empty() const {return 0.0;};
  double replace(RNA_Alphabet_Profile a,double down, RNA_Alphabet_Profile b, double over) const
  {
    TRACE(DBG_ALGEBRA,"rep","inside!!!");

    if(a.p[ALPHA_PRO_BASEPAIR]>0 && b.p[ALPHA_PRO_BASEPAIR]>0)
      {
	// pair replacement
	return a.p[ALPHA_PRO_BASEPAIR]*b.p[ALPHA_PRO_BASEPAIR]*m_s.m_bp_rep_score +
	       down+over;
      }
    else
      {
	if(a.p[ALPHA_PRO_BASE]>0 && b.p[ALPHA_PRO_BASE]>0)
	  {
		double s=0;

	    // base replacement  
		for(int i=ALPHA_PRO_BASE_A;i<=ALPHA_PRO_BASE_U;i++)
			for(int j=ALPHA_PRO_BASE_A;j<=ALPHA_PRO_BASE_U;j++)
				s+= i==j ? a.p[i]*b.p[j]*m_s.m_b_match_score : a.p[i]*b.p[j]*m_s.m_b_rep_score;

		if(s==0) // no sequence information
			s=a.p[ALPHA_PRO_BASE]*b.p[ALPHA_PRO_BASE]*m_s.m_b_rep_score;

		return s+down+over;
	  }
	else
	  {
	    // undefined operation (replace base by basepair ??)
	    return DBL_POS/4;
	  }	  
      }	     
  };

  double del(RNA_Alphabet_Profile a,double down, double over) const
  {
    if(a.p[ALPHA_PRO_BASEPAIR]>0)
      return a.p[ALPHA_PRO_BASEPAIR]*m_s.m_bp_del_score+down+over;
    else
      return a.p[ALPHA_PRO_BASE]*m_s.m_b_del_score+down+over;
  };

  double insert(double down,RNA_Alphabet_Profile b,double over) const
  {
    if(b.p[ALPHA_PRO_BASEPAIR]>0)
      return b.p[ALPHA_PRO_BASEPAIR]*m_s.m_bp_del_score+down+over;
    else
      return b.p[ALPHA_PRO_BASE]*m_s.m_b_del_score+down+over;
  };

  double choice(double a, double  b) const
  {
	  return min(a,b);
  };

  double worst_score() const
    {
      return DBL_POS;
    };

  DoubleDistProfileAlgebra(const Score &s)
    : m_s(s) {};
};

/* ****************************************** */
/*             SZAlgebra Classes              */
/* ****************************************** */

class IntSimiSZAlgebra : public IntScoreSZAlgebraType
{
 private:
  Score m_s;

 public:
  int empty() const {return 0;};	
  int replace(RNA_Alphabet a,int down, RNA_Alphabet b) const
  {
    if(a==ALPHA_BASEPAIR && b == ALPHA_BASEPAIR)
      return m_s.m_bp_rep_score+down;
    else
      {
	if(a==ALPHA_BASEPAIR || b==ALPHA_BASEPAIR)
	  return INT_MIN/4;
	else
	  {
	    if(a==b)
	      return m_s.m_b_match_score+down;
	    else
	      return m_s.m_b_rep_score+down;
	  }
      }	     
  };

  int del(RNA_Alphabet a,int down) const
  {
    if(a==ALPHA_BASEPAIR)
      return m_s.m_bp_del_score+down;
    else
      return m_s.m_b_del_score+down;
  };

  int insert(int down,RNA_Alphabet b) const
  {
    if(b==ALPHA_BASEPAIR)
      return m_s.m_bp_del_score+down;
    else
      return m_s.m_b_del_score+down;
  };

  int choice(int a, int  b) const
  {
	  return max(a,b);
  };

  int worst_score() const
  {
    return INT_MIN;
  };

  IntSimiSZAlgebra(const Score &s)
    : m_s(s) {};
};


class IntDistSZAlgebra : public IntScoreSZAlgebraType
{
 private:
  Score m_s;

 public:
  int empty() const {return 0;};
  int replace(RNA_Alphabet a,int down, RNA_Alphabet b) const
  {
    if(a==ALPHA_BASEPAIR && b == ALPHA_BASEPAIR)
      return m_s.m_bp_rep_score+down;
    else
      {
	if(a==ALPHA_BASEPAIR || b==ALPHA_BASEPAIR)
	  return INT_MAX/4;
	else
	  {
	    if(a==b)
	      return m_s.m_b_match_score+down;
	    else
	      return m_s.m_b_rep_score+down;
	  }
      }	     
  };

  int del(RNA_Alphabet a,int down) const
  {
    if(a==ALPHA_BASEPAIR)
      return m_s.m_bp_del_score+down;
    else
      return m_s.m_b_del_score+down;
  };

  int insert(int down,RNA_Alphabet b) const
  {
    if(b==ALPHA_BASEPAIR)
      return m_s.m_bp_del_score+down;
    else
      return m_s.m_b_del_score+down;
  };

  int choice(int a, int  b) const
  {
	  return min(a,b);
  };

  int worst_score() const
  {
    return INT_MAX;
  };

  IntDistSZAlgebra(const Score &s)
    : m_s(s) {};
};

/* ****************************************** */
/*           General Algebra Classe           */
/* ****************************************** */

/** Distance algebra for forests */
class IntDist_Algebra : public IntScore_AlgebraType
{
 private:
  Score m_s;

 public:
  int empty() const {return 0;};

  int replace(RNA_Alphabet a,int down, RNA_Alphabet b, int over) const
  {
    if(a==ALPHA_BASEPAIR && b == ALPHA_BASEPAIR)
      return m_s.m_bp_rep_score+down+over;
    else
      {
	if(a==ALPHA_BASEPAIR || b==ALPHA_BASEPAIR)
	  return INT_MAX/4;
	else
	  {
	    if(a==b)
	      return m_s.m_b_match_score+down+over;
	    else
	      return m_s.m_b_rep_score+down+over;
	  }
      }	     
  };

  int del(RNA_Alphabet a,int down,int over) const
  {
    if(a==ALPHA_BASEPAIR)
      return m_s.m_bp_del_score+down+over;
    else
      return m_s.m_b_del_score+down+over;
  };

  int insert(int down,RNA_Alphabet b,int over) const
  {
    if(b==ALPHA_BASEPAIR)
      return m_s.m_bp_del_score+down+over;
    else
      return m_s.m_b_del_score+down+over;
  };

  int choice(int a, int  b) const
  {
	  return min(a,b);
  };

  int worst_score() const
  {
    return INT_MAX;
  };

  IntDist_Algebra(const Score &s)
    : m_s(s) {};
};


#endif

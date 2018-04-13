/*
  Copyright by Matthias Hoechsmann (C) 2002-2003
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

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <assert.h>

template <class T>
class Matrix
{
 private:
  long m_m;
  long m_n;
  T *m_mtrx;
  
 public:
  Matrix(long m, long n) : m_m(m), m_n(n)
    {
      m_mtrx=new T[m*n];
    }
  
  ~Matrix()
    {
      delete m_mtrx;
    }

  inline const long xDim() const
    {
      return m_m;
    }

  inline const long yDim() const
    {
      return m_n;
    }

  inline const T& getAt(long x,long y) const
    {
      assert(x<m_m && y<m_n);
      return m_mtrx[y*m_m+x];
    }

  inline void setAt(long x,long y, const T &val)
    {
      assert(x<m_m && y<m_n);
      m_mtrx[y*m_m+x]=val;
    }
};
 
#endif

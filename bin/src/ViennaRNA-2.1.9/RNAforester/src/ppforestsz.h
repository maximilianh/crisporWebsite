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

#ifndef _PPFORESTSZ_H
#define _PPFORESTSZ_H

#include <cassert>
#include <iostream>

#include "misc.h"
#include "types.h"

using namespace std;

template <class L>
class PPForestSZ
{	
 protected:
  Uint m_size;
  Uint *m_lml;         // postorder index of leftmost leaf descandant of the subtree rootet at T[i]
  L *m_lb;			   // labels of nodes in postorder
  bool *m_keyroot;	   // is node a keyroot

  void calcKeyroots();

 public:
  typedef Uint size_type;
  typedef L label_type;

  PPForestSZ() : m_size(0), m_lml(NULL), m_lb(NULL),m_keyroot(NULL) {};
  PPForestSZ(Uint nrOfNodes);
  PPForestSZ(const PPForestSZ<L> &ppf);
  ~PPForestSZ();

  inline size_type size() const {return m_size;};
  inline Uint lml(Ulong i) const {return m_lml[i];};
  inline L label(Ulong i) const {return m_lb[i];};
  inline bool keyroot(Ulong i) const {return m_keyroot[i];};

  Uint numKeyroots() const
  {
	  Uint count=0;
	  for(Uint i=0;i<m_size;i++)
	  {
		  if(m_keyroot[i])
			  count++;
	  }
	  return count;
  }

  void printParameters(const string &name) const
    {
		cout << name.c_str() << "# size: " << size() << endl;
		cout << name.c_str() << "# keyroots: " << numKeyroots() << endl;
    }
};

#endif


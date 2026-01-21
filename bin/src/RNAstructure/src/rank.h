/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Authors: Chris Connett and Andrew Yohn, 2006
 */

#ifndef RANK_H
#define RANK_H

#include <vector>

#include "workslice.h"
#include "workunit.h"

class Rank {
private:
  std::vector<workunit> workunits;
  std::vector<workunit>::iterator assigned;
  bool allWorkunitsAdded;
  int worksliceSize;
  
public:
  Rank();
  
  void addWorkunit(workunit wu);

  workslice getWorkslice();
};

#endif

/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Authors: Chris Connett and Andrew Yohn, 2006
 */

#include "rank.h"

#include <algorithm>
#include <cassert>
#include <vector>

#include "workslice.h"
#include "workunit.h"

#define MIN_WORKSLICE_SIZE 64UL
#define WORKSLICES_PER_RANK 256UL

Rank::Rank()
  : allWorkunitsAdded(false) {
}

void Rank::addWorkunit(workunit wu) {
  assert(!allWorkunitsAdded);
  workunits.push_back(wu);
}

workslice Rank::getWorkslice() {
  if (!allWorkunitsAdded) {
    assigned = workunits.begin();
    allWorkunitsAdded = true;
    worksliceSize = std::max((long long unsigned int) MIN_WORKSLICE_SIZE,
                             (long long unsigned int) workunits.size() / WORKSLICES_PER_RANK);
  }

  workslice slice;

  slice.begin = assigned;
  slice.end = std::min(assigned + worksliceSize, workunits.end());
  assigned = slice.end;
  
  return slice;
}

/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Authors: Chris Connett and Andrew Yohn, 2006
 */

#ifndef WORKSLICE_H
#define WORKSLICE_H

#include <vector>

#include "workunit.h"

typedef struct {
  std::vector<workunit>::iterator begin;
  std::vector<workunit>::iterator end;
} workslice;

#endif

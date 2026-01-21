/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Authors: Chris Connett and Andrew Yohn, 2006
 */

#ifndef RANKPRODUCER_H
#define RANKPRODUCER_H

#include <pthread.h>

#include "rankmanager.h"

typedef struct {
  RankManager *manager;
  short int *lowend, *highend;

  int n1, n2;
  
  int minLoopSize;
} rankproducerargs;

int getNumInternalRanks(int n1, int n2);
int getNumExternalRanks(int n1, int n2);

void *internalrankproducer(void *arg);
void *externalrankproducer(void *arg);

#endif

/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Authors: Chris Connett and Andrew Yohn, 2006
 */

#ifndef RANKMANAGER_H
#define RANKMANAGER_H

#include <queue>
#include <pthread.h>

#include "observable.h"
#include "rank.h"
#include "TProgressDialog.h"
#include "workslice.h"

class RankManager : public Observable {
private:
  std::queue<Rank*> ranks;
  int numRanks;
  int ranksCompleted;
  int ranksGenerated;
  pthread_mutex_t lock;
  pthread_cond_t rankavailable;
  pthread_cond_t rankcomplete;
  pthread_cond_t producer_toggle;
  
  int checkInCounter;
  int numConsumers;
  
public:
  RankManager(int _numRanks, int _numConsumers);
  ~RankManager();

  void addRank(Rank *r);
  workslice getWorkslice();
};

#endif

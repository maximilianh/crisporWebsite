/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Authors: Chris Connett and Andrew Yohn, 2006
 */

#include "rankmanager.h"

#include <queue>
#include <pthread.h>

#include "observable.h"
#include "rank.h"
#include "TProgressDialog.h"
#include "workslice.h"

#define MIN_RANK_QUEUE_SIZE 5
#define MAX_RANK_QUEUE_SIZE 25

using namespace std;

RankManager::RankManager(int _numRanks, int _numConsumers)
  : numRanks(_numRanks),
    ranksCompleted(0),
    ranksGenerated(0),
    checkInCounter(0),
    numConsumers(_numConsumers) {
  pthread_mutex_init(&lock, NULL);
  pthread_cond_init(&rankavailable, NULL);
  pthread_cond_init(&rankcomplete, NULL);
  pthread_cond_init(&producer_toggle, NULL);
}

RankManager::~RankManager() {
  pthread_mutex_destroy(&lock);
  pthread_cond_destroy(&rankavailable);
  pthread_cond_destroy(&rankcomplete);
  pthread_cond_destroy(&producer_toggle);
}

void RankManager::addRank(Rank *r) {
  pthread_mutex_lock(&lock);

  // Add the rank and let the consumers know there is a new rank
  // available.
  ranks.push(r);
  ranksGenerated++;
  pthread_cond_broadcast(&rankavailable);

  // Check if the rank queue has sufficient extra ranks, and wait here
  // if it does.  A consumer will signal us when they need more.
  // Also, don't wait if all the ranks are generated, because the
  // producer will just terminate.
  if (ranks.size() >= MAX_RANK_QUEUE_SIZE && ranksGenerated < numRanks ) {
    pthread_cond_wait(&producer_toggle, &lock);
  }
  
  pthread_mutex_unlock(&lock);
}

workslice RankManager::getWorkslice() {
  pthread_mutex_lock(&lock);

  workslice slice;
  slice.end = slice.begin; // Ensure that the default values will
                           // compare equal - equal iterators in the
                           // workslice indicate an empty workslice
  bool assignmentmade = false;
  
  if (anyCanceled()) {
    //quit processing ranks because the operation has been canceled.
    while (!ranks.empty()) ranks.pop(); // clear the ranks queue
    ranksCompleted = numRanks;          // mark all ranks completed
  }

  while (!assignmentmade) {
    if (!ranks.empty()) {
      // There is work to do in the current rank OR we are about to
      // move on to a new rank and need to synchronize the threads
      // (IOW, there exists a rank)

      // Try to get work on the current rank
      slice = ranks.front()->getWorkslice();

      if (slice.begin != slice.end) {
        // We were assigned work on the current rank
        assignmentmade = true;
      } else {
        // We are about to move on to a new rank and need to
        // synchronize the threads
        
        checkInCounter++; // One more thread checked in
        if (checkInCounter == numConsumers) {
          // All threads checked in - this thread is the last
          ranksCompleted++; // Rank completed
          delete ranks.front();
          ranks.pop(); // Pop the finished rank
          //printf("rank complete, %d of %d ranks completed, and %d in the queue\n",
          //       ranksCompleted, numRanks, ranks.size());
          notifyObservers();
          
          // Reset the check in counter for the next rank
          checkInCounter = 0;

          // Notify all the other waiting threads
          pthread_cond_broadcast(&rankcomplete);

          // Check if ranks remaining in the queue are fewer than our
          // acceptable buffer size (and that there are still ranks
          // left to be generated)
          if (ranks.size() < MIN_RANK_QUEUE_SIZE && ranksGenerated < numRanks) {
            // We (the one last thread) wake the producer to produce
            // more ranks.  It will stop automatically when it has
            // generated a sufficient number of ranks.

            pthread_cond_signal(&producer_toggle);
          }
        } else {
          // Wait for other threads to check in
          pthread_cond_wait(&rankcomplete, &lock);
        }
        // Assign no work here, all will get work the next time around
        // the while loop
      }

    } else if (ranksCompleted == numRanks) {
      // Finished all ranks for which this RankManager is responsible
      
      assignmentmade = true;
      // Return the default empty workslice
      
    } else { //if (ranks.empty() && ranksGenerated < numRanks)
      // No work on current rank, but there are more ranks to be
      // generated, so wait for more work to be generated - will be
      // notified by rank producer when more work has been produced
      pthread_cond_wait(&rankavailable, &lock);
    }
  }
  
  pthread_mutex_unlock(&lock);
  
  return slice;
}

// void RankManager::updateProgressBar() {
//   RankManager::totalRanksCompleted++;
  
//   if (RankManager::progressBar != NULL) {
//     RankManager::progressBar->update((100 * RankManager::totalRanksCompleted) /
//                                      RankManager::totalNumRanks);
//   }
// }

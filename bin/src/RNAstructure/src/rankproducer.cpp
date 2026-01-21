/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Authors: Chris Connett and Andrew Yohn, 2006
 */

#include "rankproducer.h"

#include <algorithm>
#include <cmath>

#include "rank.h"
#include "rankmanager.h"

using namespace std;

int getNumInternalRanks(int n1, int n2) {
  return n1 + n2 - 5;
}

int getNumExternalRanks(int n1, int n2) {
  return n1 + 2 * n2 - 4;
}

void *internalrankproducer(void *arg) {
  rankproducerargs *rpa = static_cast<rankproducerargs*>(arg);
  RankManager *manager = rpa->manager;
  int n1 = rpa->n1;
  int n2 = rpa->n2;
  int minLoopSize = rpa->minLoopSize;

  short int *lowend = rpa->lowend;
  short int *highend = rpa->highend;
  
  int numRanks = getNumInternalRanks(n1, n2);

  int r;

  short i, j, k, l;
  short imax, jmax, kmax;

  for (r = 0; r < numRanks; r++) {
    Rank *newRank = new Rank();
    
    imax = std::min(numRanks - r, n1 - 2);
    for (i = 1; i <= imax; i++) {

      jmax = std::min(r + i + 2, n1);
      for (j = std::max(i + 2, i + 2 + r  - (n2 - 2)); j <= jmax; j++) {

        kmax = n2 - 2 - r + j - i - 2;
        for (k = 1; k <= kmax; k++) {

          l = 3 + k - 1 + r + i - j + 2;
          if (j - i > minLoopSize &&
              l - k > minLoopSize &&
              lowend[i] <= k && k <= highend[i] &&
              lowend[j] <= l && l <= highend[j]) {

            workunit wu = {i, j, k, l};
  
            newRank->addWorkunit(wu);
          }
        }
      }
    }

    manager->addRank(newRank);
  }

  return NULL;
}

void *externalrankproducer(void *arg) {
  rankproducerargs *rpa = static_cast<rankproducerargs*>(arg);
  RankManager *manager = rpa->manager;
  int n1 = rpa->n1;
  int n2 = rpa->n2;

  short int *lowend = rpa->lowend;
  short int *highend = rpa->highend;
  
  int numRanks = getNumExternalRanks(n1, n2);
  int r, majorD, minorD;

  short i, j, k, l;

  int numMajorDiagonals = n1 - 1;

  for (r = 0; r < numRanks; r++) {
    Rank *newRank = new Rank();

    for (majorD = 0; majorD < numMajorDiagonals; majorD++) {
      minorD = r - majorD;

      if (0 <= minorD && minorD <= 2 * n2 - 2) {
        i = n1 - majorD;
        j = n1 + 1;

        while (i <= n1) {
          k = max(n2 - minorD, 1);
          l = max(2 + minorD, n2 + 1);

          while (k <= n2 && l <= 2 * n2 - 1) {
            if (lowend[i] <= k && k <= highend[i] &&
                lowend[j] <= l && l <= highend[j]) {
              
              workunit wu = {i, j, k, l};
              
              //fprintf(stderr, "('v', %d, %d, %d, %d),\n", i, j, k, l);
              newRank->addWorkunit(wu);
            }

            k++;
            l++;
          }

          i++;
          j++;
        }
      }
    }

    manager->addRank(newRank);
  }

  return NULL;
}

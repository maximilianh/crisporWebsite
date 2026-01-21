/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Authors: Chris Connett and Andrew Yohn, 2006
 */

#ifndef OBSERVABLE_H
#define OBSERVABLE_H

#include <list>

#include "observer.h"

class Observable {
private:
  std::list<Observer *> observers;

public:
  Observable();
  virtual ~Observable() = 0;

  void subscribe(Observer *o);
  void unsubscribe(Observer *o);

protected:
  void notifyObservers();
  bool anyCanceled();
};

#endif

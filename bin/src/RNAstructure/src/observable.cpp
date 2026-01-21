/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Authors: Chris Connett and Andrew Yohn, 2006
 */

#include "observable.h"

#include <algorithm>
#include <list>

#include "observer.h"

using namespace std;

Observable::Observable() {
}

Observable::~Observable() {
}

void Observable::subscribe(Observer *o) {
  if(find(observers.begin(), observers.end(), o) == observers.end()) {
    observers.push_back(o);
  }
}

void Observable::unsubscribe(Observer *o) {
  list<Observer *>::iterator i =
    find(observers.begin(), observers.end(), o);
  if(i != observers.end()) {
    observers.erase(i);
  }
}

void Observable::notifyObservers() {
  list<Observer *>::iterator i = observers.begin();
  list<Observer *>::iterator end = observers.end();
  while(i != end) {
    (*(i++))->notify();
  }
}

bool Observable::anyCanceled() {
  list<Observer *>::iterator i = observers.begin();
  list<Observer *>::iterator end = observers.end();
  while(i != end) {
    if ((*(i++))->canceled()) return true;
  }
  return false;
}
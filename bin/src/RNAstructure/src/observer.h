/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Authors: Chris Connett and Andrew Yohn, 2006
 */

#ifndef OBSERVER_H
#define OBSERVER_H

class Observer {
public:
  Observer();
  virtual ~Observer() = 0;
  virtual void notify() = 0;
  virtual bool canceled() { return false; }
};

#endif

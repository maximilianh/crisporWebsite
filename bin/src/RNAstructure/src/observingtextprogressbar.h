/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Authors: Chris Connett and Andrew Yohn, 2006
 * Modified by Richard Watson, 2016
 */

#ifndef OBSERVINGTEXTPROGRESSBAR_H
#define OBSERVINGTEXTPROGRESSBAR_H

#include "TProgressDialog.h"
#include "observer.h"

class ObservingTextProgressBar : public Observer {
private:
  int value;
  int max;
  ProgressHandler *progress;
  
public:
    // 2016-06-11 RMW - instead of inheriting from TProgressDialog, ObservingTextProgressBar now delegates to a TProgressDialog.
	  // This allows the inner progress dialog to be specified by the caller, so it can be different for e.g. Text vs GUI progressbars.
	/**
	 *  pass NULL for innerProgress to create a Silent ObservingTextProgressBar
	 */
	ObservingTextProgressBar(ProgressHandler *innerProgress, int _max = 100);
  void setMaximum(int _max);
  void notify();
  bool canceled();
};

#endif

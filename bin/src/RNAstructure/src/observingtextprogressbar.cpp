/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Authors: Chris Connett and Andrew Yohn, 2006
 */

#include "observingtextprogressbar.h"
#include "observer.h"
#include "TProgressDialog.h"

ObservingTextProgressBar::ObservingTextProgressBar(ProgressHandler *innerProgress, int _max)
  : progress(innerProgress),
    value(0),
    max(_max)
{
}

void ObservingTextProgressBar::setMaximum(int _max) {
  max = _max;
}

void ObservingTextProgressBar::notify() {
	if (progress != 0) progress->update((100 * ++value) / max);
}

bool ObservingTextProgressBar::canceled() {
	return progress&&progress->canceled();
}
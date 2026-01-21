/*
 * RNAstructure
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Authors: 
 *   Chris Connett and Andrew Yohn (2006) -- authors of text version.
 *   Jessica Reuter (2009) -- SWIG version
 *   Richard Watson (2016) -- Combined Text and GUI versions and added cancellation.
 */

#include "TProgressDialog.h"

using namespace std;

const char TProgressDialog::spinchars[] = {'/', '-', '\\', '|'};

ProgressHandler::ProgressHandler() : percentProgress(0), isCanceled(false) { }

// Was the operaton cancelled?
bool ProgressHandler::canceled() const {
  return isCanceled;
}

// Get the current progress percentage
void ProgressHandler::cancel() {
  isCanceled = true;
}

// Get the current progress percentage
int  ProgressHandler::progress() const {
  return percentProgress;
}

// Reset the progress and isCanceled status
void  ProgressHandler::reset() {
  percentProgress = 0; 
  isCanceled = false;
}

// Update the value of the monitor.
void ProgressHandler::update(int percent) { percentProgress=percent; }

//---------- TProgressDialog --------------//

// Constructor
TProgressDialog::TProgressDialog(ostream * const outputStream) : output(outputStream) {
  spinstate=0; 
}

// Update the value of the monitor.
void TProgressDialog::update(int percent) {
  percentProgress = percent;
  //cout << endl << "TProgressDialog Updated: " << percent << endl;
  // output to command-line
  if (output != NULL) {
    std::ostream &out = *output;
    out << "\r";
    out.width(3);
    out << percent << "% [";
    
    for (int i = 0; i < 100; i += 2) {
      if (i <= percent) {
        out << "=";
      } else {
        out << " ";
      }
    }

    out << "] ";
    if (percent < 100) {
      out << spinchars[spinstate] << "                     ";
      // Exactly 21 spaces puts cursor at last character of an
      // 80-character wide terminal
    } else {
      out << " \n";
    }
    out.flush();
    ++spinstate %= 4;
  }
}

//----------------- PartialProgress --------------------------------//

PartialProgress::PartialProgress(ProgressHandler* clientProgress) : client(clientProgress) { workComplete=0; workInNextStep=100; }

void PartialProgress::update( int percent ) { 
    //cout << endl << "PartialProgress Updated: " << percent << " --> " << (workComplete+workInNextStep/100.0f*percent) << endl;
    if (client!=NULL) client->update(workComplete+workInNextStep/100.0f*percent);
    ProgressHandler::update(percent); // sets progressPercent
}

bool PartialProgress::canceled() const { return client==NULL?ProgressHandler::canceled():client->canceled(); }

void PartialProgress::cancel() { if (client==NULL) ProgressHandler::cancel(); else client->cancel(); } //just forward to client

void PartialProgress::setClient(ProgressHandler* newClient)  { client = newClient; } // set a new client

void PartialProgress::setWorkComplete(const int percentWorkCompleted) { 
    workComplete=percentWorkCompleted;
    if (client!=NULL) client->update(percentWorkCompleted); // update the display showing the current amount of work complete.
}

void PartialProgress::setNextStep(const int percentWorkInNextStep) { 
  workInNextStep=percentWorkInNextStep;
}

void PartialProgress::stepComplete() {
  setWorkComplete(workComplete+workInNextStep); // add the work that was specified in the previous step. (or 0 if this is the first step).
  workInNextStep=0;
}

void PartialProgress::advanceToNextStep(const int percentWorkInNextStep) { 
  stepComplete(); // add the work that was specified in the previous step and update the display.
  setNextStep(percentWorkInNextStep);
}
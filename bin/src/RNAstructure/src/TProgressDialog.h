/*
 * RNAstructure
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Authors: 
 *	 Chris Connett and Andrew Yohn (2006) -- authors of text version.
 *   Jessica Reuter (2009) -- SWIG version
 *   Richard Watson (2016) -- Combined Text and GUI versions and added cancellation.
 *   Richard Watson (2017) -- Added ProgressHandler, SimpleProgressHandler, and PartialProgress
 */

#ifndef TEXTPROGRESSBAR_H
#define TEXTPROGRESSBAR_H

#include <iostream>

#ifdef _WINDOWS_GUI // if referenced by the MFC Windows C++ GUI, use a different TProgressDialog
  #include "../Windows_interface_deprecated/TProgressDialog.h"
#else


//!  ProgressHandler is a base class for implementing objects 
//!  that handle displaying progress updates to the user.
//!  By itself ProgressHandler does nothing but store the current progress percent
//!  as well as a boolean to indicate whether or not the operation has been canceled.
//!  Derived classes (such as TProgressDialog) can override the update() method to 
//!  display progress updates to the user.
//!  This class is abstract -- it cannot be created directly (due to the protected constructor),
//!  to prevent developers from accidentally creating a "non-functional" ProgressHandler by mistake 
//!  for a Text-interface program.
//!  Text-interface programs should use the TProgressDialog subclass.
//!  Library clients (e.g. Java or python, javascript etc) should use 
//!  the SimpleProgressHandler class.
//!  Both text-interface and library clients can also use the PartialProgress class which
//!  acts as a wrapper over an inner "client" ProgressHandler and providers a higher level 
//!  of control over the progress is output to the user.
class ProgressHandler /* abstract */  { 
protected:
  // cannot create this class directly. Instead create a TProgressDialog, SimpleProgressHandler, or PartialProgress.
  ProgressHandler(); 
public:
  virtual ~ProgressHandler(){}

 /*
   * Name:        update
   * Description: Updates the progress percent
   * Argument:
   *     1.   The new progress percent value.
   */
  virtual void update( int percent );

  // Update the value of the monitor with the amount of work complete along with the total amount of work required.
  // This is a convenience method that forwards the calculated percent complete to `update`.
  inline void update(double workComplete, double totalWork) { updateFraction(workComplete/totalWork); }

  // Update the value of the monitor with the fraction (not percent) of work complete.
  // fractionComplete should be number from 0 to 1.
  // This is a convenience method that forwards the calculated percent complete to `update`.
  inline void updateFraction(double fractionComplete) { update(int(100*fractionComplete)); }

  /*
   * Name:        canceled
   * Description: returns true if the operation was canceled by the user.
   */
  virtual bool canceled() const;

    /*
   * Name:        cancel
   * Description: Signal to the current calculation that it should gracefully abort, cleanup, and exit.
   */
  virtual void cancel();

  /*
   * Name:        progress
   * Description: returns the current progress percent.
   */
  virtual int progress() const;

    /*
   * Name:        reset
   * Description: Resets the progress and canceled status.
   */
  virtual void reset();
 
 protected:
  // percent of total progress completed
  int percentProgress;
  bool isCanceled;
};

//! Basic ProgressHandler with no actual output abilities on its own.
//! All it does is keep track of the progress value that is set by
//! the working code calling update(int percent);
//! It can also be used to cancel the working code (if the latter supports it.)
class SimpleProgressHandler : public ProgressHandler { };

//!  TProgressDialog displays progress updates in Text-interface programs by printing to stdout 
//!  (or another std::ostream if specified in the constructor.)
class TProgressDialog : public ProgressHandler {
 public:
  // Public constructor and method

  /*
   * Name:        Constructor
   * Description: Constructs the TProgressDialog proxy object.
   * Argument:
   *          1. A reference to a ProgressMonitor
   */
  TProgressDialog(std::ostream *outputStream = &std::cout);

  /*
   * Name:        update
   * Description: Updates the progress percent
   * Argument:
   *     1.   The new progress percent value.
   */
  void update( int percent );

//  /*
//   * Name:        canceled
//   * Description: returns true if the operation was canceled by the user.
//   */
//  virtual bool canceled() const;
//
//    /*
//   * Name:        cancel
//   * Description: Signal to the current calculation that it should gracefully abort, cleanup, and exit.
//   */
//  virtual void cancel();
//
//  /*
//   * Name:        progress
//   * Description: returns the current progress percent.
//   */
//  virtual int progress() const;
//
//    /*
//   * Name:        reset
//   * Description: Resets the progress and canceled status.
//   */
//  virtual void reset();

 private:
  // Private variable: percent of total progress completed
  //int percentProgress;
  //bool isCanceled;
  static const char spinchars[];
  int spinstate;
  std::ostream * const output;
};


//!  PartialProgress allows better control over the progress that is displayed to the user by acting as
//!  a proxy or wrapper for an "client" progress object.
//! 
//!  Many functions (e.g. RNA::PartitionFunction) update a progress object from 0% to 100% in 
//!  the assumption that they represent the full amount of work done by a program.
//!  But programs that make multiple calls to such functions need to be able to display 
//!  accurate progress to the user.   
//!
//!  For example Rsample calls PartitionFunction, then Stochastic, then PartitionFunction again.
//!  Ordinarilly this would mean the progress shown to the user would go from 0% to 100% and then
//!  jump back to 0% to repeat two more times.  The PartialProgress class allows the client code 
//!  (e.g. in the Rsample function) to make the progress bar move from 0% to 40% for the first call to 
//!  PartitionFunction. Then from 40% to 60% for Stochastic, and finally from 60% to 100% for the  
//!  final PartitionFunction.
//!
//!  It does this by allowing the client code to set "workComplete" (the amount of work already completed) 
//!  and "workInNextStep" (the percentage of work that is to be completed in the upcoming/ongoing step).
//!  The "true" client progress (shown to the user) is calculated as:
//!  clientProgress =  workComplete + workInNextStep/100 * ongoingProgress 
//!  (where ongoingProgress is set by the inner function call --e.g. PartitionFunction when calling 
//!  ProgressHandler::update)
class PartialProgress : public ProgressHandler {
  public:
    PartialProgress(ProgressHandler* clientProgress);
    
    // client progress is calculcated as offset+scale*progressPercent
    int workComplete;   // represents percent of work already completed
    int workInNextStep; // represents percent of work to be completed by the upcoming/ongoing step.

  /*
   * Name:        update
   * Description: Updates the progress percent
   * Argument:
   *     1.   The new progress percent value.
   */
   void update( int percent );

  /*
   * Name:        canceled
   * Description: returns true if the operation was canceled by the user.
   */
  bool canceled() const;

    /*
   * Name:        cancel
   * Description: Signal to the current calculation that it should gracefully abort, cleanup, and exit.
   */
   void cancel();

   //! Call this to set the percentage of work that is to be completed in the next step of the program.
   void setNextStep(const int percentWorkInNextStep);

   //! Explicitly set the total (true) amount of work already completed.
   //! (This is called implicitly by advanceStep and stepComplete.)
   //! \param percentWorkCompleted The total true amount of work already completed.
   void setWorkComplete(const int percentWorkCompleted);

   //! Add all work completed in the previous step to the total 
   //! work complete and then update the client progress.
   //! (This is called implicitly by advanceStep.)
   void stepComplete();

   //! This calls stepComplete() to add all work completed in the previous step and then
   //! calls setNextStep to set the percentage of work that is to be completed in the next step of the program.
   //! \param percentWorkInNextStep The percentage of work that is to be completed in the next step of the program.
   void advanceToNextStep(const int percentWorkInNextStep);

   //! Set a new inner-ProgressHandler.
   void setClient(ProgressHandler* client);

   private:
     ProgressHandler* client; // the ProgressHandler passed in by client code.
};

#endif // ifdef _WINDOWS_GUI (else)
#endif // TEXTPROGRESSBAR_H

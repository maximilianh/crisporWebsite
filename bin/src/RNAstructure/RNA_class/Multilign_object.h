#ifndef MULTILIGN_H
#define MULTILIGN_H

#include <utility> //needed for library type pair
#include <vector>
#include <string>
#include "Dynalign_object.h"
#include "../src/TProgressDialog.h"

using namespace std;

//! Multilign_object Class.
/*!
    The Multilign_object class provides an entry point for the Multilign algorithm.
*/

class Multilign_object {

private:
  typedef vector<string>::size_type vs_index;
  typedef vector<vector<string> >::iterator vvs_it;
  typedef vector<vector<string> >::const_iterator vvs_cit;

public:
  //!Default constructor:
  Multilign_object();

  //!Constructor:
  //!\param inputlist is a vector of vectors of strings  storing the name of the filenames. Currently, inputList is an matrix of 4 columns: col 1 is the input seq filename; col 2 is the output ct filename; col 3 is the input constraint filename; col 4 is the input SHAPE filename. Empty string is not allowed for Col 1 and Col 2; If no SHAPE or folding constraints are given, Col 3 and Col 4 for the correponding sequence are empty strings.
  //!\param isrna is a bool indicating the sequences are RNA or DNA.  The default of true indicates RNA.
  //!\param progress is a TProgressDialog for reporting progress of the calculation to the user.  The default value of NULL means that no communication is provided.
  Multilign_object(const vector<vector<string> > &inputlist, const bool isrna=true, ProgressHandler *progress = NULL);

  Multilign_object(const bool Multifind ,const string &outputmultifind, const vector<string> &ctfiles, ProgressHandler *progress = NULL, const bool isrna=true);

  //Destructor
  ~Multilign_object();



  /// @brief    count the number of basepairs with the lowest free energies below the percent of the miminal free energy. Dsv file is used for the counting. By default, the first dsv file in the progressive dynalign calculations is used, i.e. i = 0 and j = 0.
  ///
  /// @param    i is an int value indicating which one in the iteration.
  /// @param    j is an int value indicating which iteration.
  /// @param    percent is threshold of double value in percentage.
  /// @return   an int of the number of basepairs counted.
  int CountBP(const int i=0, const int j=0, const double percent = 0.8) const;


  //! The core function doing dynalign calculation and templating
  //! In case of error, the function returns a non-zero that can be parsed by GetErrorMessage() or GetErrorMessageString().
  //! \param Dsv is a boolean value indicating to output pairwise dsv files or not. It has to be set to true currently.
  //! \param Ali is a boolean value indicating to output pairwise ali files or not. It has to be set to true currently.
  //! \param maxtrace is the maximum number of common structures to be determined.
  //! \param bpwin the the base pair window parameter, where 0 allows the structures to have similar pairs and larger windows make the structures more diverse.
  //! \param awin is the alignment window parameter, where 0 allows the alignments to be similar and larger values make the alignments more diverse.
  //! \param percent is the maximum percent difference in total folding free energy change above the lowest for suboptimal common structures.
  //! \param imaxseparation is the maximum separation between aligned nucleotides.  Values >= 0 are the traditional parameter, those below zero trigger the HMM alignment method, which is now prefered.
  //! \param gap is the cost of adding gap nucleotides in the alignment in kcal/mol.
  //! \param singleinsert is whether single basepair inserts are allowed in one sequence vs the other.
  //! \param singlefold_subopt_percent is the maximum % difference of folding energy above the lowest free energy structure for pairs in single sequence folding that will be allowed in the dynalign calculation.
  //! \param local is whether Dynalign is being run in local (true) or global mode (false).
  //! \param numProcessors is the number of processors to use for the calculation.  This requires a compilation for SMP.
  //! \return an int that indicates an error code (0 = no error, non-zero = error occurred).
  int ProgressiveMultilign(
          const short int numProcessors=1,
          const bool Dsv=1, const bool Ali=1,
          const short int maxtrace=750,
          const short int bpwin=2, const short int awin=1,
          const short int percent=20,
          const short int imaxseparation=-99,
          const float gap=0.4,
#ifdef DYNALIGN_II
		const float slope = 0.1, const float intercept = 0.5,
		const int max_elongation=5,
#else
		const bool singleinsert=true,
#endif
          const short int singlefold_subopt_percent=30,
          const bool local=false);


  int MultiTempMultilign();


  /// @brief    calculate and output multiple alignment
  ///
  /// @param    allali is the output filename of multiple alignment
  ///
  /// @return   an int of error code.
  int WriteAlignment(const string allali = "all.ali") const;


  //******************************************************
  //Functions to return error information
  //******************************************************

  //!Return an error code, where a return of zero is no error.

  //!  This function returns and error flag that is generated during construction by RNA(const char &filename, const int type, const bool IsRNA=true) or from CalculateFreeEnergy().
  //!    An error of zero is always no error.  Other codes are errors and a c-string can be fetched for the error with GetErrorMessage().
  //!\return An integer that provides the error code.
  int GetErrorCode() const {return ErrorCode;}


  //!  Return error messages based on code from GetErrorCode and other error codes.
  //!    0 = no error
  //!    1000 = Error associated with sequence 1 or with a procedure, function will get message from sequence 1 (the inherited RNA class).
  //!    2000 = Error associated with sequence 2, function will get message from sequence 2 (the RNA2 class).
  //!    3000 = Errors with each sequence, function will get messages from each.
  //!\param error is the integer error code provided by GetErrorCode().
  //!\return A string that provides an error message or from other functions that return integer error codes.
  string GetErrorMessage(const int error) const;
  
  inline const string GetErrorDetails() { return ""; } // stub to support ErrorChecker template. Implement later if useful.

    //!	Reset the underlying RNA objects internal error code, after an error is handled
    void ResetError();


  //******************************************************
  // accessors and mutators
  //******************************************************

  //This function is use to set maxpairs.
  //! \param maxpairs is int value defining how the MaxPairs will be set. By default it is set to be -1, meaning the average length of all the sequences.
  //! \return an errorcode.
  int SetMaxPairs(const int maxpairs = -1);


  /// @brief    get the value of MaxPairs
  ///
  /// @return   the value of MaxPairs
  int GetMaxPairs() const;


  /// @brief    get the average length of the input sequences
  ///
  /// @return   the average length of the input sequences.
  int AverageLength() const;

  /// @brief    set the value of iterations
  ///
  /// @param    it is an value of int assigned to iterations. By default it is set to 2.
  ///
  /// @return   an errorcode
  int SetIterations(const int it = 2);

  /// @brief    get the value of iterations
  ///
  /// @return   the value of iterations.
  int GetIterations() const;


  /// @brief    set the value of MaxDsv/maxdsvchange
  ///
  /// @param    maxdsvchange is a value of float assigned to MaxDsv. By default it is set to 1.
  ///
  /// @return   an errorcode
  int SetMaxDsv(const float maxdsvchange = 1);

  /// @brief    get the value of MaxDsv/maxdsvchange
  ///
  /// @return   the value of MaxDsv/maxdsvchange.
  float GetMaxDsv() const;

  /// @brief    get the sequence number
  ///
  /// @return   the number of input sequences
  int GetSequenceNumber() const;

  /// @brief    set the Index Sequence for Multilign calculation.
  ///
  /// @param    indexSeq is a size_t value indicating which sequence is the index sequence; by default it is the 1st one.
  ///
  /// @return   a int value of ErrorCode
  int SetIndexSeq(size_t indexSeq = 1);


  /// @brief    an overloaded function accepting a string as its parameter.
  ///
  /// @param    seqname is the seq filename that will be set as the index sequence.
  ///
  /// @return   an int value of ErrorCode
  int SetIndexSeq(const string seqname);


  /// @brief    return the filename of the index seq.
  ///
  /// @return   a string of index seq filename.
  string GetIndexSeq() const;

  /// @brief    randomize the order of inputList.
  void Randomize();

  /// @brief    add one entry into inputList.
  ///
  /// @param    seq is a string value of sequence filename to be appended1
  /// @param    ct is a string value of corresponding ct filename
  /// @param    constraint is a string value of corresponding constraint filename. By default, it is empty, meaning no folding constraint exists
  /// @param    shape is string value of corresponding SHAPE filename. By default, it is empty, meaning no SHAPE exists.
  /// @return   a is int value of ErrorCode
  int AddOneInput(const string seq, const string ct, const string constraint = "", const string shape = "");

  /// @brief    remove one entry from inputList.
  ///
  /// @param    seq is a string value of sequence filename of which the entry in inputList will be removed
  ///
  /// @return   a int value of ErrorCode
  int RemoveOneInput(const string seq);

  /// @brief    set the slope parameter for SHAPE
  ///
  /// @param    slope is a double value assigned to SHAPESlope. By default, it is set to 1.8.
  void SetSHAPESlope(const double slope = 1.8);


  /// @brief    get the SHAPESlope
  ///
  /// @return   a SHAPESlope of double value.
  double GetSHAPESlope() const;


  /// @brief    set the intercept parameter for SHAPE.
  ///
  /// @param    intercept is a double value assigned to SHAPEIntercept. By default, it is set to -0.6.
  void SetSHAPEIntercept(const double intercept = -0.6);


  /// @brief    get the SHAPEIntercept.
  ///
  /// @return   SHAPEIntercept of double value.
  double GetSHAPEIntercept() const;


  /// @brief    set the temperature to fold the sequences.
  ///
  /// @param    temp is a double value of temperature; by default it is set to 310.15K
  void SetTemperature(const double temp = 310.15);


  /// @brief    get the temperature to fold the sequences
  ///
/// @return   a double value of the set temperature.
  double GetTemperature() const;


  /// @brief    set the flag isRNA to be true or false. By default it is true. When it is true, RNA nearest neighbor parameters are used.
  ///
  /// @param    isrna
//  void SetNucType(const bool isrna = true);


  /// @brief    get the type fo nucleic acid
  ///
  /// @return   return true when it is of RNA prediction; otherwise, false.
//  bool GetNucType() const;


  /// @brief    delete intermediate pairwise dsv and aout files
  ///
  /// @return   an int value of error code.
  int CleanupIntermediateFiles() const;

  /// @brief    Provide a TProgressDialog for following calculation progress.
  ///
  /// @param    Progress is a pointer to TProgressDialog
  void SetProgress(ProgressHandler *Progress = NULL);

  /// @brief    Provide a means to stop using a TProgressDialog by assigning NULL to progress pointer.
  void StopProgress();

  /// @brief    get the progress
  ///
  /// @return   the pointer to TProgressDialog
  ProgressHandler* GetProgress() const;

  ////////////////// The following functions are used for Diagnostic purpose only ////////////
  /// @brief    For diagnostic purpose only. Output the input sequence, ct, constraints, and SHAPE filenames to stdout.

  //! Generally not needed, but for debugging input.
  void GetInputFilenames();

  /// @brief    For diagnostic purpose only. Output the paired sequence filenames to stdout.
  void GetPairs();
  vector<float> get_energies(){return energies;}
  vector<float> get_dGIndex(){return dGIndex;}
  vector<vector<string> > get_pair_alignments(){return pair_alignments;}

protected:
  int ErrorCode;
  vector<string> input_alignment;
  vector<string> input_sequences;
  vector<string> ct_files;
  string output_multifind;
  vector <float> energies, dGIndex;
  mutable vector<vector<string> > pair_alignments;
 
private:
  ///////////////////// private functions ///////////////////////////////
  //! Pair sequences for dynalign calculation
  //! \return the value of ErrorCode
  int PairSeq1();
  int PairMultifindSeq1();
  /// @brief   This function check the legality of the input filenames and prepare the parameters for the multilign calculations. It should be called before multilign calcultions at least once and whenever something related to seq/ct changes, e.g. SetIndexSeq, AddOneInput, RemoveOneInput, Randomize, etc.
  ///
  /// @return   the int value of ErrorCode.
  int PrepInput();
  int PrepMultifindInput();

  /// @brief    move the element pointed by middle before the first-pointed element.
  ///
  /// @param    first is an vector<vector<string> >::iterator
  /// @param    middle is the same type.
  void ToHead(vvs_it first, vvs_it middle);


  /// @brief    name all the dsv files
  ///
  /// @return   an int value of errorcode
  int NameDsvFiles();
  int NameMultifindDsvFiles();

  /// @brief    name all the ali files
  ///
  /// @return   an int value of errorcode
  int NameAliFiles();
  int NameMultifindAliFiles();

  ///////////////////// private data members //////////////////////////////
  ProgressHandler *progress;
  // Currently, inputList is an matrix of 4 columns:
  // col 1 is the input seq filename
  // col 2 is the output ct filename
  // col 3 is the input constraint filename. If no constraint exists, set it to empty string at current row.
  // col 4 is the input SHAPE filename. If no SHAPE exists, set it to empty string at current row.
  vector <vector<string> > inputList; // the list of vectors each containing seq, ct, constraint, SHAPE filenames
  vector <pair<vs_index, vs_index> > seqPair; // a list of pairs of index for inputList
  //bool Dsv;// generate .dsv files or not (yes by default)
  //bool Ali;// generate .ali files or not (not by default)
  string **dsvFiles;
  string **aliFiles;
  /// The following are the parameters for Multilign calculations.
  int maxPairs;
  float maxDsv;
  int iterations;
  double SHAPESlope;
  double SHAPEIntercept;
  Dynalign_object *dynobj;
  Thermodynamics thermo;
};

#endif

/*
 * A program that verifies the format of several types of files used by the 
 * RNAstructure sofware package.
 *
 * (c) 2017-2018 Mathews Lab, University of Rochester Medical Center.
 * Written by Richard M. Watson
 */

#ifndef VALIDATE_INTERFACE_H
#define VALIDATE_INTERFACE_H

#include "../RNA_class/RNA.h"
#include "../src/ParseCommandLine.h"
#include "../src/RNAFileType.h"

#define SKIP_THERMO true // whether or not to read in full thermodynamic data tables (if false only the specification file will be read)

class validate_Interface {
 public:
  // Public constructor and methods.

  /*
   * Name:        Constructor.
   * Description: Initializes all private variables.
   */
  validate_Interface();

  /*
   * Name:        parse
   * Description: Parses command line arguments to determine what options are
   *              required for a particular calculation.
   * Arguments:
   *     1.   The number of command line arguments.
   *     2.   The command line arguments themselves.
   * Returns:
   *     True if parsing completed without errors, false if not.
   */
  bool parse( int argc, char** argv );

  /*
   * Name:        run
   * Description: Run calculations.
   */
  bool run();

 private:
  // Private variables.

  // The type of file to validate
  // This can be numeric (1=CT, 2=SEQ/FASTA, 3=PFS, 4=SAV, 5=DotBracket) 
  //    or text (one of: ct, seq|fasta, pfs, sav, dot|dbn|braket) 
  //    or "auto" (the default) which uses the file extension to determine the type.
  string fileType; // (default: SEQ)
  
  RNAFileType knownType; // fileType is parsed to determine a known RNAFileType enum constant.

  // Input and output file names.
  string file;

  int maxLength; // maximum sequence length (default: INT_MAX)
  string alphabetName; // the nucleic acid alphabet to use (default: rna)
  bool allowUnknownBases; // whether errors should be suppressed for unknown base names (default: false)

  // calculate the known type from the user-supplied value.
  RNAFileType parseFileType(const string &file, const string &type);

  string errorMessage;

  bool quiet;
  bool show_info;

   // tracks count of sequences for which info is added (if show_info is true)
  int info_seq_count;

  void addInfo(const string& sequence);
  void addInfo(RNA& rna);
  bool validateOligoList(); 
  bool validateChemFile();
  bool validateFoldConstraintFile();
  bool validateMultipleSequenceFile();
  bool validateRNAFile();
  
};

#endif /* VALIDATE_INTERFACE_H */

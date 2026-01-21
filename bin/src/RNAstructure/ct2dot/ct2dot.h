/*
 * A program that converts a CT file to a dot bracket file.
 *
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#ifndef CT2DOT_INTERFACE_H
#define CT2DOT_INTERFACE_H

#include "../RNA_class/RNA.h"
#include "../src/ParseCommandLine.h"

class ct2dot_Interface {
 public:
  // Public constructor and methods.

  /*
   * Name:        Constructor.
   * Description: Initializes all private variables.
   */
  ct2dot_Interface();

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
  int run();

 private:
  // Private variables.

  // Description of the calculation type.
  string calcType;

  // Input and output file names.
  string ctFile;                 // The input CT file.
  string bracketFile;            // The output dot bracket file.

  // The number of the structure, one-indexed, in the CT file to convert.
  int number;
  bool quiet;
  DotBracketFormat format;
};

#endif /* CT2DOT_INTERFACE_H */

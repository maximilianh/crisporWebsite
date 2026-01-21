/*
 * A header file for a class that takes two CT files, then compares their pairings and outputs the data graphically.
 *
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter.
 */

#ifndef STRUCTURE_COMPARED_IMAGE_HANDLER_H
#define STRUCTURE_COMPARED_IMAGE_HANDLER_H

#include <algorithm>
#include <iomanip>
#include <iterator>

#include "../src/StructureImageHandler.h"
#include "../src/score.h"

// A small namespace that holds comparison specific drawing constants.
namespace ComparisonConstants {

	// Postscript additional text constants.
	const double ADDED_GLYPH_SIZE_PS = 4;
	const int ADDED_TEXT_SIZE_PS = 8;
	const string ADDITIONAL_TEXT_OPEN_PS = "gsave /sfm2 { findfont exch makefont setfont } bind def [8 0 0 -8 0 0] /Courier-Bold sfm2";
	const string ADDITIONAL_TEXT_CLOSE_PS = "grestore";

	// SVG additional text constants.
	const double ADDED_GLYPH_SIZE_SVG = 5.5;
	const int ADDED_TEXT_SIZE_SVG = 10;
	const string ADDITIONAL_TEXT_OPEN_SVG = "<g font-family=\"monospace\" font-size=\"10\">";
	const string ADDITIONAL_TEXT_CLOSE_SVG = "</g>";

	// Standardized labels.
	const string LABEL_ACCEPTED = "Accepted:";
	const string LABEL_ACCEPTED_DESC = "Accepted Structure:";
	const string LABEL_ACCEPTED_FILE = "Accepted Structure File Name:";
	const string LABEL_PAIRS = "Pairs:";
	const string LABEL_PAIRS_PSEUDOKNOTTED = "Pseudoknotted Pairs:";
	const string LABEL_PREDICTED = "Predicted:";
	const string LABEL_PREDICTED_DESC = "Predicted Structure:";
	const string LABEL_PREDICTED_FILE = "Predicted Structure File Name:";
};

using namespace ComparisonConstants;

// The StructureComparedImageHandler class.
class StructureComparedImageHandler : public StructureImageHandler {
 public:
  // Public constructor and methods.

  /*
   * Name:        Constructor.
   * Description: Initializes private variables.
   */
  StructureComparedImageHandler();

  /*
   * Name:        addComparisonData
   * Description: Write additional information to show comparison statistics.
   * Arguments:
   *     1. predicted
   *        The predicted structure(s) file name.
   *     2. accepted
   *        The accepted structure file name.
   *     3. number
   *        The number of the predicted structure to compare to the accepted structure.
   *     4. exact
   *        True if only exact pairs are allowed for scoring, false if not.
   *     5. isSVG
   *        True if the image written is SVG, false if not.
   *     6. filenames
   *        True if file names should be shown with descriptions, false if not.
   */
  void addComparisonData( string predicted, string accepted, int number, bool exact, bool isSVG, bool filenames );

  /*
   * Name:        getNumPredictedStructures
   * Description: Get the number of predicted structures compared.
   * Returns:
   *     The number of predicted structures.
   */
  int getNumPredictedStructures();

  /*
   * Name:        overlayStructures
   * Description: Merge structure data into a single combined structure.
   */
  void overlayStructures();

  /*
   * Name:        readAcceptedStructure
   * Description: Read data for the accepted structure, without annotation.
   * Arguments:
   *     1. accepted
   *        The accepted structure file name.
   * Returns:
   *     A string showing the completion status.
   */
  string readAcceptedStructure( string accepted );

  /*
   * Name:        readAcceptedStructureProbability
   * Description: Read data for the accepted structure, with probability annotation.
   * Arguments:
   *     1. accepted
   *        The accepted structure file name.
   *     2. probabilityFile
   *        The probability annotation file.
   *     3. text
   *        True if the annotation file is a dot plot text file, false if not.
   * Returns:
   *     A string showing the completion status.
   */
  string readAcceptedStructureProbability( string accepted, string probabilityFile, bool text = false );

  /*
   * Name:        readPredictedStructure
   * Description: Read pairs for the predicted structure.
   * Arguments:
   *     1. predicted
   *        The predicted structure(s) file name.
   *     2. number
   *        The number of the predicted structure to compare to the accepted structure.
   * Returns:
   *     A string showing the completion status.
   */
  string readPredictedStructure( string predicted, int number );

  /*
   * Name:        readPredictedStructureProbability
   * Description: Read data for the predicted structure, with probability annotation.
   * Arguments:
   *     1. predicted
   *        The predicted structure file name.
   *     2. number
   *        The number of the predicted structure to compare to the accepted structure.
   *     3. probabilityFile
   *        The probability annotation file.
   *     4. text
   *        True if the annotation file is a dot plot text file, false if not.
   * Returns:
   *     A string showing the completion status.
   */
  string readPredictedStructureProbability( string predicted, int number, string probabilityFile, bool text = false );

  /*
   * Name:        readPredictedStructureSHAPE
   * Description: Read data for the predicted structure, with SHAPE annotation.
   * Arguments:
   *     1. predicted
   *        The predicted structure file name.
   *     2. number
   *        The number of the predicted structure to compare to the accepted structure.
   *     3. SHAPEFile
   *        The SHAPE annotation file.
   * Returns:
   *     A string showing the completion status.
   */
  string readPredictedStructureSHAPE( string predicted, int number, string SHAPEFile );

  /*
   * Name:        setColorScheme
   * Description: Set the color scheme used to distinguish between pairs.
   * Arguments:
   *     1. alternative
   *        True if the alternative color scheme is used, false if not.
   *        Default is false.
   */
  void setColorScheme( bool alternative );

  /*
   * Name:        validateFiles
   * Description: Check whether the predicted and accepted files can be properly compared.
   * Arguments:
   *     1. predicted
   *        The predicted structure(s) file name.
   *     2. accepted
   *        The accepted structure file name.
   *     3. number
   *        The number of the predicted structure to compare to the accepted structure.
   * Returns:
   *     A string showing the completion status.
   */
  string validateFiles( string predicted, string accepted, int number );

 private:
  // Private variables.

  // The vector that holds a particular structure's accepted pairs.
  vector<string> acceptedPairs;

  // The vector holding the color scheme.
  // Index 0: Color of common pairs
  // Index 1: Color of pairs in predicted structure only
  // Index 2: Color of pairs in accepted structure only
  vector<string> colorScheme;

  // The strings that hold structure descriptions.
  string descriptionAccepted;
  string descriptionPredicted;

  // The vector that holds a StructureComparedImageHandler specific legend.
  vector<string> legendCompared;

  // The number of predicted structures compared.
  int numPredictedStructures;

  // The vector that holds a particular structure's predicted pairs.
  vector<string> predictedPairs;
};

#endif /* STRUCTURE_COMPARED_IMAGE_HANDLER_H */

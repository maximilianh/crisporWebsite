//! RNAstructureBackendCalculator class.
//  (c) 2010 Mathews Lab, University of Rochester Medical Center.
//  Written by Jessica S. Reuter.
/*!
  The RNAstructureBackendCalculator class handles all calculations and data
  access for the Java GUI interface.

  Note that this class should NEVER, EVER be used in any other context except
  as a proxy for the RNAstructure GUI. Note also that error checking within
  this class is not bulletproof, because the RNAstructure GUI handles some of
  it, and blocks access to certain methods if errors occur. NEVER use this
  class for anything other than an RNAstructure GUI proxy, and when using this
  class or making additions, ALWAYS make sure that any further execution in
  this class is blocked, either by the GUI or by this class, at the first sign
  of an error.

  This header file is broken up into sections, listed as follows:
  1.  File Declaration Opening Statements.
  2.  Include Statements.
  3.  Namespace Usage and Main Class Declarations.
  4.  Constructor and Destructor.
  5.  AllSub Module.
  6.  Bifold Module.
  7.  Bipartition Module.
  8.  Dynalign Module.
  9.  Efn2 Module.
  10. Fold Module.
  11. MaxExpect Module.
  12. Multilign Module.
  13. OligoScreen Module.
  14. OligoWalk Module.
  15. Partition Module.
  16. ProbKnot Module.
  17. Refold Dynalign Module.
  18. Refold Single Structure Module.
  19. Remove Pseudoknots Module.
  20. Sequence Display Module.
  21. Stochastic Module.
  22. TurboFold Module.
  23. Constraints Mutators and Accessors.
  24. Utility Methods.
  25. Instance Variables.
  26. File Declaration Closing Statement.
 */

//////////////////////////////////////////////////////////////////////////////
// File Declaration Opening Statements.
//////////////////////////////////////////////////////////////////////////////

#ifndef RNASTRUCTURE_BACKEND_CALCULATOR_H
#define RNASTRUCTURE_BACKEND_CALCULATOR_H

//////////////////////////////////////////////////////////////////////////////
// Include Statements.
//////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>

#include "../../RNA_class/Dynalign_object.h"
#include "../../RNA_class/HybridRNA.h"
#include "../../RNA_class/Multilign_object.h"
#include "../../RNA_class/Oligowalk_object.h"
#include "../../RNA_class/RNA.h"
#include "../../src/ErrorChecker.h"
#include "../../src/TProgressDialog.h"
#include "../../src/version.h"
#include "../../src/TurboFold_object.h"


//////////////////////////////////////////////////////////////////////////////
// Namespace Usage and Main Class Declaration.
//////////////////////////////////////////////////////////////////////////////

using namespace std;

class RNAstructureBackendCalculator {

//////////////////////////////////////////////////////////////////////////////
// Constructor and Destructor.
//////////////////////////////////////////////////////////////////////////////
 public:

	//  Name:
	//! Constructor.
	//
	//  Description:
	//! Initializes all possible data structures, so they can be used whenever
	//! a need for them arises.
	RNAstructureBackendCalculator();

	//  Name:
	//! Destructor.
	//
	//  Description:
	//! Deletes any data structures that were used during a particular backend
	//! calculation or data access.
	~RNAstructureBackendCalculator();

//////////////////////////////////////////////////////////////////////////////
// AllSub Module.
//////////////////////////////////////////////////////////////////////////////
 public:

	//////////////////////////////////////////////////////////////////////////////
	//  Name:
 	//! loadDoubleRNA
 	//
 	//  Description:
 	//! Build a HybridRNA data structure for modules that take two sequences as input 
 	//! (e.g. Bifold, DuplexFold, and AccessFold calculations).
	//
 	//  Parameters:
	//! \param file1   The first file initializing this structure.
	//! \param file2   The second file initializing this structure.
	//! \param isRNA   Whether the nucs are RNA (true) or DNA (false)
	//! \param fileType   The type of file indicated by file1 and file2.  See RNAInputType for details.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string loadSingleRNA( string file1, const bool isRNA, const RNAInputType fileType=FILE_SEQ);


	//////////////////////////////////////////////////////////////////////////////
	//  Name:
 	//! loadDoubleRNA
 	//
 	//  Description:
 	//! Build a HybridRNA data structure for modules that take two sequences as input 
 	//! (e.g. Bifold, DuplexFold, and AccessFold calculations).
	//
 	//  Parameters:
	//! \param file1   The first file initializing this structure.
	//! \param file2   The second file initializing this structure.
	//! \param isRNA   Whether the nucs are RNA (true) or DNA (false)
	//! \param fileType   The type of file indicated by file1 and file2.  See RNAInputType for details.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string loadDoubleRNA( string file1, string file2, const bool isRNA, const RNAInputType fileType=FILE_SEQ);


 	//  Name:
 	//! buildAllSubDataStructure
 	//
 	//  Description:
 	//! Build a data structure for AllSub calculations.
 	//
 	//  Parameters:
	//! \param file    The file initializing this structure.
	//! \param isRNA   Whether the nucs are RNA (true) or DNA (false)
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string buildAllSubDataStructure( string file, bool isRNA);

	//  Name:
	//! getSuboptimalAbsoluteDiff
	//
	//  Description:
	//! Get the suboptimal structures absolute energy difference, based on the
	//! sequence length.
	//
	//  Returns:
	//! \return The absolute energy difference.
	double getSuboptimalAbsoluteDiff();

	//  Name:
	//! getSuboptimalPercentDiff
	//
	//  Description:
	//! Get the suboptimal structures maximum percent energy difference, based
	//! on the sequence length.
	//
	//  Returns:
	//! \return The maximum percent energy difference.
	float getSuboptimalPercentDiff();

	//  Name:
	//! runAllSub
	//
	//  Description:
	//! Run AllSub calculations.
	//
	//  Parameters:
	//! \param ct         The output CT file that is written to.
	//! \param percent    The maximum percent energy difference.
	//! \param absolute   The maximum absolute energy difference.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string runAllSub( string ct, float percent, double absolute );

 public:
 

	//////////////////////////////////////////////////////////////////////////////
	// AccessFold Module.
	//////////////////////////////////////////////////////////////////////////////
	//  Name:
 	//! buildAccessFoldDataStructure
 	//
 	//  Description:
 	//! Build a data structure for AccessFold calculations.
 	//
 	//  Parameters:
	//! \param file1   The first file initializing this structure.
	//! \param file2   The second file initializing this structure.
	//! \param isRNA   Whether the nucs are RNA (true) or DNA (false)
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string buildAccessFoldDataStructure( string file1, string file2, bool isRNA );

	//  Name:
	//! runAccessFold
	//
	//  Description
	//! Run bimolecular folding.
	//
	//  Parameters:
	//! \param ctFile           The output CT file that is written to.
	//! \param percent          The maximum percent energy difference.
	//! \param maxStructures    The maximum number of structures.
	//! \param gamma            The scaling factor for accessibility information.
	//! \param windowSize       The window size.
	//! \param saveFile         Boolean flag, true if a save file should be
	//!                         generated, false if not.
	//  Returns:
	//! \return A string that provides the completion status.
	string runAccessFold(string ctFile, float percent, int maxStructures, double gamma, int windowSize, bool saveFile);

	//////////////////////////////////////////////////////////////////////////////
	// DuplexFold Module.
	//////////////////////////////////////////////////////////////////////////////
	//  Name:
 	//! buildDuplexFoldDataStructure
 	//
 	//  Description:
 	//! Build a data structure for DuplexFold calculations.
 	//
 	//  Parameters:
	//! \param file1   The first file initializing this structure.
	//! \param file2   The second file initializing this structure.
	//! \param isRNA   Whether the nucs are RNA (true) or DNA (false)
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string buildDuplexFoldDataStructure( string file1, string file2, bool isRNA );

	//  Name:
	//! runDuplexFold
	//
	//  Description
	//! Run bimolecular folding.
	//
	//  Parameters:
	//! \param ctFile           The output CT file that is written to.
	//! \param percent          The maximum percent energy difference.
	//! \param maxStructures    The maximum number of structures.
	//! \param windowSize       The window size.
	//! \param saveFile         Boolean flag, true if a save file should be
	//!                         generated, false if not.
	//  Returns:
	//! \return A string that provides the completion status.
	string runDuplexFold(string ctFile, float percent, int maxStructures, int windowSize, bool saveFile);




	//////////////////////////////////////////////////////////////////////////////
    // Bifold Module.
    //////////////////////////////////////////////////////////////////////////////

 	//  Name:
 	//! buildBifoldDataStructure
 	//
 	//  Description:
 	//! Build a data structure for bifold calculations.
 	//
 	//  Parameters:
	//! \param file1   The first file initializing this structure.
	//! \param file2   The second file initializing this structure.
	//! \param isRNA   Whether the nucs are RNA (true) or DNA (false)
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string buildBifoldDataStructure( string file1, string file2, bool isRNA );

	//  Name:
	//! runBifold
	//
	//  Description
	//! Run bimolecular folding.
	//
	//  Parameters:
	//! \param ctFile           The output CT file that is written to.
	//! \param percent          The maximum percent energy difference.
	//! \param maxStructures    The maximum number of structures.
	//! \param windowSize       The window size.
	//! \param saveFile         Boolean flag, true if a save file should be
	//!                         generated, false if not.
	//! \param intraForbidden   Boolean flag, true if intramolecular pairs are
	//!                         forbidden, false if they're allowed
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string runBifold(string ctFile, float percent, int maxStructures, int windowSize,
		bool saveFile, bool intraForbidden );

//////////////////////////////////////////////////////////////////////////////
// Bipartition Module.
//////////////////////////////////////////////////////////////////////////////
 public:
 
 	//  Name:
 	//! buildBipartitionDataStructure
 	//
 	//  Description:
 	//! Build a data structure for bipartition calculations.
 	//
 	//  Parameters:
	//! \param file1   The first file initializing this structure.
	//! \param file2   The second file initializing this structure.
	//! \param isRNA   Whether the nucs are RNA (true) or DNA (false)
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string buildBipartitionDataStructure(
		string file1, string file2, bool isRNA );

	//  Name:
	//! runBipartition
	//
	//  Description:
	//! Run the bimolecular partition function.
	//
	//  Parameters:
	//! \param pfsFile   The output partition function save file.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string runBipartition( string pfsFile );

//////////////////////////////////////////////////////////////////////////////
// Dynalign Module.
//////////////////////////////////////////////////////////////////////////////
 public:
 
 	//  Name:
 	//! buildDynalignDataStructure
 	//
 	//  Description:
 	//! Build a data structure for Dynalign calculations.
 	//
 	//  Parameters:
	//! \param file1   The first file initializing this structure.
	//! \param file2   The second file initializing this structure.
	//! \param isRNA   Whether the nucs are RNA (true) or DNA (false).
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string buildDynalignDataStructure(
		string file1, string file2, bool isRNA );

	//  Name:
	//! clearDynalignAlignmentConstraints
	//
	//  Description:
	//! Clear all alignment constraints.
	void clearDynalignAlignmentConstraints();

	//  Name:
	//! getDynalignAlignmentConstraints
	//
	//  Description:
	//! Get the Dynalign alignment constraints as an HTML string.
	//
	//  Returns:
	//! The Dynalign alignment constraints string.
	string getDynalignAlignmentConstraints();

	//  Name:
	//! getDynalignAlignmentWindowSize
	//
	//  Description:
	//! Get the Dynalign alignment window size, based on the sequence length.
	//
	//  Returns:
	//! \return The alignment window size.
	int getDynalignAlignmentWindowSize();

	//  Name:
	//! getDynalignStructureWindowSize
	//
	//  Description:
	//! Get the Dynalign structure window size, based on the sequence length.
	//
	//  Returns:
	//! \return The structure window size.
	int getDynalignStructureWindowSize();

	//  Name:
	//! readDynalignAlignmentConstraintsFile
	//
	//  Description:
	//! Read a Dynalign alignment constraints text file.
	//
	//  Parameters:
	//! \param file   The alignment constraints file to read.
	//
	//  Returns:
	//! \return A string that provides the completion status
	string readDynalignAlignmentConstraintsFile( string file );


	int GetErrorCode();
	string GetFullErrorMessage();

	// ProgressHandler& GetProgress();
	// void SetProgress(ProgressHandler &progress, const bool deleteWhenFinished = false);

	//  Name:
	//! runDynalign
	//
	//  Description:
	//! Run Dynalign calculations.
	//
	//  Parameters:
	//! \param ctFile1      The output CT file sequence 1 is written to.
	//! \param ctFile2      The output CT file sequence 2 is written to.
	//! \param saveFile     The file save data is written to.
	//! \param alignFile    The file alignment data is written to.
	//! \param percent      The maximum percent energy difference.
	//! \param structures   The maximum number of structures.
	//! \param windowStr    The structure window size.
	//! \param windowAli    The alignment window size.
	//! \param gap          The gap penalty.
	//! \param isInsert     Whether base pair inserts are allowed or not.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	#ifdef DYNALIGN_II
	string runDynalign(
		string ctFile1, string ctFile2, string saveFile, string alignFile,
		double percent, int structures, int windowStr, int windowAli,
		float gap, float slope, float intercept, int max_elongation);
	#else
	string runDynalign(
		string ctFile1, string ctFile2, string saveFile, string alignFile,
		double percent, int structures, int windowStr, int windowAli,
		float gap, bool isInsert );
	#endif
	//  Name:
	//! setDynalignAlignmentConstraint
	//
	//  Description:
	//! Set a nucleotide in sequence 1 aligned with another nucleotide in
	//! sequence 2.
	//
	//  Parameters:
	//! \param nuc1   The nucleotide in sequence 1.
	//! \param nuc2   The nucleotide in sequence 2.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string setDynalignAlignmentConstraint( int nuc1, int nuc2 );

	//  Name:
	//! writeDynalignAlignmentConstraintsFile
	//
	//  Description:
	//! Write a Dynalign alignment constraints text file.
	//
	//  Parameters:
	//! \param file   The name of the Dynalign constraints file to write.
	void writeDynalignAlignmentConstraintsFile( string file );

//////////////////////////////////////////////////////////////////////////////
// Efn2 Module.
//////////////////////////////////////////////////////////////////////////////
 public:
 
 	//  Name:
 	//! buildEfn2DataStructure
 	//
 	//  Description:
 	//! Build a data structure for efn2 calculations.
 	//
 	//  Parameters:
	//! \param file    The file initializing this structure.
	//! \param isRNA   Whether the nucs are RNA (true) or DNA (false)
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string buildEfn2DataStructure( string file, bool isRNA );

	//  Name:
	//! runEfn2
	//
	//  Description:
	//! Run efn2 calculations.
	//
	//  Parameters:
	//! \param outFile        The output energy file that is written to.
	//! \param writeDetails   Boolean flag specifying whether a thermodynamic
	//!                       details file (true) or a simple list file
	//!                       (false) is written.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string runEfn2( string outFile, bool writeDetails );

//////////////////////////////////////////////////////////////////////////////
// Fold Module.
//////////////////////////////////////////////////////////////////////////////
 public:
 
 	//  Name:
 	//! buildFoldDataStructure
 	//
 	//  Description:
 	//! Build a data structure for Fold calculations.
 	//
 	//  Parameters:
	//! \param file    The file initializing this structure.
	//! \param isRNA   Whether the nucs are RNA (true) or DNA (false)
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string buildFoldDataStructure( string file, bool isRNA );

	//  Name:
	//! getFoldWindowSize
	//
	//  Description:
	//! Get the window size used for folding, based on the sequence length.
	//
	//  Returns:
	//! \return The refold window size.
	int getFoldWindowSize();

	//  Name:
	//! runFold
	//
	//  Description:
	//! Run folding calculations on a previously folded save file.
	//
	//  Parameters:
	//! \param ctFile       The output CT file that is written to.
	//! \param percent      The maximum percent energy difference.
	//! \param structures   The maximum number of structures.
	//! \param window       The window size.
	//! \param save         Whether a save file should be written or not.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string runFold(
		string ctFile, int percent, int structures, int window, bool save );

//////////////////////////////////////////////////////////////////////////////
// MaxExpect Module.
//////////////////////////////////////////////////////////////////////////////
 public:
 
 	//  Name:
 	//! buildMaxExpectDataStructure
 	//
 	//  Description:
 	//! Build a data structure for MaxExpect calculations.
 	//
 	//  Parameters:
	//! \param file   The file initializing this structure.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string buildMaxExpectDataStructure( string file );

	//  Name:
	//! runMaxExpect
	//
	//  Description
	//! Run maximum expected accuracy.
	//
	//  Parameters:
	//! \param ctFile       The output CT file that is written to.
	//! \param score        The maximum score difference.
	//! \param structures   The maximum number of structures.
	//! \param window       The window size.
	//! \param gamma        The gamma value.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string runMaxExpect(
		string ctFile, double score, int structures, int window,
		double gamma );

//////////////////////////////////////////////////////////////////////////////
// Multilign Module.
//////////////////////////////////////////////////////////////////////////////
 public:

	//  Name:
	//! activateMultilign
	//
	//  Description:
	//! Set the Multilign module active.
	void activateMultilign();

	//  Name:
	//! addMultilignPair
	//
	//  Description:
	//! Add a sequence/CT tuple to the Multilign files vector.
	//
	//  Parameters:
	//! \param seqFile   The sequence file.
	//! \param ctFile    The ct file in the .
	void addMultilignTuple( string seqFile, string ctFile);

	//  Name:
	//! deleteMultilignPair
	//
	//  Description:
	//! Delete a sequence/CT tuple from the Multilign files vector.
	//
	//  Parameters:
	//! \param index   The index of the tuple in the Multilign files vector.
	void deleteMultilignTuple( unsigned int index );

	//  Name:
	//! getMultilignCT
	//
	//  Description:
	//! Get the specified CT output file name used in a Multilign calculation.
	//
	//  Parameters:
	//! \param index   The index of the CT in the Multilign files vector,
	//!                one-indexed.
	string getMultilignCT( int index );

	//  Name:
	//! getMultilignMaxPairs
	//
	//  Description:
	//! Get the max pairs calculated by an underlying Multilign object.
	//
	//  Returns:
	//! \return The max pairs in the underlying Multilign object.
	int getMultilignMaxPairs();

	//  Name:
	//! getMultilignSequenceSetData
	//
	//  Description:
	//! Get the Multilign sequence set formatted as a string to be shown in
	//! the GUI input window.
	//
	//  Returns:
	//! \return The sequence set formatted as a string.
	string getMultilignSequenceSetData();

	//  Name:
	//! getNumMultilignSequences
	//
	//  Description:
	//! Get the number of Multilign sequences.
	//
	//  Returns:
	//! \return The number of Multilign sequences.
	int getNumMultilignSequences();

	//  Name:
	//! runMultilign
	//
	//  Description:
	//! Run a Multilign calculation.
	//
	//  Parameters:
	//! \param percent         The maximum percent energy difference.
	//! \param maxStructures   The naximum number of structures to predict.
	//! \param bpWindow        The base pair window.
	//! \param alignWindow     The alignment window.
	//! \param gap             The gap penalty.
	//! \param insert          Whether base pair inserts are allowed.
	//! \param maxDsvChange    The maximum dsv change.
	//! \param maxPairs        The maximum number of pairs allowed.
	//! \param cycles          The number of cycles Multilign goes through.
	//! \param alignFile       The multiple alignment file.
	//! \param saveFiles       Whether save and alignment files persist.
	//! \param isRNA           Whether the nucs are RNA (true) or DNA (false)
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string runMultilign(
		int percent, int maxStructures, int bpWindow, int alignWindow,
		double gap, bool insert, double maxDsvChange, int maxPairs,
		int cycles, string alignFile, bool saveFiles, bool isRNA );

//////////////////////////////////////////////////////////////////////////////
// OligoScreen Module.
//////////////////////////////////////////////////////////////////////////////
 public:
 
 	//  Name:
 	//! buildOligoScreenDataStructure
 	//
 	//  Description:
 	//! Build a data structure for OligoScreen calculations.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string buildOligoScreenDataStructure();

	//  Name:
	//! runOligoScreen
	//
	//  Description:
	//! Run an OligoScreen calculation.
	//
	//  Parameters:
	//! \param in      The input oligo list file.
	//! \param out     The output oligo report file.
	//! \param isRNA   Whether the nucleotides are RNA (true) or DNA (false).
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string runOligoScreen( string in, string out, bool isRNA );

//////////////////////////////////////////////////////////////////////////////
// OligoWalk Module.
//////////////////////////////////////////////////////////////////////////////
 public:

 	//  Name:
 	//! buildOligoWalkDataStructure
 	//
 	//  Description:
 	//! Build a data structure for OligoWalk calculations.
 	//
 	//  Parameters:
	//! \param file    The file initializing this structure.
	//! \param isRNA   Whether the nucs are RNA (true) or DNA (false).
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string buildOligoWalkDataStructure( string file, bool isRNA );

	//  Name:
	//! canFoldOligoOligo
	//
	//  Description:
	//! Get whether a particular oligo can be folded as bimolecular.
	//
	//  Parameters:
	//! \param index   The index on the target, one-indexed.
	//
	//  Returns:
	//! \return True if the oligo can be folded, false if not.
	bool canFoldOligoOligo( int index );

	//  Name:
	//! canFoldOligoSelf
	//
	//  Description:
	//! Get whether a particular oligo can be folded as unimolecular.
	//
	//  Parameters:
	//! \param index   The index on the target, one-indexed.
	//
	//  Returns:
	//! \return True if the oligo can be folded, false if not.
	bool canFoldOligoSelf( int index );

	//  Name:
	//! determineOligoMaximum
	//
	//  Description:
	//! Determine the maximum number of oligos allowed in the calculation.
	//
	//  Parameters:
	//! \param length   The length of the oligos.
	//
	//  Returns:
	//! \return The number of oligos.
	int determineOligoMaximum( int length );

	//  Name:
	//! foldOligo
	//
	//  Description:
	//! Fold an oligo, unimolecular or bimolecular.
	//
	//  Parameters:
	//! \param sequence      The oligo sequence.
	//! \param index         The index on the target, one-indexed.
	//! \param bimolecular   True if two oligos should be folded together,
	//!                      false if not.
	//! \param isRNA         True if the oligo is RNA, false if not.
	//! \param file          The CT file written after folding.
	void foldOligo(
		string sequence, int index, bool bimolecular, bool isRNA, string file );

	//  Name:
	//! getAllOligoData
	//
	//  Description:
	//! Get the data for all oligos as one big string.
	//
	//  Parameters:
	//! \param height   The height of the graph where oligo bars are displayed.
	//
	//  Returns:
	//! \return The oligo data.
	string getAllOligoData( int height );

	//  Name:
	//! getDisplayedOligo
	//
	//  Description:
	//! Get the oligo that should be displayed for an index on the target.
	//
	//  Parameters:
	//! \param index   The index on the target, one-indexed.
	//
	//  Returns:
	//! \return The oligo.
	string getDisplayedOligo( int index );

	//  Name:
	//! getGraphRegionBegin
	//
	//  Description:
	//! Get the index where oligos start binding.
	//
	//  Returns:
	//! \return The index, one-indexed;
	int getGraphRegionBegin();

	//  Name:
	//! getGraphRegionEnd
	//
	//  Description:
	//! Get the index where oligos stop binding.
	//
	//  Returns:
	//! \return The index, one-indexed;
	int getGraphRegionEnd();

	//  Name:
	//! getMostStableOligo
	//
	//  Description:
	//! Get the oligo number with the most stable overall delta G.
	//
	//  Returns:
	//! \return The most stable oligo, one-indexed.
	int getMostStableOligo();

	//  Name:
	//! getOligoLabelData
	//
	//  Description:
	//! Get all label data for a particular oligo, in HTML.
	//
	//  Parameters:
	//! \param index   The oligo, one-indexed.
	//
	//  Returns:
	//! \return The oligo data string, delimited by semicolons.
	string getOligoLabelData( int index );

	//  Name:
	//! getOligoTargetLength
	//
	//  Description:
	//! Get the oligo target length.
	//
	//  Returns:
	//! \return   The oligo target length.
	int getOligoTargetLength();

	//  Name:
	//! getOligoTargetSequence
	//
	//  Description:
	//! Get the oligo target sequence.
	//
	//  Returns:
	//! \return   The oligo target sequence.
	string getOligoTargetSequence();

	//  Name:
	//! runOligoWalk
	//
	//  Description:
	//! Run OligoWalk calculations.
	//
	//  Parameters:
	//! \param report       The output report file that is written to.
	//! \param mode         The calculation mode.
	//! \param chemistry    The nucleic acid chemistry type.
	//! \param suboptimal   Whether the calculation includes suboptimal
	//!                     structures or not.
	//! \param length       The oligo length.
	//! \param amount       The amount, in units, of oligo in concentration.
	//! \param unit         The oligo concentration unit.
	//! \param start        The beginning index of the target region.
	//! \param stop         The end index of the target region.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string runOligoWalk(
		string report, int mode, string chemistry, bool suboptimal,
		int length, int amount, string unit, int start, int stop );

//////////////////////////////////////////////////////////////////////////////
// Partition Module.
//////////////////////////////////////////////////////////////////////////////
 public:
 
 	//  Name:
 	//! buildPartitionDataStructure
 	//
 	//  Description:
 	//! Build a data structure for partition calculations.
 	//
 	//  Parameters:
	//! \param file    The first file initializing this structure.
	//! \param isRNA   Whether the nucs are RNA (true) or DNA (false)
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string buildPartitionDataStructure( string file, bool isRNA );

	//  Name:
	//! runPartition
	//
	//  Description:
	//! Run the partition function.
	//
	//  Parameters:
	//! \param pfsFile   The output partition function save file.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string runPartition( string pfsFile );

//////////////////////////////////////////////////////////////////////////////
// ProbKnot Module.
//////////////////////////////////////////////////////////////////////////////
 public:
 
 	//  Name:
 	//! buildProbKnotDataStructure
 	//
 	//  Description:
 	//! Build a data structure for ProbKnot calculations.
 	//
 	//  Parameters:
	//! \param file   The file initializing this structure.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string buildProbKnotDataStructure( string file );

	//  Name:
	//! runPseudoknotPrediction
	//
	//  Description:
	//! Predict pseudoknots from a data structure.
	//
	//  Parameters:
	//! \param ctFile       The output CT structure file that is written to.
	//! \param iterations   The number of iterations a prediction uses.
	//! \param helix        The minimum helix length allowed for a pseudoknot.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string runPseudoknotPrediction(
		string ctFile, int iterations, int helix );

//////////////////////////////////////////////////////////////////////////////
// Refold Dynalign Module.
//////////////////////////////////////////////////////////////////////////////
 public:
 
 	//  Name:
 	//! buildRefoldDynalignDataStructure
 	//
 	//  Description:
 	//! Build a data structure for Dynalign refolding calculations.
 	//
 	//  Parameters:
	//! \param file   The file initializing this structure.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string buildRefoldDynalignDataStructure( string file );

	//  Name:
	//! runDynalignRefold
	//
	//  Description:
	//! Run a Dynalign refolding calculation.
	//
	//  Parameters:
	//! \param ctFile1      The first output CT file that is written to.
	//! \param ctFile2      The second output CT file that is written to.
	//! \param alignFile    The first output CT file that is written to.
	//! \param percent      The maximum percent energy difference.
	//! \param structures   The maximum number of structures.
	//! \param windowStr    The structure window size.
	//! \param windowAli    The alignment window size.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string runDynalignRefold(
		string ctFile1, string ctFile2, string alignFile, int percent,
		int structures, int windowStr, int windowAli );

//////////////////////////////////////////////////////////////////////////////
// Refold Single Structure Module.
//////////////////////////////////////////////////////////////////////////////
 public:
 
  	//  Name:
 	//! buildRefoldSingleDataStructure
 	//
 	//  Description:
 	//! Build a data structure for single structure refolding calculations.
 	//
 	//  Parameters:
	//! \param file   The file initializing this structure.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string buildRefoldSingleDataStructure( string file );

	//  Name:
	//! getRefoldWindowSize
	//
	//  Description:
	//! Get the window size used for refolding, based on the sequence length.
	//
	//  Returns:
	//! \return The refold window size.
	int getRefoldWindowSize();

	//  Name:
	//! runRefold
	//
	//  Description:
	//! Run folding calculations on a previously folded save file.
	//
	//  Parameters:
	//! \param ctFile       The output CT file that is written to.
	//! \param percent      The maximum percent energy difference.
	//! \param structures   The maximum number of structures.
	//! \param window       The window size.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string runRefold(
		string ctFile, int percent, int structures, int window );

//////////////////////////////////////////////////////////////////////////////
// Remove Pseudoknots Module.
//////////////////////////////////////////////////////////////////////////////
 public:
 
 	//  Name:
 	//! buildRemovePseudoknotsDataStructure
 	//
 	//  Description:
 	//! Build a data structure for pseudoknot removal calculations.
 	//
 	//  Parameters:
	//! \param file    The file initializing this structure.
	//! \param isRNA   Whether the nucs are RNA (true) or DNA (false).
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string buildRemovePseudoknotsDataStructure( string file, bool isRNA );

	//  Name:
	//! runPseudoknotRemoval
	//
	//  Description:
	//! Remove pseudoknots from a data structure.
	//
	//  Parameters:
	//! \param ctFile     The output CT structure file that is written to.
	//! \param minimize   Whether free energy should be minimized (true)
	//!                   or not (false ).
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string runPseudoknotRemoval( string ctFile, bool minimize );

//////////////////////////////////////////////////////////////////////////////
// Sequence Display Module.
//////////////////////////////////////////////////////////////////////////////
 public:

	//  Name:
	//! getSequenceComments
	//
	//  Description:
	//! Get the sequence's comment.
	//
	//  Returns:
	//! \return The sequence comment.
	string getSequenceComment();

	//  Name:
	//! getSequenceData
	//
	//  Description:
	//! Get the raw sequence.
	//
	//  Returns:
	//! \return The sequence data.
	string getSequenceData();

	//  Name:
	//! getSequenceTitle
	//
	//  Description:
	//! Get the sequence's title.
	//
	//  Returns:
	//! \return The sequence title.
	string getSequenceTitle();

	//  Name:
	//! readSequenceData
	//
	//  Description:
	//! Read data from a specific sequence.
	//
	//  Parameters:
	//! \param file   The file to read sequence data from.
	void readSequenceData( string file );

	//  Name:
	//! setSequenceComment
	//
	//  Description:
	//! Set the sequence's comment.
	//
	//  Parameters:
	//! \param comment   The sequence comment.
	void setSequenceComment( string comment );

	//  Name:
	//! setSequenceData
	//
	//  Description:
	//! Set the raw sequence.
	//
	//  Parameters:
	//! \param data   The sequence data.
	void setSequenceData( string data );

	//  Name:
	//! setSequenceTitle
	//
	//  Description:
	//! Set the sequence's title.
	//
	//  Parameters:
	//! \param title   The sequence title.
	void setSequenceTitle( string title );

	//  Name:
	//! writeFastaFile
	//
	//  Description:
	//! Write a FASTA file from the raw sequence data given.
	//
	//  Parameters:
	//! \param file   The name of the file to write.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string writeFastaFile( string file );

	//  Name:
	//! writeSequenceFile
	//
	//  Description:
	//! Write a sequence file from the raw sequence data given.
	//
	//  Parameters:
	//! \param file   The name of the file to write.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string writeSequenceFile( string file );

//////////////////////////////////////////////////////////////////////////////
// Stochastic Module.
//////////////////////////////////////////////////////////////////////////////
 public:
 
 	//  Name:
 	//! buildStochasticDataStructure
 	//
 	//  Description:
 	//! Build a data structure for stochastic calculations.
 	//
 	//  Parameters:
	//! \param file   The file initializing this structure.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string buildStochasticDataStructure( string file );

	//  Name:
	//! runStochastic
	//
	//  Description:
	//! Run stochastic sampling calculations.
	//
	//  Parameters:
	//! \param outFile    The output CT structure file that is written to.
	//! \param ensemble   The ensemble size.
	//! \param seed       The random seed.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string runStochastic( string outFile, int ensemble, int seed );

//////////////////////////////////////////////////////////////////////////////
// TurboFold Module.
//////////////////////////////////////////////////////////////////////////////
 public:

	//  Name:
	//! activateTurboFold
	//
	//  Description:
	//! Set the TurboFold module active.
	void activateTurboFold();

	//  Name:
	//! addTurboFoldTuple
	//
	//  Description:
	//! Add a sequence/CT tuple to the TurboFold files vector.
	//
	//  Parameters:
	//! \param seqFile   The sequence file.
	//! \param ctFile    The ct file.
	void addTurboFoldTuple( string seqFile, string ctFile, string pfsFile="");

	//  Name:
	//! deleteTurboFoldTuple
	//
	//  Description:
	//! Delete a sequence/CT tuple from the TurboFold files vector.
	//
	//  Parameters:
	//! \param index   The index of the tuple in the TurboFold files vector.
	void deleteTurboFoldTuple( unsigned int index );

	//  Name:
	//! getNumTurboFoldSequences
	//
	//  Description:
	//! Get the number of TurboFold sequences.
	//
	//  Returns:
	//! \return The number of TurboFold sequences.
	int getNumTurboFoldSequences();

	//  Name:
	//! getTurboFoldCT
	//
	//  Description:
	//! Get the CT output file associated with a TurboFold input sequence.
	//
	//  Parameters:
	//! \param index   The sequence number to get the CT for, one-indexed.
	//
	//  Returns:
	//! \return The appropriate CT file name.
	string getTurboFoldCT( int index );

	//  Name:
	//! getTurboFoldSaveFile
	//
	//  Description:
	//! Get the partition function save file associated with a TurboFold input
	//! sequence.
	//
	//  Parameters:
	//! \param index   The sequence number to get the save file for,
	//!                one-indexed.
	//
	//  Returns:
	//! \return The appropriate save file name.
	string getTurboFoldSaveFile( int index );

	//  Name:
	//! getTurboFoldSequenceSetData
	//
	//  Description:
	//! Get the TurboFold sequence set formatted as a string to be shown in
	//! the GUI input window. Note that this string doesn't include the
	//! partition function save files, they're just used behind the scenes.
	//
	//  Returns:
	//! \return The sequence set formatted as a string.
	string getTurboFoldSequenceSetData();

	//  Name:
	//! runTurboFoldMaximumExpectedAccuracy
	//
	//  Description:
	//! Run a TurboFold calculation in maximum expected accuracy mode.
	//
	//  Parameters:
	//! \param turboGamma        The TurboFold gamma value.
	//! \param turboIterations   The number of TurboFold iterations.
	//! \param percent           The maximum percent energy difference.
	//! \param structures        The maximum number of structures.
	//! \param window            The window size.
	//! \param meaGamma          The maximum expected accuracy gamma value.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string runTurboFoldMaximumExpectedAccuracy(
		double turboGamma, int turboIterations, double percent,
		int structures, int window, double meaGamma, string outAlnFile="" );

	//  Name:
	//! runTurboFoldPseudoknot
	//
	//  Description:
	//! Run a TurboFold calculation in pseudoknot mode.
	//
	//  Parameters:
	//! \param turboGamma        The TurboFold gamma value.
	//! \param turboIterations   The number of TurboFold iterations.
	//! \param pkIterations      The number of ProbKnot iterations.
	//! \param helix             The minimum helix length.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string runTurboFoldPseudoknot(
		double turboGamma, int turboIterations, int pkIterations, int helix, string outAlnFile="" );

  	//  Name:
	//! runTurboFoldThreshold
	//
	//  Description:
	//! Run a TurboFold calculation in threshold mode.
	//
	//  Parameters:
	//! \param turboGamma        The TurboFold gamma value.
	//! \param turboIterations   The number of TurboFold iterations.
	//! \param threshold         The probability cutoff threshold for pairs.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string runTurboFoldThreshold(
		double turboGamma, int turboIterations, double threshold, string outAlnFile="" );

 private:

	//  Name:
	//! runTurboFold
	//
	//  Description:
	//! Run a TurboFold calculation.
	//
	//  Parameters:
	//! \param gammaT       The TurboFold gamma value.
	//! \param turboIter    The number of TurboFold iterations.
	//! \param mode         The mode in which TurboFold runs.
	//!                     (Can be "MEA", "ProbKnot", or "Threshold")
	//! \param percent      The maximum percent energy difference.
	//!                     (Used in MEA mode only.)
	//! \param structures   The maximum number of structures.
	//!                     (Used in MEA mode only.)
	//! \param window       The window size.
	//!                     (Used in MEA mode only.)
	//! \param gammaM       The maximum expected accuracy gamma value.
	//!                     (Used in MEA mode only.)
	//! \param pkIter       The number of ProbKnot iterations.
	//!                     (Used in ProbKnot mode only.)
	//! \param helix        The minimum helix length.
	//!                     (Used in ProbKnot mode only.)
	//! \param cutoff       The probability cutoff threshold for base pairs.
	//!                     (Used in Threshold mode only.)
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string runTurboFold(
		double gammaT, int turboIter, string mode, double percent,
		int structures, int window, double gammaM, int pkIter, int helix,
		double cutoff, string outputAlnFile="");

	//void SetError(const int code, const string &details);
//////////////////////////////////////////////////////////////////////////////
// Constraints Mutators and Accessors.
//////////////////////////////////////////////////////////////////////////////
 public:

	//  Name:
	//! clearFoldingConstraints
	//
	//  Description:
	//! Remove all folding constraints from a particular strand.
	//
	//  Parameters:
	//! \param strand   The strand to remove constraints from.
	void clearFoldingConstraints( int strand );

	//  Name:
	//! getFoldingConstraints
	//
	//  Description:
	//! Get all folding constraints for a specific strand as an HTML string.
	//
	//  Parameters:
	//! \param strand   The strand to get constraints from.
	//
	//  Returns:
	//! \return   The constraints string.
	string getFoldingConstraints( int strand );

	//  Name:
	//! getMaxConstraintIndex
	//
	//  Description:
	//! Get the maximum index on which constraints can be specified.
	//
	//  Parameters:
	//! \param strand   The strand the nucleotide resides in.
	//
	//  Returns:
	//! \return The maximum constrainable index.
	int getMaxConstraintIndex( int strand );

	//  Name:
	//! getMaxLoop
	//
	//  Description:
	//! Get the maximum internal/bulge loop size.
	//
	//  Returns:
	//! \return The loop size.
	int getMaxLoop();

	//  Name:
	//! getMaxPair
	//
	//  Description:
	//! Get the maximum distance allowed between pairing nucleotides.
	//
	//  Returns:
	//! \return The maximum pairing distance.
	int getMaxPair();

	//  Name:
	//! getTemperature
	//
	//  Description:
	//! Get the temperature at which calculations occur.
	//
	//  Returns:
	//! \return The temperature.
	double getTemperature();

	//  Name:
	//! readFoldingConstraintsFile
	//
	//  Parameters:
	//! \param file     The folding constraints file to read.
	//! \param strand   The strand to read constraints into.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string readFoldingConstraintsFile( string file, int strand );

	//  Name:
	//! setCleavedNucleotide
	//
	//  Description:
	//! Set a nucleotide to be cleaved.
	//
	//  Parameters:
	//! \param nuc      The nucleotide to be cleaved.
	//! \param strand   The strand the nucleotide resides in.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string setCleavedNucleotide( int nuc, int strand );

	//  Name:
	//! setDoubleStrandedNucleotide
	//
	//  Description:
	//! Set a nucleotide to be double stranded.
	//
	//  Parameters:
	//! \param nuc      The nucleotide to be set double stranded.
	//! \param strand   The strand the nucleotide resides in.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string setDoubleStrandedNucleotide( int nuc, int strand );

	//  Name:
	//! setForcedHelix
	//
	//  Description:
	//! Set a forced helix.
	//
	//  Parameters:
	//! \param nuc1     The starting nucleotide of the helix.
	//! \param nuc2     The ending nucleotide of the helix.
	//! \param length   The length of the helix.
	//! \param strand   The strand the helix resides in.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string setForcedHelix( int nuc1, int nuc2, int length, int strand );

	//  Name:
	//! setMaxLoop
	//
	//  Description:
	//! Set the maximum internal/bulge loop size.
	//
	//  Parameters:
	//! \param newLoop   The new loop size.
	void setMaxLoop( int newLoop );

	//  Name:
	//! setMaxPair
	//
	//  Description:
	//! Set the maximum distance allowed between pairing nucleotides.
	//
	//  Parameters:
	//! \param newPair   The new maximum pairing distance.
	void setMaxPair( int newPair );

	//  Name:
	//! setModifiedNucleotide
	//
	//  Description:
	//! Set a nucleotide to be chemically modified.
	//
	//  Parameters:
	//! \param nuc      The nucleotide to be modified.
	//! \param strand   The strand the nucleotide resides in.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string setModifiedNucleotide( int nuc, int strand );

	//  Name:
	//! setProhibitedHelix
	//
	//  Description:
	//! Set a prohibited helix.
	//
	//  Parameters:
	//! \param nuc1     The starting nucleotide of the helix.
	//! \param nuc2     The ending nucleotide of the helix.
	//! \param length   The length of the helix.
	//! \param strand   The strand the helix resides in.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string setProhibitedHelix( int nuc1, int nuc2, int length, int strand );

	//  Name:
	//! setSHAPEFile
	//
	//  Description:
	//! Set the file to read SHAPE constraints from.
	//
	//  Parameters:
	//! \param file   The file to read constraints from.
	void setSHAPEFile( string file );

	//  Name:
	//! setSHAPEParam1
	//
	//  Description:
	//! Set the first SHAPE constraints parameter.
	//
	//  Parameters:
	//! \param value   The value of the parameter.
	void setSHAPEParam1( double value );

	//  Name:
	//! setSHAPEParam2
	//
	//  Description:
	//! Set the second SHAPE constraints parameter.
	//
	//  Parameters:
	//! \param value   The value of the parameter.
	void setSHAPEParam2( double value );

	//  Name:
	//! setSHAPEType
	//
	//  Description:
	//! Set the type of SHAPE constraints being used.
	//
	//  Parameters:
	//! \param isEnergy   True if pseudo-energy constraints, false if not.
	void setSHAPEType( bool isEnergy );

	//  Name:
	//! setSingleStrandedNucleotide
	//
	//  Description:
	//! Set a nucleotide to be single stranded.
	//
	//  Parameters:
	//! \param nuc      The nucleotide to be set single stranded.
	//! \param strand   The strand the nucleotide resides in.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string setSingleStrandedNucleotide( int nuc, int strand );

	//  Name:
	//! setTemperature
	//
	//  Description:
	//! Set the temperature at which calculations occur.
	//
	//  Parameters:
	//! \param newTemp   The new temperature.
	void setTemperature( double newTemp );

	//  Name:
	//! writeFoldingConstraintsFile
	//
	//  Parameters:
	//! \param file     The folding constraints file to write.
	//! \param strand   The strand to read constraints into.
	//
	//  Returns:
	//! \return A string that provides the completion status.
	string writeFoldingConstraintsFile( string file, int strand );

//////////////////////////////////////////////////////////////////////////////
// Utility Methods.
//////////////////////////////////////////////////////////////////////////////
 public:


	///////////////////////////////////////////////////////////////////////////
	// Class-wide Progress and Status functions.
	///////////////////////////////////////////////////////////////////////////
	//  Name:
	//! getProgressNumber
	//
	//  Description:
	//! Get the percent progress of a particular calculation, as an integer
	//! from 1 to 100.
	//
	//  Returns:
	//! \return The percent progress.
	int getProgressNumber();

	//  Name:
	//! getStructureType
	//
	//  Description:
	//! Get the type of structure being handled. The type result will be
	//! either "Radial", "Circular", or an error string.
	//
	//  Returns:
	//! \return A string that holds the type check result.
	string getStructureType();

	//  Name:
	//! cancelOperation
	//
	//  Description:
	//! Cancel the currently running operation. (Only when implemented in native methods...currently only FoldSingleWindow.)
	void cancelOperation();

	//  Name:
	//! wasCanceled
	//
	//  Description:
	//! Returns true if the most recent operation was canceled.
	bool wasCanceled();

	//  Name:
	//! setSaveFile
	//
	//  Description:
	//! Typically functions that can generate a saveFile (e.g. Fold)
	//! infer the name of the saveFile from the output file name.
	//! however, the user can override this by calling setSaveFile to 
	//! set the name explicitly before calling runFold etc.
	void setSaveFile(string path);

	//  Name:
	//! getSaveFile
	//
	//  Description:
	//! Returns the string set by setSaveFile or the empty string if 
	//! it has not been set.
	//
	//  Returns:
	//! \return A string that holds the current saveFile path, if set explicitly by the user.
	string getSaveFile();
	
	///////////////////////////////////////////////////////////////////////////
	// Static Utility functions 
	// (i.e. to enable features not available to Java)
	///////////////////////////////////////////////////////////////////////////

	//  Name:
	//! setEnvVar
	//
	//  Description:
	//! Set a system environment variable for the current process. This function adds new 
	//! environment variables or modifies the values of existing environment variables. 
	//! Environment variables define the environment in which a process executes (for 
	//! example, the default search path for libraries to be linked with a program).
	//
	//  Parameters:
	//! \param envstr     A string specifying the name and value of the environment variable to set. The 
	//!                   string should be formatted as "NAME=VALUE".
	//
	//  Returns:
	//! \return The return value is 0 if successful or non-zero in the case of an error.
	static int setEnvVar(const string& envstr);
	static string getEnvVar(const string& varname);
	//  Name:
	//! getVersion
	//
	//  Description:
	//! Get the official build version of the RNAstructure library (defined in src/version.h)
	//  Returns:
	//! \return The return value is a string containing the version number, e.g. "5.8"
	static string getVersion();


//////////////////////////////////////////////////////////////////////////////
// Instance Variables.
//////////////////////////////////////////////////////////////////////////////
 private:

	// Definition of a specific struct to hold constraints that aren't forced
	// immediately while doing calculations, and that aren't specific to a
	// single particular module.
	struct ConstraintsHolder {
		int maxLoop;
		int maxPair;
		bool shapeEnergy;
		string shapeFile;
		double shapeParam1;
		double shapeParam2;
		double temperature;
	} constraintsHolder;

	// A Dynalign data structure, its alignment constraints data structure,
	// and its error checker.
	Dynalign_object* dynalign;
	map<int, int> dynalignAlignments;
	ErrorChecker<Dynalign_object>* dynalignChecker;
	string saveFilePath;

	// A hybrid data structure and its error checker.
	HybridRNA* hybrid;
	ErrorChecker<HybridRNA>* hybridChecker;

	// The monitor to monitor progress where necessary.
	// RMW - removed ProgressMonitor, as it was not necessary: ProgressMonitor* monitor;
	ProgressHandler* progress;
	bool ownsProgress;

	// Definition of a specific struct to hold file tuples and an activation
	// boolean, and then two of these structs, one to hold Multilign data and
	// one to hold TurboFold data, if those modules are called.
	struct MultiSeqHandler {
		vector< vector<string> > files;
		bool isActive;
	};
	MultiSeqHandler multilign;
	MultiSeqHandler turboFold;

	// An oligo data structure and its error checker.
	// Also, a specific struct to hold certain pieces of oligo data.
	Oligowalk_object* oligo;
	ErrorChecker<Oligowalk_object>* oligoChecker;
	struct OligoDataHolder {
		string concentration;
		string deltaG;
		int graphBegin;
		int graphEnd;
		int mostStable;
		int oligoLength;
		string oligoType;
	} oligoDataHolder;

	// A single strand data structure and its error checker.
	RNA* rna;
	ErrorChecker<RNA>* rnaChecker;

	// A vector that holds information about a raw sequence.
	// Index 0: Title
	// Index 1: Comment
	// Index 2: Sequence Data
	vector<string> sequenceHandler;

	// Resets the progress dialog.
	void resetProgress();

	// This function is used to get the saveFile name that should be used
	// for functions like Fold etc. that can generate save files.
	// If saveFilePath is NOT empty (i.e. the user has set saveFilePath explicitly
	// via setSaveFile) then saveFilePath is returned as-is.
	// Otherwise, this function generates a new saveFilePath derived from 
	// outputFile by changing the file extension to ".sav". If outputFile
	// has no extension, then ".sav" is appended.
	// If savExtension is specified, that extension is used instead of "sav".
	// The string passed in for savExtension should NOT include the initial dot.
	string generateSavFile(const string outputFile, const char* const savExtension = NULL);

	// string errorDetails;
	// int ErrorCode;
};

//////////////////////////////////////////////////////////////////////////////
// File Declaration Closing Statement.
//////////////////////////////////////////////////////////////////////////////

#endif /* RNASTRUCTURE_BACKEND_CALCULATOR_H */

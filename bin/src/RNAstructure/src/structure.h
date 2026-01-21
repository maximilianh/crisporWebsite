#ifndef STRUCTURE_H
#define STRUCTURE_H

#include <string>
#include <stdlib.h>
#include <vector>
#include "defines.h"
#include "rna_library.h"
#include "CTCommentProvider.h"
#include "DotBracketFormat.h"
#include "../src/phmm/utils/xmath/log/xlog_math.h"

using namespace std;

#ifdef _WINDOWS_GUI
	#include "../Windows_interface_deprecated/platform.h"
#else
	#include "platform.h"
#endif //_WINDOWS

#ifdef EXTENDED_DOUBLE
	#include "extended_double.h" //inlcude code for extended double if needed
#endif//defined EXTENDED_DOUBLE

//This is a deprecated array size that must be removed:
#define maxforce 3000
#define DEFAULT_SHAPE_INTERCEPT -0.6  // default intercept used in SHAPE data pseudoenergy restraint generation (kcal/mol)
#define DEFAULT_SHAPE_SLOPE 1.8       // default slope used in SHAPE data pseudoenergy restraint generation (kcal/mol)

// Enable VERIFY_ARRAY_BOUNDS to force bounds-checking on arrays (for debugging and development only. Should be disabled for releases!)
//#define VERIFY_ARRAY_BOUNDS 

#ifdef VERIFY_ARRAY_BOUNDS 
	#define VERIFY_NUC_INDEX(INDEX) while (INDEX<1||INDEX>GetSequenceLength()) { cerr << "Nucleotide index out of bounds: " << INDEX << "(length " << GetSequenceLength() << ") at " << __FILE__ << ":" << __LINE__ << " base: " << GetBase(i) << endl; break; }
	#define VERIFY_NUC_INDEX_CT(INDEX,CT) while (INDEX<1||INDEX>CT->GetSequenceLength()) { cerr << "Nucleotide index out of bounds: " << INDEX << "(length " << CT->GetSequenceLength() << ") at " << __FILE__ << ":" << __LINE__ << endl; break; }
	#define VERIFY_STRUCTURE_INDEX(INDEX) while (INDEX<1||INDEX>GetNumberofStructures()) { cerr << "Structure index out of bounds: " << INDEX << "(length " << CT->GetNumberofStructures() << ") at " << __FILE__ << ":" << __LINE__ << endl; 
	#define VERIFY_STRUCTURE_INDEX_CT(INDEX,CT) while (INDEX<1||INDEX>CT->GetNumberofStructures()) { cerr << "Structure index out of bounds: " << INDEX << "(length " << CT->GetNumberofStructures() << ") at " << __FILE__ << ":" << __LINE__ << endl; break; }
#else
	// Define all the VERIFY_ARRAY_BOUNDS macros to do nothing.
	#define VERIFY_NUC_INDEX(INDEX)			
	#define VERIFY_STRUCTURE_INDEX(INDEX)
	#define VERIFY_NUC_INDEX_CT(INDEX,CT)
	#define VERIFY_STRUCTURE_INDEX_CT(INDEX,CT)
#endif

//Single structure is a wrapper for the information associated with just a single structure, i.e. pairs, energy, and labels.
struct singlestructure {

	//This function sizes the vectors to the approprate size:
	singlestructure(int sequencelength);

	//keep track of the basepairs
	vector<int> basepr;

	//keep track of the energy of the structure, if available
	int energy;

	//keep a string, from a ct file or sequence file with the sequence desription
	string ctlabel;

};

//! Constants used for the 'modifier' parameter for pseudo-energy-related functions
enum RestraintType { RESTRAINT_SHAPE, RESTRAINT_SHAPE_DIFF, RESTRAINT_SHAPE_AC,	RESTRAINT_SHAPE_GU, RESTRAINT_DMS, RESTRAINT_CMCT, RESTRAINT_DMSNT };

//! structure Class.
/*!
	The structure class provides a class for handling sequences and structures.
*/
//////////////////////////////////////////////////////////////////////
class structure //this structure contains all the info for a structure
{
	public:
		//******************************
		//Constructor:
		//******************************

		//!Constructor.
		//!	\param sructures is an int that specifies how many structures should be anticipated.  This sets up an initial memory allocation, but this can expand as needed.
		structure(int structures = maxstructures + 1);

		//!Destructor.
		~structure();

		//*********************************
		//Get and receive sequence and structure information:
		//*********************************

		//! Get the label for structure numer structurenumber.

		//! \param structurenumber is an int that gives the structure number being indexed, this is one indexed.
		//! \return a string that gives the label.
		string GetCtLabel(int structurenumber) const;

		//! Get the energy for structure number structurenumber.

		//! This function requires that an energy calculation has been performed.  It does not do an eenrgy calculation.
		//! \param structurenumber is an int that gives the structure number being indexed, this is one indexed.
		//! \return an int that gives the energy in kcal/mol*conversionfactor, where conversionfactor is set in /src/defines.h.
		int GetEnergy(int structurenumber) const;
		
		//! Get the number of structures stored.

		//! \return An integer that is the number of structures encoded.
		int GetNumberofStructures() const;

		//! Get the pairing partner for i in structure structurenumber.

		//! \param i is the nucleotide index, which is one-indexed.
		//! \param structurenumber is the structure number, which is one-indexed.
		//! \return The pairing partner, as a nucleotide index.
		int GetPair(int i, int structurenumber=1) const;

		//! Get the base number at position i  (i.e. the index of the base in the alphabet).
		inline int GetBase(int i) const { return numseq[i]; }

		//! Get the label associated with the sequence.

		//! \return The string read from the sequence file.
		string GetSequenceLabel() const;

		//! Get the length of the sequence.

		//! \return An integer that is the sequence length.
		inline int GetSequenceLength() const {
			//Return the value of numofbases:
			return numofbases;
		}

		const char* GetSequence() const;
		
		//! Remove the pair at index i for structure number structurenumber.

		//! If i is paired to j, the pairing for j is also removed.
		//! \param i is the nucleotide for whic pairing is to be removed.  This is one-indexed.
		//! \param structurenumber is an int that provides from which structure the pair should be removed.  This is one-indexed.
		void RemovePair(int i, int structurenumber=1);

		//! Determine if the structure has one or more Pseudoknots (crossing bonds).
		bool HasPseudoknots(int structurenumber=1) const ;

		//! This function determines which base pairs in a structure are "pseudoknots" by finding the 
		//! the largest subset of non-crossing bonds.
		//! The caller can retrieve the list of "normal" bonds (normalPairs) or the list of 
		//! pseudoknots (pseudoknotPairs) or both.
		//!
		//! This function does NOT modify the structure itself.
		//!
		//! Note that the resulting list of pseudoknots may itself have 2nd-order pseudoknots (crossing bonds).
		//! Example: currentPairs     (input)  = [ 6, 5, 7, 8, 2, 1, 3, 4 ] == { 1:6  2:5  3:7  4:8 }
		//!          normalPairs      (output) = [ 6, 5, 0, 0, 2, 1, 0, 0 ] == { 1:6  2:5           }
		//!          pseudoknotPairs  (output) = [ 0, 0, 7, 8, 0, 0, 3, 4 ] == {           3:7  4:8 }
		//!            (note that the two pairs in the example's pseudoknot group are mutually crossing)
		//! /param pseudoknotPairs  A pointer to a result vector that should be filled with base pairing information for pseudoknot pairs 
		//!             (i.e. any pair that is NOT in the optimal set of non-crossing pairs. 
		//!             This pointer can be NULL in which case it is ignored.
		//! /param normalPairs  A pointer to a result vector that should be filled with base pairing information for all pairs that are in
		//!             the largest subset of non-crossing pairs. These represent the "normal", NON-pseudoknot basepairs. 
		//!             This can be NULL in which case it is ignored.
		void FindPseudoknots(const int structurenumber, vector<int> *pseudoknotPairs=NULL, vector<int> *normalPairs=NULL) const;

		//! Fills the results output vector with the "pseudoknot rank" of each base pair in this structure.
		//! This method is sequence and energy agnostic, and it does NOT modify the structure itself.
		//! 
		//! In general, a structure contains multiple base pairs some of which may "cross" each other.  
		//! It is possible to separate all pairs into distinct groups such that in each group no pair 
		//! crosses any other. We can then rank the groups by the total number of pairs in each.
		//!
		//! Thus the "pseudoknot rank" is hereby defined as: 
		//!     0 for unpaired bases
		//!     1 for bases with pairs in the largest non-crossing pair group (i.e. the normal, "NON-pseudoknot" pairs)
		//!     2 for bases with pairs in the 2nd-largest non-crossing pair group (i.e.  first-order pseudoknots)
		//!     3 for bases with pairs in the 3nd-largest non-crossing pair group (i.e.  second-order pseudoknots)
		//!     ... and so on.
		//!  Thus any pair with rank of 1 is what we consider a "normal", NON-pseudoknot pair, while 
		//!  any pair with a rank greater than 1 indicates a "pseudoknot" pair.
		//!  Importantly, note that this is NOT the same as the "crossing count" (the number of pairs in 
		//!  the structure that cross a given pair).
		//!  Two pairs could have the same crossing count, but have different rank values.
		//!  
		//!  Example:
		//!   Consider a structure with the following base pairs:  1:11, 2:10, 3:9, 4:13, 5:12, 6:14 and 7:8
		//!           ┌──14   Pairs 1, 2, and 3 do not cross each other, but they are crossed pairs 4, 5, and 6.
		//!       ┌───│──13   Pairs 4 and 5 do not cross each other, but they are crossed by pair 6.
		//!       │ ┌─│──12   Pair 7 does not cross any other pairs.
		//!    1──│─│─│──11   So the largest group of mutually non-crossing pairs (rank = 1) is { 1, 2, 3, and 7 }.
		//!    2──│─│─│──10   The second largest group (rank = 2) is formed by pairs { 4, and 5 }.
		//!    3──│─│─│───9   The third largest group (rank = 3) is formed by pair { 6 }.
		//!    4──┘ │ │        
		//!    5────┘ │       The input pairs vector and the corresponding getPseudoknotPairRanks results vector for this 
		//!    6──────┘       example would be:       pairs   = [ (0), 11, 10, 9, 13, 12, 14, 8, 7, 3, 2, 1, 5, 4, 6 ]
		//!    7──────────8                           results = [ (0),  1,  1, 1,  2,  2,  3, 1, 1, 1, 1, 1, 2, 2, 3 ]
		//!
		//! /param pairs (intput) A reference to a vector containing base pairing information in the format used by the singlestructure class 
		//!              (i.e. if a basepair exists between i and j, then pairs[i]==j and pairs[j]==i ).
		//! /param results (output) A reference to a vector that will be filled with the pair rank of each corresponding pair in `pairs`.
		//!              (i.e. if a basepair exists between i and j in pairs, then results[i] and results[j] will both be set to the 
		//!              rank of the basepair between i and j.
		void GetPseudoknotRanks(vector<int> &results, const int structurenumber=1) const;

		//! Remove all basepairs that represent PseudoKnots (i.e. all bonds with a "crossing level" of 2 or higher).
		//! This method is sequence and energy agnostic. 
		//! \param structurenumber The index of the structure to break.
		//! \param brokenPairs A pointer to a vector of integers that will be filled with the 
		//!                    base-pairing information for bases that were broken. This vector has the same format as
		//!                    the singlestructure::basepr vector.
		//!					   If this parameter is NULL, no information about broken pairs will be returned.
		void BreakPseudoknots(int structurenumber=1, vector<int> *brokenPairs=NULL);

		//! Set the label for a structure, using a string.

		//! \param label is a string that will be strored.
		//! \param structurenumber is the index to which structure will hold the label.  This is one-indexed.
		void SetCtLabel(const string &label, const int structurenumber);
		

		//! Set the label for a structure,using a pointer to cstring.

		//! \param label is a pointer to char that provides the label.
		//! \param structurenumber is the index to which structure will hold the label.
		void SetCtLabel(const char *label, const int structurenumber);

		//! Set the energy for structure numer structurenumber.

		//! \param structurenumber is an int that gives the structure number being indexed, this is one indexed.
		//! \param energy is an int that sets the energy in kcal/mol*conversionfactor, where conversionfactor is set in /src/defines.h.
		void SetEnergy(int structurenumber, int energy);

		
		//Set the pairing partner for i in structure structurenumber.

		//! This function sets nucleotide i paired to j in structure number structurenumber.
		//! \param i is the nucleotide index of the first pairing partner, which is one-indexed.
		//! \param j is the nucleotide index of the second pairing partner.
		//! \param structurenumber is the structure number, which is one-indexed.
		void SetPair(int i, int j, int structurenumber=1);

		//! Set the label from a sequence, using a string.

		//! \param label is a string that will be stored.
		void SetSequenceLabel(const string& label);

		//! base pairing cutoff is float between .1 to 1 decide minimum threshold of base pairing
		void SetBasePairingCutoff(float i);

		//! base pairing cutoff is float between .1 to 1 decide minimum threshold of base pairing
		float GetBasePairingCutoff();

		//! check if positions i and j can base pair in ct
		bool can_pair(int i, int j, structure* ct);

		//Check if base pairing between i, j is flanked by non-pairing nucleotides (not gaps)
		bool can_pair_isolated(int i, int j, structure* ct);

		//! Set the sequence for this structure.
		//! This includes the following operations:
		//!   - Allocate space for the bases.
		//!   - Verify each base exists in the alphabet.
		//!   - Setup the arrays nucs, numseq, and hnumber.
		//!   - Check to see if any nucleotide needs to be single-stranded and call AddSingle for them.
		//! The sequence can contain whitespace, which is ignored.
		//! If an error occurs (such as the data-table has not yet been loaded) the return value 
		//! will be an RNA error code (i.e. it corresponds to one of those defined in RNA::GetErrorMessage)
		//! Additionally lastErrorDetails may be set (which can be queried with GetErrorDetails())
		//!
		//! \param sequence The nucleotide sequence for the structure.  This should contain only valid bases. No whitespace is allowed.
		int SetSequence(const string& sequence);
		int SetSequenceFast(const string& sequence);
		int SetSequenceFast(const string& sequence, const short* numseq1);
		int SetTmpSequence(const string& sequence);
		int FillNoOfGapsMatrix();
		int FillGappedtoUngappedMap();
		int SetSequenceWithoutGaps(const string& sequence);
		//********************************
		// Get and Set constraint information
		//*********************************
		

		//! Allocate a bool array with information about whether any given pair is allowed.

		//! This function must be called after reading a sequence and before any template information is read or written.
		//! This mechanism is orthogonal to the functions AddForbiddenPair, GetForbiddenPair5, and GetForbiddenPair3.  It exists for the convenience of coding functiopns that need to forbid a large number of pairs.
		//! The memory use is cleaned up in the destructor.
		void allocatetem();


		//! Add a nucleotide to the list of those that must pair.

		//! \param i is an int that indicates the nucleotide position, one indexed. 
		void AddDouble(int i);

		//! Add a pair of nucleotides to the list of those not allowed to form.

		//! \param i is an int that indicates the 5' nucleotide position, one indexed.
		//! \param j is an int that indicates the 3' nucleotide position, one indexed.
		void AddForbiddenPair(int i, int j);

		//! Add a nucleotide to the list of Us in GU pairs.

		//! Note that there is no error checking.  This nucleotide must be a U.
		//! \param i is an int that indicates the nucleotide position, one indexed.
		void AddGUPair(int i);

		//! Add a nucleotide to the list of those accessible to traditional chemical modification.

		//! These nucleotides can only be unpaired, at the end of a helix, in a GU pair, or adjacent to a GU pair.
		//! \param i is an int that indicates the nucleotide position, one indexed.
		void AddModified(int i);

		//! Add a pair of nucleotides to the list of those that must form.

		//! Note that there is no error checking.  This should be an allowed base pair.
		//! \param i is an int that indicates the 5' nucleotide position, one indexed.
		//! \param j is an int that indicates the 3' nucleotide position, one indexed.
		void AddPair(int i, int j);

		//! Add a nucleotide to the list of those not able to pair.

		//! \param i is an int that indicates the nucleotide position, one indexed. 
		void AddSingle(int i);

		//! Add a domain constraint that requires all nucleotides n, where i <= n <= j, to only pair between i and j
		void AddDomain(int i, int j);

		//!	Indicate if pairing distance is limited for structrure prediction methods.


		//! \return A bool that is true if the pairing distance has a limit.
		inline bool DistanceLimited() {

			return limitdistance;

		}

		//! Get a nucleotide that must be base paired.

		//! \return An int that gives the nucleotide position, one indexed.
		//! \param i is an int that gives the constraint number, which should be between 0 and GetNumberofDoubles()-1, inclusive.
		int GetDouble(int i);

		//! Get a nucleotide that must not be in a specific pair.

		//! \return An int that gives the nucleotide position for the 5' partner, one indexed.
		//! \param i is an int that gives the constraint number, which should be between 0 and GetNumberofForbiddenPairs()-1, inclusive.
		int GetForbiddenPair5(int i);

		//! Get a nucleotide that must not be in a specific pair.

		//! \return An int that gives the nucleotide position for the 3' partner, one indexed.
		//! \param i is an int that gives the constraint number, which should be between 0 and GetNumberofForbiddenPairs()-1, inclusive.
		int GetForbiddenPair3(int i);

		//! Get a nucleotide that must be a U in a GU pair.

		//! \return An int that gives the nucleotide position, one indexed.
		//! \param i is an int that gives the constraint number, which should be between 0 and GetNumberofGU()-1, inclusive.
		int GetGUpair(int i);

		//! Get a nucleotide that is accessible to chemical modification.

		//! These nucleotides can only be unpaired, at the end of a helix, in a GU pair, or adjacent to a GU pair.
		//! \return An int that gives the nucleotide position, one indexed.
		//! \param i is an int that gives the constraint number, which should be between 0 and GetNumberofMofidied()-1, inclusive.
		int GetModified(int i);

		//! Get the number of nucleotides forced to be double-stranded.

		//! \return An int that is the number of nucleotides constrained to be double stranded.
		int GetNumberofDoubles();

		//! Get the number of pairs that are forbidden.

		//! \return An int that is the number of forbidden pairs.
		int GetNumberofForbiddenPairs();
		
		//! Get the number of Us forced to be in GU pairs.

		//! \return An int that is the number of nucleotides constrained to be in GU pairs.
		int GetNumberofGU();
		
		//! Get the number of nucleotides that are accessible to chemical modification.

		//! \return An int that is the number of nucleotides constrained to be chemically modified.
		int GetNumberofModified();
		
		//! Get the number of nucleotides forced to be single-stranded.

		//! \return An int that is the number of nucleotides constrained to be single stranded.
		int GetNumberofSingles();
		
		//! Get the number of pairs that are constrained to occur.

		//! \return An int that is the number of forced pairs.
		int GetNumberofPairs();

		//! \For increasing speed of AlignmentFold
		//! \This is used by multi_can_pair in AlignmentFold
		//! \Fills pairing matrix based on bp_cutoff
		void FillPairingMatrix();

		//! \Creates 1-D array of length of alignment
		//! \element i is true if atleast 50% sequences contain gaps
		void FillGapMatrix();


		//! \return An int that is the number of folding domains.
		int GetNumberofDomains();
		
		//! Get a nucleotide that must be in a specific pair.

		//! \return An int that gives the nucleotide position for the 5' partner, one indexed.
		//! \param i is an int that gives the constraint number, which should be between 0 and GetNumberofPairs()-1, inclusive.
		int GetPair5(int i);

		//! Get a nucleotide that must be in a specific pair.

		//! \return An int that gives the nucleotide position for the 3' partner, one indexed.
		//! \param i is an int that gives the constraint number, which should be between 0 and GetNumberofPairs()-1, inclusive.
		int GetPair3(int i);

		//!	Provide the maximum distance between nucleotides that can pair.

		//! \return An int that is the maximum distance in sequence for nucleotides that can pair.
		inline int GetPairingDistanceLimit() {

			return maxdistance;

		}

		//! Get a nucleotide that must be single stranded.

		//! \return An int that gives the nucleotide position, one indexed.
		//! \param i is an int that gives the constraint number, which should be between 0 and GetNumberofSingles()-1, inclusive.
		int GetSingle(int i);
		
		int GetNumberofSequences();

        
        //! Remove the specified number of single-stranded constraints (pop off from end)
        void RemoveSingleStrandConstraints(int number);
        


		//! Returns a nucleotide that defines the 5' end of a domain
		int GetDomain5(int i);
		
		
		//! Returns a nucleotide that defines the 3' end of a domain
		int GetDomain3(int i);
		

		//! Reset, i.e. remove, all constraints
		void RemoveConstraints();


		//! Set a maximum pairing distance between nucleotides that can pair.  

		//! For nucleotides i and j, if |j-i|>=ct->maxdistance, the nucleotides cannot pair.
		//! This also sets a bool so that DistanceLimited() will return true.
		//! \param maxdistance is an int that will be the maximum distance.
		void SetPairingDistance(int maxdistance);


		//*********************************
		//Functions for disk I/O
		//*********************************

		//! Write a ct file to disk.
		//! \param ctoutfile is a pointer to a c-string that provides a filename.
		//! \param append is a bool that indicates if these structures should be appended to the end of the file.  
		//!            The default, false, is to overwrite any existing file.
		//! \param commentProvider a reference to an object that returns application-defined comments for each structure 
		//!		in a CT file. (e.g. "ENERGY = ..." ) 
		//!         This parameter allows a program to customize when comments are written and how they appear. 
		//!     For example many programs write "ENERGY = ..." but only when GetEnergy() != 0.
		//!     However, a program like MaxExpect might write "SCORE = ..." (even if GetEnergy() returns zero).
		//!     Similarly, Design might want to add comments like "NED = 0.02" which would require a different source of information than
		//!     GetEnergy() -- or at the very least, different formatting/precision. 
		//!         To provide custom comments, pass an in an object derived from CTCommentProvider that has overridden 
		//!     the getComment function. To disable comments, pass in CTComments::None. 
		//!     The default is CTComments::Energy, which provides legacy "ENERGY = ..." comments.
		int ctout(const char * const ctoutfile, const bool append=false, CTCommentProvider &commentProvider=CTComments::Energy) const;

		//! Write a dot-bracket file.
		//! Pseudoknots are encoded using alternate symbol pairs (parentheses, square brackets, curly brances, 
		//! angled brackets and then upper/lower case letters (e.g. 'A' opens, 'a' closes).
		//! \param filename is a const char pointer to a Null-terminated cstring that provides a filename.
		//! \param structurenumber the 1-based index of the structure to write. If this is -1 (the default) all structures are written.
		//! \param format One of the DotBracketFormat enum values that specify the format of dot-bracket 
		//!        files that contain multiple structures (i.e. structurenumber==-1).
		//!        This determines whether the structure title and sequence are written once at the top followed by 
		//!        multiple lines of brackets or the structure title and sequence are written out individually for each structure etc.
		//! \param commentProvider Allows customization of structure labels. See ctout parameter information for commentProvider.
		//! \return Returns 0 if successful. Otherwise returns an error code corresponding to those defined in RNA::GetErrorMessage(int code).
		int writedotbracket(const char * const filename, const int structurenumber=-1, 
			const DotBracketFormat format=DBN_FMT_MULTI_TITLE,
			CTCommentProvider &commentProvider=CTComments::Energy, 
			const bool append = false) const;

		//! Open a CT File.
		//! This opens a ct file and stores all the information in this instance of structure.
		//! A non-zero return indicates and error in reading the file.
		//! \param ctfile is a pointer to a Null-terminated cstring that gives the filename, including any necessary path information.
		//! \return The return value is 0 on success and non-zero on error. The error code returned corresponds to those defined in RNA::GetErrorMessage.
		int openct(const char *ctfile);

		//! This opens a dot-bracket (DBN) file and stores all the information in this instance of structure.
		//! A non-zero return indicates and error in reading the file.
	    //! Note: This can parse dot-bracket files with or without pseudoknots and automatically distinguishes 
	    //!       between  the following DBN formats (see DotBracketFormat for details): 
		//!       DBN_FMT_SINGLE_TITLE,  DBN_FMT_SIDE_TITLES, or DBN_FMT_MULTI_TITLE. 
		//!       This function cannot parse the format DBN_FMT_MULTI_TITLE_AND_SEQ because a structure class can only contain a single sequence.
		//! \param bracketFile is a pointer to a Null-terminated cstring that gives the filename, including any necessary path information.
		//! \return The return value is 0 on success and non-zero on error. The error code returned corresponds to those defined in RNA::GetErrorMessage.
		int opendbn(const char *bracketFile);

		//! Write the sequence to a file. The output format can be SEQ, FASTA, or plain text.
		//! \param seqfile is a const char pointer to a cstring that gives the filename, including any path information.
		//! \param seqFileType indicates the type of sequence file to write: 0=Plain Text (No label or delimiters etc), 1=SEQ, 2=FASTA (default)
		//! \param append if true, the file will be appended to instead of overwritten if it exists.
		//! \return An int that indicates an error state: 1 on no error and 0 on error.
		int writeseq (const char *seqfile, int seqFileType = 2, bool append = false);

		//! Open a sequence file.
		//! This function works on .seq and FASTA files as well as plain-text sequences.
		//! \param seqfile is a const char pointer to a cstring that gives the filename, including any path information.
		//! \return An int that indicates an error state: 1 on no error and 0 on error.
		int openseq (const char *seqfile);
		
		//! Open a sequence file.
		//! This function works on .seq and FASTA files as well as plain-text sequences.
		//! Importantly, this function differs from the legacy openseq in its return codes. openseq returns 1 on 
		//! success and 0 on error, while this function returns 0 on success and an RNA error code on error.
		//! 
		//! \param seqfile is a const char pointer to a cstring that gives the filename, including any path information.
		//! \return An int that indicates an error state: 0 on no error or a number indicating a more detailed error code. See RNA::GetErrorMessage(int)
		int openseqx (const char *seqfile);

		//*******************************
		//Functions that act on whole structures
		//*******************************

		//! Add another empty structure to the list of singlestructures.
		// DHM: Remember if this is structure 1 to set the label from some sequence label.!
		void AddStructure();

		//! Remove all pairs from a structure, i.e. make it a clean slate for a new set of pairs
		
		//! \param struturenumber is an index to the structure to be cleaned.
		void CleanStructure(int structurenumber);

		

		//! Remove the last structure.
		void RemoveLastStructure();
		
		//! Remove all structures.
		void RemoveAllStructures();

		//! Remove the structure at structurenumber.
		
		//! If the last structure is being removed, it is more efficient to use RemoveLastStructure();
		//! \param structurenumber is an int that is the index to which structure should be removed.  This is one indexed.
		void RemoveStructure(int structurenumber);

		//********************************
		//Functions for accessing the datatable pointer or opening the datatables
		//********************************

		//! Get the datatable pointer.

		//! This function returns the pointer to the datatable, data.
		//! \return The pointer to the underlying datatable.
		datatable *GetThermodynamicDataTable();

		//! Set the datatable pointer.

		//! This function sets the pointer to the datatable, data.
		//! \param DataTablePointer is the pointer to a datatable.
		void SetThermodynamicDataTable(datatable *DataTablePointer);

		//! Returns true if the data property has been set to a valid datatable and its alphabet has been successfully read.
		bool IsAlphabetLoaded();
		//! Returns true if the data property has been set to a valid datatable and the thermodyanamic tables have been read.
		bool IsThermoDataLoaded();

		//********************************
		//Additional functions
		//********************************

		//! Sort structures by energy.

		//! This function sorts structures in energy from lowest to highest.
		//! It is important that the structure energies be present by structure prediction or by efn2.
		void sort();


		//! Find problem in the set of structures.

		//! This function is for debugging.  It checks each pair in each structure to look for inconsistencies. 
		//! \return true when an inconstency with pairing is found and false otherwise.
		bool ProblemwithStructures();

		//! This function reads a SHAPE reactivity datafile and parse the data into single-stranded amd chemical modification constraints.
		//! This function is largely depracated by the pseudo-free energy approach.  It is still available for experimentation.
		//! \return 0 on success or an error code compatible with RNA::GetErrorMessage (for example 1=file not found, 2=could not open file)
		int ReadSHAPE(const char *filename, float SingleStrandThreshold, float ModificationThreshold);//Read SHAPE reactivity data from a file

		//void ReadSHAPE(const char *filename, bool calculate=true);//Read SHAPE reactivity data from a file
		//! This function reads a SHAPE reactivity datafile and saves the data for a linear penalty.
		//! \param calculatePseudoEnergies (default true) indicate whether these data are being read for folding.  
		//!        (false means the raw values need to be stored.)
		//! \param modifier One of the RestraintType enum values to indicate which type of restraint is to be calculated (e.g. SHAPE, DMS, CMCT, etc)
		//! \return 0 on success or an error code compatible with RNA::GetErrorMessage (for example 1=file not found, 2=could not open file)
		int ReadSHAPE(const char *filename, RestraintType modifier=RESTRAINT_SHAPE, bool calculatePseudoEnergies=true);//Read SHAPE reactivity data from a file

		//! Read offset files
		//! This function must be called ofter reading SHAPE files, if read, because the same infrastructure is used.
		//! Either filename can be NULL, in which case that offset is not recorded.
		//! \param SSOffset provides free energies to add to nucleotides if they are single-stranded
		//! \param DSOffset provides free energies to add to nucleotides if they are double-stranded
		//! \return 0 on success or an error code compatible with RNA::GetErrorMessage (for example 1=file not found, 2=could not open file)
		int ReadOffset(const char *SSOffset, const char *DSOffset);//Read Free Energy Offset Files.

		//! This function reads an experimental pair bonus file, similar to SHAPE, but just straightforward
		//! application as kcal bonuses.  As with SHAPE, bonus is applied at 0x, 1x, and 2x for
		//!  single stranded, edge base pairs, and internal base pairs.
		//! \return 0 on success or an error code compatible with RNA::GetErrorMessage (for example 1=file not found, 2=could not open file)		
		int ReadExperimentalPairBonus(const char *filename, double const experimentalOffset = 0.0, double const experimentalScaling = 1.0 );

		//! Write SHAPE and SHAPEss parameters out to file, exactly as they are currently stored.
		//! Currently this is just for debugging purposes.
		//! \return 0 on success or 2 if the file cannot be opened for writing.
		//! \param printHeaders If true, the sequence label and headers "SHAPE" and "SHAPEss" 
		//!        will be written on separate lines following a hash (#) symbol.
		int WriteSHAPE(const string &outfile, bool printHeaders = true);


		//! Deletes the SHAPE, SHAPEss, and SHAPEss_region arrays if the have been allocated.
		//! Also sets the `shaped` field to false and sets the SHAPE, SHAPEss, and SHAPEss_region 
		//! pointers to NULL. This can be called any time to remove SHAPE data.
		void DeleteSHAPE();

		//! Allocates the SHAPE and SHAPEss arrays if they haven't already been created.
		//! Also initializes SHAPEss_region, a triangular 2-D array that stores ss SHAPE energies for loops.
		void AllocateSHAPE();

		//! Return an exact copy of SHAPE data (e.g. for backup-and-restore purposes). 
		//! The return value is a pointer to a dynamically allocated array of doubles.
		//! \param includeSHAPEss Indicates whether SHAPEss data should be included. 
		//!        If true, the array will contain both SHAPE and SHAPEss values and will
		//!          have a size of 4*numbases+2 (with SHAPE data from arr[0] to arr[2N] and 
		//!          SHAPEss from arr[2N+1] to arr[4N+1]).
		//! 	   If false, the array will only contain SHAPE data and will have a size of 
		//!          2*numofbases+1.
		double* CopySHAPE(const bool includeSHAPEss);

		//! Load data from an array of doubles into the SHAPE array (e.g. for backup-and-restore purposes).
		//! \param shapeArray The array of values to load into the SHAPE (and optionally SHAPEss) arrays.
		//!    The data is copied -- this function does NOT take ownership of the passed-in array.
		//! \param includeSHAPEss Indicates whether or not SHAPEss data is included in the passed-in array. 
		//!          If true, the array must contain both SHAPE and SHAPEss values and must
		//!            have a size of at least 4*numbases+2 (with SHAPE data from arr[0] to arr[2N] 
		//!            and SHAPEss from arr[2N+1] to arr[4N+1])
		//!          If false, the array need only contain SHAPE data and must have a size of 
		//!            at least 2*numofbases+1.
		//! If the passed-in array is NULL, this structure's SHAPE data will be deleted.
		void LoadSHAPE(const double* shapeArray, const bool includeSHAPEss);

		double **constant;//constant is used to hold an array of equilibrium constants.  In partition function calculations, 
					//the equilibrium constant is multiplied by constant[j][i] when the i-j pair is formed. 
					//NOTE: The use of constant is NOT orthogonal to using chemical modification data.  They cannot
					//both be used at once.
		void allocateconstant();//Function to allocate memory for constant array.
		bool SHAPEFileRead;


		//!	Returns extended details about the last error. (e.g. error messages produced during file read operations that are otherwise lost.)
		const string& GetErrorDetails();


		//Is a nucleotide index a certain type of nucleotide?:
		//return true if i is type c or false otherwise
		// inline bool IsNuc(int i, char c) {
			
		// 	if (std::find(data->alphabet[numseq[i]].begin(), data->alphabet[numseq[i]].end(), c) != data->alphabet[numseq[i]].end())
		// 		return true;
		// 	else return false;

		// }
		inline bool IsNuc(int i, char c) 
		{
			VERIFY_NUC_INDEX(i);
			return (std::find(data->alphabet[numseq[i]].begin(), data->alphabet[numseq[i]].end(), c) != data->alphabet[numseq[i]].end());
		}
		inline bool IsNuc(int i, char c, datatable *data) {
			VERIFY_NUC_INDEX(i);
			return (std::find(data->alphabet[numseq[i]].begin(), data->alphabet[numseq[i]].end(), c) != data->alphabet[numseq[i]].end());
		}

		void RemoveEnergyLabels(const char* customLabel=NULL);

		//! Generates constraint matrix that can be passed to partition-cuda.
		//! Currently only supports forbidden pairs and specific forced base pairs.
		//! Constraints that force a base to be double stranded, but without a 
		//!    specific pairing partner are ignored.
		int *generate_constraint_matrix();
		
		//! Sets the default behavior regarding non-critical warnings. 
		//! 1 (ON)  -- Write warnings to stdout.
		//! 2 (ERR) -- Write warnings to STDERR for better error detection in scripts etc..
		//! 0 (OFF) -- Suppress warnings (not recommended unless the user is aware of 
		//!            problematic data and does not need the clutter of warning output.
		//! The default value of this variable is taken from the environment variable
		//! RNA_WARNINGS which should be set to "OFF", "ON", or "ERR" or 
		//! the numeric equivalents "0", "1", or "2". 
		static int ShowWarnings;

		//! SumShapeRepeats affects how multiple datapoints for the same nucleboase are handled 
		//! in SHAPE data. True (the default) causes ReadSHAPE to add multiple datapoints for use
		//! in the "resample with replacement" technique, but setting the environment variable 
		//! AVG_SHAPE_REPEATS to 1 causes multiple values to be averaged instead.
		static bool SumShapeRepeats; // default true.

		//static const double kT;


#ifdef SWIG
		// this hides the following definitions from SWIG when generating proxy-code for Java etc.
		private:
#endif // SWIG

		//! Returns a reference to an output stream appropriate for warnings (according to 
		//! the structure::ShowWarnings setting.)
		ostream& cwarn(); 

		string sequencelabel;//a label that was read from disk along with a sequence

		
		short int* numseq;
		int *hnumber;
		
		int inter[3],allocatedstructures;
		char *nucs;
		bool intermolecular,allocated,templated,stacking;
		bool **tem;//tem stores template information as to whether a pair is allowed

		void allocate(int size = maxbases);
		void deallocate();
		void allocatestructure(int structures);
		
		// Stores base pairing information calculated using ct->bp_cutoff for AlignmmentFold
		vector<vector<bool>> pairing_matrix;
		
		// Stores information about the presence of gaps at each position of the alignment.
		vector<bool> gap_matrix;

		// Stores no. of gaps between positions i and j in a given single sequence.
		vector<vector<int> > no_of_gaps_matrix;

		// maps position of bases in gapped sequence to ungapped sequence. 
		vector<int> gapped_to_ungapped_map;

		// to temporarily store nucs, numseq etc for multi_erg3 
		structure* tmp_ct_hp;

		// to temporarily store nucs, numseq etc for multi_erg2 
		structure* tmp_ct_i;

		// string to store loop without gaps
		string loop_str;
		

		short int min_gu, min_g_or_u;//NMR-derived constraint variables
		short int neighbors[maxforce][maxneighborlength],nneighbors;//also NMR-derived index this from zero in both dimensions
		//regional NMR constraints:
		short int nregion,rmin_gu[maxregions],rmin_g_or_u[maxregions];
		short int rneighbors[maxregions][maxforce][maxneighborlength],rnneighbors[maxregions],start[maxregions],stop[maxregions];
		//microarray type constraints:
		short int nmicroarray,microstart[maxregions],microstop[maxregions],microunpair[maxregions];
		bool *fcedbl;//pointer to a 2-D array used in Dynalign to track nucleotides that must be double-stranded
		
		
		double *SHAPE;//double array to contain SHAPE data -- values less than -500 are ignored
		double **EX;// double array that contains experimental bonuses/penalties
		bool shaped;//keeps track of whether SHAPE data was loaded
		bool experimentalPairBonusExists;//keeps track of whether experimental bonus data was loaded
		bool ssoffset;//keeps track of whether a single stranded offset was read from disk
		double SHAPEslope,SHAPEintercept;//values of slope and intercept for SHAPE data modification of pairing stability
		//SINGLE STRANDED SHAPE ENERGY VARIABLES AND FUNCTIONS
		double *SHAPEss; //short int array that contains SHAPE data for single-stranded segments
		double SHAPEslope_ss, SHAPEintercept_ss; //values of the slope and intercept for SHAPE data modifying single stranded loop stability
		short int **SHAPEss_region;  //2-d short int array containing energy values for hairpin loop combinations
		int SHAPEss_calc(int index_i, int index_j);  //Returns pseudoenergy term for a hairpin loop using single stranded SHAPE data
		short int SHAPEss_give_value(int index);  //Returns the single stranded SHAPE pseudo energy for a given nucleotide
		double CalculatePseudoEnergy(const double data, const RestraintType modifier, const double, const double, const int ntcode, const bool);
		double Gammadist(const double data, const double shape, const double loc, const double scale);
		double Potential(const double data, const std::vector< std::vector<double> > &params, const double kT, const int ntcode = 1);
        
		void ReadProbabilisticPotentialParams();//Read chemical modifier distributions from file
		
        //Parameters for distributions
		std::vector< std::vector<double> > SHAPE_params;
		std::vector< std::vector<double> > DMS_params;
		std::vector< std::vector<double> > DMS_paramsnt;
        std::vector< std::vector<double> > CMCT_params;



		bool distsread;//keep track if the distribution files have been read from disk.

		// energyLabelWriter is a pointer to a function that creates an energy label to be inserted into structure 
		// titles when writing CT and dot-bracket files
		const char* const (*energyLabelWriter)(const int structurenumber);

		//! store sequences from RNA multiple sequence alignment into vector of structures
		vector<structure> sequences_from_alignment;

		//! store indices of nucleotides sequence alignment without gaps into vector of structures
		// Store in structure class, a 2D vector that encode the alignment :
		// rows would be sequences, columns would be an index in the sequences
		// e.g.:
		// 1 2 3 4 5 0 6 7 8
		// 0 1 2 0 0 3 4 5 6
		// Then, a hairpin(or bulge or internal loop) closed by column i and j could refer to this table to know the sequence position.
		// Fill this table after reading the FASTA alignment.
		// vector<vector<int>> sequences_index;

		//! stores sequence without gaps from multiple sequence alignment
		structure* sequence_without_gaps;



		//! base pairing cutoff, default .5
		float bp_cutoff;

		//! stores number of sequences in sequences_from_alignment vector, default 1
		int number_of_sequences;

		//! returns pointer to individual sequence i in vector<structure> 
		structure* get_individual_sequence(int i);

		//! opens alignment file in FASTA format and populates vector<structure> 
		int open_alignment(const char* alignment_file);

	private:
		
		//!	Set extended details about the last error. (e.g. error messages produced during file read operations that are otherwise lost.)
		void SetErrorDetails(const string &details); // , bool outputToCErr = true

		//! Fills SHAPEss_region, a 2-d array with pseudo energy terms for loops from i-j using 
		//! single-stranded (ss) SHAPE parameters or offsets.
		//! The SHAPEss and SHAPEss_region arrays must already be allocated (by calling AllocateSHAPE)
		//! and SHAPEss must already contain ss SHAPE pseudo-energies.
		//! This is called from e.g. ReadSHAPE and ReadOffset.
		void FillSHAPEssRegions();

		int numofbases;//number of nucleotides in sequence
		bool limitdistance;//toggle to indicate that there is a limit on the maximum distance between nucs in base pairs
		int maxdistance;//maximum distance between nucs in base pairs
		
		vector<singlestructure> arrayofstructures;//This holds an array of structures, i.e. base pairing information and comments
			
		//variables for holding folding constraints:
		vector<int> doublestranded; //nucleotides that must be double stranded
		vector<int> singlestranded; //nucleotides that must be single stranded
		vector<int> GUpair; //Us in GU pairs
		vector<int>	modified; //nucleotides accessible to tradictional chemical modification agents
		vector<int> pair5; //5' partner in forced pair
		vector<int> pair3; //3' partner in forced pair
		vector<int> forbid5; //5' partner in a forbidden pair
		vector<int> forbid3; //3' partner in a forbidden pair
		vector<int> domains5; //domain definitions
		vector<int> domains3; //domain definitions

		string lastErrorDetails;
		datatable *data;//store the thermodynamic data, and access the alphabet of nucleotides

		
};


//char *tobase (int i);//convert a numeric value for a base to the familiar
								//character


//void tonum(char *base,structure *ct,int count); //converts base to a numeric

#ifndef SWIG // The following should be excluded from SWIG code generation

int ecompare(const void *i, const void *j);

integersize ergcoaxflushbases(int i, int j, int ip, int jp, datatable *data);

//this function calculates flush coaxial stacking
//it requires sequence in i,j,ip, and jp
integersize ergcoaxinterbases1(int i, int j, int ip, int jp, int k, int l, datatable *data); 


//this funtion calculates an intervening mismatch coaxial stack
integersize ergcoaxinterbases2(int i, int j, int ip, int jp, int k, int l, datatable *data); 

//this funtion calculates an intervening mismatch coaxial stack

integersize ergcoaxflushbases(int i, int j, int ip, int jp, structure *ct, datatable *data);
integersize multi_ergcoaxflushbases(int i, int j, int ip, int jp, structure *ct, datatable *data);

integersize ergcoaxinterbases1(int i, int j, int ip, int jp, structure *ct, datatable *data);
integersize multi_ergcoaxinterbases1(int i, int j, int ip, int jp, structure *ct, datatable *data);

integersize ergcoaxinterbases2(int i, int j, int ip, int jp, structure *ct, datatable *data);
integersize multi_ergcoaxinterbases2(int i, int j, int ip, int jp, structure *ct, datatable *data);

int decon1(int x,int alphabetsize);//used by ergmulti to find a nucleotide from a base pair
int decon2(int x, int alphabetsize);//used by ergmulti to find a nucleotide from a base pair
integersize ergmulti(int st, int ip, structure *ct, datatable *data, bool simplemb);
//calculate the multi branch loop free energy for a loop starting at nuc ip
//	in structure number st of ct
integersize ergexterior(int st, structure *ct, datatable *data, int start=1, int stop=0);
//calculate the exterior loop free energy in structure number ip

// Calculate covariation energy term for column i and j in multiple sequence alignment
// Bij = Cij - phi_1 * qij where phi_1 = 1
integersize multi_erg_covar(int i, int j, structure* ct);

//calculates energy of stacked base pairs
integersize erg1(int i,int j,int ip,int jp,structure *ct,datatable *data);		
integersize multi_erg1(int i, int j, int ip, int jp, structure* ct, datatable* data);

//calculates energy of a bulge/internal loop
integersize erg2(int i,int j,int ip,int jp,structure *ct,datatable *data,char a, char b);		
integersize multi_erg2(int i, int j, int ip, int jp, structure* ct, datatable* data, char a, char b);

// calculates the energy of an interior part of an internal loop (includes asymmetry)
integersize erg2in(int i,int j,int ip,int jp,structure *ct, datatable *data,char a, char b);
//calculates the energy of an exterior part of an internal loop (only has length and terminal stack components)		
integersize erg2ex(int i,int j,int size,structure *ct, datatable *data);
	

//calculate the energy of a hairpin loop:
integersize erg3(int i,int j,structure *ct,datatable *data,char dbl);
integersize multi_erg3(int i, int j, structure* ct, datatable* data, char dbl);

// structure* remove_gaps(int i, int j, structure* ct);
structure* remove_gaps(int i, int j, structure* ct_i, structure* tmp_ct);

//structure* remove_gaps(int i, int j, int ip, int jp, structure* ct);
structure* remove_gaps(int i, int j, int ip, int jp, structure* ct_i, structure* temp_ct_2);

int no_of_gaps(int i, int j, structure* ct);

//calculate the energy of a dangling end:
integersize erg4(int i, int j, int ip, int jp, structure *ct, datatable *data, bool lfce);
integersize multi_erg4(int i, int j, int ip, int jp, structure* ct, datatable* data, bool lfce);

//calculates energy of a dangling base
//erg4 without parameter usage counting
integersize erg4_nc(int i, int j, int ip, int jp, structure *ct, datatable *data, bool lfce);

inline integersize SHAPEend(int i, structure *ct);//calculate the SHAPE pseudo energy for a single
									//paired nucleotide

//this function calculates whether a terminal pair i,j requires the end penalty
inline integersize penalty(int i,int j,structure* ct, datatable *data) {
	integersize energy;
	VERIFY_NUC_INDEX_CT(i,ct);
	VERIFY_NUC_INDEX_CT(j,ct);
	if (data->AUappliestoGU&&(ct->IsNuc(i,'U')||ct->IsNuc(j,'U'))) energy= data->auend;
	else if (!data->AUappliestoGU&& (ct->IsNuc(i, 'A') || ct->IsNuc(j, 'A'))) energy = data->auend;
	else energy= 0;//no end penalty
	return (energy/*+SHAPEend(i,ct)+SHAPEend(j,ct)*/);
}

//this function calculates whether a terminal pair i, j requires the end penalty for multiple sequence alignment
inline integersize multi_penalty(int i, int j, structure* ct, datatable* data)
{
	integersize energy = 0;
	int N = ct->number_of_sequences;
	for (int m = 0; m < N; m++)
	{	
		structure* seq_m = ct->get_individual_sequence(m);
		VERIFY_NUC_INDEX_CT(i, seq_m);
		VERIFY_NUC_INDEX_CT(j, seq_m);
		if (seq_m->IsNuc(i, 'U') || seq_m->IsNuc(j, 'U'))
		{
				energy = energy + seq_m->GetThermodynamicDataTable()->auend;
		}
	}
	return energy;
}


inline integersize penalty_nc(int i, int j, structure* ct, datatable *data) 
{
	integersize energy;

	#ifdef COUNTING
	int tmp_count = data->auend.get;
	#endif

	energy = penalty(i, j, ct, data);

	#ifdef COUNTING	
	data->auend.get = tmp_count;
	#endif

	return energy;
}

//returns the free energy of coaxial stacking of i-j onto ip-jp
integersize ergcoax(int i, int j, int ip, int jp, int k, structure *ct, datatable *data);
	
//When adding the included and excluded fragments for a structure with pair i-j, the SHAPE free energy needs correction
//so that it is not counted twice.

//Write the structure class-specific items in a save file
void writestructuresave(ofstream *out, structure *ct);

//Read the structure class-specific items in a save file
void openstructuresave(ifstream *out, structure *ct);

//! If the label contains "ENERGY = <NUMBER>" at the start of the string, it will be removed,
//! leaving the remaining text intact. Leading whitespace is also removed.
//! A custom label can be specified instead of "ENERGY", in which case, the full text to be removed 
//! is "<customLabel> = <NUMBER>".
//! /param text A reference to the string that contains the structure label. It will be modified in-place.
//! /param customLabel A c-string that specifies what label to search for (and remove) from the start of the text. 
//!            The default is "ENERGY". 
void eraseEnergyLabel(string &text, const char*const customLabel="ENERGY");

#endif // Place functions that should be available to client languages (e.g. Java, python etc) BELOW this line.

//! Determine if the structure has one or more Pseudoknots (crossing bonds).
bool hasPseudoknots(const vector<int> &pairs);

//! Given a vector containing base-pairing information (currentPairs), this function determines 
//! which pairs are "pseudoknots" by finding the the largest subset of non-crossing bonds.
//! The caller can retrieve the list of "normal" bonds (normalPairs) or the list of pseudoknots (pseudoknotPairs)
//! or both.
//! This function does NOT modify the input vector (but it IS safe to pass a pointer to the input vector in as 
//!   as either pseudoknotPairs or normalPairs in which case it will be written to, thus altering the structure.)
//! Note that the list of pseudoknots may itself have 2nd-order pseudoknots (i.e. crossing bonds).
//! Example: currentPairs     (input)  = [ 6, 5, 7, 8, 2, 1, 3, 4 ] == { 1:6  2:5  3:7  4:8 }
//!          normalPairs      (output) = [ 6, 5, 0, 0, 2, 1, 0, 0 ] == { 1:6  2:5           }
//!          pseudoknotPairs  (output) = [ 0, 0, 7, 8, 0, 0, 3, 4 ] == {           3:7  4:8 }
//!            (note that the two pairs in the example's pseudoknot group are mutually crossing)
//! /param currentPairs A reference to a vector containing base pairing information in the format used by the singlestructure class 
//!             (i.e. if a basepair exists between i and j, then pairs[i]==j and pairs[j]==i ).
//! /param pseudoknotPairs  A pointer to a result vector that should be filled with base pairing information for pseudoknot pairs 
//!             (i.e. any pair that is NOT in the optimal set of non-crossing pairs. 
//!             This pointer can be NULL in which case it is ignored.
//! /param normalPairs  A pointer to a result vector that should be filled with base pairing information for all pairs that are in
//!             the largest subset of non-crossing pairs. These represent the "normal", NON-pseudoknot basepairs. 
//!             This can be NULL in which case it is ignored.
void findPseudoknots(const vector<int> &currentPairs, vector<int> *pseudoknotPairs = NULL, vector<int> *normalPairs = NULL);

#endif //STRUCTURE_H

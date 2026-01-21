
//Use the precompiler to make sure the class definition is not included more than once.
#if !defined(DYNALIGNOBJECT_H)
#define DYNALIGNOBJECT_H


#include "TwoRNA.h"
#include "../src/dynalign.h"
#include "../src/dynalignheap.h"

//! Dynalign_object Class.
/*!
	The Dynalign_object class provides an entry point for the Dynalign algorithm. 
	The class is inherited from the TwoRNA class, which itself contains two instances to the class RNA.
*/


//Note the stylized comments provide facility for automatic documentation via doxygen.

class Dynalign_object:  public TwoRNA {


	public:

		//******************************************************
		//Constructors
		//******************************************************

		//!Constructor 
		//!This is a default constructor that calls the TwoRNA default constructor.
		Dynalign_object();



		//!This constuctor is available to reach the TwoRNA constructors, from which this class is inherited.
		//!This constructor uses two cstrings to provide sequences.  IsRNA is true for RNA folding and flase for DNA.
		//! Constructor 

		//! This constuctor is available to reach the TwoRNA constructors, from which this class is inherited.
		//! This constructor uses two cstrings to provide sequences.  IsRNA is true for RNA folding and false for DNA.
				//!	Input sequences should contain A,C,G,T,U,a,c,g,t,u,x,X.
		//!	Capitalization makes no difference.
		//!	T=t=u=U.  If IsRNA is true, the backbone is RNA, so U is assumed.  If IsRNA is false, the backbone is DNA, so T is assumed.
		//!	x=X= nucleotide that neither stacks nor pairs.
		//!	For now, any unknown nuc is considered 'X'.
		//! Both sequences are passed to underlying RNA classes for each sequence.
		//!	\param sequence1 is a NULL terminated c string for sequence 1.
		//!	\param sequence2 is a NULL terminated c string for sequence 2.
		//!	\param IsRNA is a bool that indicates whether these sequences are RNA or DNA.  true=RNA.  false=DNA.  Default is true.  Both sequences must have the same backbone.
		Dynalign_object(const char sequence1[], const char sequence2[], const bool IsRNA=true); 

		//!Constructor

		//!This constuctor is available to reach the TwoRNA constructors, from which this class is inherited.
		//!This constructor uses two ctsirngs as filenames, accompanied by integers to set the file type.  IsRNA is true for RNA folding and false for DNA.  
		//!	The existing files, specified by filenames, can either be a ct file, a sequence, or an RNAstructure save file. 
		//!	Therefore, the user provides a flag for the file type: 
		//!		type = 1 => .ct file, type = 2 => .seq file, type = 3 => partition function save (.pfs) file, type = 4 => folding save file (.sav).
		//! The file opening is performed by the constructors for the RNA classes that underlie each sequence.
		//!	This constructor generates internal error codes that can be accessed by GetErrorCode() after the constructor is called.  0 = no error.
		//! The errorcode can be resolved to a c string using GetErrorMessage.		
		//!	Note that the contructor needs to be explicitly told, via IsRNA, what the backbone is because files do not store this information.
		//! Note also that save files explicitly store the thermodynamic parameters, therefore changing the backbone type as compaared to the original calculation will not change structure predictions.
		//! \param filename1 is a null terminated c string and refers to sequence 1.
		//! \param filename2 is a null terminated c string and refers to sequence 2.
		//! \param type1 is an integer that indicates the file type for sequence 1.
		//! \param type2 is an integer that indicates the file type for sequence 2.
		//!	\param IsRNA is a bool that indicates whether these sequences are RNA or DNA.  true=RNA.  false=DNA.  Default is true.  Only one backbone is allowed for both sequences.
		Dynalign_object(const char filename1[], const RNAInputType type1, const char filename2[], const RNAInputType type2, const bool IsRNA=true);


		//!Constructor
		//!	Constructor that copies thermodynamic parameter tables from an existing Thermodynamics (or RNA etc) instance.
		//! The file opening is performed by the constructors for the RNA classes that underlie each sequence.
		//!	This constructor generates internal error codes that can be accessed by GetErrorCode() after the constructor is called.  0 = no error.
		//! The errorcode can be resolved to a c string using GetErrorMessage.		
		//! \param filename1 is a null terminated c string and refers to sequence 1.
		//! \param filename2 is a null terminated c string and refers to sequence 2.
		//! \param type1 is an integer that indicates the file type for sequence 1.
		//! \param type2 is an integer that indicates the file type for sequence 2.
		//! \param thermo1 is a pointer to the Thermodynamics object to copy. The internal datatables object is stored by reference, 
		//!        so subsequence changes made to it WILL affect this object.
		Dynalign_object(const char filename1[], const RNAInputType type1, const char filename2[], const RNAInputType type2, const Thermodynamics *thermo1);


		//!Constructor

		//!This constructor allows the user to read a dynalign save file (.dsv) to get base pairing information.
		//!	This constructor generates internal error codes that can be accessed by GetErrorCode() after the constructor is called.  0 = no error.
		//! The errorcode can be resolved to a c string using GetErrorMessage.
		//! \param filename[] is a cstring that indicates the filename of a .dsv file.
		Dynalign_object(const char filename[]);


		//!Constructor
		//!This constructor is used to perform Dynaligh refolding.
		//!This does not allow any changes in constraints, but does allow the creation of different set of suboptimal structures.
		//!	This constructor generates internal error codes that can be accessed by GetErrorCode() after the constructor is called.  0 = no error.
		//! The errorcode can be resolved to a c string using GetErrorMessage.
		//!\param filename is the name of a Dynalign save file name (.dsv).
		//!\param maxtrace is the maximum number of common structures to be determined.  The recommended default is 20.
		//!\param bpwin the the base pair window parameter, where 0 allows the structures to have similar pairs and larger windows make the structures more diverse.  The recommended default is 5.
		//!	\param awin is the alignment window parameter, where 0 allows the alignments to be similar and larger values make the alignments more diverse.  The recommended default is 1.
		//!	\param percent is the maximum percent difference in total folding free energy change above the lowest for suboptimal common structures.  The recommended default is 20.
		Dynalign_object(const char* filename, const short maxtrace, const short bpwin, const short awin, const short percent);

		//******************************************************
		//Functions that predict common secondary structure 
		//******************************************************

		//! Predict the lowest free energy structure common to two sequences and suboptimal solutions with the Dynalign algorithm.

		//!	In case of error, the function returns a non-zero that can be parsed by GetErrorMessage() or GetErrorMessageString().
		//!	\param maxtrace is the maximum number of common structures to be determined.  The defaults is 20.
		//!	\param bpwin the the base pair window parameter, where 0 allows the structures to have similar pairs and larger windows make the structures more diverse. The default is 5.
		//!	\param awin is the alignment window parameter, where 0 allows the alignments to be similar and larger values make the alignments more diverse.  The default is 1.
		//!	\param percent is the maximum percent difference in total folding free energy change above the lowest for suboptimal common structures.  The defaults is 20.
		//!	\param imaxseparation is the maximum separation between aligned nucleotides.  Values >= 0 are the traditional parameter, those below zero trigger the HMM alignment method, which is now prefered.
		//!	\param gap is the cost of adding gap nucleotides in the alignment in kcal/mol.
		//! \param singleinsert is whether single basepair inserts are allowed in one sequence vs the other.
		//! \param savefile is c-string with the name of a dynalign savefile (*.dsv) to be created.
		//! \param optimalonly can be used to turn on a calculation of only the energy (when true) and not the structures.
		//! \param singlefold_subopt_percent is the maximum % difference of folding energy above the lowest free energy structure for pairs in single sequence folding that will be allowed in the dynalign calculation.
		//! \param local is whether Dynalign is being run in local (true) or global mode (false).
		//! \param numProcessors is the number of processors to use for the calculation.  This requires a compilation for SMP.
		//! \param maxpairs is under development for multiple sequence folding.  Use -1 (default) for now.
		//! \return An int that indicates an error code (0 = no error, non-zero = error occurred).
#ifdef DYNALIGN_II
	int Dynalign(const short int maxtrace=20, 
                     const short int bpwin=5, const short int awin=1, const short int percent=20, const short int imaxseparation=-99, const float slope = 0.1, const float intercept = 0.5, const float gap=0.4, const int max_elongation = 5, 
			const char savefile[]=NULL, const bool optimalonly=false, const short int singlefold_subopt_percent=30, const bool local=false, 
			const short int numProcessors=1, const int maxpairs=-1);
#else
		int Dynalign(const short int maxtrace=20, 
			const short int bpwin=5, const short int awin=1, const short int percent=20, const short int imaxseparation=-99, const float gap=0.4, const bool singleinsert=true, 
			const char savefile[]=NULL, const bool optimalonly=false, const short int singlefold_subopt_percent=30, const bool local=false, 
			const short int numProcessors=1, const int maxpairs=-1);
#endif

		//******************************************************
		//Functions that write output to disk
		//******************************************************

		//! Write the alignment to disk.

		//! This function should be called after loading a dynalign save file or after a dynalign calculation has been performed.
		//! This function generates no error flag.  Nothing can go wrong... 
		//!	\param filename is the file to which the alignment should be written.
		void WriteAlignment(const char filename[]);


		//******************************************************
		//Functions that specify or report alignment restraints 
		//******************************************************
		
		//!Force an alignment during a Dynalign calculation).

		//!Nucleotide i from sequence 1 will be aligned to nucleotide k in sequence 2 in subsequent Dynalign calculation.  
		//!The function returns 0 with no error and a non-zero otherwise that can be parsed by GetErrorMessage() or GetErrorMessageString().
		//!\param i is the index of nucleotide from sequence 1.
		//!\param k is the index of nucleotide from sequence 2.
		//!\return An integer that indicates an error code (0 = no error, 100 = nucleotide i out of range, 101 = nucleotide k out of range). 
		int ForceAlignment(const int i, const int k);

		//!Get an alignment constraint.
		
		//!\param i is the nucleotide number.
		//!\param seq is the sequence (1 or 2) from which i is derived.
		//!\return An integer that indicates the nucleotide to which i is forced to be aligned, where 0 indicates no alignment.
		int GetForcedAlignment(const int i, const int seq);

		//! Read alignment constraints from disk.
		
		//!The file format is:
		//!	i1 k1
		//! i2 k2
		//! -1 -1
		//!Where each line gives a aligned pair (i from sequence 1 and k from sequence 2).  The file terminates with -1 -1 to indicate the file end.
		//!The function returns 0 with no error and a non-zero otherwise that can be parsed by GetErrorMessage() or GetErrorMessageString().
		//!\param filename is a c string that is the file name to be read.
		//!\return An integer that indicates an error code (0 = no error, 102 = file not found, 103 = error reading constraint file). 
		int ReadAlignmentConstraints(const char filename[]);



		//!Read a ct file to determine what pairs will be allowed for sequence 1 in a subsequent dynalign calculation.

		//!This results in all pairs but those in the ct being disallowed.
		//!\param ctfilename is the name of the ct file to be read to provide the template.
		//!\return An integer that indicates an error code (0=no error, 104=file not found, 105=template is already specified)
		int Templatefromct(const char ctfilename[]);


		
		//Read a dynalign save file to determine what pairs will be allowed for sequence 1 in a subsequent dynalign calculation.

		//!This reads a dsv file and only allows pairs with folding free energy change between the lowest and lowest + maxdsvchange
		//!in a subsequent dynalign calculation.
		//!\param dsvfilename is the name of the ct file to be read to provide the template.
		//!\param maxdsvchange in a float that gives a percent difference in free energy above the lowest free energy change.
		//!\return An integer that indicates an error code (0=no error, 106=file not found, 105=template is already specified)
		int Templatefromdsv(const char dsvfilename[], const float maxdsvchange);

		


		//********************************************************************************
		//Functions that report results of folding.
		//********************************************************************************
		
		//!Report the best energy for pair i-j from sequence number sequence (1 or 2).

		//!This function reports the lowest ffolding free energy for any pairs between i-j in sequence number sequence (1 or 2).
		//!This requires a search over all possible pairs in the second sequence.
		//!NOTE:  This function ONLY works after reading a Dynalign save file (.dsv) using the constructor.
		//!This is because the Dynalign energies are not normally stored after calling Dynalign.
		//!This function generates internal error codes that can be accessed by GetErrorCode() after the constructor is called.  0 = no error, 107 = Data not available, 108 = nucleotide out of range.
		//!The errorcode can be resolved to a c string using GetErrorMessage.		
		//!\param sequence is an integer indicating the sequence # (must be 1 or 2).
		//!\param i is the 5' nucleotide in a pair.
		//!\param j is the 3' nucleotide in a pair.
		//!\return A double that gives an energy in kcal/mol.
		double GetBestPairEnergy(const int sequence,const int i, const int j);


		//!Report the lowest total free energy change from a Dynalign calculation.

		//!NOTE:  This function ONLY works after reading a Dynalign save file (.dsv) using the constructor.
		//!This is because the Dynalign energies are not normally stored after calling Dynalign.
		//!This function generates internal error codes that can be accessed by GetErrorCode() after the constructor is called.  0 = no error, 107 = Data not available, 108 = nucleotide out of range.
		//!The errorcode can be resolved to a c string using GetErrorMessage.
		//!\return a double that gives an energy in kcal/mol.
		double GetLowestEnergy();
	

		//******************************************************
		//Functions to return error information
		//******************************************************

		//Return an error code, where a return of zero is no error.

		//	This function returns and error flag that is generated during construction by RNA(const char &filename, const int type, const bool IsRNA=true) or from CalculateFreeEnergy().
		//		An error of zero is always no error.  Other codes are errors and a c-string can be fetched for the error with GetErrorMessage().
		//\return An integer that provides the error code.
		//int GetErrorCode();

		//!	Return error messages based on code from GetErrorCode and other error codes.		

		//!		0 = no error
		//!		100-999 = Error associated with Dynalign, to be handled here.
		//!		>=1000 = Errors for underlying sequence, get message from TwoRNA base class.
		//!	Current errors handled here are:
		//!		100 "Nucleotide from sequence 1 is out of range.\n";
		//!		101 "Nucleotide from sequence 2 is out of range.\n";
		//!		102 "Alignment constraint file not found.\n";
		//!		103 "Error reading alignment constraint file.\n";
		//!		104 "CT file not found.\n";
		//!		105	"A template has already been specified; only one is allowed.\n";
		//!		106	"DSV file not found.\n";
		//!		107 "Data not available to calculate energy.\n"
		//!		108 "Nucleotide out of range.\n";
		//!		109 "Value of maxpairs is too large to be achievable.\n"
		//!		110 "Error reading thermodynamic parameters."
		//!\param error is the integer error code provided by GetErrorCode() or by a call to a function that returns an error.
		//!\return A pointer to a c string that provides an error message or from other functions that return integer error codes.
		const char* GetErrorMessage(const int error);

		//********************************************************************************
		//Functions that specify or report constraints on folding.
		//These constraints only affect subsequent calls to structure prediction routines.
		//********************************************************************************

		//******************************************************************
		//Functions that provide a connection to TProgressDialog, for following calculation progress:
		//******************************************************************
		
		//!Provide a TProgressDialog for following calculation progress.
		//!A TProgressDialog class has a public function void update(int percent) that indicates the progress of a long calculation.
		//!\param Progress is a TProgressDialog class.
		void SetProgress(ProgressHandler& Progress);


		//!Provide a means to stop using a TProgressDialog.
		//!StopProgress tells the RNA class to no longer follow progress.  This should be called if the TProgressDialog is deleted, so that this class does not make reference to it.
		void StopProgress();

		//******************************************************
		//Destructor
		//******************************************************
		~Dynalign_object();	

	private:
		
		//Common initializations to both constructors.
		void CommonConstructor();

		//Allocate space in the forcealign array for storing alignment constraints.
		void AllocateForceAlign();

		//store either a ct or dsv filename to use as a template
		void storetemplatefilename(const char *name);

		//variables for the backend
		short **align;//keep track of the alignment
		short **forcealign;//keep track of any forced alignments
		bool dsv_templated, ct_templated;
		char *templatefilename;
		float MAXDSV;
		int modificationflag;


		//arrays for storing Dynalign information
		dynalignarray *w,*vmod;
		varray *v;
		wendarray *w5,*w3;
		short *lowend,*highend;
		datatable *data;
		short gap;
		short lowest;

		int Maxtrace;//The maximum # of tracebacks in a Dynalign calculation.

		bool savefileread;//a bool to indicate that a savefile was read.

		double ***array;//a 3-D array to hold to store pair energies
	


};

#endif //DYNALIGNOBJECT_H defined

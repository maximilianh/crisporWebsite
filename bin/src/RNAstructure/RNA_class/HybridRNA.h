
#if !defined(HYBRIDRNA_H)
#define HYBRIDRNA_H


#include "TwoRNA.h"

//! HybridRNA Class.
/*!
	The HybridRNA class provides an entry point for all the bimolecular structure prediction routines of RNAstructure. 
	The class is inherited from the RNA class and contains an instance of TwoRNA, which itself contains two instances to the class RNA.
*/

//Use the precompiler to make sure the class definition is not included more than once.

//Note the stylized comments provide facility for automatic documentation via doxygen.

class HybridRNA:  public RNA {


	public:

		//******************************************************
		//Constructors
		//******************************************************

		//!Constructors

		//!The constuctors are available to reach the TwoRNA constructors, from which this class is built.
		//!Please see the TwoRNA constructor entries for documentation.
		
		//!	TwoRNA Constructor that copies thermodynamic parameter tables from an existing Thermodynamics (or RNA etc) instance.
		//! The file opening is performed by the constructors for the RNA classes that underlie each sequence.
		//!	This constructor generates internal error codes that can be accessed by GetErrorCode() after the constructor is called.  0 = no error.
		//! The errorcode can be resolved to a c string using GetErrorMessage.		
		//! \param filename1 is a null terminated c string and refers to sequence 1.
		//! \param filename2 is a null terminated c string and refers to sequence 2.
		//! \param type1 is an integer that indicates the file type for sequence 1.
		//! \param type2 is an integer that indicates the file type for sequence 2.
		//! \param copyThermo is a pointer to the Thermodynamics object to copy. The internal datatables object is stored by reference, 
		//!        so subsequent changes made to it WILL affect this object.
		HybridRNA(const char filename1[], const int type1, const char filename2[], const int type2, Thermodynamics* copyThermo);
		
		//! Constructor - user provides a sequences as c strings.

		//!	Input sequences should contain A,C,G,T,U,a,c,g,t,u,x,X.
		//!	Capitalization makes no difference.
		//!	T=t=u=U.  If IsRNA is true, the backbone is RNA, so U is assumed.  If IsRNA is false, the backbone is DNA, so T is assumed.
		//!	x=X= nucleotide that neither stacks nor pairs.
		//!	For now, any unknown nuc is considered 'X'.
		//! Both sequences are passed to underlying RNA classes for each sequence.
		//!	\param sequence1 is a NULL terminated c string for sequence 1.
		//!	\param sequence2 is a NULL terminated c string for sequence 2.
		//!	\param IsRNA is a bool that indicates whether these sequences are RNA or DNA.  true=RNA.  false=DNA.  Default is true.  Both sequences must have the same backbone.
		HybridRNA(const char sequence1[], const char sequence2[], const bool IsRNA=true); 
		

		//! Constructor - user provides a filenames for existing files as a c string.

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
		HybridRNA(const char filename1[], const int type1, const char filename2[], const int type2, const bool IsRNA=true);


		//!	TwoRNA Constructor that copies thermodynamic parameter tables from an existing Thermodynamics (or RNA etc) instance.
		//! The file opening is performed by the constructors for the RNA classes that underlie each sequence.
		//!	This constructor generates internal error codes that can be accessed by GetErrorCode() after the constructor is called.  0 = no error.
		//! The errorcode can be resolved to a c string using GetErrorMessage.		
		//! \param filename1 is a null terminated c string and refers to sequence 1.
		//! \param filename2 is a null terminated c string and refers to sequence 2.
		//! \param type1 is an integer that indicates the file type for sequence 1.
		//! \param type2 is an integer that indicates the file type for sequence 2.
		HybridRNA(const char filename1[], const int type1, const char filename2[], const int type2, const char *alphabet);

		//******************************************************
		//Functions that predict bimolecular secondary structure 
		//******************************************************

		
		//! Predict the lowest free energy secondary structure for two interacting strands and generate suboptimal structures.
		//! This method does not allow intramolecular pairs.  It considers accessibility with a heuristic that uses the partition function.

		//! If the temperature has not been specified using the RNA base class SetTemperature and no free energies have been calculated, the
		//!		thermodynamic parameters have not been read and therefore they will be read by this function call.  The 
		//!		parameter files should be located in the directory specified by the environment variable $DATAPATH of the pwd.
		//!	In case of error, the function returns a non-zero that can be parsed by GetErrorMessage() or GetErrorMessageString().
		//! \param gamma is a scaling factor that weights accessibility.  The defaults is 0.4
		//!	\param percent is the maximum % difference in free energy in suboptimal structures from the lowest free energy structure.  The default is 50.
		//!	\param maximumstructures is the maximum number of suboptimal structures to generate.  The fefault is 20.
		//!	\param window is a parameter that specifies how different the suboptimal structures should be from each other (0=no restriction and larger integers require structures to be more different).  The default is 0.
		//!	\param maxinternalloopsize is the maximum number of unpaired nucleotides in bulge and internal loops.  This is used to accelerate the prediction speed.  The default is 30.
		//! \return An int that indicates an error code (0 = no error, 5 = error reading thermodynamic parameter files, 14 = traceback error).
		int AccessFold(const double gamma=0.4, const float percent=50, const int maximumstructures=20, const int window=0, const int maxinternalloopsize = 30);
		
		//! Predict the lowest free energy secondary structure and generate suboptimal structures using a heuristic.

		//! This function predicts the lowest free energy structure and suboptimal structures.
		//! If the temperature has not been specified using the RNA base class SetTemperature and no free energies have been calculated, the
		//!		thermodynamic parameters have not been read and therefore they will be read by this function call.  The 
		//!		parameter files should be located in the directory specified by the environment variable $DATAPATH of the pwd.
		//!	In case of error, the function returns a non-zero that can be parsed by GetErrorMessage() or GetErrorMessageString().
		//!	\param percent is the maximum % difference in free energy in suboptimal structures from the lowest free energy structure.  The default is 10.
		//!	\param maximumstructures is the maximum number of suboptimal structures to generate.  The defaults is 20.
		//!	\param window is a parameter that specifies how different the suboptimal structures should be from each other (0=no restriction and larger integers require structures to be more different).  The default is 0.
		//!	\param savefile is c string containing a file path and name for a savefile (.sav)that can be used to generate energy dot plots and to refold the secondary structure using different suboptimal structure parameters.  The default is "", which results in no save file written.
		//!	\param maxinternalloopsize is the maximum number of unpaired nucleotides in bulge and internal loops.  This is used to accelerate the prediction speed.  The default is 30.
		//! \return An int that indicates an error code (0 = no error, 5 = error reading thermodynamic parameter files, 14 = traceback error).
		int FoldBimolecular(const float percent=10, const int maximumstructures=20, const int window=0, const char savefile[]="", const int maxinternalloopsize = 30);

		//! Predict the lowest free energy secondary structure for two strands that cannot form intramolecular pairs and generate suboptimal structures using a heuristic.

		//! This function predicts the lowest free energy bimolecular structure and suboptimal structures.
		//!	This function does not allow any folding constraints.
		//! If the temperature has not been specified using SetTemperature and no free energies have been calculated, the
		//!		thermodynamic parameters have not been read and therefore they will be read by this function call.  The 
		//!		parameter files should be located in the directory specified by the environment variable $DATAPATH of the pwd.
		//!	In case of error, the function returns a non-zero that can be parsed by GetErrorMessage() or GetErrorMessageString().
		//!	\param percent is the maximum % difference in free energy in suboptimal structures from the lowest free energy structure.  The default is 40.
		//!	\param maximumstructures is the maximum number of suboptimal structures to generate.  The default is 10.
		//!	\param window is a parameter that specifies how different the suboptimal structures should be from each other (0=no restriction and larger integers require structures to be more different).  The default is 0.
		//!	\param maxinternalloopsize is the maximum number of unpaired nucleotides in bulge and internal loops.  This is used to accelerate the prediction speed.  The default is 30.
		//! \return An int that indicates an error code (0 = no error, 5 = error reading thermodynamic parameter files, 14 = traceback error).
		int FoldDuplex(const float percent=40, const int maximumstructures=10, const int window=0, const int maxinternalloopsize = 30);

		
		//! Predict the bimolecular partition function for a sequence (with no intramolecular pairs).

		//! This function must be called to predict base pair probabilities, perform stochastic traceback, or for maximizing expected accuracy.
		//! This predicts the partition function without intramolecular pairs.
		//! If the temperature has not been specified using SetTemperature and no free energies have been calculated, the
		//!		thermodynamic parameters have not been read and therefore they will be read by this function call.  The 
		//!		parameter files should be located in the directory specified by the environment variable $DATAPATH of the pwd.
		//!	In case of error, the function returns a non-zero that can be parsed by GetErrorMessage() or GetErrorMessageString().
		//!\param savefile is a c string that contains the path and filename for creating a save file.  This defaults to "", which indicates no file is to be written.
		//! \return An int that indicates an error code (0 = no error, 5 = error reading thermodynamic parameter files).
		int PartitionFunctionBimolecular(const char savefile[]=""); 

		//******************************************************************
		//Functions that provide access to underlying RNA classes from squences, an instance of TwoRNA:
		//These functions give access to the input sequences.
		//******************************************************************
		
		//!Access the underlying RNA class from an instance of TwoRNA.
		//!This is provided for use with two sequence methods.
		//!Generally, there is no need for end users to use this function.

		//!\return A pointer to the underlying RNA class for sequence 1.
		RNA *GetRNA1();

		//!Access the underlying RNA class from an instance of TwoRNA.
		//!This is provided for use with two sequence methods.
		//!Generally, there is no need for end users to use this function.

		//!\return A pointer to the underlying RNA class for sequence 2.
		RNA *GetRNA2();
		

		//******************************************************
		//Functions to return error information
		//******************************************************

		//!Return an error code, where a return of zero is no error.

		//!	This function returns and error flag that is generated during construction by RNA(const char &filename, const int type, const bool IsRNA=true) or from CalculateFreeEnergy().
		//!		An error of zero is always no error.  Other codes are errors and a c-string can be fetched for the error with GetErrorMessage().
		//!\return An integer that provides the error code.
		int GetErrorCode();

		//!	Return error messages based on code from GetErrorCode and other error codes.		

		//!		0 = no error
		//!		<100 = Error to be fetched from RNA base class.
		//!		100-999 = Error associated with bimolecular folding, to be handled here.
		//!		>=1000 = Errors for underlying sequence, get message from TwoRNA base class.
		//!\param error is the integer error code provided by GetErrorCode().
		//!\return A pointer to a c string that provides an error message or from other functions that return integer error codes.
		const char* GetErrorMessage(const int error);

		//********************************************************************************
		//Functions that specify or report constraints on folding.
		//These constraints only affect subsequent calls to structure prediction routines.
		//********************************************************************************

		//! Get whether intramolecular pairs are allowed.
		
		//!\return A bool that indicates whether intramolecular pairs are forbidden (true = forbidden, false = not).
		bool GetForbidIntramolecular();

		//! Set whether intramolecular pairs are allowed.
		
		//! If true is passed to this function, intramolecular pairs will be forbidden in FoldBimolecular. 
		//!\param forbid is a bool that indicates whether intramolecular pairs are forbid.
		void SetForbidIntramolecular(const bool forbid);

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
		~HybridRNA();	

	private:
		//add another structure class to contain the predicted bimolecular structures
		//RNA *RNADuplex;
		//bool RNADuplexallocated;

		//a function to encapsulate common constructor functions.
		void commonconstructor(const char fileNameOrsequence1[], const int fileType1, const char fileNameOrsequence2[], const int fileType2);

		//a flag to indicate if unimoleculr pairs are forbidden
		bool forbidunimolecular;

		//This function makes a composite sequence in the RNA::ct with each of the sequences in TwoRNA
		void SetupBimolecular();

		//Use a two RNA class to store the two sequences.
		TwoRNA *sequences;




};

#endif //HYBRID_RNA defined

#if !defined(OLIGOWALKOBJECT_H)
#define OLIGOWALKOBJECT_H


#include "RNA.h"
#include "../src/intermolecular.h"


//! Oligowalk_object Class.
/*!
	The Oligowalk_class inhereits from RNA and provides the OligoWalk functionality.  
	Additionally, it provides OligoScreen.
*/


//Note the stylized comments provide facility for automatic documentation via doxygen.

class Oligowalk_object:  public RNA {


	public:

		//******************************************************
		//Constructors
		//******************************************************


		//!Constructor - user provides a sequence as a c string.

		//! This constructor simply calls the underlying RNA(const char sequence[], const bool IsRNA=true) constructor and sets IsRNA=true.
		//!	Input sequence should contain A,C,G,T,U,a,c,g,t,u,x,X.
		//!	Capitalization makes no difference.
		//!	T=t=u=U.  U is assumed by OligoWalk.
		//!	x=X= nucleotide that neither stacks nor pairs.
		//!	For now, any unknown nuc is considered 'X'.
		//! Note that sequences will subsequently be indexed starting at 1 (like a biologist), so that the 0th position in the sequence array will be nucleotide 1.
		//!	\param sequence is a NULL terminated c string.
		Oligowalk_object(const char sequence[]); 
		
		//!Constructor - user provides a filename for existing file as a c string.

		//! This constructor simply calls the underlying RNA(const char filename[], const int type, const bool IsRNA=true) constructor and sets IsRNA=true.
		//!	The existing file, specified by filename, can either be a ct file, a sequence, or an RNAstructure save file. 
		//!	Therefore, the user provides a flag for the file: 
		//!		type = 1 => .ct file, type = 2 => .seq file, type = 3 => partition function save (.pfs) file, type = 4 => folding save file (.sav).
		//! For OligoWalk calculations, types 1 and 2 are the relevant types.
		//!	This constructor generates internal error codes that can be accessed by GetErrorCode() after the constructor is called.  0 = no error.
		//! The errorcode can be resolved to a c string using GetErrorMessage.		
		//! Note also that save files explicitly store the thermodynamic parameters, therefore changing the backbone type as compaared to the original calculation will not change structure predictions.
		//! \param filename is null terminated c string.
		//! \param type is an integer that indicates the file type.
		Oligowalk_object(const char filename[], const int type);

		//! Default Constructor - user provides nothing.

		//! This constructor calls the underlying RNA(IsRNA=true) constructor.
		//! This basic constructor is provided for performing OligoScreen calculations, which do not need an input sequence.
		//! This constructor needs to be called with RNA=true for Oligos to be RNA and RNA=false for Oligos to be DNA.
		//! \param IsRNA is a bool where true=RNA and false=DNA.
		Oligowalk_object(const bool IsRNA=true);

		//******************************************************
		//Function to perform and support the OligoWalk calculation.
		//******************************************************

		//!Perform an OligoWalk calculation.
		//! Note that this can only be performed once.
		//! \param oligo_length is an int that gives the length of the oligonucleotides.
		//! \param isDNA is a bool that indicates the oligonucleotide chemistry, where true = DNA and false = RNA.
		//! \param option is an int that gives the calculation type, where option 1 = break local target structure to bind oligo, option 2 = refold target RNA after oligo binding, and option 3 = no target structure considered.
		//! \param oligo_concentration is a double that indicates the oligonucleotide concentration in M.
		//! \param usesub is an int that indicates whether suboptimal structures are to be used, where 0 = none and 3 = use heuristic method.
		//! \param start is an int that indicates the starting location of the walk.
		//! \param stop is an int that indicates the ending location of the walk.
		//! \return An integer that indicates an error code that can be parsed by GetErrorMessage() or GetErrorMessageString(), 0 = no error.
		int Oligowalk(const int oligo_length, const bool isDNA, const int option, const double oligo_concentration, const int usesub, const int start, const int stop); 

		//! Get the breaking target DG for a given nucleotide.
		//! This can only be called after performing a valid OligoWalk calculation.
		//!	In case of error, the function returns a free energy change of zero.
		//!		Note!: That a free energy change of zero is also a valid folding free energy change.
		//!	Errors will also generate an internal error code, accessible with GetErrorCode().
		//! \param index is an int that specifies the 5' end of an oligonucleotide binding site on the target.
		//! \return A double that is the free energy change in kcal/mol.
		double GetBreakTargetDG(const int index);

		//! Get the duplex DG for a given nucleotide.
		//! This can only be called after performing a valid OligoWalk calculation.
		//!	In case of error, the function returns a free energy change of zero.
		//!		Note!: That a free energy change of zero is also a valid folding free energy change.
		//!	Errors will also generate an internal error code, accessible with GetErrorCode().
		//! \param index is an int that specifies the 5' end of an oligonucleotide binding site on the target.
		//! \return A double that is the free energy change in kcal/mol.
		double GetDuplexDG(const int index);

		//! Get the bimolecular oligo-oligo DG for a given nucleotide.
		//! This can only be called after performing a valid OligoWalk calculation.
		//!	In case of error, the function returns a free energy change of zero.
		//!		Note!: That a free energy change of zero is also a valid folding free energy change.
		//!	Errors will also generate an internal error code, accessible with GetErrorCode().
		//! \param index is an int that specifies the 5' end of an oligonucleotide binding site on the target.
		//! \return A double that is the free energy change in kcal/mol.
		double GetOligoOligoDG(const int index);

		//! Get the oligo-self DG for a given nucleotide.
		//! This can only be called after performing a valid OligoWalk calculation.
		//!	In case of error, the function returns a free energy change of zero.
		//!		Note!: That a free energy change of zero is also a valid folding free energy change.
		//!	Errors will also generate an internal error code, accessible with GetErrorCode().
		//! \param index is an int that specifies the 5' end of an oligonucleotide binding site on the target.
		//! \return A double that is the free energy change in kcal/mol.
		double GetOligoSelfDG(const int index);

		//! Get the overall DG for a given nucleotide.
		//! This can only be called after performing a valid OligoWalk calculation.
		//!	In case of error, the function returns a free energy change of zero.
		//!		Note!: That a free energy change of zero is also a valid folding free energy change.
		//!	Errors will also generate an internal error code, accessible with GetErrorCode().
		//! \param index is an int that specifies the 5' end of an oligonucleotide binding site on the target.
		//! \return A double that is the free energy change in kcal/mol.
		double GetOverallDG(const int index);

		//! Get the Tm for a given nucleotide.
		//! This can only be called after performing a valid OligoWalk calculation.
		//!	In case of error, the function returns a free energy change of zero.
		//!		Note!: That a free energy change of zero is also a valid folding free energy change.
		//!	Errors will also generate an internal error code, accessible with GetErrorCode().
		//! \param index is an int that specifies the 5' end of an oligonucleotide binding site on the target.
		//! \return A double that is the Tm in degrees C.
		double GetTm(const int index);

		//!	Write a report for an OligoWalk calculation.
		//! This must be called after performing a valid OligoWalk calculation.
		//! \param outputfilename is a c-string that provides a filename to which the report will be written.
		//! \param oligo_length is an int that gives the length of the oligonucleotides.
		//! \param isDNA is a bool that indicates the oligonucleotide chemistry, where true = DNA and false = RNA.
		//! \param option is an int that gives the calculation type, where option 1 = break local target structure to bind oligo, option 2 = refold target RNA after oligo binding, and option 3 = no target structure considered.
		//! \param oligo_concentration is a double that indicates the oligonucleotide concentration in M.
		//! \param usesub is an int that indicates whether suboptimal structures are to be used, where 0 = none and 3 = use heuristic method.
		//! \param start is an int that indicates the starting location of the walk.
		//! \param stop is an int that indicates the ending location of the walk.
		//! \return An integer that indicates an error code that can be parsed by GetErrorMessage() or GetErrorMessageString(), 0 = no error.
		int WriteReport(const char outputfilename[], const int oligo_length, const bool isDNA, const int option, const double oligo_concentration, const int usesub, const int start, const int stop);


		//******************************************************
		//Function to perform OligoScreen calculation.
		//******************************************************

		//! This function runs OligoScreen.

		//! Read a list of oligonucleotides in infilename and output thermodynamic characteristics in outfilename.
		//! The backbone type is set when the constructor is called.  Note that oligoscreen has no target sequence, so the default constructor OligoWalk(bool IsRNA) should be used.
		//! \param infilename is a c-string that indicates the filename to be read.
		//! \param outfilename is a c-string that indicates the name of an output file to be written with report.
		//! \return An int that indicates an error code that can be parsed by GetErrorMessage() or GetErrorMessageString(), 0 = no error.
		int OligoScreen(const char infilename[], const char outfilename[]);

		//******************************************************
		//Functions to perform error reporting.
		//******************************************************

		//!	Return error messages based on code from GetErrorCode and other error codes.		
		//!		100 = no OligoWalk data present
		//!		101 = OligoWalk has already been performed
		//!		other codes are handled by the RNA base class function GetErrorMessage.
		//!\param error is the integer error code provided by GetErrorCode()/
		//!\return A pointer to a c string that provides an error message or from other functions that return integer error codes.
		const char* GetErrorMessage(const int error);

		//******************************************************
		//Destructor:
		//******************************************************

		//! This is the destructor.
		~Oligowalk_object();



	private:

		//Common code for all constructors
		void CommonConstructor();

		//Tables needed for storing OligoWalk results
		int **table,**numofsubstructures;
		siPREFILTER *prefilter;

		//Keep track of oligo_length to facilitate deletion of table and numofsubstructures arrays
		int length;


};

#endif //!defined(DYNALIGNOBJECT_H)
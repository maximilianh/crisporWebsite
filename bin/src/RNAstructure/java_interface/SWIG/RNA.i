%module RNABackend;

#define COMPILE_SMP
#define _JAVA_GUI

%{
#include "../../RNA_class/RNA.h"
// #include "../../src/StructureTools.h"
#include "../../RNA_class/TwoRNA.h"
#include "../../RNA_class/HybridRNA.h"
#include "../../RNA_class/Dynalign_object.h"
%}
%include "std_string.i"
%include "std_vector.i"

namespace std {
   %template(CIntVector) vector<int>;
   %template(CDblVector) vector<double>;
};

//  ****************************************************************************************************
//  Ignore Overloads (for C++ functions with default arguments) 
//  ****************************************************************************************************
//  This section adds %ignore and %rename directives to stop SWIG from generating interface code for the
//  overloads of C++ functions with default arguments. Directives use pattern-matching to determine which
//  elements of an interface they should affect.
//  The %ignore directive before a functional form that INCLUDES default arguments causes ALL overloads 
//  of the function that where those arguments are not supplied to be ignored (including the overload in 
//  which all arguments are specified).
//     e.g.  %ignore MyFunc(arg1=default, arg2=default); // ignore MyFunc(arg1, arg2), MyFunc(arg1), and MyFunc()
//  In contrast, without the default arguments, %ignore causes only that specific overload to be ignored:
//     e.g.  %ignore MyFunc(arg1);      // ignore only MyFunc(arg1), but not MyFunc(arg1,arg2) or MyFunc()
//  The %rename("%s") directive essentially UN-ignores matching overloads:
//     e.g.  %rename("%s")  MyFunc(arg1, arg2);       // UN-ignore a specific form of the function

// The next line includes auto-generated instructions to ignore overloads of functions with default arguments.
   %include "remove-default-args.i"

// If there are any overloads that SHOULD be kept, list them below using the %rename syntax shown above.

	//// -- In general, keep these forms that only allow structurenumber to be omitted. They are probably useful overloads.
	//// RNA::SpecifyPair(const int i, const int j, const int structurenumber = 1);
	//// RNA::RemovePairs(const int structurenumber = 1);
	//// RNA::RemoveBasePair(const int i, const int structurenumber =1);
	//// RNA::CalculateFreeEnergy(const int structurenumber = 1, const bool UseSimpleMBLoopRules = false);

//  ****************************************************************************************************
//  Ignore functions that are redundant due to target language type-mapping
//  ****************************************************************************************************
//  The Java datatype mapping for both `char*` and `string` is java.lang.String, so both `string MyFunc()` 
//  and `char* MyFunc()`  are translated to Java as `String MyFunc()` so they are redundant.
//  SWIG catches this (if the functions have the same name), but explicitly ignoring one of the functions
//  avoids a warning. In cases where the functions have different names, SWIG will hapilly generate both,
//  but this results in excessive/redundant functions.
	%ignore structure::SetCtLabel(char *,int);  // exact duplicate of SetCtLabel(string,int) 
	%ignore RNA::GetErrorMessageString;         // unnecessary duplicate of char* GetErrorMessage
	%ignore operator[];
	%ignore PartialProgress::workComplete;
	
//  ****************************************************************************************************
//  Ignore unnecessary functions with pointer parameters to avoid SWIGTYPE_p_XXX types
//  ****************************************************************************************************
	%ignore structure::constant;  		// avoid SWIGTYPE_p_p_double
	%ignore singlestructure::basepr;  	// avoid SWIGTYPE_p_vector_int
	%ignore TProgressDialog::TProgressDialog(std::ostream *outputStream); 	// avoid the ostream object
	%ignore eraseEnergyLabel(string &label);
	%ignore RNA::Rsample; // avoid RsampleData and vector<double>
	%ignore structure::SetCtLabel(char const *,int const); //duplicate of structure::SetCtLabel(std::string const &,int const)
        %ignore structure::LoadSHAPE;
        %ignore structure::CopySHAPE;

	// Ignore the CTCommentProvider parameter for WriteCt, WriteDotBracket, ctout, and writedotbracket
	%ignore RNA::WriteCt(const char filename[], bool append=false, CTCommentProvider &commentProvider=CTComments::Energy) const;
	%ignore RNA::WriteDotBracket(const char filename[], const int structurenumber=-1, const DotBracketFormat format=DBN_FMT_MULTI_TITLE, CTCommentProvider &commentProvider=CTComments::Energy) const;
	%rename(%s) RNA::WriteCt(const char filename[]) const;
	%rename(%s) RNA::WriteDotBracket(const char filename[], const int structurenumber=-1, const DotBracketFormat format=DBN_FMT_MULTI_TITLE) const;
	%ignore structure::ctout(const char * const ctoutfile, const bool append=false, CTCommentProvider &commentProvider=CTComments::Energy) const;
	%ignore structure::writedotbracket(const char * const filename, const int structurenumber=-1, const DotBracketFormat format=DBN_FMT_MULTI_TITLE, CTCommentProvider &commentProvider=CTComments::Energy, const bool append = false) const;
	%rename(%s) structure::ctout(const char * const ctoutfile) const;		
	%rename(%s) structure::writedotbracket(const char * const filename, const int structurenumber=-1, const DotBracketFormat format=DBN_FMT_MULTI_TITLE) const;
	%rename(GetMyErrorMessage) HybridRNA::GetErrorMessage(const int error);

	%rename("%s") BreakPseudoknot(const bool minimum_energy=true, const int structurenumber, const bool useFastMethod = true); // keep the last parameter optional
	%ignore structure::BreakPseudoknots(int structurenumber=1, vector<int> *brokenPairs=NULL);
	%rename("%s") structure::BreakPseudoknots(int structurenumber);
	%rename("%s") findPseudoknots(const vector<int> &currentPairs, vector<int> *pseudoknotPairs = NULL, vector<int> *normalPairs= NULL); // keep all parameters optional.

//  ****************************************************************************************************
//  Ignore unnecessary global functions/consts (which would normally be put in a *Proxy.java file).
//  ****************************************************************************************************
	%ignore ::tonumi(char *base);
	%ignore ::tobase (int i);
	%ignore ::tonum(char *base,structure *ct,int count);
	%ignore ::ecompare(const void *i, const void *j);
	%ignore ::swap;
	%ignore ::maxforce;
	%ignore ::TOLERANCE;	

	%ignore ergcoaxflushbases;
	%ignore ergcoaxinterbases1;
	%ignore ergcoaxinterbases2;
	%ignore ergcoax;
	%ignore decon1;
	%ignore decon2;
	%ignore ergmulti;
	%ignore ergexterior;
	%ignore erg1;
	%ignore erg2;
	%ignore erg2in;
	%ignore erg2ex;
	%ignore erg3;
	%ignore erg4;
	%ignore penalty;
	
	%ignore openstructuresave;
	%ignore writestructuresave;
	%ignore SHAPEend;


//  ****************************************************************************************************
//  Provide empty class/struct definitions (to generate MyClass instead of SWIGTYPE_p_MyClass)
//  ****************************************************************************************************
	struct datatable {};  // provide an empty definition to generate a datatable.java type instead of SWIGTYPE_p_datatable.

//  ****************************************************************************************************
//  Allow some c++ functions to return Java arrays
//  ****************************************************************************************************
	%include "arrays_java.i"
	//%apply int[] { int* arr };
	%apply int[] { intArr };
	%apply double[] { double* arr };
	%apply float[] { float* arr };

//  ****************************************************************************************************
//  Include c++ enums as Java enums (instead of classes).
//  ****************************************************************************************************
	%include "enums.swg"

//  ****************************************************************************************************
//  Finally, after all above directives, include the headers for which Java proxies should be generated.
//  ****************************************************************************************************

//    %include "../../src/RNAFileType.h"
	%include "../../src/DotBracketFormat.h"
    %include "../../src/TProgressDialog.h"
    %include "../../src/structure.h"
    %include "../../RNA_class/thermodynamics.h"
    %include "../../RNA_class/RNA.h"
    // %include "../../src/StructureTools.h"
	%include "../../RNA_class/TwoRNA.h"
	%include "../../RNA_class/HybridRNA.h"
	%include "../../RNA_class/Dynalign_object.h"

//  ****************************	
//  Extend some classes to provide better compatibility
//  ******************************	
	%extend singlestructure {
		int getSequenceLength() { return $self->basepr.size()-1; }
		void setSequenceLength(int length) { $self->basepr.resize(length+1); $self->energy = 0; }
		int getBasePair(int basepos) { return $self->basepr[basepos]; }
		void setBasePair(int basepos, int pairpos) { $self->basepr[basepos] = pairpos; }
	};
	%extend RNA {
		int getInts(intArr arr, int size) {
			if (size >= 1000)
				for (int i = 0; i < size; i++) {
					arr[i] = i;
				}
			return 1000;
		}
		void initInts() {
			gotInts = new int[1000];
			RNA_getInts($self, gotInts, 1000);
		}
		int getIntAt(int i) { 
			return gotInts[i]; 
		}
	};
	%runtime  %{
		typedef int *intArr;
		intArr gotInts;
	%}
	

	const char* getDataPath(const char* const defaultPath = DATAPATH_DEFAULT);
    void setDataPath(const char* const path);

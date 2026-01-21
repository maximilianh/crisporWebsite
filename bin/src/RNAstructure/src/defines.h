// defines.h
#ifndef DEFINES_H
#define DEFINES_H

#define maxfil 350    //maximum length of file names
#define INFINITE_ENERGY 14000  //an arbitrary value given to infinity
#define maxtloop 200 //maximum tetraloops allowed (info read from tloop, triloop, and hexaloop)
#define DYNALIGN_INFINITY 14000  //an arbitrary value given to infinity
#define maxstructures 1010 //maximum number of structures in ct file -- deprecated
#define maxbases 3000   //maximum number of bases in a structure -- deprecated
#define ctheaderlength 250 //maximum length of string containing info on sequence
#define ga_bonus -10 //the value of a bonus for the "almost coaxial stacking"
								//case in efn2
#define amax 400 //this is a maximum line length for void linout (below)
#define col 80  //this is the number of columns in an output file
#define numlen 8  //maximum digits in a number
#define maxforce 3000 //maximum number of bases that can be forced single or paired 
#define maxneighborlength 25 //maximum length of paired neighbors from NMR 
#define maxregions 10//maximum number of regions for NMR constraints
#define maxgu 100 //maximum number of u's in gu pair

#define minloop 3 //The minimum substructure with a basepair possible
#define maxsequencelinelength 4000
#define RT_37C 0.61633 //gas constant times 37 deg. C  (previously named "rt", but this could easily cause problems e.g. "int rt = 3;" would be an error.)
#define TEMP_37C 310.15 // The temperature 37 degrees celsius converted to units Kelvin
#define Rgas 1.987213//gas constant (to correct for the energy units)
#define RKC 1.987213e-3 //gas constant in Kcal/mol
#define SINGLE 1
#define PAIR 2
#define NOPAIR 4
#define DUBLE 8 //corruption of DOUBLE because Visual Studio uses DOUBLE
#define INTER 16
#define scalingdefinition 0.6 //0.6  //per nuc scale in partition function
#define PFMAX 1e300  //maximum size of storage variable in partition function
#define	PFMIN 1e-300 //minimum size of storage variable

// Defines that determine partition function scaling.  The default behavior is to calculate in the explicit log scale.

// When PF_LINEAR is defined, the partition function computes in the linear scale using a per-nucleotide scaling factor
// #define PF_LINEAR

// When PF_LOG_CLASS is defined, the partition function computes using the log_double class.  The calculations are implicitly 
// log scale but treated like they are in the linear scale without a scaling factor
// #define PF_LOG_CLASS //use log_double class

// Defines for how calculations are performed in the explicit log scale
#define USE_XLOG_ZERO
#define USE_LOGP1_LOOKUP_SUM
#define LOG_LOOKUP_ENTRIES 10000

#define EPSILON 1e-300 //estimate of machine precision
#define SCALEUP 1.05 //use to scale up partition function scaling 
#define SCALEDOWN 0.95 //used to scale down the partition function scaling
						//Note that scaleup and scaledown need to be safe to raise
						//to the power of the sequence length.

#define PERCENTDOTS 30 //Used by Dynalign  The maximum percent difference
						//in folding free energy compared to the lowest 
						//free energy structure (determined by single
						//sequence structure prediction) for pairs
						//that will be allowed in a Dynalign calculation
						//by default.


//The SIMPLEMBLOOP affects the behavior of the partition function AND stochastic traceback codes
//#define SIMPLEMBLOOP  //use this line to use an energy model that only includes 3' dangling ends
						//on each helix in multibranch and exterior loops
#undef SIMPLEMBLOOP  //use this line to use an energy model that explicitly considers all dangles 
						//and/or coaxial stacks in helices


//The dynamic programming algorithms require integer math:
//The following three definitions are for tenths precision in parameters
	//allowing short integers for the free energy array..
	//DIGITS is used to format the energy output.
#define conversionfactor 10
#define conversionprecision 1 // this should always be defined as LOG10(conversionfactor). It is the units of precision that an energy has.  (e.g. it should be 2 for conversionfactor=100 or 0 for conversionfactor=1 )
// #define integersize short
#define integersize int
#define DIGITS "%.1f"
//The following three definitions are for hundreths precision in parameters
	//requiring full integer arrays for the the free energy arrays.
// #define conversionfactor 100
// #define conversionprecision 2 // this should always be defined as LOG10(conversionfactor). It is the units of precision that an energy has.  (e.g. it should be 2 for conversionfactor=100 or 0 for conversionfactor=1 )
// #define integersize int
// #define DIGITS "%.2f"

// hairpin with less than 3 nucleotides 
// for use in AlignmentFold
#define hairpinlt3 conversionfactor*8
#define hairpinlt3k -12.98 // for AlignmentPartition, equal to boltzman(hairpinlt3, 310.15)
#define safiversion 7//version is the save file version for single sequences
#define pfsaveversion 9//this is the version of save file format for partition functions
#define dynalignsavefileversion 1//version of the save files used by Dynalign
#define T37inK 310.15//37 degrees C in Kelvin
#define DH_DS_CORR 0.9996

// DATAPATH variables.  Note: Use getDataPath() and setDataPath() from rna_library.h
#define DATAPATH_ENV_VAR "DATAPATH"
#define DATAPATH_DEFAULT "."

// ####### Utility / Debugging Macros #######
#define TO_STRING_MACRO(x) #x //helper for expansion of the TO_STRING macro.
#define TO_STRING(x) TO_STRING_MACRO(x) // converts any macro to a c-string of its expanded value. For exmaple: TO_STRING(conversionfactor) => "10"

//Macros for cout-style debugging. Requires "#include <iostream>"  Usage: CDBG("Some debug message or value")
#define CDBG_STREAM std::cout  // can be redefined locally to switch to cerr etc.
#define CDBG(MSG) CDBG_STREAM << "DEBUG:\t" << MSG << "\t(" << __FILE__ << ":" << __LINE__ << ")" << std::endl;
#define CDBG_LINE CDBG("")

#endif // DEFINES_H

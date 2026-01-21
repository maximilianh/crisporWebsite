// defines.h

#if !defined(DEFINES_H)
#define DEFINES_H




#define maxfil 350    //maximum length of file names
#define infinity 14000  //an arbitrary value given to infinity
#define maxtloop 100 //maximum tetraloops allowed (info read from tloop, triloop, and hexaloop)
#define maxstructures 5 //maximum number of structures in ct file -- deprecated
#define maxbases 3000   //maximum number of bases in a structure -- deprecated
#define ctheaderlength 125 //maximum length of string containing info on sequence
#define ga_bonus -10 //the value of a bonus for the "almost coaxial stacking"/case in efn2
#define amax 400 //this is a maximum line length for void linout (below)
#define col 80  //this is the number of columns in an output file
#define numlen 8  //maximum digits in a number
#define maxforce 1000 //maximum number of bases that can be forced single or paired 
#define maxneighborlength 25 //maximum length of paired neighbors from NMR 
#define maxregions 10//maximum number of regions for NMR constraints
#define maxgu 100 //maximum number of u's in gu pair

#define minloop 3 //The minimum substructure with a basepair possible
#define maxsequencelinelength 4000
#define rt 0.61633 //gas constant times 37 deg. C 
#define R 1.987213//gas constant (to correct for the energy units)
#define RKC 1.987213e-3 //gas constant in Kcal/mol
#define SINGLE 1
#define PAIR 2
#define NOPAIR 4
#define DOUBLE 8
#define INTER 16
#define scalingdefinition 0.6  //per nuc scale in partition function
#define PFPRECISION double


//The dynamic programming algorithms require integer math:
//The following two definitions are for tenths precision in parameters
	//allowing short integers for the free energy array..
#define conversionfactor 10
#define conversionprecision 1 // this should always be defined as LOG10(conversionfactor). It is the units of precision that an energy has.  (e.g. it should be 2 for conversionfactor=100 or 0 for conversionfactor=1 )
#define integersize short
//The following two definitions are for hundreths precision in parameters
	//requiring full integer arrays for the the free energy arrays.
//#define conversionfactor 100
//#define conversionprecision 2 // this should always be defined as LOG10(conversionfactor). It is the units of precision that an energy has.  (e.g. it should be 2 for conversionfactor=100 or 0 for conversionfactor=1 )
//#define integersize int


#endif

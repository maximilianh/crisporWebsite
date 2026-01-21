
#if !defined(ALGORITHM_H)
#define ALGORITHM_H

#include "DynProgArray.h"
#include "forceclass.h"
#include "dotarray.h"
#include "rna_library.h"
#include "structure.h"
#include "TProgressDialog.h"

//***********************************Structures:

/////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////

void de_allocate (int **v,int i);//deallocates memory for a 2d array
void de_allocate (bool **v,int i);//alternative form of de_allocate
void de_allocate (short int **v,int i);//alternative form of de_allocate



//**********************************prototypes:



void getout (char *energyfile);//get name of file to output
										//	energy info

//! Calculate the energy of one or all sub-structure(s) in the specified structure object.
//! Optionally write thermodynamic details to the specified output file.
//! \param data A pointer to the thermodynamic datatables to use for the calculation.
//! \param structnum The 1-based index of the structure to calculate energy for.  Zero means calulate the energy for all structures.
//! \param simplemb Use the simple multibranch-loop algorithm.
//! \param outputfilename A c-string specifying the path to an output file to which thermodynamic details should be written. If NULL, no output file will be written.
//! \return Returns an error code if outputfilename was specified (not NULL) but it could not be opened or written to. Otherwise returns 0 (success).
int efn2(datatable *data,structure *ct, int structnum, bool simplemb, const char *outputfilename);

//! Calculate the energy of one or all sub-structure(s) in the specified structure object.
//! Optionally write thermodynamic details to the specified output stream.
//! \param data A pointer to the thermodynamic datatables to use for the calculation.
//! \param structnum The 1-based index of the structure to calculate energy for.  Zero means calulate the energy for all structures.
//! \param simplemb Use the simple multibranch-loop algorithm.
//! \param output A pointer to an opened output stream to which thermodynamic details should be written. If NULL (the default), no output will be written.
void efn2(datatable *data,structure *ct, int structnum = 0, bool simplemb = false, ostream *output=NULL);//energy calculator

double ergexteriordiff(datatable *data,structure *ct, int structnum, bool simplemb, int min_index, int max_index);

void energyout(structure *ct,char *enrgyfile);

//dynamic programming algorithm for secondary structure prediction by free energy minimization
	//this is the dynamic folding algorithm of Zuker
         //cntrl6 = #tracebacks
         //cntrl8 = percent sort
         //cntrl9 = window
		//ProgressHandler is an interface for returning the progress of the calculation (TProgressDialog is a sub class of ProgressHandler suitable for Text output)
		//Savfile is for creating a file with arrays and parameters for refolding with different 
			//suboptimal tracebacks
		//quickenergy indicates whether to determine the lowest free energy for the sequence without a structure
		//quickstructure is a bool that will generate only the lowest free energy structure.  No savefiles can generated. 
		//maxinter is the maximum number of unpaired nucleotides allowed in an internal loop
		//allow_isolated determines whether the heuristic to prevent isolated base pairs is used.  The default is false, i.e. filter isolated base pairs from forming.
	//This returns an error code, where zero is no error and non-zero indicates a traceback error.
int dynamic (structure *ct,datatable *data,int cntrl6,int cntrl8,int cntrl9,
			 ProgressHandler* update=0, bool quickenergy = false, char* savfile = 0, int maxinter = 30, bool quickstructure = false, bool simple_iloops = true, bool disablecoax=false, bool allow_isolated=false);


void fill(structure *ct, DynProgArray<integersize> &v, DynProgArray<integersize> &w, DynProgArray<integersize> &wmb, forceclass &fce, int &vmin,bool *lfce, bool *mod,
          integersize *w5, integersize *w3, bool quickenergy,
          datatable *data, DynProgArray<integersize> *w2, DynProgArray<integersize> *wmb2, DynProgArray<integersize> *we,ProgressHandler* update = 0, int maxinter = 30, bool quickstructure = false, bool simple_iloops = true, bool disablecoax=false,bool allow_isolated=false);


//The fill step of the dynamic programming algorithm for free energy minimization:
void fill(structure *ct, DynProgArray<integersize> &v, DynProgArray<integersize> &w, DynProgArray<integersize> &wmb, forceclass &fce, int &vmin,bool *lfce, bool *mod,
		  integersize *w5, integersize *w3, bool qickenergy,
		  datatable *data, DynProgArray<integersize> *w2, DynProgArray<integersize> *wmb2, ProgressHandler* update=0, int maxinter = 30, bool quickstructure = false, bool simple_iloops = true, bool disablecoax=false, bool allow_isolated = false);

//this overloaded dynamic function is used by NAPSS program to generate a special format dotplot
void dynamic (structure *ct,datatable* data,int cntrl6, int cntrl8,int cntrl9,
              DynProgArray<integersize> *v, DynProgArray<integersize> *vmb/*tracks MB loops*/, DynProgArray<integersize> *vext/*tracks exterior loops*/,
              ProgressHandler* update=0, bool quickenergy = false, char* savefile = 0, int maxinter = 30, bool quickstructure = false, bool simple_iloops = true, bool disablecoax=false, bool allow_isolated = false);
//this overloaded fill function is used to NAPSS program to generate a special format dotplot
void fill(structure *ct, DynProgArray<integersize> &v, DynProgArray<integersize> &vmb, DynProgArray<integersize> &vext, DynProgArray<integersize> &w, DynProgArray<integersize> &wmb, forceclass &fce, 
          int &vmin, bool *lfce, bool *mod,integersize *w5, integersize *w3, bool quickenergy,
          datatable *data, DynProgArray<integersize> *w2, DynProgArray<integersize> *wmb2, ProgressHandler* update=0, int maxinter = 30, bool quickstructure = false, bool simple_iloops = true, bool disablecoax=false, bool allow_isolated = false);

void errmsg(int err,int err1);//function for outputting info in case of an error
void update (int i);//function informs user of progress of fill algorithm

//filter is used to choose structures to output after efn2
//	this can make the number of structures more reasonable for inspection
//	it takes a structure, ct, which also contains the final output,
//	percent sort, maximum number of structures, and window size
void filter(structure* ct, int percent, int max, int window);

//Use the fill information to generate a set of suboptimal structures using the mfold heuristic.
	//This returns an error code, where zero is no error and non-zero indicates a traceback error.
int traceback(structure *ct, datatable *data, DynProgArray<integersize> *v, DynProgArray<integersize> *w, DynProgArray<integersize> *wmb, DynProgArray<integersize> *w2,DynProgArray<integersize> *wmb2, integersize *w3, integersize *w5, forceclass *fce,
	bool *lfce,integersize vmin, int cntrl6, int cntrl8, int cntrl9, bool *mod);

//this function is used to calculate the values of all the dots in a dot plot
void dotefn2(structure *ct, datatable *data, DynProgArray<integersize> *v, DynProgArray<integersize> *w, DynProgArray<integersize> *w2,
	int *w3, int *w5, short int **fce, bool *lfce,int vmin,dotarray *dots,
   ProgressHandler* PD = 0);
void calcpnum(dotarray *dots, int *pnum, int increment, short int numofbases,
	ProgressHandler *PD = 0);
void savefile(int i, std::ofstream* sav);//this function is used to make a save file
											//after the fill algorithm
short int readfile(std::ifstream *read);//this function is used to read save files
void savedot(dotarray *dots,structure *ct, char *filename); //save dot plot info
void readdot(dotarray *dots, structure *ct, char *filename);//read a dot plot file
void dpalign(dotarray *dots1,dotarray *dots2,structure* ct1,structure *ct2,short int *align);
short int getbestdot(dotarray *dots1,dotarray *dots2, structure* ct1,
	structure *ct2, short int i, short int j);//return the best dot for base i
   //in dots1 and j in dots2
//dpalign will align two dot plots and store the info in the array align
void energydump (structure *ct, DynProgArray<integersize> *v, datatable *data, int n,char *filename, int i, int j);
//energydump will spit out the composite free energies for a traceback
void energydump (structure *ct, datatable *data,DynProgArray<integersize> *v, int n,char *filename);
//energydump2 will spit out the composite free energies for a traceback -- with
//the au penalty associated with the correct entity
int checknp(bool lfce1,bool lfce2); 
//this function is used by the fill and trace to check whether nucleotides 
//contained in a dangling end are forced double stranded

void opensav(char* filename, structure* ct, int cntrl6, int cntrl8,int cntrl9);//opens a save file with information filled by
   									//fill algorithm

#ifdef DYNALIGN_II
int trace(structure *ct, datatable *data, int ii, int ji,
          DynProgArray<integersize> *v, DynProgArray<integersize> *w, DynProgArray<integersize> *wmb, DynProgArray<integersize> *w2, DynProgArray<integersize> *wmb2, 
          bool *lfce, forceclass *fce, integersize *w3, integersize *w5,bool *mod,DynProgArray<integersize> *we = NULL,integersize energy = 0,short open = 0, short pair = 0, bool quickstructure = false);
#else
int trace(structure *ct, datatable *data, int ii, int ji,
		DynProgArray<integersize> *v, DynProgArray<integersize> *w, DynProgArray<integersize> *wmb, DynProgArray<integersize> *w2, DynProgArray<integersize> *wmb2, 
          bool *lfce, forceclass *fce, integersize *w3, integersize *w5,bool *mod,bool quickstructure=false);
#endif
void readsav(const char *filename, structure *ct, DynProgArray<integersize> *w2, DynProgArray<integersize> *wmb2, 
			 integersize *w5, integersize *w3, bool *lfce, bool *mod, datatable *data,
			 DynProgArray<integersize> *v, DynProgArray<integersize> *w, DynProgArray<integersize> *wmb, forceclass *fce, int *vmin);

void writehelixfile(char *filename,structure *ct,int StructureNumber);//write a file with the helices

//reads the save file

//outputs a structure in cct format, which is shortened format compared to ct.
void cctout(structure *ct, char *filename);

bool multi_can_pair(int i, int j, structure *ct);
integersize multi_eparam(int i, structure* ct);
integersize multi_tstkm(int a, int b, int c, int d, structure* ct);
integersize multi_tstack(int a, int b, int c, int d, structure* ct);
integersize covar(int i, int j, structure* ct);
#endif

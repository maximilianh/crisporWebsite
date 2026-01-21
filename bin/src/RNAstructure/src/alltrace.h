#if !defined(ALLTRACE_H)
#define ALLTRACE_H

#include "TProgressDialog.h"
#include "defines.h"
#include "structure.h"
#include "algorithm.h"
#include "rna_library.h"
#include "stackclass.h"

////////////////////////////////////////////////////////////////////////
//DynProgArray encapsulates the large 2-d arrays of w and v, used by the dynamic
//	algorithm
class atDynProgArray {
   int Size;

   public:
   	
      int k;
      integersize **dg;
      integersize infinite;
	  bool allocated;
      

      //the constructor allocates the space needed by the arrays
   		atDynProgArray(int size);
		atDynProgArray();
		void allocate(int size);

      //the destructor deallocates the space used
      ~atDynProgArray();

      //f is an integer function that references the correct element of the array
   	inline integersize &f(int i, int j);
};

//f is an integer function that references the correct element of the array
inline integersize &atDynProgArray::f(int i, int j) {
      	
   if (i>j) {
        return infinite;
    }
   else if (i>Size) return f(i-Size,j-Size);//dg[i-Size][j-Size];
   else return dg[i][j-i];
         
}

//Function to generate all suboptimal structures for sequence in structure ct and store them in structure ct.
	//data is a pointer to datatable, where the thermodynamic parameters are stored.
	//percentdeta is the maximum percent difference from the lowest free energy structure.
	//absolutedelta is the maximum folding free energy change (in kcal/mol * conversionfactor, as defined in defines.h)
	//update is a TProgressDialog, used to track progress.
	//save is the name of a savefile, which generates save files that can be used by realltrace
	//NoMBLoop = whether multibranch loops are allowed, wehere true indicates NO multibranch loops
void alltrace(structure* ct,datatable* data, short percentdelta, short absolutedelta, ProgressHandler* update, char* save, bool NoMBLoop=false);
void readalltrace(char *filename, structure *ct, 
			 short *w5,  
			 atDynProgArray *v, atDynProgArray *w, atDynProgArray *wmb, atDynProgArray *wmbl, atDynProgArray *wl, atDynProgArray *wcoax,
			 atDynProgArray *w2, atDynProgArray *wmb2, forceclass *fce, bool *lfce, bool *mod, datatable *data);

void realltrace(char *savefilename, structure *ct, short percentdelta, short absolutedelta, char *ctname = NULL);


#define startingsize 500  //maximum number of structure fragments to start in alltracestructurestack (below)
//#define startingrefinementstacksize 25


//a stack to keep track of partially refined structures
	//contains a stackclass to keep track of where refinements need to occur 
class alltracestructurestack {

	public:
		short **basepairs;
		int maximumsize;
		int current; //current location in stack
		alltracestructurestack(short size, int sizeofstack=startingsize);
		stackclass *refinementstack;//a stack to keep track of where refinement is occuring
		~alltracestructurestack();
		void addpair(short i, short j, int index);
		void push();
		void push(integersize totalenergy, bool topair,short pairi, short pairj, bool tostack, short i, short j, short open, integersize energy, short pair, bool topair2,short pairi2, short pairj2, bool tostack2, short i2, short j2, short open2, integersize energy2, short pair2,bool tostack3=false, short i3=0, short j3=0, short open3=0, integersize energy3=0, short pair3=0);
		short numberofnucs; //length of the sequence
		void allocatearrays();
		void deletearrays();
		void pull();
		void pushtorefinement(short a, short b, short c, integersize d, short e);
		bool pullfromrefinement(short *a, short *b, short *ct, integersize *d, short *e);
		short *energy;
		integersize peekatenergy();
		void placeenergy(short energy);
		
		void flushbullpen();

		bool refined;
		bool bullpentopair,bullpentopair2;
		bool bullpentostack,bullpentostack2,bullpentostack3;
		short bullpeni,bullpenj,bullpenopen,bullpenpair,bullpenpairi,bullpenpairj;
		integersize bullpenenergy,bullpentotalenergy;
		short bullpeni2,bullpenj2,bullpenopen2,bullpenpair2,bullpenpairi2,bullpenpairj2;
		integersize bullpenenergy2;
		short bullpeni3,bullpenj3,bullpenopen3,bullpenpair3;
		integersize bullpenenergy3;
		short readpair(short i);
		void stackup(int index);

		//The following is infrastructure to keep track of stacked nucleotides:
		short stack1[2],stack2[2];
		void nstack(short i, short j, short k=0, short l=0);
		short **stacks;
		short readstacking(short i);

		#if defined (pfdebugmode) 
		ofstream out;

		#endif

};

#endif //!defined ALLTRACE_H

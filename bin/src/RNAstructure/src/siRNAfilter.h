#ifndef SIPREFILTER_H
#define SIPREFILTER_H

#include "rna_library.h"
#include "algorithm.h"
//=======================================================================
//siPREFILTER is a class to calculate the pre-filter score of functional siRNA
//Those siRNA with certain score will be chosed as functional candidate for folding target
class siPREFILTER {
		
	  int j;
	  double stack[5][5];
	  double end[5];
public:
      int useit;
	  bool usescore;
	  int *score;
	  float *melt;
	  double DG5,DG3;
	  double *enddiff;
	  datatable *data,*dhdata;

	  siPREFILTER(datatable& DATA,datatable& dhDATA,int useprefilter,bool scoreit,int size,bool isdna=false);
			//size is the total number of oligos
	        //the constructor allocates the space needed by the arrays
      ~siPREFILTER();
	        //the destructor deallocates the space used
	  void count(structure *ct,int i,int test=0);//i is the number of the oligo
};

#endif


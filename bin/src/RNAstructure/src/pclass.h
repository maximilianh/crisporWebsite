#if !defined(PCLASS_H)
#define PCLASS_H


/*======================================================================
	Two files (pclass.h and pclass.cpp) are made based on pfunction.h,
pfunction.cpp from the old version of RNAStructure.

	A new class ,Pclass, would be used instead of pfunction(...).
An inherited class from Pclass, OligoPclass, was also added to reused the 
partition function array (v,w,w5,...) for OligoWalk program.

														----Nov., 2005
															John(Zhi Lu)
  =======================================================================*/
//#pragma once 

#include <math.h>
#include "rna_library.h"
#include "algorithm.h"
#include "pfunction.h"
#include "DynProgArray.h"
#include "../src/phmm/utils/xmath/log/xlog_math.h"










void readpfsave(char *filename, structure *ct,PFPRECISION *w5, PFPRECISION *w3, 
				DynProgArray<PFPRECISION> *v, DynProgArray<PFPRECISION> *w, DynProgArray<PFPRECISION> *wmb, 
				DynProgArray<PFPRECISION> *wl, DynProgArray<PFPRECISION> *wmbl, DynProgArray<PFPRECISION> *wcoax,
				forceclass *fce, PFPRECISION *scaling, bool *mod, bool *lfce, 
				pfdatatable *data);
// double calculateprobability(int i, int j, DynProgArray<PFPRECISION> *v, 
// 								 PFPRECISION *w5, structure *ct, pfdatatable *data,
// 								 bool *lfce, bool *mod, PFPRECISION scaling,
// 								 forceclass *fce);
//inline void rescale(int i, int j,structure *ct, pfdatatable *data, DynProgArray<PFPRECISION> *v, DynProgArray<PFPRECISION> *w, DynProgArray<PFPRECISION> *wl, DynProgArray<PFPRECISION> *wcoax,
//			 DynProgArray<PFPRECISION> *wmb,DynProgArray<PFPRECISION> *wmbl, PFPRECISION *w5, PFPRECISION *w3, PFPRECISION **wca, PFPRECISION **curE, PFPRECISION **prevE, PFPRECISION rescalefactor); 
		//function to rescale all arrays when partition function calculation is headed out of bounds

int pfshape(structure *ct, PFPRECISION  temp);//scale the shape energy into partitionfunction

//=======================================================================
/*
Pclass would be used instead of pfunction(). 
This will fill the array (v,w,...) of partition function calculation with
O(N3) in calculation time, with limited or unlimited internal loop sizes.
The old algorithm of O(N4) is also avaiable by calling oldfill,oldpartition()
*/
class Pclass {

protected:	
	int i,j,h,d,maxj,lowi;
	int dp,ip,jp,ii,jj,jpf,jf,bl,ll;
	int k,l,m,n,o,p;
	int before,after;
	
	bool *lfce,*mod;//[maxbases+1][maxbases+1];
	PFPRECISION e;
	PFPRECISION twoscaling,rarray;
	
	forceclass *fce;	
	
	vector< vector<bool> > inc;
	

public:
	int number;	
	pfdatatable* data;
	structure* ct;
	PFPRECISION *w5,*w3,**wca,**curE,**prevE,**tempE;
	DynProgArray<PFPRECISION> *w,*v,*wmb,*wl,*wmbl,*wcoax;
	
	Pclass(structure *CT,pfdatatable *DATA);
	~Pclass();
	
	inline void limitdist();//This function handles the case where base pairs are not
					//not allowed to form between nucs more distant
					//than ct->maxdistance
	void store(char *save);//This function will save the arrays in a binary file

	//-------------------------------------------------------------------------------------------------
	//old fill routine with O(N4) for unlimited internal loops
	void fillw3();//This is a inline function to calc w3 in the old fill routine (O(N4) in time)
	void oldfill();//This is a old fill routine,O(N4) for unlimited internall loop
	void oldpartition(bool quickQ,PFPRECISION *Q=NULL,ProgressHandler* update=NULL,char *save=NULL);
					//This is partition function using oldfill()
	
	//--------------------------------------------------------------------------------------------------
	//improved fill routine with O(N3) for unlimited internal loops
	void fillw5();//This is a inline function to calc w5 in the improved fill routine (O(N3) in time)
	void interprefill();//This function prefill curE and prevE for internal loops's energy for the improved algorithm
	void fill();//This function fill the array for a partition function of the sequence
	void partition( bool quickQ,PFPRECISION *Q=NULL,ProgressHandler* update=NULL,char *save=NULL);
					//If quickQ == true, return the partition function value in Q
					//(Q is scaled in case overflowed, not a real value)
					//otherwise, save the partial partition functions
					//in the datafile named save
					//If updates on progress are unwanted, set update=NULL
};



//=======================================================================
/*
OligoPclass is inherited from Pclass
with additional functions for reusing the partition function array in the
OligoWalk refolding.

*/
class OligoPclass:public Pclass {

public:
	PFPRECISION *copyw5,**copywca;
	DynProgArray<PFPRECISION> *copyw,*copyv,*copywmb,*copywl,*copywmbl,*copywcoax;
			//These arrays and classes store the informations to be copied in refolding
	OligoPclass(structure *CT, pfdatatable *DATA);//initiate the copy arrays
	~OligoPclass();

	//----------------------------------------------------------------------------------------------
	//folding the whole length of target sequence	
	void partition4refill(PFPRECISION *Q,char *save=NULL);
			//normal partition function plus a copy of the arrays of v,w,w5 ...
	void refill(structure *CT,PFPRECISION *Qc, int start, int stop,PFPRECISION &rescaleinrefill,char *save=NULL);
			//recalculate the partionfunction for constrained sequence
			//assisted with copying the stored information in copyw,copyv ...	
			
	//----------------------------------------------------------------------------------------------
	//refolding at different region site on target
	void scanfill(structure *CT,PFPRECISION *Q,int reverse=0,char *save=NULL);//refold different region on target
	void scanconstrain(structure *CT,PFPRECISION *Qc,int start, int stop,PFPRECISION &rescaleinrefill,char *save=NULL);
																//refold the region on target with constrain
	
	void reset4oligo(structure *CT);//reset the arrays to be refilled again for other sequence
};


#endif


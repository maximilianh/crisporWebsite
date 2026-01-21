// OligoObject.cpp: implementation of the COligoObject class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "RNAstructure.h"
#include "OligoObject.h"
#include "../src/defines.h"
#include "../src/algorithm.h"
#include "globals.h"
#include "../src/rna_library.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

COligoObject::COligoObject(char *Datapath,char *Startpath,bool *Isdna, bool *Istargetdna, int *Option, 
		 int *Length,double *C,
		int *Usesub)
{

	strcpy(datapath,Datapath);
	strcpy(startpath,Startpath);
	isdna = *Isdna;
	istargetdna = *Istargetdna;
	option = *Option;
	
	length = *Length;
	c = *C;
	usesub = *Usesub;
	allocated = false;
	
	siRNA = false;//default is not siRNA design
	T = (float) 310.15;
	

}

COligoObject::~COligoObject()
{

	int i;


	if (allocated) {

		
		for (i = 0; i < ct.GetSequenceLength() - length + 2; i++) {
   			delete[] table[i];
			delete[] numofsubstructures[i];
		}
		delete[] table;
		delete[] numofsubstructures;

		
		
		if (isdna) {
			delete ddata;
			
			//delete hybriddata;
		
		
		}
		delete helixstack;

		if (siRNA) delete[] mask;

	}

}


bool COligoObject::AllocateTable() {
	int i,j,k,l;
	char loop2[maxfil], stackf[maxfil], tstackh[maxfil], tstacki[maxfil],
   		tloop[maxfil], miscloop[maxfil], danglef[maxfil], int22[maxfil],
		int21[maxfil], triloop[maxfil], coax[maxfil], tstackcoax[maxfil],
		tstackm[maxfil], coaxstack[maxfil], tstack[maxfil],
		int11[maxfil],hexaloop[maxfil],tstacki23[maxfil],tstacki1n[maxfil];
	rddata *enthalpyhybrid;//Needed for DNA oligos at temps other than 37 degrees C.

	allocated=true;

	//check if stop needs to be reduced for the length of the oligos
	if (stop>(ct.GetSequenceLength()-length+1)) stop = ct.GetSequenceLength()-length+1;



	table = new int*[ct.GetSequenceLength() - length + 2];

	for (i = 0; i < ct.GetSequenceLength() - length + 2; i++) {
		if (siRNA) table[i] = new int[7];
   		else table[i] = new int[6];
	}

	//allocate memory of number of suboptimal structures
	numofsubstructures= new int*[ct.GetSequenceLength() - length +2];
	for (i = 0; i < ct.GetSequenceLength() - length + 2; i++)	{
		numofsubstructures[i]= new int [2];
		numofsubstructures[i][0]=0;
		numofsubstructures[i][1]=0;
	}

	if (isdna) {
		ddata = new datatable();
		hybriddata = new rddata;


	}
	helixstack = new thermo(datapath);
	


	getdat (loop2,stackf,tstackh,tstacki,tloop,miscloop,danglef,
   	int22,int21,coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,hexaloop,tstacki23, tstacki1n, datapath,true);

	if (opendat (loop2,stackf,tstackh,tstacki,tloop,miscloop,danglef,int22,int21,
   	coax,tstackcoax,coaxstack,tstack, tstackm, triloop, int11, hexaloop, tstacki23, tstacki1n, &data)==0) {



      
	  return false;

	}

	getenthalpydat(loop2, stackf, tstackh, tstacki,
			tloop, miscloop, danglef, int22,
			int21,coax, tstackcoax,
			coaxstack, tstack, tstackm, triloop,
			int11, hexaloop, tstacki23, tstacki1n, datapath, true);//true = is an RNA

	

	if (opendat(loop2, stackf, tstackh, tstacki,
			tloop, miscloop, danglef, int22,
			int21,coax, tstackcoax,
			coaxstack, tstack, tstackm, triloop,
			int11,hexaloop,tstacki23, tstacki1n, &dhdata)==0) {
	
				
				return false;

	}

	if (T<310||T>311) {

		//change the temperature from 310.15
		if (newtemp(true,&data)==0) {
			//if newtemp returned zero, pass a warning to the user
			return false;

		}

	}




	if (isdna) {
		

		getdat (loop2,stackf,tstackh,tstacki,tloop,miscloop,danglef,
   	      	int22,int21,coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,hexaloop,
			tstacki23, tstacki1n,datapath,false);

		if (opendat (loop2,stackf,tstackh,tstacki,tloop,miscloop,danglef,int22,int21,
   	      	coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,hexaloop,tstacki23,
			tstacki1n,ddata)==0) {
       	
		
			
			return false;
		}

		if (T<310||T>311) {

			denthalpy = new datatable();

			getenthalpydat(loop2, stackf, tstackh, tstacki,
				tloop, miscloop, danglef, int22,
				int21,coax, tstackcoax,
				coaxstack, tstack, tstackm, triloop,
				int11, hexaloop, tstacki23, tstacki1n, datapath, false);

	

			if (opendat(loop2, stackf, tstackh, tstacki,
					tloop, miscloop, danglef, int22,
					int21,coax, tstackcoax,
					coaxstack, tstack, tstackm, triloop,
					int11,hexaloop,tstacki23, tstacki1n, denthalpy)==0) {
			
						
						return false;

			}

			//change the temperature from 310.15
			if (newtemp(false,ddata)==0) {
				//if newtemp returned zero, pass a warning to the user
				return false;

			}

			delete denthalpy;

		}

		strcpy(stackf,datapath);
		strcat(stackf,"\\");
		strcat (stackf,"stackdr.dat");
		if (readrd (hybriddata,stackf)==0) {
      			
			
			return false;

		}

		if (T<310||T>311) {
			//The temperature has changed.
			//Read the enthalpy data into a rddata.
			
			strcpy(stackf,datapath);
			strcat(stackf,"\\");
			strcat (stackf,"stackdr.dh");
			enthalpyhybrid = new rddata;
			if (readrd (enthalpyhybrid,stackf)==0) {
      					
				return false;

			}

			for (i=0;i<5;i++) {
				for (j=0;j<5;j++) {
					for (k=0;k<5;k++) {
						for (l=0;l<5;l++) {
							hybriddata->stack[i][j][k][l]=Tscale(T,hybriddata->stack[i][j][k][l],enthalpyhybrid->stack[i][j][k][l]);
						}
					}
				}
			}
			hybriddata->init=Tscale(T,hybriddata->init,enthalpyhybrid->init);
			delete enthalpyhybrid;

		}

		strcpy(helixstack->DH,datapath);
		strcat(helixstack->DH,"\\");
		strcat(helixstack->DH,"stackdr.dh");

		strcpy(helixstack->DS,datapath);
		strcat(helixstack->DS,"\\");
		strcat(helixstack->DS,"stackdr.ds");

		strcpy(helixstack->HELIX,datapath);
		strcat(helixstack->HELIX,"\\");
		strcat(helixstack->HELIX,"helixdr.dat");
		
	}
	
	if (helixstack->read()==0) {
      	
		return false;
	}

	if (siRNA) {
		//mask will store whether an oligo meets siRNA design criteria
		mask = new bool [ct.GetSequenceLength() - length + 2];

	}
	else mask = NULL;
	

	return true;

}



int COligoObject::newtemp(bool isrna, datatable *freeenergy) {
	//The temperature of folding was changed from 37 degrees C, alter the data files:
	

	if (isrna) dG_T(T,*freeenergy,dhdata,*freeenergy);//isrna==true -> RNA, so use dhdata for the enthalpy
	else dG_T(T,*freeenergy,*denthalpy,*freeenergy);//isrna==false -> DNA, so use denthalpy for the enthalpy


	return 1;

}






	


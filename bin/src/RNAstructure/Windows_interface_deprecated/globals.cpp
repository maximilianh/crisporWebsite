//This file contains the global interface functions 

#include "stdafx.h"
#include "error.h"

void errmsg(int a,int b) 
{

		
	CError error;
	error.DoModal();

}

/////////////////////////////////////////////////////////////////////////////
// CRNAstructureApp message handlers




/*	Function getdat

	Function gets the names of data files to open

*/

void getdat(char *loop, char *stackf, char *tstackh, char *tstacki,
	char *tloop, char *miscloop, char *danglef, char *int22,
	char *int21,char *coax, char *tstackcoax,
    char *coaxstack, char *tstack, char *tstackm, char *triloop,
    char *int11, char *hexaloop, char *tstacki23, char *tstacki1n,
	char *datapath, bool isRNA)

{

strcpy (loop,datapath);

strcpy (stackf,datapath);

strcpy (tstackh,datapath);

strcpy (tstacki,datapath);

strcpy (tloop,datapath);

strcpy (miscloop,datapath);

strcpy (danglef,datapath);

strcpy (int22,datapath);

strcpy (int21,datapath);

strcpy (triloop,datapath);

strcpy (coax,datapath);

strcpy (tstackcoax,datapath);

strcpy (coaxstack,datapath);

strcpy (tstack,datapath);

strcpy (tstackm,datapath);

strcpy (int11,datapath);

strcpy (hexaloop,datapath);

strcpy (tstacki23,datapath);

strcpy (tstacki1n,datapath);



strcat (loop,"\\");
strcat (stackf,"\\");
strcat (tstackh,"\\");
strcat (tstacki,"\\");
strcat (tloop,"\\");
strcat (miscloop,"\\");
strcat (danglef,"\\");
strcat (int22,"\\");
strcat (int21,"\\");
strcat (triloop,"\\");
strcat (coax,"\\");
strcat (tstackcoax,"\\");
strcat (coaxstack,"\\");
strcat (tstack,"\\");
strcat (tstackm,"\\");
strcat (int11,"\\");
strcat (hexaloop,"\\");
strcat (tstacki23,"\\");
strcat (tstacki1n,"\\");

if (!isRNA) {
	//load DNA parameters


	strcat (loop,"dna");
	strcat (stackf,"dna");
	strcat (tstackh,"dna");
	strcat (tstacki,"dna");
	strcat (tloop,"dna");
	strcat (miscloop,"dna");
	strcat (danglef,"dna");
	strcat (int22,"dna");
	strcat (int21,"dna");
	strcat (triloop,"dna");
	strcat (coax,"dna");
	strcat (tstackcoax,"dna");
	strcat (coaxstack,"dna");
	strcat (tstack,"dna");
	strcat (tstackm,"dna");
	strcat (int11,"dna");
	strcat (hexaloop,"dna");
	strcat (tstacki23,"dna");
	strcat (tstacki1n,"dna");

}

strcat (loop,"loop.dat");
strcat (stackf,"stack.dat");
strcat (tstackh,"tstackh.dat");
strcat (tstacki,"tstacki.dat");
strcat (tloop,"tloop.dat");
strcat (miscloop,"miscloop.dat");
strcat (danglef,"dangle.dat");
strcat (int22,"int22.dat");
strcat (int21,"int21.dat");
strcat (triloop,"triloop.dat");
strcat (coax,"coaxial.dat");
strcat (tstackcoax,"tstackcoax.dat");
strcat (coaxstack,"coaxstack.dat");
strcat (tstack,"tstack.dat");
strcat (tstackm,"tstackm.dat");
strcat (int11,"int11.dat");
strcat (hexaloop,"hexaloop.dat");
strcat (tstacki23,"tstacki23.dat");
strcat (tstacki1n,"tstacki1n.dat");
}

/*	Function getdat

	Function gets the names of enthalpy data files to open

*/

void getenthalpydat(char *loop, char *stackf, char *tstackh, char *tstacki,
	char *tloop, char *miscloop, char *danglef, char *int22,
	char *int21,char *coax, char *tstackcoax,
    char *coaxstack, char *tstack, char *tstackm, char *triloop,
    char *int11, char *hexaloop, char *tstacki23, char *tstacki1n,
	char *datapath, bool isRNA)

{

strcpy (loop,datapath);

strcpy (stackf,datapath);

strcpy (tstackh,datapath);

strcpy (tstacki,datapath);

strcpy (tloop,datapath);

strcpy (miscloop,datapath);

strcpy (danglef,datapath);

strcpy (int22,datapath);

strcpy (int21,datapath);

strcpy (triloop,datapath);

strcpy (coax,datapath);

strcpy (tstackcoax,datapath);

strcpy (coaxstack,datapath);

strcpy (tstack,datapath);

strcpy (tstackm,datapath);

strcpy (int11,datapath);

strcpy (hexaloop,datapath);

strcpy (tstacki23,datapath);

strcpy (tstacki1n,datapath);



strcat (loop,"\\");
strcat (stackf,"\\");
strcat (tstackh,"\\");
strcat (tstacki,"\\");
strcat (tloop,"\\");
strcat (miscloop,"\\");
strcat (danglef,"\\");
strcat (int22,"\\");
strcat (int21,"\\");
strcat (triloop,"\\");
strcat (coax,"\\");
strcat (tstackcoax,"\\");
strcat (coaxstack,"\\");
strcat (tstack,"\\");
strcat (tstackm,"\\");
strcat (int11,"\\");
strcat (hexaloop,"\\");
strcat (tstacki23,"\\");
strcat (tstacki1n,"\\");

if (!isRNA) {
	//load DNA parameters


	strcat (loop,"dna");
	strcat (stackf,"dna");
	strcat (tstackh,"dna");
	strcat (tstacki,"dna");
	strcat (tloop,"dna");
	strcat (miscloop,"dna");
	strcat (danglef,"dna");
	strcat (int22,"dna");
	strcat (int21,"dna");
	strcat (triloop,"dna");
	strcat (coax,"dna");
	strcat (tstackcoax,"dna");
	strcat (coaxstack,"dna");
	strcat (tstack,"dna");
	strcat (tstackm,"dna");
	strcat (int11,"dna");
	strcat (hexaloop,"dna");
	strcat (tstacki23,"dna");
	strcat (tstacki1n,"dna");

}

strcat (loop,"loop.dh");
strcat (stackf,"stack.dh");
strcat (tstackh,"tstackh.dh");
strcat (tstacki,"tstacki.dh");
strcat (tloop,"tloop.dh");
strcat (miscloop,"miscloop.dh");
strcat (danglef,"dangle.dh");
strcat (int22,"int22.dh");
strcat (int21,"int21.dh");
strcat (triloop,"triloop.dh");
strcat (coax,"coaxial.dh");
strcat (tstackcoax,"tstackcoax.dh");
strcat (coaxstack,"coaxstack.dh");
strcat (tstack,"tstack.dh");
strcat (tstackm,"tstackm.dh");
strcat (int11,"int11.dh");
strcat (hexaloop,"hexaloop.dh");
strcat (tstacki23,"tstacki23.dh");
strcat (tstacki1n,"tstacki1n.dh");
}


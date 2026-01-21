//This file contains the global interface functions 

#include <iostream>
#include <cstdlib>
#include <cstring>
using namespace std;

/////////////////////////////////////////////////////////////////////////////
// CRNAstructureApp message handlers




/*	Function getdat

	Function gets the names of data files to open

*/

void getdat(char *fspec, char *loop, char *stackf, char *tstackh, char *tstacki,
	char *tloop, char *miscloop, char *danglef, char *int22,
	char *int21,char *coax, char *tstackcoax,
    char *coaxstack, char *tstack, char *tstackm, char *triloop,
    char *int11, char *hexaloop, char *tstacki23, char *tstacki1n,
	const char *datapath, bool isRNA, bool isEnthalpy)

{

//copy the path to each name
  strcpy(fspec,datapath);
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
  


  if( !isRNA) {
	  //these are dna parameters and so they need to start with "dna"
	  strcat (fspec,"dna.");
	strcat (loop,"dna.");
	strcat (stackf,"dna.");
	  strcat (tstackh,"dna.");
	  strcat (tstacki,"dna.");
	  strcat (tloop,"dna.");
	  strcat (miscloop,"dna.");
	  strcat (danglef,"dna.");
	  strcat (int22,"dna.");
	  strcat (int21,"dna.");
	  strcat (triloop,"dna.");
	  strcat (coax,"dna.");
	  strcat (tstackcoax,"dna.");
	  strcat (coaxstack,"dna.");
	  strcat (tstack,"dna.");
	  strcat (tstackm,"dna.");
	  strcat (int11,"dna.");
	  strcat (hexaloop,"dna.");
	  strcat (tstacki23,"dna.");
	  strcat (tstacki1n,"dna.");  
	 
  }
  else { //Right now, this means that the files are for RNA
	  //these are dna parameters and so they need to start with "dna"
	  strcat (fspec,"rna.");
	strcat (loop,"rna.");
	strcat (stackf,"rna.");
	  strcat (tstackh,"rna.");
	  strcat (tstacki,"rna.");
	  strcat (tloop,"rna.");
	  strcat (miscloop,"rna.");
	  strcat (danglef,"rna.");
	  strcat (int22,"rna.");
	  strcat (int21,"rna.");
	  strcat (triloop,"rna.");
	  strcat (coax,"rna.");
	  strcat (tstackcoax,"rna.");
	  strcat (coaxstack,"rna.");
	  strcat (tstack,"rna.");
	  strcat (tstackm,"rna.");
	  strcat (int11,"rna.");
	  strcat (hexaloop,"rna.");
	  strcat (tstacki23,"rna.");
	  strcat (tstacki1n,"rna.");  
	 
  }
	 
	//append the actual file name
  strcat(fspec, "specification.dat");//include full extension here, it is the same for DG and DH
  strcat (loop,"loop.");
  strcat (stackf,"stack.");
  strcat (tstackh,"tstackh.");
  strcat (tstacki,"tstacki.");
  strcat (tloop,"tloop.");
  strcat (miscloop,"miscloop.");
  strcat (danglef,"dangle.");
  strcat (int22,"int22.");
  strcat (int21,"int21.");
  strcat (triloop,"triloop.");
  strcat (coax,"coaxial.");
  strcat (tstackcoax,"tstackcoax.");
  strcat (coaxstack,"coaxstack.");
  strcat (tstack,"tstack.");
  strcat (tstackm,"tstackm.");
  strcat (int11,"int11.");
  strcat (hexaloop,"hexaloop.");
  strcat (tstacki23,"tstacki23.");
  strcat (tstacki1n,"tstacki1n.");
  
  if (isEnthalpy) {
	  //thse are enthalpy partameters so they need to end in .dh
	strcat (loop,"dh");
	strcat (stackf,"dh");
	  strcat (tstackh,"dh");
	  strcat (tstacki,"dh");
	  strcat (tloop,"dh");
	  strcat (miscloop,"dh");
	  strcat (danglef,"dh");
	  strcat (int22,"dh");
	  strcat (int21,"dh");
	  strcat (triloop,"dh");
	  strcat (coax,"dh");
	  strcat (tstackcoax,"dh");
	  strcat (coaxstack,"dh");
	  strcat (tstack,"dh");
	  strcat (tstackm,"dh");
	  strcat (int11,"dh");
	  strcat (hexaloop,"dh");
	  strcat (tstacki23,"dh");
	  strcat (tstacki1n,"dh");  
  }
  else {
	  //these are free energy parameters and the files end in .dat
	strcat (loop,"dg");
	strcat (stackf,"dg");
	  strcat (tstackh,"dg");
	  strcat (tstacki,"dg");
	  strcat (tloop,"dg");
	  strcat (miscloop,"dg");
	  strcat (danglef,"dg");
	  strcat (int22,"dg");
	  strcat (int21,"dg");
	  strcat (triloop,"dg");
	  strcat (coax,"dg");
	  strcat (tstackcoax,"dg");
	  strcat (coaxstack,"dg");
	  strcat (tstack,"dg");
	  strcat (tstackm,"dg");
	  strcat (int11,"dg");
	  strcat (hexaloop,"dg");
	  strcat (tstacki23,"dg");
	  strcat (tstacki1n,"dg");  


  }
}

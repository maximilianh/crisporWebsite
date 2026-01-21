// napss_dot a single RNA sequence.  Programmed by DHM, 8/17/05

//This has been modified to output dot plot information 11/26/07 in conjunction with a change to algorithm.cpp

#include "stdafx.h"
#include "platform.h"
#include "rna_library.h"
#include "structure.h"
#include "algorithm.h"
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
using namespace std;


void getdat(char *loop, char *stackf, char *tstackh, char *tstacki,
	char *tloop, char *miscloop, char *danglef, char *int22,
	char *int21,char *coax, char *tstackcoax,
    char *coaxstack, char *tstack, char *tstackm, char *triloop,
    char *int11, char *hexaloop, char *tstacki23, char *tstacki1n,
	char *datapath, bool isRNA);

int main(int argc, char* argv[])
{

	if (argc!=3) {
		
		cerr << "Usage: napss_dot sequence_in dotplot_out\n";
		exit(0);

	}
	char loop[maxfil],stackf[maxfil],tstackh[maxfil],tstacki[maxfil],
		tloop[maxfil],miscloop[maxfil],danglef[maxfil],int22[maxfil],
		int21[maxfil],coax[maxfil],tstackcoax[maxfil],
		coaxstack[maxfil],tstack[maxfil],tstackm[maxfil],triloop[maxfil],int11[maxfil],hexaloop[maxfil],
		tstacki23[maxfil], tstacki1n[maxfil],datapath[maxfil];const char *pointer;

	int maxtracebacks,percent,windowsize;
	
   
	structure ct;
	datatable data;

	int ir,i;

	maxtracebacks = 0;
	percent = 0;
	windowsize = 0;

	pointer = getDataPath();
	if (pointer!=NULL) strcpy(datapath,pointer);
	else strcpy(datapath,"");

//  strcpy(datapath,"../data/");

	//open the data files -- must reside with executable
	//open the thermodynamic data tables
	getdat (loop, stackf, tstackh, tstacki,tloop, miscloop, danglef, int22,
          int21,coax, tstackcoax,coaxstack, tstack, tstackm, triloop,
          int11, hexaloop, tstacki23, tstacki1n, datapath, true);//the true indicates RNA parameters
	if (opendat (loop,stackf,tstackh,tstacki,tloop,miscloop,danglef,int22,int21,
   		coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,hexaloop,tstacki23, tstacki1n,&data)==0) {
      	
		cout << "A data file was lost\n";
	}


	//open the sequence
	openseq(&ct, argv[1]);

	//now do the structure prediction and output the dotplot:
	dynamic(&ct,&data,maxtracebacks,percent,windowsize,0,false,argv[2]);

	return 0;
}

void errmsg(int err,int erri) {

if (err==30) {
	cout << "End Reached at traceback #"<<erri<<"\n";
   exit(1);
}
if (err==100) {
	cout << "error # "<<erri;
   exit(1);
}
switch (err) {
	case 1:
   	cout << "Could not allocate enough memory";
      break;
   case 2:
   	cout << "Too many possible base pairs";
      break;
   case 3:
   	cout << "Too many helixes in multibranch loop";
   case 4:
   	cout << "Too many structures in CT file";
   default:
   	cout << "Unknown error";
}
cin >> err;
exit(1);
return;

}

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

 
  if( !isRNA) strcat( datapath,"dna");
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



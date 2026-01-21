/************************************************************************
* partition_inter.cpp:							*
* Usage: foo.o  [lisfile] [extension]					*
* input: .lis .seq .dat .dh						*
* output: .pfs								*
* Calculate partition funtion with N3 instead of N4 for  unlimited	*
* internal loop; using algorithm_inter.cpp file instead of algorithm.cpp*
* based on partition_tem.cpp:						*
* Calculate partition function at certain temperature with 		*
* dG_tem.cpp file:							*
* dg = dh - ds*T= dh - (dh-dg37)*T/310.15				*   	
* temperatures are read from lis files                                  *   
* dh are read from .dh files and stored in Class dhdata			*
* dg37 are read from .dat files and stored in Class data		*
* dg are calculated and are stored in Class dg				*
*************************************************************************
* Created:	2003;			Modified:	Dec.2005	*
* Copyright: Zhi John Lu & David Mathews				*	
*************************************************************************/

#include "structure.cpp"
#include "rna_library.cpp"
#include "algorithm.cpp"
#include "pfunction_inter.cpp"
#include "dG_tem.cpp"
#include "../src/phmm/utils/xmath/log/xlog_math.h"

#define maxsequences 2000
#define path "../archive/"
//#define path ""


void getdat(char *loop, char *stackf, char *tstackh, char *tstacki,
		char *tloop, char *miscloop, char *danglef, char *int22,
      char *int21,char *coax, char *tstackcoax,
      char *coaxstack, char *tstack,char *tstackm,char *triloop, char *int11, char *hexaloop,
	  char *tstacki23, char *tstacki1n);

void getdh(char *loop, char *stackf, char *tstackh, char *tstacki,
                char *tloop, char *miscloop, char *danglef, char *int22,
		      char *int21,char *coax, char *tstackcoax,
		            char *coaxstack, char *tstack,char *tstackm,char *triloop, char *int11, char *hexaloop,
			              char *tstacki23, char *tstacki1n);

int main(int argc, char *argv[]) {


	char loop[maxfil],stackf[maxfil],tstackh[maxfil],tstacki[maxfil],
		tloop[maxfil],miscloop[maxfil],danglef[maxfil],int22[maxfil],
      int21[maxfil],coax[maxfil],tstackcoax[maxfil],
      coaxstack[maxfil],tstack[maxfil],tstackm[maxfil],triloop[maxfil],int11[maxfil],hexaloop[maxfil],
	  tstacki23[maxfil], tstacki1n[maxfil];
   
   datatable *data,*dhdata,*dg;
   structure *ct;
   //TProgressDialog PD;

 	char list[maxfil],outfile[maxfil],label[maxfil],prefix[maxsequences][maxfil];
   char sequence[maxfil],dsv[maxfil];
   
   
   register int i,ir;
   float temperature[maxsequences],T[maxsequences];// degree for temperature and K for T
   PFPRECISION pfTem[maxsequences];//input temperature for pftable
   data=new datatable();
   dhdata=new datatable();
   dg=new datatable();
   TProgressDialog PD;
   
   pfdatatable *pfdata;

   if (argc!=3) {
   	cout << "Usage: rnatotal [lisfile] [extension]\n";



   	cout << "Enter the name of a .lis file: \n";
   	cin >> list;


   	cout << "Enter an extension for this set of foldings:  \n";
   	cin >> label;
   }
   else {
   	  strcpy(list,argv[1]);
      strcpy(label,argv[2]);
   }


	//getinfo (&cntrl6,&cntrl8,&cntrl9);

   ifstream inf;
	inf.open(list);



   //read the list of structures to be folded from a list file:
	for (ir=1;ir<=maxsequences;ir++) {
		
		inf >> prefix[ir];
		inf >> temperature[ir];
		if(inf.eof()) break;
		T[ir]=temperature[ir]+273.15;
		cout << "ir = "<<ir<<"ct = "<<prefix[ir]<<"temp =  " << temperature[ir] << "T = " <<T[ir]<<"\n"<<flush;
		
	}

	ir--;

   getdat (loop,stackf,tstackh,tstacki,tloop,miscloop,danglef,
   	int22,int21,coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,hexaloop,tstacki23, tstacki1n);
	if (opendat (loop,stackf,tstackh,tstacki,tloop,miscloop,danglef,int22,int21,
   	coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,hexaloop,tstacki23, tstacki1n,data)==0) {
      	cout << "A data file was lost";
         cin >> i;
   }
   getdh (loop,stackf,tstackh,tstacki,tloop,miscloop,danglef,
 
           int22,int21,coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,hexaloop,tstacki23, tstacki1n);
	if (opendat (loop,stackf,tstackh,tstacki,tloop,miscloop,danglef,int22,int21,
        coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,hexaloop,tstacki23, tstacki1n,dhdata)==0) {
	cout << "A dhdata file was lost";
	 cin >> i;
   }



   for (i=1;i<=ir;i++) {//fold and score every structure


   	ct = new structure;
      

      strcpy(sequence,path);
      strcat(sequence,prefix[i]);
      strcat(sequence,".seq");
    cout << i << "  "<<sequence<<"\n"<<flush;

     
      strcpy(dsv,prefix[i]);
      strcat(dsv,label);
      strcat(dsv,".pfs");
      

      openseq (ct,sequence);


	dG_T(T[i], *data, *dhdata, *dg);
	pfTem[i]=T[i];
	pfdata = new pfdatatable(dg,scalingdefinition,pfTem[i]);
	pfunction(ct,pfdata, &PD, dsv);

    

   	delete ct;
        delete pfdata;
	}

   delete data;
   delete dhdata;
   delete dg;
   return 0;

}





/*	Function getdat

	Function gets the names of data files to open

*/

void getdat(char *loop, char *stackf, char *tstackh, char *tstacki,
		char *tloop, char *miscloop, char *danglef, char *int22,
      char *int21,char *coax, char *tstackcoax,
      char *coaxstack, char *tstack, char *tstackm,char *triloop, char *int11,
	  char *hexaloop,char *tstacki23, char *tstacki1n )

{
strcpy (loop,"loop.dat");
strcpy (stackf,"stack.dat");
strcpy (tstackh,"tstackh.dat");
strcpy (tstacki,"tstacki.dat");
strcpy (tloop,"tloop.dat");
strcpy (miscloop,"miscloop.dat");
strcpy (danglef,"dangle.dat");
strcpy (int22,"int22.dat");
strcpy (int21,"int21.dat");
strcpy (coax,"coaxial.dat");
strcpy (tstackcoax,"tstackcoax.dat");
strcpy (coaxstack,"coaxstack.dat");
strcpy (tstack,"tstack.dat");
strcpy (tstackm,"tstackm.dat");
strcpy (triloop,"triloop.dat");
strcpy (int11,"int11.dat");
strcpy (hexaloop,"hexaloop.dat");
strcpy (tstacki23,"tstacki23.dat");
strcpy (tstacki1n,"tstacki1n.dat");
}

/*      Function getdh      
Function gets the names of dhdata files to open */

void getdh(char *loop, char *stackf, char *tstackh, char *tstacki,char *tloop, char *miscloop, char *danglef, char *int22,char *int21,char *coax, char *tstackcoax,char *coaxstack, char *tstack, char *tstackm,char *triloop, char *int11,char *hexaloop,char *tstacki23, char *tstacki1n )
{
strcpy (loop,"loop.dh");
strcpy (stackf,"stack.dh");
strcpy (tstackh,"tstackh.dh");
strcpy (tstacki,"tstacki.dh");
strcpy (tloop,"tloop.dh");
strcpy (miscloop,"miscloop.dh");
strcpy (danglef,"dangle.dh");
strcpy (int22,"int22.dh");
strcpy (int21,"int21.dh");
strcpy (coax,"coaxial.dh");
strcpy (tstackcoax,"tstackcoax.dh");
strcpy (coaxstack,"coaxstack.dh");
strcpy (tstack,"tstack.dh");
strcpy (tstackm,"tstackm.dh");
strcpy (triloop,"triloop.dh");
strcpy (int11,"int11.dh");
strcpy (hexaloop,"hexaloop.dh");
strcpy (tstacki23,"tstacki23.dh");
strcpy (tstacki1n,"tstacki1n.dh");
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

void update(int i) {

//	cout<< i<<"\n"<<flush;
}






# include <math.h>
//# include <stdlib.h>
#include <cstring>
#include <cstdlib>
#include "algorithm.h"
#include "thermo.h"

using namespace std;

thermo::thermo(const string& path) {
   short int i,j,k,l;
	//initailize the values in the tables
   for (i=0;i<5;i++) {
   	for (j=0;j<5;j++) {
      	for (k=0;k<5;k++) {
         	for (l=0;l<5;l++) {

      			dh[i][j][k][l] = 0;
         		ds[i][j][k][l] = 0;
            }
         }
      }
   }
   if (path.empty()) {
   	//set default names for the data files
   	DH="stack.dh";
   	DS="stack.ds";
      HELIX="helix.dat";
   } else {
      DH=path+"/stack.dh";
      DS=path+"/stack.ds";
      HELIX=path+"/helix.dat";
   }
};


//Throws an exception (string error message) if the datatables could not be read.
int thermo::read() {
   ifstream dhf(DH.c_str()),dsf(DS.c_str()),hf(HELIX.c_str());
   char lineoftext[100];
   short int i,j,k,l,count;

   //check that all the files exist with a C i/o function
	if ( !(dhf.good()&&dsf.good()&&hf.good()) ) return 0;

/* Read info from stack.dh */
//add to the stack table the case where X (represented as 0) is looked up:

for (count=1;count<=42;count++) dhf >> lineoftext;//get past text in file
for (i=1;i<=4;i++) {
	for (count=1;count<=60;count++) dhf >> lineoftext;
	for (k=1;k<=4;k++) {
		for (j=1;j<=4;j++) {
			for (l=1;l<=4;l++) {
					dhf >> lineoftext;
					if (strcmp(lineoftext,".")){
						dh[i][j][k][l] =(short) floor(conversionfactor*(atof(lineoftext))+.5);
					}
					else dh[i][j][k][l] = INFINITE_ENERGY;


			}

		}
	}
}

/* Read info from stack.ds */
//add to the stack table the case where X (represented as 0) is looked up:

for (count=1;count<=42;count++) dsf >> lineoftext;//get past text in file
for (i=1;i<=4;i++) {
	for (count=1;count<=60;count++) dsf >> lineoftext;
	for (k=1;k<=4;k++) {
		for (j=1;j<=4;j++) {
			for (l=1;l<=4;l++) {
					dsf >> lineoftext;
					if (strcmp(lineoftext,".")){
						ds[i][j][k][l] =(short) floor(conversionfactor*(atof(lineoftext))+.5);
					}
					else ds[i][j][k][l] = INFINITE_ENERGY;


			}

		}
	}
}

//read info from helix.dat
for (count = 1;count <=4;count++) hf >> lineoftext;
hf >> lineoftext;
dhi = (short) floor(conversionfactor*(atof(lineoftext))+.5);
for (count = 1;count <=4;count++) hf >> lineoftext;
hf >> lineoftext;
dsi = (short) floor(conversionfactor*(atof(lineoftext))+.5);
for (count = 1;count <=5;count++) hf >> lineoftext;
hf >> lineoftext;
dss = (short) floor(conversionfactor*(atof(lineoftext))+.5);
for (count = 1;count <=5;count++) hf >> lineoftext;
hf >> lineoftext;
dha = (short) floor(conversionfactor*(atof(lineoftext))+.5);
for (count = 1;count <=5;count++) hf >> lineoftext;
hf >> lineoftext;
dsa = (short) floor(conversionfactor*(atof(lineoftext))+.5);

return 1;

}

//Mix&match, copyright 1996 by David Mathews
//	Programmed July 1996, based on the FORTRAN Program Selector,
//	also by David Mathews (written August 1995)


//Program reads a CT file with multiple sub-optimal structures
//Outputs a CT file with one structure, with the lowest free energy possible
//	by recombination of all subdomains



# include <math.h>



# define maxmod 500
//#define maxtloop 200
//#define infinity 16000
//#define maxstructures 1010 //maximum number of structures in ct file
//#define maxbases 3000   //maximum number of bases in a structure
//#define ctheaderlength 100 //maximum length of string containing info on sequence
//#define maxnopair 300 //maximum number of bases that can be forced single


#include "interface.h"
#include "structure.h"
#include "algorithm.h"
#include "rna_library.h"
using namespace std;


#define maxdepth 50
struct stack {

	int i[maxdepth],j[maxdepth],p;
   stack() {
    	p = 0;
   }
   void put(int ip,int jp) {
    	p++;
      i[p] = ip;
      j[p] = jp;
   }
   int get(int *ip,int *jp) {
      if (p==0) return 0;

    	(*ip) = i[p];
      (*jp) = j[p];
      p--;
      return 1;
   }



};




//void efn2(datatable *data,structure *ct); //found in algorithm.cpp



//twoefn2 is a copy of efn2, taylored for mix&match to fill an array


void twoefn2(datatable *data,structure *ct,int i, int j, short int **energy,int count)
{
int k,sum,ip,jp,m,n;
//stackstruct stack;
//bool flag;
//short int n5,n3;

/*	stack = a place to keep track of where efn2 is located in a structure
	inter = indicates whether there is an intermolecular interaction involved in
   	a multi-branch loop
*/

//stack.sp = 0;  //set stack counter






   ct->energy[count]=0;
	//push(&stack,1,ct->numofbases,1,0);//put the whole structure onto the stack

	//subroutine://loop starts here (I chose the goto statement to speed operation)
	//pull(&stack,&i,&j,&open,&null,&stz);//take a substructure off the stack

   //cout << "pulled "<<i<<"  energy = "<<ct->energy[count]<<"\n";
   //cin >> k;






 		while(i<j) {
			while (ct->basepr[count][i+1]==j-1) {//are i,j and i+1,j-1 stacked?
				ct->energy[count] = ct->energy[count] + erg1(i,j,i+1,j-1,ct,data);
				i++;
				j--;
            energy[i][j] = ct->energy[count];
			}

			sum = 0;
			k = i + 1;
				// now efn2 is past the paired region, so define
				//		 the intervening non-paired segment

			while (k<j) {
				if (ct->basepr[count][k]>k)	{
					sum++;
					ip = k;
					k = ct->basepr[count][k] + 1;
					jp = k-1;
				}
				else if (ct->basepr[count][k]==0) k++;
			}

			if (sum==1) {//bulge/internal loop
				ct->energy[count] = ct->energy[count] +
            	erg2(i,j,ip,jp,ct,data,0,0);

            if ((ct->basepr[count][i+1]==0)&&(ct->basepr[count][j-1]==0)) {
             	energy[i+1][j-1] = energy[i][j]+data->tstkm[ct->numseq[i]]
               	[ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] +
                  2*data->eparam[6] + data->eparam[5] + data->eparam[10];

               for (m=i+2;m<ip;m++) {
               	for (n=j-2;n>jp;n--) {
                   	energy[m][n] = energy[i+1][j-1] +
                     	(m-i + j - 2 - n)*data->eparam[6];
                  }
               }

            }
            if (ct->basepr[count][i+1]==0) {
             	energy[i+1][j] = energy[i][j] + erg4(i,j,i+1,1,ct,data,false) +
               	data->eparam[6] + data->eparam[5] + data->eparam[10];

               for (m=i+2;m<ip;m++) {
                	energy[m][j] = energy[i+1][j] +
                  	data->eparam[6]*(m-i-1);
               }

            }
            if (ct->basepr[count][j-1]==0) {
             	energy[i][j-1] = energy[i][j] + erg4(i,j,j-1,2,ct,data,false) +
               	data->eparam[6] + data->eparam[5]+data->eparam[10];

               for (m=j-2;m>jp;m--) {
                	energy[i][m] = energy[i][j-1] + (j-m-1)*data->eparam[6];
               }
            }


				i = ip;
				j = jp;
            energy[i][j] = ct->energy[count];
			}
			else {//multi-branch loop or hairpin loop
          	if ((ct->basepr[count][i+1]==0)&&(ct->basepr[count][j-1]==0)) {
             	energy[i+1][j-1] = energy[i][j]+data->tstkm[ct->numseq[i]]
               	[ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] +
                  2*data->eparam[6] + data->eparam[5]+data->eparam[10];

               for (m=i+2;m<j-1;m++) {
               	for (n=j-2;n>m;n--) {
                   	energy[m][n] = energy[i+1][j-1] +
                     	(m-i + j - 2 - n)*data->eparam[6];
                  }
               }

            }
            if (ct->basepr[count][i+1]==0) {
             	energy[i+1][j] = energy[i][j] + erg4(i,j,i+1,1,ct,data,false) +
               	data->eparam[6] + data->eparam[5]+data->eparam[10];

               for (m=i+2;m<j-1;m++) {
                	energy[m][j] = energy[i+1][j] +
                  	data->eparam[6]*(m-i-1);
               }

            }
            if (ct->basepr[count][j-1]==0) {
             	energy[i][j-1] = energy[i][j] + erg4(i,j,j-1,2,ct,data,false) +
               	data->eparam[6] + data->eparam[5]+data->eparam[10];

               for (m=j-2;m>i+1;m--) {
                	energy[i][m] = energy[i][j-1] + (j-m-1)*data->eparam[6];
               }
            }



         	i = j+1;//this will break the loop
         }

      }

      return;


}



bool compat(structure *ct,int data,int n) {//this checks whether a structure
															//is compatable with the
                                             //modification data

//returns true if the structure is compatable with the modification data


//ct is a structure containing sequence information and base pairing
//			patterns for multiple sub-optimal structures
//data is a single modification location (ie: position of a modified base)
//n is the number of structure that is being checked in ct

if (ct->basepr[n][data]==0) return true;//the base isn't paired
else if (ct->basepr[n][data+1]==0) return true;//the neigboring base isn't paired
else if (ct->basepr[n][data-1]==0) return true;//the neighboring base isn't paired
else if (ct->basepr[n][data+1]!=(ct->basepr[n][data]-1)) return true;
			//the base is at the end of a helix because the next bases in the
         	//helix aren't paired to each other
else if (ct->basepr[n][data-1]!=(ct->basepr[n][data]+1)) return true;
			//the base is at the end of a helix because the next bases in the
         	//helix aren't paired to each other
else if ((ct->numseq[data]==3)&&(ct->numseq[ct->basepr[n][data]]==4))
	return true;	//the base is a G paired to a U
else if ((ct->numseq[data]==4)&&(ct->numseq[ct->basepr[n][data]]==3))
	return true; 	//base is a U paired to a G
else if ((ct->numseq[data+1]==3)&&(ct->numseq[ct->basepr[n][data+1]]==4))
	return true;	//the neighboring pair is a G-U
else if ((ct->numseq[data+1]==4)&&(ct->numseq[ct->basepr[n][data+1]]==3))
	return true; 	//the neighboring pair is a G-U
else if ((ct->numseq[data-1]==3)&&(ct->numseq[ct->basepr[n][data-1]]==4))
	return true;	//the neigboring pair is a G-U
else if ((ct->numseq[data-1]==4)&&(ct->numseq[ct->basepr[n][data-1]]==3))
	return true; 	//the neighboring pair is a G-U
else return false;	//the base is in a helix, therefore the structure is
								//incomapatble with the modification observed at data
}


int externalcompat(structure *ct,int i,int ip,int j,int jp,int k,int numofmod,int *mod){
	int a;

   //for (a=i;a<ip;a++) {
   //	if ((ct->basepr[k][a]<i)&&(ct->basepr[k][a]>0)) return 0;
   //   else if ((ct->basepr[k][a]>ip)&&(ct->basepr[k][a]<jp)) return 0;
   //   else if (ct->basepr[k][a]>j) return 0;
   //}
   //for (a=jp+1;a<=j;a++) {
   //	if ((ct->basepr[k][a]<i)&&(ct->basepr[k][a]>0)) return 0;
   //   else if ((ct->basepr[k][a]>ip)&&(ct->basepr[k][a]<jp)) return 0;
   //   else if (ct->basepr[k][a]>j) return 0;
   //}
   for (a=1;a<=numofmod;a++) {
    	if (mod[a] == i) return 1;
      else if (mod[a] == ip) return 1;
      else if (mod[a] == jp) return 2;
      else if (mod[a] == j) return 2;
      else if ((mod[a]>i)&&(mod[a]<ip)) {
			if (!compat(ct,mod[a],k)) return 1;
      }
      else if ((mod[a]<j)&&(mod[a]>jp)) {
      	if (!compat(ct,mod[a],k)) return 2;
      }
   }
   return 3;

}

void mixmatch (structure* ct,structure* ctout,datatable *data, char *datafile,
	ProgressHandler *update) {

//This is the main algorithm. It receives:
//	ct - a structure containing structural information for many sub-optimal
//				structures
//	ctout - a structure that will receive the structural information for the
//				minimum free energy compatable structure
//	data - a structure that contains all the free energy information that is
//				required by efn2

short int **subenergy;
//int **subct;
int i,j,n,mod[maxmod],numofmod,k,size,pos,ip,jp,l,ext;
structure substructure;
bool nopair;
short int **subcompat;
substructure.numofstructures=1;
stack st;


for (i=1;i<=ct->numofbases;i++) ctout->basepr[1][i] = 0;

substructure.allocate(ct->numofbases);

//ofstream out("dump");


//**subenergy = array that contains the lowest free energy found for
//		subsequence i to j from all suboptimal structures
//**subct = array that contains the number of the lowest free energy suboptimal
//					structure found on subsequence i,j
//**subcompat = array that indicates whether a specific suboptimal structure
//					is compatable with the chemical modification data from
//					subsequence i to j
//arrays are diagonal so subct[i][j] is replaced with subenergy[j][i]



//read the modification data from the file datafile
ifstream in;
in.open(datafile);
for (i=1;i<=maxmod&&!(in.eof());i++) {
	in>>mod[i];
}

numofmod = i-2;



//dynamically alocate enough space in arrays **subenergy,**subct,and **subcompat

	subenergy=new short int *[ct->numofbases+3];
	subcompat=new short int *[ct->numofbases+3];
   //subct= new int *[ct->numofbases+3];
   for (i=0;i<=(ct->numofbases+2);i++) {
   	subenergy[i] = new short int [ct->numofbases+3];
		subcompat[i] = new short int [ct->numofbases+3];
      //subct[i] = new int [ct->numofbases+3];
   }

   





//initialize values in subenergy and subct
for (i=0;i<=ct->numofbases+2;i++) {
	for (j=i;j<=ct->numofbases+2;j++) {
   	subenergy[i][j]=INFINITE_ENERGY;
      subenergy[j][i] = 0;
   }
}



for (n=1;n<=ct->numofstructures;n++) {//this is a loop over all suboptimal
													//			structures
      update->update(int((100*double(n))/(double(4*ct->numofbases))));
		for (i=1;i<=ct->numofbases;i++) {//initialize the values in subcompat
		    for (j=i;j<=ct->numofbases;j++) {
				subcompat[i][j] = true;
		    }
		}

      for (k=1;k<=numofmod;k++) {//check the compatability of structure
      									//	n with each modified base
         	if (!compat(ct,mod[k],n)) {//if a structure is incompatable with
            									//	a modification data point, make this
                                       //	reflected in subcompat:
					for (i=mod[k]+1;i<=ct->numofbases;i++) {
						subcompat[mod[k]][i]=false;
					}
					for (i=mod[k]-1;i>=1;i--) {
						for (j=mod[k];j<=ct->numofbases;j++) {
							subcompat[i][j] = false;
						}
					}
				}
      }
   //let the user know the progress
//	cout <<"Working on suboptimal structure # "<<n<<"\n"<<flush;


	for (i=1;i<=(ct->numofbases-2);i++) {

			for(j=ct->numofbases;j>=i+2;j--) {

      		//check structure n from base i to base j for compatability

				if (!subcompat[i][j]) goto sub2;//if incompatable, exit one loop


            //Now make sure there are base pairs and that they are
            //	alowed, ie: any pair is to a base between i and j, inclusive

            nopair = true;

            for (k=i;k<=j;k++) {
            	if ((ct->basepr[n][k]>j)) goto sub2;//incompatable pair
               else if ((ct->basepr[n][k]!=0)&&(ct->basepr[n][k]<i)) goto sub2;
               else if ((ct->basepr[n][k]!=0)) nopair = false;
            }

            if (nopair) goto sub1; //no pairs, so exit two loops

            //Now structure n from i to j is an allowed structure

            //make a substructure from i to j

            substructure.numofbases=j-i+1;
            for (k=1;k<=substructure.numofbases;k++) {
            	if (ct->basepr[n][i+k-1]!=0)
            		substructure.basepr[1][k]=(ct->basepr[n][i+k-1]-i+1);
               else substructure.basepr[1][k] = 0;
               substructure.numseq[k]=ct->numseq[i+k-1];
            }

            //calculate the energy

            efn2(data,&substructure);


            //if this is the most favorable substructure yet found, record it

            if (substructure.energy[1]<subenergy[i][j]) {
            	subenergy[i][j] = substructure.energy[1];
               subenergy[j][i] = n;
            }
	 		sub2:
			continue;
         }

      sub1:
		continue;
      }



}


//for (j=3;j<=ct->numofbases;j++) {
//	for (i=1;i<=j;i++) {
//   	out << "i= "<<i<<"  j= "<<j<<"  subenergy = "<<subenergy[i][j]
//      	<<"  subct = "<<subct[i][j]<<"\n";
//   }
//}



//The tables are full, search for the best structure

//cout <<"tables are full\n"<<flush;

//Now find the optimal re-construction of the sub-domains
//		For a subsequence from i to j not used, set subct[i][j] to zero
for (size=3;size<=ct->numofbases;size++) {
	//update the user on progress
	//cout << "size ="<<size<<"\n";
   update->update(int((100*double(size+ct->numofbases))/
   	(double(2*ct->numofbases))));
	for (pos=0;pos<=ct->numofbases-size;pos++) {

   	i = ct->numofbases-size-pos+1;
      j = ct->numofbases-pos;



      if (subenergy[i][j]>subenergy[i+1][j]) {
      	subenergy[i][j]=subenergy[i+1][j];
         subenergy[j][i] = 0;

         //out << "i+1,j "<<i<<","<<j<<"\n";

      }
      if	(subenergy[i][j]>subenergy[i][j-1]) {
      	subenergy[i][j] = subenergy[i][j-1];
         subenergy[j][i] = 0;

         //out << "i,j-1 "<<i<<","<<j<<"\n";

      }
      if (subenergy[i][j]>subenergy[i+2][j]) {
      	subenergy[i][j] = subenergy[i+2][j];
         subenergy[j][i] = 0;

         //out << "i+2,j "<<i<<","<<j<<"\n";

      }
      if (subenergy[i][j]>subenergy[i][j-2]) {
      	subenergy[i][j]=subenergy[i][j-2];
         subenergy[j][i]= 0;

         //out << "i,j-2 "<<i<<","<<j<<"\n";

      }
      if (size>5) {
      	for (k=i+2;k<=j-3;k++) {
         	if ((subenergy[k+1][j]+subenergy[i][k])<subenergy[i][j]) {
            	subenergy[i][j] = subenergy[k+1][j]+subenergy[i][k];
               subenergy[j][i]=0;

               //out << "i,j;k "<<i<<","<<j<<" "<<k<<"\n";


            }
         }
      }
      if (size>12) {
      for (k=1;k<=ct->numofstructures;k++) {
         //check each structure:
         if (ct->basepr[k][i] == j) {
         	for (ip=i;ip<j-5;ip++) {
         		for (jp=j;jp>ip+5;jp--) {
          			subcompat[ip][jp] = INFINITE_ENERGY;
            	}
         	}
         	twoefn2(data,ct,i,j,subcompat,k);
         	for (ip=i;ip<j-5;ip++) {
         		for (jp=j;jp>ip+5;jp--) {
            		ext = externalcompat(ct,i,ip,j,jp,k,numofmod,mod);

               	if (ext==1) goto jumpi;//incompatability between i and ip
               	else if (ext==2) goto jumpj;//incompatability between j and jp
                  else if (ext == 3) {
                  	if ((subcompat[ip][jp]+subenergy[ip+1][jp-1])<
                  		(subenergy[i][j])) {
                     	subenergy[i][j] = subcompat[ip][jp]
                        	+subenergy[ip+1][jp-1];
                     	subenergy[j][i] = -k;
                 	}



                  }



            	}
            	jumpj:
				continue;
         	}
         	jumpi:
			continue;
         }
       	/*for (ip=i;ip<j-5;ip++) {
         	for (jp=j;jp>ip+5;jp--) {

            	if (ct->basepr[k][ip]==jp) {
                	//there is an external domain possible
                  ext = externalcompat(ct,i,ip,j,jp,k,numofmod,mod);

                  if (ext==1) {
                  	goto jumpi;//there is an incompatability between i - ip
                  }
                  else if (ext==2) {
                   	goto jumpj;//there is an incompatability between jp - j
                  }
                  else if (ext == 3) {
                  substructure.numofbases = ip - i + j - jp + 5;
                  for (l=i;l<=ip;l++) {
                   	substructure.numseq[l-i+1] = ct->numseq[l];
                     if ((ct->basepr[k][l]>=i)&&(ct->basepr[k][l]<=ip)) {
                     	substructure.basepr[1][l-i+1] = ct->basepr[k][l] - i+1;
                     }
                     else if ((ct->basepr[k][l]>=jp)&&(ct->basepr[k][l]<=j)) {
                      	substructure.basepr[1][l-i+1] = ct->basepr[k][l] -jp+5+ip-i;
                     }
                     else substructure.basepr[1][l-i+1] = 0;
                  }
                  substructure.numseq[ip-i+2] = 5;
                  substructure.numseq[ip-i+3] = 5;
                  substructure.numseq[ip-i+4] = 5;

                  substructure.basepr[1][ip-i+2] = 0;
                  substructure.basepr[1][ip-i+3] = 0;
                  substructure.basepr[1][ip-i+4] = 0;

                  substructure.intermolecular = true;

                  substructure.inter[0] = ip-i+2;
                  substructure.inter[1] = ip-i+3;
                  substructure.inter[2] = ip-i+4;

                  for (l=jp;l<=j;l++) {
                   	substructure.numseq[ip-i+5+l-jp] = ct->numseq[l];
                     if ((ct->basepr[k][l]>=i)&&(ct->basepr[k][l]<=ip)) {
                     	substructure.basepr[1][ip-i+5+l-jp] = ct->basepr[k][l] - i + 1;
                     }
                     else if ((ct->basepr[k][l]>=jp)&&(ct->basepr[k][l]<=j)) {
                      	substructure.basepr[1][ip-i+5+l-jp] = ct->basepr[k][l] -jp+5+ip-i;
                     }
                     else substructure.basepr[1][ip-i+5+l-jp] = 0;
                  }

                  //now substructure is complete
                  efn2(data,&substructure);

                  if ((substructure.energy[1]-34+subenergy[ip+1][jp-1])<
                  	(subenergy[i][j])) {
                     subenergy[i][j] = substructure.energy[1]-34+subenergy[ip+1][jp-1];
                     subct[i][j] = -k;
                 	}
                  }
               }

            }
            jumpj:
         }
         jumpi: */
      }

      }

   }
}
//cout << "Lowest free energy found\n"<<flush;

//for (j=3;j<=ct->numofbases;j++) {
//	for (i=1;i<=j;i++) {
//   	out << "i= "<<i<<"  j= "<<j<<"  subenergy = "<<subenergy[i][j]
//      	<<"  subct = "<<subct[i][j]<<"\n";
//   }
//}


//i = 1;
//j = ct->numofbases;

st.put(1,ct->numofbases);
//Now build a structure that uses the subdomains in the lowest free energy
//	structure found above
while (st.get(&i,&j)) {

//	cout << "i = "<<i<<"\n  j = "<<j<<"\n";

   if (subenergy[j][i]!=0) {
   	if (subenergy[j][i]>0) {
   		for (k=i;k<=j;k++) {
      		ctout->basepr[1][k]=ct->basepr[subenergy[j][i]][k];
      	}
      	//i = j + 1;
      	//j = ct->numofbases;
      	//if (st.get(&i,&j) == 0) {
         //	i = j+1;
         //   j = ct->numofbases;
         //}

      }
      else {//best structure from i to j has an exterior domain:
      	k = -subenergy[j][i];
         	for (ip=i;ip<j-5;ip++) {
         		for (jp=j;jp>ip+5;jp--) {
          			subcompat[ip][jp] = INFINITE_ENERGY;
            	}
         	}
         	twoefn2(data,ct,i,j,subcompat,k);
         	for (ip=i;ip<j-5;ip++) {
         		for (jp=j;jp>ip+5;jp--) {
            		ext = externalcompat(ct,i,ip,j,jp,k,numofmod,mod);

               	//if (ext==1) goto jumpi;//incompatability between i and ip
               	if (ext==2) goto jumpj2;//incompatability between j and jp
                  else if (ext == 3) {
                  	if ((subcompat[ip][jp]+subenergy[ip+1][jp-1])==
                  		(subenergy[i][j])) {

              				goto jump;

                 		}



                  }



            	}
            	jumpj2:
				continue;
         	}


      	/*for (ip=i;ip<j-3;ip++) {
       		for (jp=j;jp>ip+3;jp--) {

            	if (ct->basepr[k][ip] == jp) {
               substructure.numofbases = ip - i + j - jp + 5;
               if (externalcompat(ct,i,ip,j,jp,k,numofmod,mod)) {
                  for (l=i;l<=ip;l++) {
                   	substructure.numseq[l-i+1] = ct->numseq[l];
                     if ((ct->basepr[k][l]>=i)&&(ct->basepr[k][l]<=ip)) {
                     	substructure.basepr[1][l-i+1] = ct->basepr[k][l] - i + 1;
                     }
                     else if ((ct->basepr[k][l]>=jp)&&(ct->basepr[k][l]<=j)) {
                      	substructure.basepr[1][l-i+1] = ct->basepr[k][l] -jp+5+ip-i;
                     }
                     else substructure.basepr[1][l-i+1] = 0;
                  }
                  substructure.numseq[ip-i+2] = 5;
                  substructure.numseq[ip-i+3] = 5;
                  substructure.numseq[ip-i+4] = 5;

                  substructure.intermolecular = true;

                  substructure.inter[0] = ip-i+2;
                  substructure.inter[1] = ip-i+3;
                  substructure.inter[2] = ip-i+4;

                  substructure.basepr[1][ip-i+2] = 0;
                  substructure.basepr[1][ip-i+3] = 0;
                  substructure.basepr[1][ip-i+4] = 0;

                  for (l=jp;l<=j;l++) {
                   	substructure.numseq[ip-i+5+l-jp] = ct->numseq[l];
                     if ((ct->basepr[k][l]>=i)&&(ct->basepr[k][l]<=ip)) {
                     	substructure.basepr[1][ip-i+5+l-jp] = ct->basepr[k][l] - i + 1;
                     }
                     else if ((ct->basepr[k][l]>=jp)&&(ct->basepr[k][l]<=j)) {
                      	substructure.basepr[1][ip-i+5+l-jp] = ct->basepr[k][l] -jp+5+ip-i;
                     }
                     else substructure.basepr[1][ip-i+5+l-jp] = 0;
                  }

                  //now substructure is complete
                  efn2(data,&substructure);

                  if ((substructure.energy[1]+subenergy[ip+1][jp-1]-34)==
                  	(subenergy[i][j]))
                     //the correct ip and jp have been found

                     goto jump;

               }
               }

            }
         }
         jump:*/

         jump:
         for (l=i;l<=ip;l++) {
          	ctout->basepr[1][l] = ct->basepr[k][l];
         }
         for (l=jp;l<=j;l++) {
            ctout->basepr[1][l] = ct->basepr[k][l];

         }

         st.put(ip+1,jp-1);
         i = ip+1;
         j = jp-1;



      }




   }
   else if (subenergy[i][j-1]==subenergy[i][j]) {
   	j--;
      st.put(i,j);

   }
   else if (subenergy[i][j-2]==subenergy[i][j]) {
   	j=j-2;
      st.put(i,j);

   }
   else if (subenergy[i+1][j]==subenergy[i][j]) {
   	ctout->basepr[1][j] = 0;
      i++;
      st.put(i,j);

   }
   else if (subenergy[i+2][j]==subenergy[i][j]) {
   	ctout->basepr[1][j] = 0;
      i=i+2;
      st.put(i,j);

   }
   else {
   	for (k=i+2;k<=j-2;k++) {
      	if (subenergy[i][j]==(subenergy[i][k]+subenergy[k+1][j])) {
         	st.put(i,k);
            st.put(k+1,j);
            break;
         }
      }
   }
}



//record information in ctout:
for (i=1;i<=ct->numofbases;i++) {
	ctout->numseq[i] = ct->numseq[i];
   
}
ctout->numofstructures=1;
ctout->numofbases=ct->numofbases;
ctout->energy[1]=subenergy[1][ct->numofbases];
strcpy(ctout->ctlabel[1],"mix&match result\n");



//clean up memory use:
for (i=0;i<ct->numofbases+2;i++) {
   	//delete[] subct[i];
	delete[] subenergy[i];
   delete[] subcompat[i];
}
//delete subct;
delete[] subenergy;
delete[] subcompat;

}

//Bimol provides bimolecular secondary structure prediction where no intramolecular pairing is allowed.
//Written by Laura DiChiacchio in 2008 and modernized by David Mathews in 2014.

#include "algorithm.h"
#include "rna_library.h"
#include "structure.h"
#include "stackstruct.h"
#include "stackclass.h"
#include "platform.h"
#include "defines.h"
#include "../RNA_class/RNA.h"
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <chrono>

using namespace std;


void bimoltracebackV(int i, int j, int jct3, int k, int N1, int N2, int maxloop, short **V, structure *ct3, datatable *data);
void bimoltracebackVp(int i, int j, int jct3, int k, int N1, int N2, int maxloop, short **Vp, structure *ct3, datatable *data);

void accessfoldtracebackV(int i, int j, int jct3, int k, int N1, int N2, int maxloop, double **V, double *G, structure *ct3, datatable *data);
void accessfoldtracebackVp(int i, int j, int jct3, int k, int N1, int N2, int maxloop, double **Vp, double *G, structure *ct3, datatable *data);

void bimol(structure *ct1, structure *ct2, structure *ct3, int maxloop, int maxtracebacks, int percent, int windowsize, 
           datatable *data) {
  

  int N1 = ct1->GetSequenceLength();
  int N2 = ct2->GetSequenceLength();
  int maxseq = ct1->GetSequenceLength() + ct2->GetSequenceLength() + 3;

  int i,j;

    
  int a, b, k, k1, k2, k3;
  int vmin, crit, num, numbp;
  int q, up, ir, c, cur;
  int iret, jret, sort;
  int *energy, *heapi, *heapj; 
  
  ct3->allocate(maxseq);
  
  //Set the location of the thermodynamic data tables
  ct3->SetThermodynamicDataTable(data);

  vector< vector<bool> > inc = data->pairing;
  
  //ct3 must contain a sequence that concatonates the sequences of ct1 and ct2  
  string label;
  label = ct1->GetSequenceLabel();
  
  //Make sure there is no newline at the end of this line.
	if (label.length()>0) if (label[label.size()-1]=='\n') {
		label.erase(label.size()-1,1);		
	}

  label+="_";
  label+= ct2->GetSequenceLabel();
  ct3->SetSequenceLabel(label);

  
  
  for (i=1; i<=N1; i++) {
    ct3->numseq[i] = ct1-> numseq[i];
    ct3->nucs[i] = ct1->nucs[i];
    ct3->hnumber[i] = i;
  }
  
  //add intermolecular linker
  for (i=N1+1; i<(N1+4); i++) {

	  //ToDo: Could guard this code to make sure the linker exists in this dataset.

	//use the first linker entry (likely the only linker entry)
	ct3->numseq[i]=data->basetonum(data->linker[0]);
    ct3->nucs[i] = data->linker[0];
    ct3->hnumber[i] = 0;
  }
  
  for (i=N1+4, j=1; i<(N1+N2+4); i++, j++) {
    ct3->numseq[i] = ct2->numseq[j];
    ct3->nucs[i] = ct2->nucs[j];
    ct3->hnumber[i] = j;
  }
        
  //dynamic allocation of memory for V[i][j], the mfe from 1 to i in sequence 1 and j to N2 in sequence 2, where i is paired to j.
  short **V;
  V = new short *[N1+1];
  for (int k=0; k<=N1; k++) {
    V[k] = new short [N2+1];
  }
   
  //dynamic allocation of memory for Vp[i][j], the mfe from i to N1 in sequence 1 and 1 to j in sequence 2, where i is paired to j.
  short **Vp;
  Vp = new short *[N1+1];
  for (int k = 0; k <= N1; k++) {
    Vp[k] = new short[N2+1];
  }
  
  //main loops over i and j to fill the V array
  for (i=N1; i >= 1; i--) {
#ifdef SMP
#pragma omp parallel for
#endif
    for (j=1; j <= N2; j++) {
      int jct3 = N1+3+j;
      int ip, jp, jpct3;
      int big, biggest, end;
      //initialize V[i][j]
      V[i][j] = INFINITE_ENERGY;
      
	  //first, check if this is a valid pair
      if (inc[ct3->numseq[i]][ct3->numseq[jct3]]) {
        
        //check case for a terminal mismatch on a starting pair
        if (i < N1 && j > 1) {
          V[i][j] = min(V[i][j], data->tstack[ct3->numseq[i]][ct3->numseq[jct3]][ct3->numseq[i+1]][ct3->numseq[jct3-1]] +   
                        penalty(i, jct3, ct3, data));
        } 
        //i+1 dangles if nothing can stack on j
        else if (i < N1 && j == 1) {
          V[i][j] = min(V[i][j], erg4(i, jct3, i+1, 1, ct3, data, false) + penalty(i, jct3, ct3, data));
        }  
        //j-1 dangles if nothing can stack on i
        else if (i == N1 && j > 1) {
          V[i][j] = min(V[i][j], erg4(i, jct3, jct3-1, 2, ct3, data, false) + penalty(i, jct3, ct3, data));
        }
		//here nothing can stack because this pair involves terminal nucleotides
        else (V[i][j] = min(V[i][j], penalty(i, jct3, ct3, data)));
        
		//now loop over all possible closing pairs at the far end of a base pair stack, bulge loop, or interior loop
        biggest = N1-i+j-3; 
        biggest = min(biggest, maxloop);
        
        for (int a = 0; a <= biggest; a++) {
          end = i+a+1;
          big = min(end, N1);
          for (ip = max(i+1, i-j+a+3); ip <= big; ip++) {
            jp = j+ip-i-a-2;
            jpct3 = jp+N1+3;
             
			//make sure this pair at the far end is valid
            if (inc[ct3->numseq[ip]][ct3->numseq[jpct3]]) {
              
              //base pair stack
			  if (a==0){//ip == i+1 && jp == j-1) {
                V[i][j] = min(V[i][j], erg1(i, jct3, ip, jpct3, ct3, data)+V[ip][jp]);
              }
              
              //either internal loop or bulge
              else V[i][j] = min(V[i][j], erg2(i, jct3, ip, jpct3, ct3, data, 0, 0)+V[ip][jp]);
              
              
            }
          }//end loop over ip
        
        }//end loop over a
      }//end check on valid pair
    }//end loop over j
  }//end loop over i
    
  
  //main loops over i and j for Vp, V calculations in reverse direction 
  for (i=1; i <= N1; i++) {
#ifdef SMP
#pragma omp parallel for
#endif
    for (j = N2; j >= 1; j--) {
      int jct3 = maxseq-N2+j;
      int ip, jp, jpct3;
      int big, biggest, end;
	  //initialize Vp
      Vp[i][j] = INFINITE_ENERGY; 
      
	  //Make sure this is a valid pair
      if (inc[ct3->numseq[i]][ct3->numseq[jct3]]) {
    
        //consider a terminal mm
        if (i > 1 && j < N2) {
          Vp[i][j] = min(Vp[i][j], data->tstack[ct3->numseq[jct3]][ct3->numseq[i]][ct3->numseq[jct3+1]][ct3->numseq[i-1]] + 
                         penalty(jct3, i, ct3, data));
        }   
        //i-1 dangles if nothing can stack on j
        else if (i > 1 && j == N2) {
          Vp[i][j] = min(Vp[i][j], erg4(jct3, i, i-1, 2, ct3, data, false) + penalty(jct3, i, ct3, data));
        }
        //j+1 dangles if nothing can stack on i
        else if (i == 1 && j < N2) {
          Vp[i][j] = min(Vp[i][j], erg4(jct3, i, jct3+1, 1, ct3, data, false) + penalty(jct3, i, ct3, data));
        }     
        //terminal bp, so no stacks are allowed
        else (Vp[i][j] = min(Vp[i][j], penalty(jct3, i, ct3, data)));
        

		//Consider an internal loop, bulge loop, or helical stack closed at the far end by another pair
        biggest = N2-j+i-3;
        biggest = min(biggest, maxloop);
     
        for (int a = 0; a <= biggest; a++) {
          end = j+a+1;
          big = min(end, N2); 
          for (jp = max(j+1, j-i+a+3), jpct3 = jp+N1+3; jp <= big; jp++, jpct3++) {
            ip = i+jp-j-a-2;
            
			//check that the pair at the far end is valid
            if (inc[ct3->numseq[ip]][ct3->numseq[jpct3]]) {
            
              //bp stack
			  if (a==0){//ip == i-1 && jp == j+1) {
                Vp[i][j] = min(Vp[i][j], erg1(jct3, i, jpct3, ip, ct3, data)+Vp[ip][jp]);
              }
              
              //internal loop or bulge
              else Vp[i][j] = min(Vp[i][j], erg2(jct3, i, jpct3, ip, ct3, data, 0, 0)+Vp[ip][jp]);
            }//end valid pair at far end
        
          }//end loop over ip
        }//end loop over a
      }//end valid pair
      
    }//end loop over j
      
  }//end loop over i

   
  //begin traceback
  //for (k = 1; k <= maxtracebacks; k++) {

    //initialize vmin
    vmin = INFINITE_ENERGY;
    
    
    //find lowest dG, set vmin = min dG
    for (i = N1; i >= 1; i--) {
      for (j = 1; j <= N2; j++) {
        
        if ((V[i][j] + Vp[i][j] + data->init) < vmin) {
          vmin = V[i][j] + Vp[i][j] + data->init;
        }
      }
    }

	//This code previously prevented tracebacks for structures with positive free energy change.
	//Now, just trap and forbid tracebacks for very high folding free energy changes.  These, for example, 
	//	are structures that have no viable pairs.
	//Return the mfe as the unpaired structure, however.
	if (vmin >= 0) {
		//Add an empty structure to ct:
		ct3->AddStructure();
		//Set the energy to zero, a random coil:
		ct3->SetEnergy(1, 0);
		if (vmin >= INFINITE_ENERGY / 2) return;//no viable structure found other than the empty structure, i.e. random coil
	}
    
    //mark will keep track of bps that have been seen before
    bool **mark;
    mark = new bool *[maxseq+1];
    for (i = 0; i <= maxseq; i++) {
      mark[i] = new bool [maxseq+1]; 
    }
    
    //initialize mark
    for (int count = 1; count <= maxseq; count++) {
      for (int count2 = 1; count2 <= maxseq; count2++) {
        mark[count][count2] = false;
      }
    }
    
    //set critical value for vmin - suboptimal structures within percentage difference of lowest dG set by user
	crit = (short)(abs(vmin) * (float(percent) / 100.0));
	crit = crit + vmin;


    
    sort = 90000;
    
    energy = new int[sort+1];
    heapi = new int[sort+1];
    heapj = new int[sort+1];
    
    num = 0;
    for (i = N1; i >= 1; i--) {
      for (j = 1; j <= N2; j++) {
        if (num == sort) {
          
          //allocate more space for heapi, heapj, energy because it is about to overflow
          delete[] heapi;
          delete[] heapj;
          delete[] energy;
          sort = 10*sort;
          heapi = new int[sort+1];
          heapj = new int[sort+1];
          energy = new int[sort+1];
          i = N1;
          j = 1;
          num = 0;
        
        }
      
        //add critical values to heap
        if (V[i][j] + Vp[i][j] + data->init <= crit) {
      
        num++;
        heapi[num] = i;
        heapj[num] = j;
        energy[num] = V[i][j] + Vp[i][j] + data->init;
        
        }
      }
    }
    
    //make a heap
    for (q = 2; q <= num; q++) {
      
      cur = q;
      up = cur/2;
      
      while ((energy[cur] < energy[up] && up >= 1)) {
        
        swap(heapi[cur], heapi[up]);
        swap(heapj[cur], heapj[up]);
        swap(energy[cur], energy[up]);
        cur = cur/2;
        up = up/2;
        
      }
    }
    
    //sort the heap
    for (ir = num-1; ir >= 1; ir--) {
      swap(heapi[ir+1], heapi[1]);
      swap(heapj[ir+1], heapj[1]);
      swap(energy[ir+1], energy[1]);
      
      up = 1;
      c = 2;
      
      while (c <= ir) {
        
        if (c != ir) {
          if (energy[c+1] < energy[c]) c++;
        }
        if (energy[c] < energy[up]) {
          swap(heapi[c], heapi[up]);
          swap(heapj[c], heapj[up]);
          swap(energy[c], energy[up]);
          up = c;
          c = 2*c;
          
        }
        else c = ir+1;
      }
    }
    
	//ct3->numofstructures = 0; //initialize the number of structures
	bool failedprevious = false;//track whether a structure was traced, but not different enough from previous structures

    //for max number of structures input by user
    for (int cntr = num; cntr > 0 && (((ct3->GetNumberofStructures()) < maxtracebacks)||(failedprevious)); cntr--) {
      
      //find next unmarked bp
      if (!mark[heapi[cntr]][N1+3+heapj[cntr]]) {
        
        iret = heapi[cntr];
        jret = heapj[cntr];

		if (failedprevious) ct3->CleanStructure(ct3->GetNumberofStructures()); //last structure should not be kept, it should have pairs scrubbed
		else ct3->AddStructure();//Need a new structure to contain this traceback
		failedprevious=false;
                
        int jct3 = jret + N1 + 3;



		
        
        //perform traceback
        bimoltracebackV(iret, jret, jct3, ct3->GetNumberofStructures(), N1, N2, maxloop, V, ct3, data);
        bimoltracebackVp(iret, jret, jct3, ct3->GetNumberofStructures(), N1, N2, maxloop, Vp, ct3, data);
      
        ct3->SetEnergy(ct3->GetNumberofStructures(),energy[cntr]);
        
        numbp = 0;

		//check if this structure is different enough from previous structures using the idea of a window
        for (k1 = 1; k1 <= maxseq; k1++) {
          
          if (k1 < (ct3->GetPair(k1,ct3->GetNumberofStructures()))) {
            if (!(mark[k1][ct3->GetPair(k1,ct3->GetNumberofStructures())])) numbp++;
          }
        }
        for (k1 = 1; k1 <= maxseq; k1++) {
          if (k1 < (ct3->GetPair(k1,ct3->GetNumberofStructures()))) {
            mark[k1][ct3->GetPair(k1,ct3->GetNumberofStructures())] = true;
            if (windowsize > 0) {
              for (k2 = max(1, k1 - windowsize); k2 <= min(maxseq, k1 + windowsize); k2++) {
                for (k3 = max(k2, (ct3->GetPair(k1,ct3->GetNumberofStructures()) - windowsize)); k3 <= min(maxseq, (ct3->GetPair(k1,ct3->GetNumberofStructures()) + windowsize)); k3++) {
                      
                      mark[k2][k3] = true;
                      
                    }
              }
            }
          }
        }
        
        if (numbp <= windowsize && ct3->GetNumberofStructures() > 1) {
          failedprevious=true;
        }
        else {


			//place the structure name (from ctlabel[1]) into each structure
			ct3->SetCtLabel(ct3->GetSequenceLabel(), ct3->GetNumberofStructures());
        } 
        
      }
    }
	if (failedprevious) ct3->RemoveLastStructure();//last structure should not be kept
        
    de_allocate (mark, maxseq+1);
    delete[] energy;
    delete[] heapi;
    delete[] heapj;
    
  //}
    
    for (int l = 0; l <= N1; l++) {
      delete[] V[l];
      delete[] Vp[l];
    }
    
    delete[] V;
    delete[] Vp;
    
}



void bimoltracebackV(int i, int j, int jct3, int k, int N1, int N2, int maxloop, short **V, structure *ct3, datatable *data) {
  
  bool done, found;
  int jpct3, ip, jp;
  int biggest, big, end; 
  vector< vector<bool> > inc = data->pairing;
  
  //done will switch to true when the end of this traceback is complete
  done = false; 

  while (!done) {

    //register the current base pair and find the next, if there are any
     ct3->SetPair(i,jct3,k);
     
      
	 
	 
      
	 if (i!=N1&&j!=1) {
		 if (V[i][j] == data->tstack[ct3->numseq[i]][ct3->numseq[jct3]][ct3->numseq[i+1]][ct3->numseq[jct3-1]] + penalty(i, jct3, ct3, data)) {
          
			done = true;
		 }
          
     }
	 else if (i!=N1) {
		  
		 if (V[i][j] == erg4(i, jct3, i+1, 1, ct3, data, false) + penalty(i, jct3, ct3, data)) {
          
			done = true;
          
		}
	 }
	 else if (j!=1){
		 if (V[i][j] == erg4(i, jct3, jct3-1, 2, ct3, data, false) + penalty(i, jct3, ct3, data)) {
        
			done = true;
        
		 }
	 }
	 else {
		 if (V[i][j] == penalty(i, jct3, ct3, data)) {
        
			done = true;
		 }
        
	 }

	 if (!done) {
		//look for the next pair because we are not yet done tracing back  

	    //found will change to true when the next pair is found
        found = false;
        
        biggest = N1-i+j-3; 
        biggest = min(biggest, maxloop);
        
        for (int a = 0; a <= biggest && !found; a++) {
          end = i+a+1;
          big = min(end, N1);
          for (ip = max(i+1, i-j+a+3); ip <= big  && !found ; ip++) {
            jp = j+ip-i-a-2;
            jpct3 = jp+N1+3;
            
            if (ip == i+1 && jp == j-1 && V[i][j] == erg1(i, jct3, ip, jpct3, ct3, data)+V[ip][jp] && 
                 inc[ct3->numseq[ip]][ct3->numseq[jpct3]]) {
              
              i++;
              j--;
              jct3--;
              found = true;
                         
            }
            else if (V[i][j] == erg2(i, jct3, ip, jpct3, ct3, data, 0, 0)+V[ip][jp] &&  
                     inc[ct3->numseq[ip]][ct3->numseq[jpct3]]) {
              
              i = ip;
              j = jp;
              jct3 = jpct3;
              found = true;
            }
          }
          
        }
        
        if (found == false) {
          found = true;
          cerr << "Error in tracebackV at "<< i << " " << jct3 << V[i][j] << "\n";
          done = true;//throw in the towel because an error was found 
        }
      }//end if !done
	}//end while !done
}
      
void bimoltracebackVp(int i, int j, int jct3, int k, int N1, int N2, int maxloop, short **Vp, structure *ct3, datatable *data) { 
  
  bool done, found;
  int jpct3, ip, jp;
  int biggest, big, end;
  vector< vector<bool> > inc = data->pairing;
  
  //done will switch to true when this traceback is complete
  done = false;
  
  while (!done) {
    
    //record the current base pair and then go look for the next
    ct3->SetPair(i,jct3,k);
    
	if (i!=1&&j!=N2) {
		if(Vp[i][j] == data->tstack[ct3->numseq[jct3]][ct3->numseq[i]][ct3->numseq[jct3+1]][ct3->numseq[i-1]] + penalty(jct3, i, ct3, data)) {
      
			done = true;
		}
      
    }
	else if (i!=1) {
		if (Vp[i][j] == erg4(jct3, i, i-1, 2, ct3, data, false) + penalty(jct3, i, ct3, data)) {
        
			done = true;
		}
        
    }
	else if (j!=N2) {
		if (Vp[i][j] == erg4(jct3, i, jct3+1, 1, ct3, data, false) + penalty(jct3, i, ct3, data)) {
      
			done = true;
		}
      
    }
	else {
		if (Vp[i][j] == penalty(jct3, i, ct3, data)) {
    
			done = true;
		}
    }
	if (!done) {
      
      found = false;
      

      biggest = N2-j+i-3;
      biggest = min(biggest, maxloop);
     
      for (int a = 0; a <= biggest && !found; a++) {
        end = j+a+1;
        big = min(end, N2); 
        for (jp = max(j+1, j-i+a+3), jpct3 = jp+N1+3; jp <= big && !found; jp++, jpct3++) {
          ip = i+jp-j-a-2;
          
          if (ip == i-1 && jp == j+1 && Vp[i][j] == erg1(jct3, i, jpct3, ip, ct3, data)+Vp[ip][jp] && 
              inc[ct3->numseq[ip]][ct3->numseq[jpct3]]) {
            
            i--;
            j++;
            jct3++;
            found = true;
            
          }
            
          else if (Vp[i][j] == erg2(jct3, i, jpct3, ip, ct3, data, 0, 0)+Vp[ip][jp] && 
                   inc[ct3->numseq[ip]][ct3->numseq[jpct3]]) {
            
            i = ip;
            j = jp;
            jct3 = jpct3;
            found = true;
            
          }
        }
      }
      
      if (found == false) {
        found = true;
        cerr << "Error in tracebackVp at " << i << " " << jct3 << " " << Vp[i][j] << "\n";
        done = true; 
      } 
    
    }//end if !done
      
  }//end while !done
}


void accessfold(structure *ct1, structure *ct2, structure *ct3, int maxloop, int maxtracebacks, int percent, int windowsize, 
           datatable *data, double gamma, bool IsRNA, double temperature) {
  
  int i, j, ip, jp, jct3, jpct3;
  int big, biggest, end;
  int N1 = ct1->GetSequenceLength();
  int N2 = ct2->GetSequenceLength();
  int maxseq = ct1->GetSequenceLength() + ct2->GetSequenceLength() + 3;
    
  int a, b, k, k1, k2, k3;
  int num, numbp, error;
  int q, up, ir, c, cur;
  int iret, jret, sort;
  int *heapi, *heapj; 
  double vmin, crit, *energy;
  vector< vector<bool> > inc = data->pairing;

  //ct3->allocate(maxseq);
//  ct3->allocatestructure();          called during accessfold_interface
  ct3->allocate(maxseq);
  
  //set the locationof the data tables in ct3:
  ct3->SetThermodynamicDataTable(data);
  
  //ct3 must contain a sequence that concatonates the sequences of ct1 and ct2  
  string label;
  label = ct1->GetSequenceLabel();

  //Make sure there is no newline at the end of this line.
	if (label[label.size()-1]=='\n') {
		label.erase(label.size()-1,1);		
	}

  label+="_";
  label+= ct2->GetSequenceLabel();
  ct3->SetSequenceLabel(label);

   
  
  for (i=1; i<=N1; i++) {
    ct3->numseq[i] = ct1-> numseq[i];
    ct3->nucs[i] = ct1->nucs[i];
    ct3->hnumber[i] = i;
  }
  
  //add intermolecular linker
  for (i=N1+1; i<(N1+4); i++) {
    //ToDo: Could guard this code to make sure the linker exists in this dataset.

	//use the first linker entry (likely the only linker entry)
	ct3->numseq[i]=data->basetonum(data->linker[0]);
    ct3->nucs[i] = data->linker[0];
    ct3->hnumber[i] = 0;
  }
  
  for (i=N1+4, j=1; i<(N1+N2+4); i++, j++) {
    ct3->numseq[i] = ct2->numseq[j];
    ct3->nucs[i] = ct2->nucs[j];
    ct3->hnumber[i] = j;
  }
        
  //dynamic allocation of memory for V[i][j], the mfe from 1 to i in sequence 1 and j to N2 in sequence 2, where i is paired to j.
  double **V;
  V = new double *[N1+1];
  for (int k=0; k<=N1; k++) {
    V[k] = new double [N2+1];
  }
   
  //dynamic allocation of memory for Vp[i][j], the mfe from i to N1 in sequence 1 and 1 to j in sequence 2, where i is paired to j.
  double **Vp;
  Vp = new double *[N1+1];
  for (int k = 0; k <= N1; k++) {
    Vp[k] = new double[N2+1];
  }

//RNA class constructor - use RNA(filename, type, IsRNA) - type = 2 for .seq, type = 3 for .pfs 6/11
	
  char *sequence;

  //Convert the sequence in ct1 to a character array to feed to RNA class:

  //Allocate the space in the array
  sequence = new char [ct1->GetSequenceLength()+1];

  for (int index=0;index<ct1->GetSequenceLength();++index) {
	  //For each nucleotide, fetch sequence
	  sequence[index]=ct1->nucs[index+1];
  }
  //null terminate the string
  sequence[ct1->GetSequenceLength()]='\0';


	RNA *rna1 = new RNA(sequence, IsRNA);
	error = rna1->GetErrorCode();
  	if (error!=0) {
    		cerr << rna1->GetErrorMessage(error);
    		delete rna1;
			delete[] sequence;
			return;
  	}

	rna1->SetTemperature(temperature);
	error = rna1->PartitionFunction();
	if (error!=0) {
    		cerr << rna1->GetErrorMessage(error);
    		delete rna1;
			delete[] sequence;
			return;
  	}
	delete[] sequence;

	//Convert the sequence in ct2 to a character array to feed to RNA class:

  //Allocate the space in the array
  sequence = new char [ct2->GetSequenceLength()+1];

  for (int index=0;index<ct2->GetSequenceLength();++index) {
	  //For each nucleotide, fetch sequence
	  sequence[index]=ct2->nucs[index+1];
  }
  //null terminate the string
  sequence[ct2->GetSequenceLength()]='\0';


	RNA *rna2 = new RNA(sequence, IsRNA);
	error = rna2->GetErrorCode();
  	if (error!=0) {
    		cerr << rna2->GetErrorMessage(error);
			delete rna1;
    		delete rna2;
			delete[] sequence;
			return;
  	}

	rna2->SetTemperature(temperature);
	error = rna2->PartitionFunction();
	if (error!=0) {
    		cerr << rna2->GetErrorMessage(error);
    		delete rna1;
			delete rna2;
			delete[] sequence;
			return;
  	}
	delete[] sequence;

// dynamic allocation of P - 'unpairing' probability in ct3 6/11

	double *P;
	P = new double[maxseq+1];

// dynamic allocation of Gp - pseudo-free energy - pairing penalty 6/11
	
	double *Gp;
	Gp = new double[maxseq+1];
	

//perform partition function calculation for each seq independently 6/11
        
       //rna1->PartitionFunction(); // comment out if using pfs
	//error = rna1->PartitionFunction();	
  	//if (error!=0) {
	//    cerr << rna1->GetErrorMessage(error);
	//    delete rna1;
	// }

    //   rna2->PartitionFunction();
	//error = rna2->PartitionFunction();	
  	//if (error!=0) {
	//   cerr << rna2->GetErrorMessage(error);
	//   delete rna2;
 	// }

	//calculate probability P for each i of being unpaired 6/11

  // auto clock = std::chrono::system_clock::now().time_since_epoch();
  // auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(clock);
  // string name = "accessfold_";
  // name += std::to_string(ms.count());
  // name += ".txt";
  // ofstream out(name.c_str());

  // {
  //   out << "RNA1 Pair Probs\n";
  //   RNA &r = (*rna1);
  //   for(int i=1;i<r.GetSequenceLength();i++)
  //     for(int j=i+1;j<=r.GetSequenceLength();j++)
  //       out << i << "\t" << j << "\t" << r.GetPairProbability(i,j) << "\n";
  // }
  // {
  //   out << "\nRNA2 Pair Probs\n";
  //   RNA &r = (*rna2);
  //   for(int i=1;i<r.GetSequenceLength();i++)
  //     for(int j=i+1;j<=r.GetSequenceLength();j++)
  //       out << i << "\t" << j << "\t" << r.GetPairProbability(i,j) << "\n";
  // }
  // out.close();
  

	for (i = 1; i <= maxseq; i++) {
		P[i] = 1;
	} //initialize P[i] for all i over both seqs

	for (i = 1; i <= N1; i++) {
		for (j = i+3; j <= N1; j++) {
			P[i] = P[i] - rna1->GetPairProbability(i,j);
		}
		for (k = i-3; k >= 1; k--) {
			P[i] = P[i] - rna1->GetPairProbability(k,i);
		}
	} //calculate P for seq1

	for (i=1; i <= N2-3; i++) {
		for (j = i+3; j <= N2; j++) {
			P[i+N1+3] = P[i+N1+3] - rna2->GetPairProbability(i,j);
		}
		for (k=i-3; k>= 1; k--) {
			P[i+N1+3] = P[i+N1+3] - rna2->GetPairProbability(k,i);
    } 
	} //calcuate P for seq2

    //calculate Gp[i] for all i over both seqs 6/11
  for (i = 1; i <= maxseq; i++) {
    if (P[i]<=0.0) P[i] = 1E-100;
    Gp[i] = gamma*conversionfactor*RT_37C*log(P[i]);
    //if (std::isnan(Gp[i])) Gp[i] = INFINITE_ENERGY; // required because P[i]
  }
  
  // Debug values in P[i] and Gp[i]
  // {ofstream out("accessfold_PGP.log"); for(i=1; i<=maxseq; i++){out<<P[i]<<"\t"<<Gp[i]<<"\n";}}
  
  //main loops over i and j to fill the V array - pseudo-G term incorporated 6/11
  for (i=N1; i >= 1; i--) {
    for (j=1, jct3 = N1+4; j <= N2; j++, jct3++) {

      //initialize V[i][j]
      V[i][j] = INFINITE_ENERGY;
      
	  //first, check if this is a valid pair
      if (inc[ct3->numseq[i]][ct3->numseq[jct3]]) {
        
        //check case for a terminal mismatch on a starting pair
        if (i < N1 && j > 1) {
          V[i][j] = min(V[i][j], data->tstack[ct3->numseq[i]][ct3->numseq[jct3]][ct3->numseq[i+1]][ct3->numseq[jct3-1]] +   
                        penalty(i, jct3, ct3, data) - Gp[i] - Gp[jct3]);
        } 
        //i+1 dangles if nothing can stack on j
        else if (i < N1 && j == 1) {
          V[i][j] = min(V[i][j], erg4(i, jct3, i+1, 1, ct3, data, false) + penalty(i, jct3, ct3, data) - Gp[i] - Gp[jct3]);
        }  
        //j-1 dangles if nothing can stack on i
        else if (i == N1 && j > 1) {
          V[i][j] = min(V[i][j], erg4(i, jct3, jct3-1, 2, ct3, data, false) + penalty(i, jct3, ct3, data) - Gp[i] - Gp[jct3]);
        }
		//here nothing can stack because this pair involves terminal nucleotides
        else (V[i][j] = penalty(i, jct3, ct3, data) - Gp[i] - Gp[jct3]);
        
		//now loop over all possible closing pairs at the far end of a base pair stack, bulge loop, or interior loop
        biggest = N1-i+j-3; 
        biggest = min(biggest, maxloop);
        
        for (int a = 0; a <= biggest; a++) {
          end = i+a+1;
          big = min(end, N1);
          for (ip = max(i+1, i-j+a+3); ip <= big; ip++) {
            jp = j+ip-i-a-2;
            jpct3 = jp+N1+3;
             
			//make sure this pair at the far end is valid
            if (inc[ct3->numseq[ip]][ct3->numseq[jpct3]]) {
              
              //base pair stack
			  if (a==0){//ip == i+1 && jp == j-1) {
                V[i][j] = min(V[i][j], erg1(i, jct3, ip, jpct3, ct3, data)+V[ip][jp] - Gp[i] - Gp[jct3]);
              }
              
              //either internal loop or bulge
              else V[i][j] = min(V[i][j], erg2(i, jct3, ip, jpct3, ct3, data, 0, 0)+V[ip][jp] - Gp[i] -Gp[jct3]);
              
              
            }
          }//end loop over ip
        
        }//end loop over a
      }//end check on valid pair
    }//end loop over j
  }//end loop over i
  
  // Debug values in V[i][j]
  // {ofstream out("accessfold_V.log"); for(i=1; i<=N1; i++){for (j=1; j<=N2; j++) { out << V[i][j] << "\n";}}}

  //main loops over i and j for Vp, V calculations in reverse direction 
  for (i=1; i <= N1; i++) {
    for (j = N2, jct3 = maxseq; j >= 1; j--, jct3--) {

	  //initialize Vp
      Vp[i][j] = INFINITE_ENERGY; 
      
	  //Make sure this is a valid pair
      if (inc[ct3->numseq[i]][ct3->numseq[jct3]]) {
    
        //consider a terminal mm
        if (i > 1 && j < N2) {
          Vp[i][j] = min(Vp[i][j], data->tstack[ct3->numseq[jct3]][ct3->numseq[i]][ct3->numseq[jct3+1]][ct3->numseq[i-1]] + 
                         penalty(jct3, i, ct3, data) - Gp[i] - Gp[jct3]);
        }   
        //i-1 dangles if nothing can stack on j
        else if (i > 1 && j == N2) {
          Vp[i][j] = min(Vp[i][j], erg4(jct3, i, i-1, 2, ct3, data, false) + penalty(jct3, i, ct3, data) - Gp[i] - Gp[jct3]);
        }
        //j+1 dangles if nothing can stack on i
        else if (i == 1 && j < N2) {
          Vp[i][j] = min(Vp[i][j], erg4(jct3, i, jct3+1, 1, ct3, data, false) + penalty(jct3, i, ct3, data) - Gp[i] - Gp[jct3]);
        }     
        //terminal bp, so no stacks are allowed
        else (Vp[i][j] = min(Vp[i][j], penalty(jct3, i, ct3, data) - Gp[i] - Gp[jct3]));
        

		//Consider an internal loop, bulge loop, or helical stack closed at the far end by another pair
        biggest = N2-j+i-3;
        biggest = min(biggest, maxloop);
     
        for (int a = 0; a <= biggest; a++) {
          end = j+a+1;
          big = min(end, N2); 
          for (jp = max(j+1, j-i+a+3), jpct3 = jp+N1+3; jp <= big; jp++, jpct3++) {
            ip = i+jp-j-a-2;
            
			//check that the pair at the far end is valid
            if (inc[ct3->numseq[ip]][ct3->numseq[jpct3]]) {
            
              //bp stack
			  if (a==0){//ip == i-1 && jp == j+1) {
                Vp[i][j] = min(Vp[i][j], erg1(jct3, i, jpct3, ip, ct3, data)+Vp[ip][jp] - Gp[i] - Gp[jct3]);
              }
              
              //internal loop or bulge
              else Vp[i][j] = erg2(jct3, i, jpct3, ip, ct3, data, 0, 0)+Vp[ip][jp] - Gp[i] - Gp[jct3];
            }//end valid pair at far end
        
          }//end loop over ip
        }//end loop over a
      }//end valid pair
      
    }//end loop over j
      
  }//end loop over i

  // Debug values in Vp[i][j]
  // {ofstream out("accessfold_VP.log"); for(i=1; i<=N1; i++){for(j=1; j<=N2; j++){out<<Vp[i][j]<<"\n";}}}
   
  //begin traceback
  //for (k = 1; k <= maxtracebacks; k++) {

    //initialize vmin
    vmin = 0;
    
    
    //find lowest dG, set vmin = min dG
    for (i = N1; i >= 1; i--) {
      for (j = 1; j <= N2; j++) {
        
        if ((double)(V[i][j] + Vp[i][j] + data->init) < vmin) {
          vmin = (double)(V[i][j] + Vp[i][j] + data->init);
        }
      }
    }
    
    //mark will keep track of bps that have been seen before
    bool **mark;
    mark = new bool *[maxseq+1];
    for (i = 0; i <= maxseq; i++) {
      mark[i] = new bool [maxseq+1]; 
    }
    
    //initialize mark
    for (int count = 1; count <= maxseq; count++) {
      for (int count2 = 1; count2 <= maxseq; count2++) {
        mark[count][count2] = false;
      }
    }
    
    //set critical value for vmin - suboptimal structures within percentage difference of lowest dG set by user
    crit = vmin - (double) ((((float) percent)/100.0)*((float) vmin));
    
    sort = 90000;
    
    energy = new double[sort+1];
    heapi = new int[sort+1];
    heapj = new int[sort+1];
    
    num = 0;
    for (i = N1; i >= 1; i--) {
      for (j = 1; j <= N2; j++) {
        if (num == sort) {
          
          //allocate more space for heapi, heapj, energy because it is about to overflow
          delete[] heapi;
          delete[] heapj;
          delete[] energy;
          sort = 10*sort;
          heapi = new int[sort+1];
          heapj = new int[sort+1];
          energy = new double[sort+1];
          i = N1;
          j = 1;
          num = 0;
        
        }
      
        //add critical values to heap
        if (V[i][j] + Vp[i][j] + data->init <= crit) {
      
        num++;
        heapi[num] = i;
        heapj[num] = j;
        energy[num] = V[i][j] + Vp[i][j] + data->init;
        
        }
      }
    }
    
    //make a heap
    for (q = 2; q <= num; q++) {
      
      cur = q;
      up = cur/2;
      
      while ((energy[cur] < energy[up] && up >= 1)) {
        
        swap(heapi[cur], heapi[up]);
        swap(heapj[cur], heapj[up]);
        swap(energy[cur], energy[up]);
        cur = cur/2;
        up = up/2;
        
      }
    }
    
    //sort the heap
    for (ir = num-1; ir >= 1; ir--) {
      swap(heapi[ir+1], heapi[1]);
      swap(heapj[ir+1], heapj[1]);
      swap(energy[ir+1], energy[1]);
      
      up = 1;
      c = 2;
      
      while (c <= ir) {
        
        if (c != ir) {
          if (energy[c+1] < energy[c]) c++;
        }
        if (energy[c] < energy[up]) {
          swap(heapi[c], heapi[up]);
          swap(heapj[c], heapj[up]);
          swap(energy[c], energy[up]);
          up = c;
          c = 2*c;
          
        }
        else c = ir+1;
      }
    }
    
	//ct3->numofstructures = 0; //initialize the number of structures

    //for max number of structures input by user
    for (int cntr = num; cntr > 0 && ((ct3->GetNumberofStructures()) < maxtracebacks); cntr--) {
      
      //find next unmarked bp
      if (!mark[heapi[cntr]][N1+3+heapj[cntr]]) {
        
        iret = heapi[cntr];
        jret = heapj[cntr];
        ct3->AddStructure();
		//if (ct3->GetNumberofStructures()==1) {
		//	ct3->SetCtLabel(ct3->GetSequenceLabel(),1);
		//}
		//ct3->numofstructures++;
        //ct3->checknumberofstructures();        
        jct3 = jret + N1 + 3;

		//reset all pairs in the current structure to zero
        //for (i = 0; i <= maxseq; i++) {
        //  ct3->basepr[ct3->numofstructures][i] = 0; 
        //}
		ct3->CleanStructure(ct3->GetNumberofStructures());
        
        //perform traceback
        accessfoldtracebackV(iret, jret, jct3, ct3->GetNumberofStructures(), N1, N2, maxloop, V, Gp, ct3, data);
        accessfoldtracebackVp(iret, jret, jct3, ct3->GetNumberofStructures(), N1, N2, maxloop, Vp, Gp, ct3, data);
      
        //ct3->energy[ct3->numofstructures] = (int) energy[cntr];
		ct3->SetEnergy(ct3->GetNumberofStructures(),(int) energy[cntr]);
        
        numbp = 0;

		//check if this structure is different enough from previous structures using the idea of a window
        for (k1 = 1; k1 <= maxseq; k1++) {
          
          if (k1 < (ct3->GetPair(k1,ct3->GetNumberofStructures()))) {
            if (!(mark[k1][ct3->GetPair(k1,ct3->GetNumberofStructures())])) numbp++;
          }
        }
        for (k1 = 1; k1 <= maxseq; k1++) {
          if (k1 < (ct3->GetPair(k1,ct3->GetNumberofStructures()))) {
            mark[k1][ct3->GetPair(k1,ct3->GetNumberofStructures())] = true;
            if (windowsize > 0) {
              for (k2 = max(1, k1 - windowsize); k2 <= min(maxseq, k1 + windowsize); k2++) {
                for (k3 = max(k2, (ct3->GetPair(k1,ct3->GetNumberofStructures() )- windowsize)); k3 <= min(maxseq, (ct3->GetPair(k1,ct3->GetNumberofStructures()  )+windowsize)); k3++) {
                      
                      mark[k2][k3] = true;
                      
                    }
              }
            }
          }
        }
        
        if (numbp <= windowsize && ct3->GetNumberofStructures() > 1) {
          ct3->RemoveLastStructure();
        }
        //else {
          
        //  if (ct3->GetNumberofStructures() != 1) //strcpy(ct3->ctlabel[ct3->numofstructures], ct3->ctlabel[1]);
        //} 
        
      }
    }
        
    de_allocate (mark, maxseq+1);
    delete[] energy;
    delete[] heapi;
    delete[] heapj;
    
  //}
    
    for (int l = 0; l <= N1; l++) {
      delete[] V[l];
      delete[] Vp[l];
    }
    
    delete[] V;
    delete[] Vp;

	delete rna1;
	delete rna2;
	delete[] P;
	delete[]Gp;
    
}



void accessfoldtracebackV(int i, int j, int jct3, int k, int N1, int N2, int maxloop, double **V, double *Gp, structure *ct3, datatable *data) {
  
  bool done, found;
  int jpct3, ip, jp;
  int biggest, big, end; 
  vector< vector<bool> > inc = data->pairing;
  
  //done will switch to true when the end of this traceback is complete
  done = false; 

  while (!done) {

    //register the current base pair and find the next, if there are any
     //ct3->basepr[k][i] = jct3;
     //ct3->basepr[k][jct3] = i;
	 ct3->SetPair(i,jct3,k);
      
	 
	 
      
	 if (i!=N1&&j!=1) {
		 if (V[i][j] == data->tstack[ct3->numseq[i]][ct3->numseq[jct3]][ct3->numseq[i+1]][ct3->numseq[jct3-1]] + penalty(i, jct3, ct3, data) - Gp[i] - Gp[jct3]) {
          
			done = true;
		 }
          
     }
	 else if (i!=N1) {
		  
		 if (V[i][j] == erg4(i, jct3, i+1, 1, ct3, data, false) + penalty(i, jct3, ct3, data) - Gp[i] - Gp[jct3]) {
          
			done = true;
          
		}
	 }
	 else if (j!=1){
		 if (V[i][j] == erg4(i, jct3, jct3-1, 2, ct3, data, false) + penalty(i, jct3, ct3, data) - Gp[i] - Gp[jct3]) {
        
			done = true;
        
		 }
	 }
	 else {
		 if (V[i][j] == penalty(i, jct3, ct3, data) - Gp[i] - Gp[jct3]) {
        
			done = true;
		 }
        
	 }

	 if (!done) {
		//look for the next pair because we are not yet done tracing back  

	    //found will change to true when the next pair is found
        found = false;
        
        biggest = N1-i+j-3; 
        biggest = min(biggest, maxloop);
        
        for (int a = 0; a <= biggest && !found; a++) {
          end = i+a+1;
          big = min(end, N1);
          for (ip = max(i+1, i-j+a+3); ip <= big; ip++) {
            jp = j+ip-i-a-2;
            jpct3 = jp+N1+3;
            
            if (ip == i+1 && jp == j-1 && V[i][j] == erg1(i, jct3, ip, jpct3, ct3, data)+V[ip][jp] - Gp[i] - Gp[jct3] && 
                 inc[ct3->numseq[ip]][ct3->numseq[jpct3]]) {
              
              i++;
              j--;
              jct3--;
              found = true;
                         
            }
            else if (V[i][j] == erg2(i, jct3, ip, jpct3, ct3, data, 0, 0)+V[ip][jp] - Gp[i] - Gp[jct3] &&  
                     inc[ct3->numseq[ip]][ct3->numseq[jpct3]]) {
              
              i = ip;
              j = jp;
              jct3 = jpct3;
              found = true;
            }
          }
          
        }
        
        if (found == false) {
          found = true;
          cerr << "Error in tracebackV at "<< i << " " << jct3 << V[i][j] << "\n";
          done = true;//throw in the towel because an error was found 
        }
      }//end if !done
	}//end while !done
}
      
void accessfoldtracebackVp(int i, int j, int jct3, int k, int N1, int N2, int maxloop, double **Vp, double *Gp, structure *ct3, datatable *data) { 
  
  bool done, found;
  int jpct3, ip, jp;
  int biggest, big, end;
  vector< vector<bool> > inc = data->pairing;
  
  //done will switch to true when this traceback is complete
  done = false;
  
  while (!done) {
    
    //record the current base pair and then go look for the next
    //ct3->basepr[k][i] = jct3;
    //ct3->basepr[k][jct3] = i;
	ct3->SetPair(i,jct3,k);
    
	if (i!=1&&j!=N2) {
		if(Vp[i][j] == data->tstack[ct3->numseq[jct3]][ct3->numseq[i]][ct3->numseq[jct3+1]][ct3->numseq[i-1]] + penalty(jct3, i, ct3, data) - Gp[i] - Gp[jct3]) {
      
			done = true;
		}
      
    }
	else if (i!=1) {
		if (Vp[i][j] == erg4(jct3, i, i-1, 2, ct3, data, false) + penalty(jct3, i, ct3, data) - Gp[i] - Gp[jct3]) {
        
			done = true;
		}
        
    }
	else if (j!=N2) {
		if (Vp[i][j] == erg4(jct3, i, jct3+1, 1, ct3, data, false) + penalty(jct3, i, ct3, data)- Gp[i] - Gp[jct3]) {
      
			done = true;
		}
      
    }
	else {
		if (Vp[i][j] == penalty(jct3, i, ct3, data) - Gp[i] - Gp[jct3]) {
    
			done = true;
		}
    }
	if (!done) {
      
      found = false;
      

      biggest = N2-j+i-3;
      biggest = min(biggest, maxloop);
     
      for (int a = 0; a <= biggest && !found; a++) {
        end = j+a+1;
        big = min(end, N2); 
        for (jp = max(j+1, j-i+a+3), jpct3 = jp+N1+3; jp <= big && !found; jp++, jpct3++) {
          ip = i+jp-j-a-2;
          
          if (ip == i-1 && jp == j+1 && Vp[i][j] == erg1(jct3, i, jpct3, ip, ct3, data)+Vp[ip][jp] - Gp[i] - Gp[jct3] && 
              inc[ct3->numseq[ip]][ct3->numseq[jpct3]]) {
            
            i--;
            j++;
            jct3++;
            found = true;
            
          }
            
          else if (Vp[i][j] == erg2(jct3, i, jpct3, ip, ct3, data, 0, 0)+Vp[ip][jp] - Gp[i] - Gp[jct3] && 
                   inc[ct3->numseq[ip]][ct3->numseq[jpct3]]) {
            
            i = ip;
            j = jp;
            jct3 = jpct3;
            found = true;
            
          }
        }
      }
      
      if (found == false) {
        found = true;
        cerr << "Error in tracebackVp at " << i << " " << jct3 << " " << Vp[i][j] << "\n";
        done = true; 
      } 
    
    }//end if !done
      
  }//end while !done
}

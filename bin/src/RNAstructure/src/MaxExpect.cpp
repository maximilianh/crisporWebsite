/*
 * include necessary classes
 */
//#include "stdafx.h" //keep it unchanged
#include "structure.h"
#include "defines.h" 

#include "MaxExpectStack.h"



#include "rna_library.h"
#include "algorithm.h"
#include "pfunction.h"
#include "MaxExpect.h"
//#include "platform.h"

#include <cstring>
#include <iostream>
using namespace std;

#define maxsort 90000 //Starting point for the maximum number of basepairs within %cntrl8
							//of the minimum free energy.

#define MIN_HP_LENGTH 2  // minimum length for a hairpin loop, set small because 
						//this code is used when breaking pseudoknots as well as predicted MEA structures
#define MIN_SIB_LENGTH 7 // minimum length for a stack, internal loop, or bulge loop
#define MIN_MBWOS_LENGTH 10 // minimum length for a multi-branch w/o the stack
#define MIN_MB_LENGTH 12 // minimum length for a multi-branch

#undef DEBUG 

const double DOUBLE_DELTA = 1e-13; //the amount the doubles can differ - an estimate of precision for this calculation


/*
 * Method for executing the traceback through v and w to find
 * the optimal structure with the highest bp probability 
 */
void traceBack(structure *ct, double **vwArray, double **bpProbArray, double gamma, int ip, int jp) {
	int i, j;
	int branchPt;
	bool foundTrace = false; // used to determine if the traceback value was found
	
	
	expectMaxStack *expectMxStk;

	// create the stack to hold the pairs to process with the probability
	expectMxStk = new expectMaxStack(ct->GetSequenceLength());

	// place the 1 and N on the stack
	expectMxStk->push(ip, jp);

	// continue processing while the stack has values
	// (don't forget the v(j,i) corresponds to w(i,j)
	while(expectMxStk->pull(&i, &j))
	{
		foundTrace = false;

		#if defined DEBUG
			printf(" Stack pull: i:j %i:%i - probability: %21.17f checking if V = W\n", i, j, vwArray[i][j]);
		#endif
		
		// check for the end of a hairpin
		if ( (j-i) < MIN_HP_LENGTH-1)
		{
			#if defined DEBUG
				cout<<" Hairpin reached\n";
			#endif
		}
		// check to see if v value is equal to w
		else if (doubleEqual(vwArray[j][i], vwArray[i][j]))
		{
			// add to the bp list set to each other
			ct->SetPair(i,j,ct->GetNumberofStructures());

			

			#if defined DEBUG
				printf("   Basepair found, pushing i+1:j-1 %i:%i with a probability of: %21.17f to the stack\n",
					i+1, j-1, vwArray[i+1][j-1]);
			#endif

			expectMxStk->push(i+1, j-1);

		} // end if vArray = bp probability total, then bp was found
		// else the bases assessed are not paired
		else
		{
			#if defined DEBUG
				cout<<" Not paired, checking 5' neighbor, then 3' neighbor, then stack/branch\n";
			#endif

			// check 5' neighbor
			if (doubleEqual(vwArray[i][j], (vwArray[i+1][j]) + vwArray[i][i]) )
			{
				foundTrace = true;
				expectMxStk->push(i+1, j);
				#if defined DEBUG
					cout<<"  5' neighbor: pushing: i+1:j "<<i+1<<":"<<j<<"\n";
				#endif
			}
			// check 3' neighbor
			else if (doubleEqual(vwArray[i][j], (vwArray[i][j-1]) + vwArray[j][j]) )
			{
				foundTrace = true;
				expectMxStk->push(i, j-1);
				#if defined DEBUG
					cout<<"  3' neighbor: pushing i:j-1 "<<i<<":"<<j-1<<"\n";
				#endif
			}
			// else must be a branch
			else
			{
				foundTrace = getStructure(i, j, vwArray[i][j], vwArray, &branchPt);

				#if defined DEBUG
					printf("  Looked for multibranch on stack: [%i][%i] and received %i\n",i,j,foundTrace);
				#endif

				if (foundTrace)
				{
					#if defined DEBUG
						printf("  Placing multibranch on stack: [%i][%i] and [%i][%i]\n",i,branchPt,branchPt+1,j);
					#endif
						
					expectMxStk->push(i, branchPt);
					expectMxStk->push(branchPt+1, j);
				} // end if trace found
			} // end else branch found

			if (!foundTrace)
				cout<<"**WARNING:  Something went wrong in non-pair\n";

		} // end else the base paris assessed are not paired

	} // end while stack has values

	delete expectMxStk;

} // end traceback

/*
 * Method for executing the traceback through v and w to find
 * the optimal structure with the highest bp probability 
 */
void traceBackExternal(structure *ct, double **vwArray, double **vwPArray, double **bpProbArray, double gamma, int ip, int jp) {
	int i, j, k;
	//int branchPt;
	bool foundTrace = false; // used to determine if the traceback value was found
	
	
	expectMaxStack *expectMxStk;

	// create the stack to hold the pairs to process with the probability
	expectMxStk = new expectMaxStack(ct->GetSequenceLength());

	
	//ip must pair to jp by design
	//determine how V was reached
	

	if (ip>1&&jp<ct->GetSequenceLength()) expectMxStk->push(ip-1, jp+1);
	else if (ip>1) traceBack(ct,vwArray,bpProbArray,gamma,1,ip-1);
	else if (jp<ct->GetSequenceLength()) traceBack(ct,vwArray,bpProbArray,gamma,jp+1,ct->GetSequenceLength());
	

	// continue processing while the stack has values
	// (don't forget the v(j,i) corresponds to w(i,j)
	while(expectMxStk->pull(&i, &j))
	{
		foundTrace = false;

		#if defined DEBUG
			printf(" Stack pull: i:j %i:%i - probability: %21.17f checking if V = W\n", i, j, vwArray[i][j]);
		#endif
		
		
		
		// check to see if v value is equal to w
		if (doubleEqual(vwPArray[j][i], vwPArray[i][j]))
		{
			// add to the bp list set to each other
			ct->SetPair(i,j,ct->GetNumberofStructures());
			

			#if defined DEBUG
				printf("   Basepair found, pushing i+1:j-1 %i:%i with a probability of: %21.17f to the stack\n",
					i+1, j-1, vwPArray[i+1][j-1]);
			#endif

			if (i>1&&j<ct->GetSequenceLength()) expectMxStk->push(i-1, j+1);
			else if (i>1) traceBack(ct,vwArray,bpProbArray,gamma,1,i-1);
			else if (j<ct->GetSequenceLength()) traceBack(ct,vwArray,bpProbArray,gamma,j+1,ct->GetSequenceLength());
			

		} // end if vwPArray = bp probability total, then bp was found
		// else the bases assessed are not paired
		else
		{
			#if defined DEBUG
				cout<<" Not paired, checking 5' neighbor, then 3' neighbor, then stack/branch\n";
			#endif

			// check 5' neighbor
			if (i>1) {
				if (doubleEqual(vwPArray[i][j], vwPArray[i-1][j] + vwArray[i][i]) )
				{
					foundTrace = true;
					expectMxStk->push(i-1, j);
					#if defined DEBUG
						cout<<"  5' neighbor: pushing: i+1:j "<<i+1<<":"<<j<<"\n";
					#endif
				}
			}
			if (j<ct->GetSequenceLength()&&!foundTrace) {
			// check 3' neighbor
				if (doubleEqual(vwPArray[i][j], (vwPArray[i][j+1]) + vwArray[j][j]) )
				{
					foundTrace = true;
					expectMxStk->push(i, j+1);
					#if defined DEBUG
						cout<<"  3' neighbor: pushing i:j-1 "<<i<<":"<<j-1<<"\n";
					#endif
				}
			}
			if (i==1&&j==ct->GetSequenceLength()&&!foundTrace) {
				if (doubleEqual(vwPArray[i][j],vwArray[j][j]+vwArray[i][i])) foundTrace=true;

			}
			// else must be a branch
			if (!foundTrace)
			{

				#if defined DEBUG
					printf("  Looked for multibranch on stack: [%i][%i] and received %i\n",i,j,foundTrace);
				#endif

				for (k=2;k< i&&!foundTrace; k++) {
					if (doubleEqual(vwPArray[i][j] ,vwArray[k][i]+vwPArray[k-1][j])) {
						foundTrace = true;
						traceBack(ct,vwArray,bpProbArray,gamma,k,i);
						expectMxStk->push(k-1,j);

					}


				}

				for (k=j+1;k<ct->GetSequenceLength()&&!foundTrace;k++) {

					if (doubleEqual(vwPArray[i][j],vwArray[j][k]+vwPArray[i][k+1])) {


						foundTrace=true;
						traceBack(ct,vwArray,bpProbArray,gamma,j,k);
						expectMxStk->push(i,k+1);
				

					}
				
				}




			} // end else branch found

			if (!foundTrace)
				cout<<"**WARNING:  Something went wrong in non-pair\n";

		} // end else the base pairs assessed are not paired

	} // end while stack has values

	delete expectMxStk;

} // end traceback



//Trace is responsible for coordinating traceback of suboptimal (and optimal) structures
void trace(structure *ct, double **vwArray, double **vwPArray, double **bpProbArray, double gamma, double maxPercent, int maxStructures, int Window) {

	bool **mark;
	register int number;
	int sort;
	int i;
	int iret,jret,numbp,count,count2,k1,k2,num;
	int *heapi,*heapj;
	int cur,c,k,j,cntr;
	double crit;


	double *energy;




	//mark keeps track of which pairs have been formed by the suboptimal routine


	sort = maxsort;
	number = ct->GetSequenceLength();


	//dynamically allocate space for mark: (This keeps track of pairs that have been predicted in structures)
	mark = new bool *[number + 1];
	for (i=0;i<=(number);i++)
	mark[i] = new bool [number + 1];
	
	for (count=1;count<=(number);count++) {
		for (count2=1;count2<=(number);count2++) {
   		mark[count][count2]=false;
	   }
	}

	
	//ct->numofstructures = 0;


	//Determine the miniumum allowable score using the best score and maxPercent
	crit= ((vwArray[1][number])*((maxPercent+DOUBLE_DELTA)/100.0));//DOUBLE_DELTAis a tolerance term
	crit =  vwArray[1][number] - crit;



	energy = new double [sort+1];
	heapi = new int [sort+1];
	heapj = new int [sort+1];

	num = 0;
	for (i=1;i<ct->GetSequenceLength();i++) {
		for (j=i+MIN_HP_LENGTH-1;j<=ct->GetSequenceLength();j++) {
			if (num==sort) {
				 //allocate more space for the heap
   				delete[] heapi;
				delete[] heapj;
				delete[] energy;
				sort = 10*sort;
				heapi = new int [sort+1];
   				heapj = new int [sort+1];
   				energy = new double [sort+1];
				i = 1;
				j = 2;
				num = 0;

		   }


			//check the best score for a structure conating the i-j pair
			//Put it in the heap if the score is good enough
			if ((vwArray[j][i]+vwPArray[j][i]-2*gamma*bpProbArray[j][i])>=crit&&bpProbArray[j][i]>DOUBLE_DELTA) {

   				num++;
   				heapi[num]=i;
   				heapj[num]=j;
				energy[num] = (vwArray[j][i]+vwPArray[j][i]-2*gamma*bpProbArray[j][i]);
	   			
			}

		}//end for j
		
	}//end for i



	//sort the base pair list:


	///////////////////////////////////

	//make a heap:

	int q,up,ir;
	for (q=2;q<=num;q++) {
		cur = q;
		up = cur/2;
		while ((energy[cur]>energy[up])&&up>=1) {
			swap(heapi[cur],heapi[up]);
			swap(heapj[cur],heapj[up]);
			swap(energy[cur],energy[up]);
			cur = cur/2;
			up = up/2;
		}
	}

	//sort the heap:

	for (ir=num-1;ir>=1;ir--) {
		swap(heapi[ir+1],heapi[1]);
		swap(heapj[ir+1],heapj[1]);
		swap(energy[ir+1],energy[1]);

		up =1 ;
		c = 2;
		while (c<=ir) {
			if (c!=ir) {
				if (energy[c+1]>energy[c]) c++;
			}
			if (energy[c]>energy[up]) {
				swap(heapi[c],heapi[up]);
				swap(heapj[c],heapj[up]);
				swap(energy[c],energy[up]);
				up=c;
				c=2*c;
			}
			else c = ir+1;
		}
	}





	//Use a bool to keep track if the last structure needed to be removed:
	bool removestructure=false;

	for (cntr=num;cntr>0&&(ct->GetNumberofStructures()<maxStructures||removestructure);cntr--) {
		//This is routine to select the region of the structure to be
	   //	folded, it allows for sub-optimal structure predictions
	   //err=0;
	   //Select the next valid unmarked basepair
		if (!mark[heapi[cntr]][heapj[cntr]]) {
   		
		   iret = heapi[cntr];
		   jret = heapj[cntr];
		   
		   if (!removestructure) ct->AddStructure();
		   

		  
   		
			//zero the BP array, if the previous structure needed to be removed:
		   if(removestructure) ct->CleanStructure(ct->GetNumberofStructures()); 

		   //perform the traceback

		   
		   ct->SetPair(iret,jret,ct->GetNumberofStructures());
		   
			traceBack(ct, vwArray, bpProbArray, gamma, iret+1, jret-1);//internal fragment
		   //traceBackExternal(ct, vwArray, vwPArray, bpProbArray, gamma, 2, 72);
			traceBackExternal(ct, vwArray, vwPArray, bpProbArray, gamma, iret, jret);//external fragment

       		ct->SetEnergy(ct->GetNumberofStructures(),energy[cntr]);
				
			//count the number of new base pairs not within window of existing
            		//base pairs
        	numbp = 0;
			for (k=1;k<=number;k++) {
            	if (k<(ct->GetPair(k,ct->GetNumberofStructures()))) {
               		if (!(mark[k][ct->GetPair(k,ct->GetNumberofStructures())])) numbp++;
				}
			}
			for (k=1;k<=(number);k++) {
            	if (k<ct->GetPair(k,ct->GetNumberofStructures())) {
               		//Mark "traced back" base pairs and also base pairs
					  //	which are within a window of cntrl9
					mark[k][ct->GetPair(k,ct->GetNumberofStructures())]=true;
					if (Window>0) {
                  		for (k1=max(1,k-Window);k1<=min(number,k+Window);k1++) {
                     		for (k2=max(k1,ct->GetPair(k,ct->GetNumberofStructures())-Window);k2<=min(number,ct->GetPair(k,ct->GetNumberofStructures())+Window);k2++) {
	                        	
								mark[k1][k2]=true;
							}
						}
					}
				}
			}
			if (numbp<=Window&&ct->GetNumberofStructures()>1) {
            		removestructure=true;

				   
			}
			else {

            	//place the structure name (from ctlabel[1]) into each structure
				ct->SetCtLabel(ct->GetSequenceLabel(),ct->GetNumberofStructures());
				
					
					

			}
         		
			 
	    }//end if !mark
   		
	}//end for over cntr


	//Remove last structure, if required
	if (removestructure) {

		ct->RemoveStructure(ct->GetNumberofStructures());

	}

	de_allocate (mark,number+1);
	delete[] energy;
	delete[] heapi;
	delete[] heapj;





}//end trace


void MaxExpectFill(structure *ct, DynProgArray<PFPRECISION> *v, PFPRECISION *w5, pfdatatable *pfdata, bool *lfce, bool *mod, forceclass *fce, double maxPercent, int maxStructures, int Window, double gamma, ProgressHandler *progress) {
	double **bpProbArray; //contains the raw bp probabilities for each bp
	double *bpSSProbArray; //contains the raw single strand probability for a base
	double **vwArray;  //contains v and w recursion values
	double **vwPArray; //the v' and w' recursion values

	double *w3Array,*w5Array;//w3Array[i] is the maximum score from nucletides i to ct->GetSequenceLength()
							//w5Array[i] is the maximum score from nucleotides 1 to i

	

	//allocate main arrays and initialize the Arrays to defaults
	bpProbArray = new double *[ct->GetSequenceLength()+1];
	bpSSProbArray = new double [ct->GetSequenceLength()+1];
	vwArray = new double *[ct->GetSequenceLength()+1];
	vwPArray = new double *[ct->GetSequenceLength()+1];

	int i, j, iSmall, jBig, Length;
	
	
	//double sumPij = 0; //holds the sum of probabilities for base pairs based on a specific i over js


	for (i=0;i<=ct->GetSequenceLength();i++) {
		bpProbArray[i] = new double [ct->GetSequenceLength()+1];
		vwArray[i] = new double [ct->GetSequenceLength()+1];
		vwPArray[i] = new double [ct->GetSequenceLength()+1];

		//ct->basepr[1][i] = 0;
		bpSSProbArray[i] = 0;

		for (j=0;j<=ct->GetSequenceLength();j++) {
			bpProbArray[i][j]=0;
			vwArray[i][j]=-0;
			vwPArray[i][j]=-0;
		}
	}

	#if defined DEBUG
		w3Array = new double [ct->GetSequenceLength()+1];
		w5Array = new double [ct->GetSequenceLength()+1];
	#else
		w3Array = NULL;
		w5Array = NULL;

	#endif


		ct->nucs[0] = ' ';

	// Recursion rules investigate for vwArray
	//    1)  if the base pair (BP) probability value is 0, skip that pair
	//    2)  hairpin turns at 5 BP
	//    3)  stack/internal/bulge pairing at 7 BPs 
	//    4)  multibranching at 12 BPs (2 hairpins and a stack)
	// Because of the rules for hairpin, the probabilities 
	//    can be taken from the bpProbArray directly until j-i > 5
	
	// Calculate the single stranded probabilities for each base
	// Pi = 1 - (for all j, sum(Pij)
	// fill in w for the diagonal for the Pi,i
	for (i=1; i<=ct->GetSequenceLength(); i++)
	{

		

		bpSSProbArray[i] = 1;

		for (j=1; j<=ct->GetSequenceLength(); j++)
		{
			if (i!=j)
			{
				if (i < j)
				{
					iSmall = i;
					jBig = j;
				}
				else
				{
					iSmall = j;
					jBig = i;
				}
					
				//subtract the paired probability
				bpSSProbArray[i] = bpSSProbArray[i] - calculateprobability(iSmall, jBig, v, w5, ct, pfdata, lfce, mod, pfdata->scaling, fce);

			} //end if

		} // end for each paired base

		vwArray[i][i] = bpSSProbArray[i];
	} // end loop over each base pair


	//Calculate the base pair probabilities to start...
	for (Length = MIN_HP_LENGTH; Length <=ct->GetSequenceLength(); Length++)
	{

		//begin populating v and w along the diagonal starting with the
		//   shortest hairpin length
		for (i=1, j=i+Length-1; j<=ct->GetSequenceLength(); i++, j++)
		{
			// check for canonical pair
			// JG 20071226 - modified to just check probability rather than canonical
			bpProbArray[j][i] = calculateprobability(i, j, v, w5, ct, pfdata, lfce, mod, pfdata->scaling, fce);

		}
	}

	
	MEAFill(ct, bpProbArray, bpSSProbArray, vwArray, vwPArray, w5Array, w3Array, &pfdata->pairing, gamma, maxPercent,progress);



	// start traceback 
	trace(ct, vwArray, vwPArray, bpProbArray, gamma, maxPercent, maxStructures, Window);

	


	

	// Deallocate memory for the MaxExpect calculation
	//Arrays with functionality in the fill step
	for (i=0;i<=ct->GetSequenceLength();i++) {
		delete[] bpProbArray[i];
	}
	delete[] bpProbArray;

	delete[] bpSSProbArray;

	for (i=0; i<=ct->GetSequenceLength(); i++) {
		delete[] vwArray[i];
		delete[] vwPArray[i];
	}
	delete[] vwArray;
	delete[] vwPArray;

	#if defined DEBUG
		delete[] w5Array;
		delete[] w3Array;
	#endif


}


//This is actual fill routine for maximum expewcted accuracy structure prediction:
//bool OnlyCanonical indicates whether only Canonical pairs should be allowed
void MEAFill(structure *ct, double **bpProbArray, double *bpSSProbArray, double **vwArray, double **vwPArray, double *w5Array, double *w3Array, vector<vector<bool> >* pairing, double gamma, double maxPercent, ProgressHandler *progress, bool OnlyCanonical) {

	vector< vector<bool> >& inc = *pairing; //just an alias

//<<<<<<< MaxExpect.cpp

//	int size = 0;
	//int kSS = 0; //internal single stranded probability loop variable
//	double valueArray[maxfil];

//	double max = 0;    //holds max values from called methods

//	int i,j,Length,k;
//=======
	int Length;
//>>>>>>> 1.19

	// Populate the V Array
	// Start a for loop the will increment the length of the BP
	//    expand the lenth of BPs using the bpLength loop
	//
	// The i axis will be "horizontal" and ASCENDING
	//    the j axis will be "vertical" and ASCENDING 
	//    such that the max value will be in the [max][1] array spot
	//
	// Populate the wArray in parallel with vArray using the vwArray
	//    but switch the indexes to use the other side of the "square"
	//
	// **WARNING**
	// **WARNING** the diagonal of vwArray is used for w of the base being single stranded
	for (Length = MIN_HP_LENGTH; Length <=ct->GetSequenceLength(); Length++)
	{

		if (((Length%10)==0)&&progress!=NULL) {
			progress->update((100*Length)/(2*ct->GetSequenceLength()));
			if (progress->canceled()) break;
		}

		//begin populating v and w along the diagonal starting with the
		//   shortest hairpin length
//<<<<<<< MaxExpect.cpp
//		for (i=1, j=i+Length-1; j<=ct->GetSequenceLength(); i++, j++)
//		{
			
//=======
#ifdef SMP
#pragma omp parallel for
#endif
		for (int j=Length; j<=ct->GetSequenceLength();j++){
			int i=j-Length+1;	
			int k, size;
			int kSS = 0; //internal single stranded probability loop variable
			double max = 0;
			double valueArray[maxfil];
//>>>>>>> 1.19

			//if (!isCanonical(ct->nucs[i], ct->nucs[j]))
			if ((!inc[ct->numseq[i]][ct->numseq[j]])&&OnlyCanonical)
			{
				vwArray[j][i] = -DOUBLE_INFINITY;
				#if defined DEBUG
					printf("v fill non-canonical v[%i][%i]\n",i,j);
				#endif
			}
			else // is canonical pair or we are allowing non-canonicals
			{
				#if defined DEBUG
					printf("v fill canonical v[%i][%i]\n",i,j);
				#endif

				// call the method to calculate the probability
				// utilize gamma to adjust for pairing weight
				// JG 20071226 - moved up:  bpProbArray[j][i] = calculateprobability(i, j, v, w5, ct, pfdata, lfce, mod, scaling, fce);
	
				//*************************************************
				//  vArray logic START
				//  vArray uses the lower left of the array
				//  v = 2*gamma*probability + w(subloop)
				//*************************************************
				//if (bpLength == 2)
				//	vwArray[j][i] = 2 * gamma * bpProbArray[j][i];
				//if (bpLength == 3)
				//	vwArray[j][i] = 2 * gamma * bpProbArray[j][i] + vwArray[i+1][i+1];
				//else
					vwArray[j][i] = 2 * gamma * bpProbArray[j][i] + vwArray[i+1][j-1];

				#if defined DEBUG
				
					printf("  V[%i][%i]\t=\t%21.17f\n",i,j,vwArray[j][i]);
					printf("    Pair ProbArray[%i][%i] is: %21.17f (*2gamma = %21.17f)\n",i,j,bpProbArray[j][i],
							bpProbArray[j][i]*2*gamma);
					printf("    SSi  ProbArray[%i]     is: %21.17f\n",i, bpSSProbArray[i]);
					printf("    SSj  ProbArray    [%i] is: %21.17f\n",j, bpSSProbArray[j]);
				#endif

			} // end else was a canonical pair
			//*************************************************
			//  vArray logic END
			//*************************************************
			
			//*************************************************
			//  wArray logic START
			//  wArray uses the upper right of the array
			//*************************************************
			// use the max of vArray, or its wArray neighbors
			// or the multibranch-stack
			size = 4;
			valueArray[0] = vwArray[j][i]; // vArray i,j value
			valueArray[1] = vwArray[i+1][j] + bpSSProbArray[i]; // 5' neighbor SS
			valueArray[2] = vwArray[i][j-1] + bpSSProbArray[j]; // 3' neighbor SS
			valueArray[3] = -DOUBLE_INFINITY; 

			if (Length >= MIN_MBWOS_LENGTH)
			{

				
				for (k=i+1; k < j; k++)
				{
					if ((vwArray[i][k] + vwArray[k+1][j])>valueArray[3]) {

						valueArray[3] = vwArray[i][k] + vwArray[k+1][j];

					}

					
				} // end for multibranch choices


			} // end if mb w/o stack check

			




			///////////////////////////////

			getMax(&max, valueArray, size);

			vwArray[i][j] = max;

			#if defined DEBUG
				printf("  W[%i][%i]\t=\t%21.17f\n",i,j,vwArray[i][j]);
			#endif

			//*************************************************
			//  wArray logic END
			//*************************************************

			i++;

		} // end for population of the v and w arrays along the diagonal
	} // end increasing BP length loop

	#if defined DEBUG
		//w3[1] and w5[ct->GetSequenceLength()] should == vwArray[1][ct->GetSequenceLength()]
		//If debugging, calculate w3 and w5 to check this is true
		//Now fill w5:
		w5Array[1] = bpSSProbArray[1];
		for (i=2;i<=ct->GetSequenceLength();i++) {
			w5Array[i] = w5Array[i-1]+bpSSProbArray[i]; //add an unpaired nucleotide

			if (w5Array[i]<vwArray[i][1]) w5Array[i] = vwArray[i][1]; //check whether a whole branch is the best score

			//Now check for bifurcations
			if (i>=MIN_MBWOS_LENGTH) {
				for (j=1;j<=i-MIN_HP_LENGTH;j++) {

					if (w5Array[i]<(w5Array[j]+vwArray[i][j+1])) w5Array[i] = w5Array[j]+vwArray[i][j+1];

				}//end loop over j

			} //end if i>=MIN_MBWOS_LENGTH

		}//end loop over i

		//now fill w3:
		w3Array[ct->GetSequenceLength()]=bpSSProbArray[ct->GetSequenceLength()];
		for (i=ct->GetSequenceLength()-1;i>=1;i--) {
			w3Array[i] = w3Array[i+1] +bpSSProbArray[i]; //add an unpaired nucleotide 

			if (w3Array[i]<vwArray[ct->GetSequenceLength()][i]) w3Array[i] = vwArray[ct->GetSequenceLength()][i]; //check whether a whole branch is the best score

			//Now check for bifurcations
			if ((ct->GetSequenceLength()-i+1)>=MIN_MBWOS_LENGTH) {
				for (j=ct->GetSequenceLength();j>=i+MIN_HP_LENGTH-1;j--) {

					if (w3Array[i]<(w3Array[j]+vwArray[j-1][i])) w3Array[i] = w3Array[j]+vwArray[j-1][i];

				}//end loop over j

			} //end if i>=MIN_MBWOS_LENGTH

		}
	#endif


	// Populate the V' Array
	
	// **WARNING**
	// **WARNING** the diagonal of vwArray is used for w of the base being single stranded
	if (progress==NULL||!progress->canceled())
	for (Length = ct->GetSequenceLength(); Length >= MIN_HP_LENGTH; Length--) 
	{
//<<<<<<< MaxExpect.cpp
//		if (((Length%10)==0)&&progress!=NULL) progress->update((200*(ct->GetSequenceLength()-Length))/(2*ct->GetSequenceLength()));
//		for (i=1, j=i+Length-1; j<=ct->GetSequenceLength(); i++, j++)
//=======
		if (((Length%10)==0)&&progress!=NULL) progress->update((200*(ct->GetSequenceLength()-Length))/(2*ct->GetSequenceLength()));
#ifdef SMP
#pragma omp parallel for
#endif
		for (int j=Length;j<=ct->GetSequenceLength();j++)
//>>>>>>> 1.19
		{
			int i=j-Length+1;	
			int k, size;
			int kSS = 0; //internal single stranded probability loop variable
			double max = 0;
			double valueArray[maxfil];
			



			//if (!isCanonical(ct->nucs[i], ct->nucs[j]))
			if ((!inc[ct->numseq[i]][ct->numseq[j]])&&OnlyCanonical)
			{
				vwPArray[j][i] = -DOUBLE_INFINITY;
				#if defined DEBUG
					printf("v fill non-canonical v[%i][%i]\n",i,j);
				#endif
			}
			else // is canonical pair
			{
				#if defined DEBUG
					printf("v fill canonical v[%i][%i]\n",i,j);
				#endif

				if (i>1&&j<ct->GetSequenceLength()) {
					vwPArray[j][i] = 2 * gamma * bpProbArray[j][i] + vwPArray[i-1][j+1];
				}
				else if (i>1) {
					vwPArray[j][i] = 2 * gamma * bpProbArray[j][i] + vwArray[1][i-1];	
				
				}
				else if (j<ct->GetSequenceLength()) {
					vwPArray[j][i] = 2 * gamma * bpProbArray[j][i] + vwArray[j+1][ct->GetSequenceLength()];	
				}
				else vwPArray[j][i] = 2 * gamma * bpProbArray[j][i];


				#if defined DEBUG
				
					printf("  V[%i][%i]\t=\t%21.17f\n",i,j,vwArray[j][i]);
					printf("    Pair ProbArray[%i][%i] is: %21.17f (*2gamma = %21.17f)\n",i,j,bpProbArray[j][i],
							bpProbArray[j][i]*2*gamma);
					printf("    SSi  ProbArray[%i]     is: %21.17f\n",i, bpSSProbArray[i]);
					printf("    SSj  ProbArray    [%i] is: %21.17f\n",j, bpSSProbArray[j]);
				#endif

			} // end else was a canonical pair
			//*************************************************
			//  vArray logic END
			//*************************************************
			
			//*************************************************
			//  wArray logic START
			//  wArray uses the upper right of the array
			//*************************************************
			// use the max of vArray, or its wArray neighbors
			// or the multibranch-stack
			size = 4;
			valueArray[0] = vwPArray[j][i]; // vArray i,j value
			if (i>1) valueArray[1] = vwPArray[i-1][j] + bpSSProbArray[i]; // 5' neighbor SS
			else valueArray[1] = -DOUBLE_INFINITY;
			if (j<ct->GetSequenceLength()) valueArray[2] = vwPArray[i][j+1] + bpSSProbArray[j];// 3' neighbor SS
			else if (i==1) valueArray[2] = bpSSProbArray[j] + bpSSProbArray[i];//case where i==1 and j==N
			else valueArray[2] = -DOUBLE_INFINITY;
			valueArray[3] = -DOUBLE_INFINITY; 

			
			for (k=2; k < i; k++)
			{
				if ((vwArray[k][i] + vwPArray[k-1][j])>valueArray[3]) {

					valueArray[3] = vwArray[k][i] + vwPArray[k-1][j];

				}
				

				
			} // end for multibranch choices -1

			for (k=j+1; k < ct->GetSequenceLength(); k++)
			{
				if ((vwArray[j][k] + vwPArray[i][k+1])>valueArray[3]) {

					valueArray[3] = vwArray[j][k] + vwPArray[i][k+1];

				}
				

				
			} // end for multibranch choices -2



			///////////////////////////////

			getMax(&max, valueArray, size);

			vwPArray[i][j] = max;

			#if defined DEBUG
				printf("  W[%i][%i]\t=\t%21.17f\n",i,j,vwArray[i][j]);
			#endif

			//*************************************************
			//  wArray logic END
			//*************************************************

			i++;

		} // end for population of the v and w arrays along the diagonal
	} // end increasing BP length loop


}


/*
 * Method for recursion through the probability array (bpProbArray)
 * Input:  The array cantaining the partition function probabilities
 *         The "structure" array
 */
void bpMatch(structure *ct, char* pfsfile,
		      double gamma, double maxPercent, int maxStructures, int Window, ProgressHandler *progress) {
	//int i, j;
	short vers;


	// get the number of bases from the save file
	ifstream sav(pfsfile,ios::binary);
	read(&sav,&(vers));//read the version of the save file
		//right now there is no infrastructure to indicate the wrong version is being read. 
		//This should be changed in the future...
	int sequencelength;
	read(&sav,&(sequencelength));
	sav.close();

	ct->allocate(sequencelength);


	// define the variables for the partition function and probability methods
	DynProgArray<PFPRECISION> *w = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
	DynProgArray<PFPRECISION> *v = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
	DynProgArray<PFPRECISION> *wmb = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
	forceclass *fce = new forceclass(ct->GetSequenceLength());
	DynProgArray<PFPRECISION> *wl = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
	DynProgArray<PFPRECISION> *wlc = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
	DynProgArray<PFPRECISION> *wcoax = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
	DynProgArray<PFPRECISION> *wmbl = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
	PFPRECISION *w5 = new PFPRECISION [ct->GetSequenceLength()+1];
	PFPRECISION *w3 = new PFPRECISION [ct->GetSequenceLength()+2];
	PFPRECISION scaling;
	bool *lfce = new bool [2*ct->GetSequenceLength()+1];
	bool *mod = new bool [2*ct->GetSequenceLength()+1];
	pfdatatable *pfdata = new pfdatatable();
	datatable *data = new datatable();
	

	// read the pfs file
	readpfsave(pfsfile,ct,w5,w3,v,w,wmb,wl,wlc,wmbl,wcoax,fce,&scaling,mod,lfce,pfdata,data);

	//Run the fill routine
	MaxExpectFill(ct, v, w5, pfdata, lfce, mod, fce, maxPercent, maxStructures, Window, gamma, progress);


	// Deallocate memory for partition function primitives
	//partition function save file primitives
	delete pfdata;
	delete data;
	delete w;
	delete v;
	delete wmb;
	delete fce;
	delete wl;
	delete wcoax;
	delete wmbl;
	delete[] w5;
	delete[] w3;
	delete[] lfce;
	delete[] mod;


	

} //end probrecursion

// This method attempts to find the substructure
// Return true if the value is found
bool getStructure(int i, int j, double branchValue,
	                 double **vwArray, int *branchPt) 
{
	// look thorough the multibranch possibilities
	for ((*branchPt) = i; *branchPt <= j; (*branchPt)++)
	{
		// don't allow the branch to have a 0 value
		if ( (0 != vwArray[i][*branchPt]) && (0 != vwArray[*branchPt+1][j]) )
		{
			if( doubleEqual(branchValue, 
		             (vwArray[i][*branchPt] + vwArray[*branchPt+1][j])))
			{
				return true;
			} // end if found branch value
		} // end if not 0 branches
	} // end for multibranch choices

	return false;

} // end getStructure



// Get max returns the maximum value
// Max assumes that 0 is the minimum value
void getMax(double *max, double *valueArray, int size)
{
	// initialize variables and pre-process the inputs
	*max = -DOUBLE_INFINITY;

	if (0 != size)
	{
		*max = valueArray[0];

		for (int i=1; i < size; i++)
			*max = (valueArray[i] > *max) ? valueArray[i] : *max;
	} // end if the array size is not 0
		
} //end getMax method

// Method for scoring the structure to the correct structure
/*void scoreStructure(structure *ct, structure *correctct,
		    int *SensitivityScore, int *SensitivityScorebps,
		    int *PPVScore, int *PPVScorebps,
		    int *testpairs, int *correctpairs)
{
	int i;

	// initialize variables
	(*testpairs)=0;
	(*correctpairs)=0;
	(*SensitivityScore) = 0;
	(*SensitivityScorebps) = 0;
	(*PPVScore) = 0;
	(*PPVScorebps) = 0;


	if (correctct->GetNumberofStructures()!=1) {
		cout << "There is more than one structure in the phylogenetic ct file.";
		cout << "\nThis is not allowed.\n";
	}

	if (correctct->GetSequenceLength()!=ct->GetSequenceLength()) {
		cout << "The two ct files have different lengths.\n";
		cout << "This is not allowed.";
	}

	//count the bases in the test and phylogenetic ct:
	for (i=1;i<=ct->GetSequenceLength();i++) {
       		if (ct->basepr[1][i]>i) (*testpairs)++; //phylogenic pairs
       		if (correctct->basepr[1][i]>i) (*correctpairs)++; //test pairs
       	}

	//Now check every pair for an exact score match or base pair slip match scorebps
	for (i=1;i<=ct->GetSequenceLength();i++) {
		// PPM match numerator (is the phylogenetic base pair in the test))
		if (ct->basepr[1][i]>i) {
                        if (correctct->basepr[1][i]==ct->basepr[1][i])
			{
                		(*PPVScore)++;
                		(*PPVScorebps)++;
				#if defined DEBUG
					cout<<"S1: i: "<<i<<" ct: "<<ct->basepr[1][i]<<" correctct: "<<correctct->basepr[1][i];
				#endif
			}
                        else if ((correctct->basepr[1][i]+1)==ct->basepr[1][i]) //i slip
			{
                                (*PPVScorebps)++;
				#if defined DEBUG
					cout<<"S2: i: "<<i<<" ct: "<<ct->basepr[1][i]<<" correctct: "<<correctct->basepr[1][i];
				#endif
			}
                        else if ((correctct->basepr[1][i]-1)==ct->basepr[1][i]) //i slip
			{
                                (*PPVScorebps)++;
				#if defined DEBUG
					cout<<"S3: i: "<<i<<" ct: "<<ct->basepr[1][i]<<" correctct: "<<correctct->basepr[1][i];
				#endif
			}
                        else if ((correctct->basepr[1][i+1])==ct->basepr[1][i]) //j slip
			{
                                (*PPVScorebps)++;
				#if defined DEBUG
					cout<<"S4: i: "<<i<<" ct: "<<ct->basepr[1][i]<<" correctct: "<<correctct->basepr[1][i];
				#endif
			}
                        else if (i-1)
                        {
				#if defined DEBUG
					cout<<"S5out: i: "<<i<<" ct: "<<ct->basepr[1][i]<<" correctct: "<<correctct->basepr[1][i];
				#endif

                                if((correctct->basepr[1][i-1])==ct->basepr[1][i]) //j slip
				{
					(*PPVScorebps)++;
					#if defined DEBUG
						cout<<"S5in: i: "<<i<<" ct: "<<ct->basepr[1][i]<<" correctct: "<<correctct->basepr[1][i];
					#endif
				}
                        }
			#if defined DEBUG
				cout<<"\n";
			#endif
		} // end PPV check

		// Sensitivity numerator (is the test pair in the phylogenetic structure)
		if (correctct->basepr[1][i]>i) {
                        if (ct->basepr[1][i]==correctct->basepr[1][i])
			{
                		(*SensitivityScore)++;
                		(*SensitivityScorebps)++;
				#if defined DEBUG
					cout<<"P1: i: "<<i<<" ct: "<<ct->basepr[1][i]<<" correctct: "<<correctct->basepr[1][i];
				#endif
			}
                        else if ((ct->basepr[1][i]+1)==correctct->basepr[1][i]) //i slip
			{
                                (*SensitivityScorebps)++;
				#if defined DEBUG
					cout<<"P2: i: "<<i<<" ct: "<<ct->basepr[1][i]<<" correctct: "<<correctct->basepr[1][i];
				#endif
			}
                        else if ((ct->basepr[1][i]-1)==correctct->basepr[1][i]) //i slip
			{
                                (*SensitivityScorebps)++;
				#if defined DEBUG
					cout<<"P3: i: "<<i<<" ct: "<<ct->basepr[1][i]<<" correctct: "<<correctct->basepr[1][i];
				#endif
			}
                        else if ((ct->basepr[1][i+1])==correctct->basepr[1][i]) //j slip
			{
                                (*SensitivityScorebps)++;
				#if defined DEBUG
					cout<<"P4: i: "<<i<<" ct: "<<ct->basepr[1][i]<<" correctct: "<<correctct->basepr[1][i];
				#endif
			}
                        else if (i-1)
                        {
				#if defined DEBUG
					cout<<"P5out: i: "<<i<<" ct: "<<ct->basepr[1][i]<<" correctct: "<<correctct->basepr[1][i];
				#endif

                                if((ct->basepr[1][i-1])==correctct->basepr[1][i]) //j slip
				{
					(*SensitivityScorebps)++;
					#if defined DEBUG
						cout<<"P5in: i: "<<i<<" ct: "<<ct->basepr[1][i]<<" correctct: "<<correctct->basepr[1][i];
					#endif
				}
                        }
			#if defined DEBUG
				cout<<"\n";
			#endif

		} // end Sensitivity check
	} // end score for loop

} // end scoreStructure*/

bool isCanonical(char i, char j)
{
	bool isCanonical = false;
	// A = 1; C = 2; G = 3; U = 4
	switch (i)
	{
		case 'A':
			isCanonical = (j =='U') ? true : false;
			break;
		case 'C':
			isCanonical = (j =='G') ? true : false;
			break;
		case 'U':
			isCanonical = (j =='A') ? true :
			         ( (j =='G') ? true : false );
			break;
		case 'G':
			isCanonical = (j =='U') ? true :
			         ( (j =='C') ? true : false );
			break;
		default:
			isCanonical = false;

	} // end switch on get pairs

	return isCanonical;
} // end isCanonical

// Function for comparing doubles
// Use a constant parameter defined variable for the precision
bool doubleEqual(double double1, double double2)
{
	//const double DOUBLE_DELTA = 0.000000000000001; //the amount the doubles can differ
	const double DOUBLE_DELTA = 1e-13; //the amount the doubles can differ

	bool doubleEqual = false;
	double doubleMin, doubleMax;

	doubleMin = double1 - double1*DOUBLE_DELTA;
	doubleMax = double1 + double1*DOUBLE_DELTA;

	#if defined DEBUG
		printf("     min %21.17f : max %21.17f : 1 %21.17f : 2 %21.17f\n", doubleMin, doubleMax, double1, double2);
	#endif

	if (double2 <= doubleMax && double2 >= doubleMin)
		doubleEqual = true;

	return doubleEqual;
		
} //end doubleEqual

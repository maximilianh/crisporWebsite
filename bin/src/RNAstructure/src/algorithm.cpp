/* 	RNA Secondary Structure Prediction, Using the Algorithm of Zuker
	C++ version by David Mathews, copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006

	Addition programming by Zhi John Lu, 2005, 2006, to revise internal loop discovery to an O(N^3) algorithm.

	Programmed for Isis Pharmaceuticals, the Turner Lab, and the Mathews Lab
	Department of Biochemistry and Biophysics, University of Rochester Medical Center

	Revised on the basis of current research:
	Mathews, Sabina, Zuker, & Turner.  1999.  JMB: 288:911-940.
	Mathews, Childs, Disney, Schroeder, Zuker, & Turner.  2004. PNAS: 101:7287-7292. 
 */

#include "structure.h"
#include "algorithm.h"
#include "stackclass.h"
#include "stackstruct.h"
#ifdef _WINDOWS_GUI
#include "../Windows_interface_deprecated/platform.h"
#else 
#include "platform.h" 
#endif //_WINDOWS
#include <cstring>
#include <cmath>
#include <cstdlib>
#ifdef _CUDA_CALC_DYNALIGN_
#include "../fold-smp/frna.h"
#include "../fold-smp/fparam.h"
#include "../partition-smp/base.h"
#include "../partition-smp/util.h"
#endif
using namespace std;


//#define maxfil 250    //maximum length of file names

#define maxtloop 200 //maximum tetraloops allowed (info read from tloop)
//#define ctheaderlength 125 //maximum length of string containing info on sequence
#define numlen 8  //maximum digits in a number
#define maxnopair 600 //maximum number of bases that can be forced single
#undef debugmode //a flag to turn off debugging features
//#define debugmode //a flag to turn on debugging features
//flags for debugging
#undef timer //#define timer //flag to indicate the code execution should be timed



#ifndef INSTRUMENTED //this has to be here or NAPSS build fails. !?
void arraydump(	DynProgArray<integersize>& v,
				DynProgArray<integersize>& w,
				DynProgArray<integersize>& wmb, 
				const integersize* w5, 
				const integersize* w3, 
				const int n, 
				const char* filename	);
#endif

//********************************functions:

/*	Function efn2

	Calculates the free energy of each structure in a structure called structure.

	Structures cannot have pseudoknots

	Structnum indicates which structure to calculate the free energy
	the default, 0, indicates "all" structures

	simplemb = true means that the energy function for mb loops is linear,
	so identical to the dynamic programming algorithms.
	The default is false, so logarithmic in the number of unpaired nucleotides.

	This neew version of efn2 is now stacking aware -- if ct->stacking==true, it
	will use the stacking information stored in ct, rather than inferring stacking
	by free energy minimization.

 */
#ifndef INSTRUMENTED //Execute if pre-compiler flag INSTRUMENTED is not defined. Instrumented is only defined when compiling NAPSS.
// If INSTRUMENTED is defined, functions 'dynamic' and 'fill' will be overloaded

// Calculate the energy for all structures and optionally write details to a file
int efn2(datatable *data, structure *ct, int structnum, bool simplemb, const char *outputfilename) 
{
	ofstream out;
	if (outputfilename!=NULL) { //open a file to contain the thermodynamic details
		out.open(outputfilename);
		// if (!out.good()) return 34; Do not return here. Even if the output fails, the energies should still be calculated because the structure itself is affected (i.e. by SetEnergy)
	} 
	efn2(data, ct, structnum, simplemb, out.good()?&out:NULL);
	return (outputfilename==NULL||out.good()) ? 0 : 34; // return 0 unless the output file could not be opened. 34 means "Error opening output file for writing" in  RNA::GetErrorMessage
}

double ergexteriordiff(datatable *data, structure *ct, int structnum, bool simplemb, int min_index, int max_index)
{
	return ((float) (ergexterior(structnum, ct, data, min_index, max_index)-ergexterior(structnum, ct, data)))/conversionfactor;
}

// Calculate the energy for all structures and optionally write details to an output stream
void efn2(datatable *data, structure *ct, int structnum, bool simplemb, ostream *out)
{	
	int i, j, k, open, null, stz, count, sum, ip, jp;
	stackstruct stack;
	forceclass fce(ct->GetSequenceLength());
	int start, stop;
	
	char temp[maxfil];
	int tempsum;
	int energy;

	/*	stack = a place to keep track of where efn2 is located in a structure
		inter = indicates whether there is an intermolecular interaction involved in
		a multi-branch loop
	 */

	stack.sp = 0;  //set stack counter

	if (ct->intermolecular) {//this indicates an intermolecular folding
		for (i=0;i<3;i++) {
			forceinterefn(ct->inter[i], ct, &fce);
		}
	}
	if (structnum!=0) {
		start = structnum;
		stop = structnum;

	}
	else {
		start = 1;
		stop = ct->GetNumberofStructures();
	}



	for (count = start; count <= stop; count++) {//one structure at a time

		if (out!=NULL) { //open an output file for debugging


			*out << "Thermodynamic details for structure # "<<count<<"\n";
		}

		energy = ergexterior(count, ct, data);
		ct->SetEnergy(count, energy);

		if (out!=NULL) {
			sprintf(temp, DIGITS, (float (energy)/conversionfactor));
			//gcvt((float (ct->energy[count]))/conversionfactor,6,temp);
			*out << "\tExterior loop = "<<temp<<"\n";	
		}

		i=1;
		while (i<ct->GetSequenceLength()) 
		{ 	
			if (ct->GetPair(i, count)!=0) 
			{
				push(&stack, i, ct->GetPair(i, count), 1, 0);
				i = ct->GetPair(i, count);
			}
			i++;
		}

		subroutine://loop starts here (I chose the goto statement to make the code more readable)
		pull(&stack, &i, &j, &open, &null, &stz);//take a substructure off the stack

		while (stz != 1)
		{
			while (ct->GetPair(i, count)==j)
			{ //are i and j paired?
				if (out!=NULL) 
				{
					tempsum=0;
				}

//				while (ct->GetPair(i+1,count)==j-1) {//are i,j and i+1,j-1 stacked?
//					energy += erg1(i,j,i+1,j-1,ct,data);
				while (ct->GetPair(i+1, count)==j-1) {//are i,j and i+1,j-1 stacked?
					ct->SetEnergy(count, ct->GetEnergy(count) + erg1(i, j, i+1, j-1, ct, data));

					if (out!=NULL) {
						tempsum = tempsum + erg1(i, j, i+1, j-1, ct, data); 
						sprintf(temp, DIGITS, (float (erg1(i, j, i+1, j-1, ct, data)))/conversionfactor);
						//gcvt((float (erg1(i,j,i+1,j-1,ct,data)))/conversionfactor,6,temp);
						*out << "\t\tStack = "<<temp<<"  for stack of "<<
							i<<"-"<<j<<"\n";
					}
					i++;
					j--;
				}

				if (out!=NULL) {
					sprintf(temp, DIGITS, (float (tempsum))/conversionfactor);
					//gcvt((float (tempsum))/conversionfactor, 6, temp);
					*out << "\tHelix total = "<<temp<<"\n";	
				}

				sum = 0;
				k = i + 1;
				// now efn2 is past the paired region,  so define
				//		 the intervening non-paired segment

				while (k<j) 
				{
					if (ct->GetPair(k, count)>k)	
					{
						sum++;
						ip = k;
						k = ct->GetPair(k, count) + 1;
						jp = k-1;
					}
					else if (ct->GetPair(k, count)==0) k++;
				}

				if (sum==0) 
				{//hairpin loop		
					ct->SetEnergy(count, ct->GetEnergy(count) + erg3(i, j, ct, data, fce.f(i, j)));
					if (out!=NULL) {
						sprintf(temp, DIGITS, (float (erg3(i, j, ct, data, fce.f(i, j))))/conversionfactor);
						//gcvt((float (erg3(i, j, ct, data, fce.f(i, j))))/conversionfactor, 6, temp);
						*out << "\tHairpin = "<<temp<<"  for closure of "<<
							i<<"-"<<j<<"\n";	
					}
					goto subroutine;
				}
				else if (sum==1) 
				{//bulge/internal loop
					ct->SetEnergy(count, ct->GetEnergy(count) + erg2(i, j, ip, jp, ct, data, fce.f(i, ip), fce.f(jp, j)));
					if (out!=NULL) 
					{

						sprintf(temp, DIGITS, (float (erg2(i, j, ip, jp, ct, data, fce.f(i, ip), fce.f(jp, j))))/conversionfactor);
						//gcvt((float (erg2(i, j, ip, jp, ct, data, fce.f(i, ip), fce.f(jp, j))))/conversionfactor, 6, temp);
						*out << "\tInternal/bulge = "<<temp<<"  for closure of "<<
							i<<"-"<<j<<"\t" << "Size = " << (ip-i+j-jp-2);
						if (	(ip-i+j-jp-2) == 1) {
							//This is a bulge loop with a single bulge,  so report the stack information
							sprintf(temp, DIGITS, ((float)data->stack[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][ct->numseq[jp]])/conversionfactor);
							//gcvt(((float)data->stack[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][ct->numseq[jp]])/conversionfactor, 6, temp);
							*out << " BP stack across bulge = " << temp;
						}
						*out << "\n";	
					}
					i = ip;
					j = jp;
				}
				else 
				{//multi-branch loop
				
//					energy += ergmulti(count,  i,  ct,  data,  simplemb);
					ct->SetEnergy(count, ct->GetEnergy(count) + 
						ergmulti(count,  i,  ct,  data,  simplemb));

					if (out!=NULL) 
					{
						sprintf(temp, DIGITS, (float (ergmulti(count,  i,  ct,  data,  simplemb)))/conversionfactor);
						//gcvt((float (ergmulti(count,  i,  ct,  data,  simplemb)))/conversionfactor, 6, temp);
						*out << "\tMultiloop = "<<temp<<"  for closure of "<<
							i<<"-"<<j<<"\n";	
					}

					//put the exiting helices on the stack:
					sum++;//total helixes = sum + 1
					i++;
					for (k=1;k<sum;k++) 
					{
						while (ct->GetPair(i, count)==0) i++;
						push (&stack, i, ct->GetPair(i, count), 1, 0);
						i = ct->GetPair(i, count)+1;
					}
					goto subroutine;

				}
			}    
		}
		//Strore the energy in the underlying structure class:
//		ct->SetEnergy(count, energy);
		if (out!=NULL) {
			sprintf(temp, DIGITS, (float (ct->GetEnergy(count)))/conversionfactor);
			*out << "\n\n\tTotal energy = "<<temp<<"\n\n\n\n";
		}
	}
	return;
}

//This function adds the pair between i and j to structure ct.  
//This pair is added to the last structure structured in ct.
//This function is careful to register the pair for the real nucleotide indexes. 
//This is important because the pair could have been discovered in an exterior fragment.
void registerbasepair(structure *ct,  short i,  short j) 
{	//i and j are paired -- put that information in structure ct
	if (j <= ct->GetSequenceLength()) 
	{
		ct->SetPair(i, j, ct->GetNumberofStructures());
	}
	else if (i > ct->GetSequenceLength()) 
	{
		i = i - ct->GetSequenceLength();
		j = j - ct->GetSequenceLength();
		ct->SetPair(i, j, ct->GetNumberofStructures());
	}
	else 
	{
		j = j - ct->GetSequenceLength();
		ct->SetPair(i, j, ct->GetNumberofStructures());
	}
	return;
}

// This function, given a pair between ii and ji, will determine base pairs
// in the lowest free energy structure that contain that pair.
// Pairs are stored in structure ct at the last structure number recorded in ct.
// This returns an error code, where zero is no error and non-zero indicates a traceback error.
// If quickstructure=true, then the function knows that the fragment to traceback is from 1 to N, and the W5 array should be traced.
#ifdef DYNALIGN_II
int trace(structure *ct, datatable *data, int ii, int ji,
          DynProgArray<integersize> *v, DynProgArray<integersize> *w, DynProgArray<integersize> *wmb, DynProgArray<integersize> *w2, DynProgArray<integersize> *wmb2, 
          bool *lfce, forceclass *fce, integersize *w3, integersize *w5, bool *mod, DynProgArray<integersize> *we, integersize energy, short open, short pair, bool quickstructure)
#else
int trace(		structure *ct, 
				datatable *data, 
				int ii, 
				int ji,
				DynProgArray<integersize> *v, 
				DynProgArray<integersize> *w, 
				DynProgArray<integersize> *wmb, 
				DynProgArray<integersize> *w2, 
				DynProgArray<integersize> *wmb2, 
				bool *lfce, 
				forceclass *fce, 
				integersize *w3, 
				integersize *w5,
				bool *mod,
				bool quickstructure	)
#endif
{	
	integersize(*ptr_to_erg1)(int i, int j, int ip, int jp, structure * ct, datatable * data);
	integersize(*ptr_to_erg2)(int i, int j, int ip, int jp, structure * ct, datatable * data, char a, char b);
	integersize(*ptr_to_erg3)(int i, int j, structure * ct, datatable * data, char dbl);
	integersize(*ptr_to_erg4)(int i, int j, int ip, int jp, structure * ct, datatable * data, bool lfce);
	integersize(*ptr_to_ergcoaxflushbases) (int i, int j, int ip, int jp, structure * ct, datatable * data);
	integersize(*ptr_to_ergcoaxinterbases1) (int i, int j, int ip, int jp, structure * ct, datatable * data);
	integersize(*ptr_to_ergcoaxinterbases2) (int i, int j, int ip, int jp, structure * ct, datatable * data);
	integersize(*ptr_to_penalty) (int i, int j, structure * ct, datatable * data);

	if (ct->GetNumberofSequences() == 1) {
		ptr_to_erg1 = erg1;
		ptr_to_erg2 = erg2;
		ptr_to_erg3 = erg3;
		ptr_to_erg4 = erg4;
		ptr_to_ergcoaxflushbases = ergcoaxflushbases;
		ptr_to_ergcoaxinterbases1 = ergcoaxinterbases1;
		ptr_to_ergcoaxinterbases2 = ergcoaxinterbases2;
		ptr_to_penalty = penalty;
	}
	else if (ct->GetNumberofSequences() > 1) {
		ptr_to_erg1 = multi_erg1;
		ptr_to_erg2 = multi_erg2;
		ptr_to_erg3 = multi_erg3;
		ptr_to_erg4 = multi_erg4;
		ptr_to_ergcoaxflushbases = multi_ergcoaxflushbases;
		ptr_to_ergcoaxinterbases1 = multi_ergcoaxinterbases1;
		ptr_to_ergcoaxinterbases2 = multi_ergcoaxinterbases2;
		ptr_to_penalty = multi_penalty;
	}

	stackclass *stack;
	short number;
	short i, j, k, a, b, d;

	#ifndef DYNALIGN_II
    short open, pair;
	integersize energy;
	#else
	#endif

	bool found;
	int tracebackerror = 0;

//	#define can_pair(i, j) (data->can_pair(i, j, ct->numseq))
//	#define can_pair(i, j) (multi_can_pair(i, j, ct))
#define can_pair(i, j) (ct->can_pair(i, j, ct))

	number = ct->GetSequenceLength();
	stack = new stackclass();

	#ifdef DYNALIGN_II
	stack->push(ii, ji, open, energy, pair);
	#else

	//initialize the stack:
	if (quickstructure) 
	{
		//Tracebeck from W5:
		stack->push(1, ct->GetSequenceLength(), 1, w5[ct->GetSequenceLength()], 0);
		ct->AddStructure();
		ct->SetEnergy(1, w5[ct->GetSequenceLength()]);
		//ct->numofstructures=1;
	}
	else stack->push(ii, ji, 0, v->f(ii, ji), 1); //normal traceback from a pair

	#endif
	
	//zero all basepairs:
	//Move this to calling function...
	//if (ji<=(number)) {
	//	for (k=ii;k<=ji;k++) ct->basepr[ct->numofstructures][k] = 0;
	//}
	//else {
	//	for (k=1;k<=(ji-(number));k++) ct->basepr[ct->numofstructures][k]=0;
	//	for (k=ii;k<=(number);k++) ct->basepr[ct->numofstructures][k]=0;
	//}


	// Add forced pairs to the underlying structure ct
	for (k = 0; k < ct->GetNumberofPairs(); k++) 
	{
		ct->SetPair(ct->GetPair5(k), ct->GetPair3(k), ct->GetNumberofStructures());
	}

	while (stack->pull(&i, &j, &open, &energy, &pair))
	{
		found = true;
		// keep taking pairs off the stack until the stack is empty:
		if (pair == 0) 
		{
			found = false;
			// a bifurcated segment has been taken off the stack:
			if (open) 
			{
				// this is an exterior loop
				// so either i == 1 or j == n
				#ifndef DYNALIGN_II
				// remove nucleotides:
				if (j == number) 
				{
					while(energy == w3[i + 1] + ct->SHAPEss_give_value(i) && i < number) 
					{
						energy -= ct->SHAPEss_give_value(i);
						i++;				
						if (i == number) break;
					}
				}

				if (i == 1) 
				{
					while(energy == w5[j-1] + ct->SHAPEss_give_value(j) && j != 1) 
					{
						energy -= ct->SHAPEss_give_value(j);
						j--;
					}
				}

				if (i != number && j !=1 ) 
				{
					if (energy == v->f(i + 1, j) + ptr_to_penalty(i + 1, j, ct, data) + ptr_to_erg4(j, i + 1, i, 2, ct, data, lfce[i])) 
					{
						i++;
						stack->push(i, j, 0, v->f(i, j), 1);
						found = true;
					}
					else if ((mod[i + 1] || mod[j]) && can_pair(i + 1, j)) 
					{
						if (energy == 
								v->f(i + 2, j - 1) 
								+ ptr_to_penalty(i + 1, j, ct, data) 
								+ ptr_to_erg4(j, i + 1, i, 2, ct, data, lfce[i]) 
								+ ptr_to_erg1(i + 1, j, i + 2, j - 1, ct, data)) 
						{
							stack->push(i + 2, j - 1, 0, v->f(i + 2, j - 1), 1);
							registerbasepair(ct, i + 1, j);
							found = true;
						}
					}

					if (!found && (energy == v->f(i, j - 1) + ptr_to_penalty(i, j - 1, ct, data) + ptr_to_erg4(j - 1, i, j, 1, ct, data, lfce[j]))) 
					{
						j--;
						stack->push(i, j, 0, v->f(i, j), 1);
						found = true;
					}
					else if ((!found) && (mod[i] || mod[j-1]) && can_pair(i, j - 1)) 
					{
						if (energy == 
								v->f(i + 1, j - 2) 
								+ ptr_to_penalty(i, j - 1, ct, data) 
								+ ptr_to_erg4(j - 1, i, j, 1, ct, data, lfce[j])
								+ ptr_to_erg1(i, j - 1, i + 1, j - 2, ct, data)) 
						{
							registerbasepair(ct, i, j - 1);
							stack->push(i + 1, j - 2, 0, v->f(i + 1, j - 2), 1);
							found = true;
						}
					}

					if (!found && (energy ==
									v->f(i + 1, j - 1)
									+ ptr_to_penalty(i + 1, j - 1, ct, data)
								//	+ data->tstack[ct->numseq[j - 1]][ct->numseq[i + 1]][ct->numseq[j]][ct->numseq[i]] 
									+ multi_tstack(j - 1, i + 1, j, i, ct)
									+ checknp(lfce[i], lfce[j]) 
									+ ct->SHAPEss_give_value(i) 
									+ ct->SHAPEss_give_value(j))) 
					{
						j--;
						i++;
						stack->push(i, j, 0, v->f(i, j), 1);
						found = true;
					}
					else if (!found && (mod[i + 1] || mod[j - 1]) && can_pair(i + 1, j - 1)) 
					{
						if (energy == 
								v->f(i + 2, j - 2) 
								+ ptr_to_penalty(i + 1, j - 1, ct, data) 
							//	+ data->tstack[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
								+ multi_tstack(j - 1, i + 1, j, i, ct)
								+ checknp(lfce[i], lfce[j])
								+ ptr_to_erg1(i + 1, j - 1, i + 2, j - 2, ct, data)
								+ ct->SHAPEss_give_value(i)
								+ ct->SHAPEss_give_value(j)
							) 
						{
							registerbasepair(ct, i + 1, j - 1);
							stack->push(i + 2, j - 2, 0, v->f(i + 2, j - 2), 1);
							found = true;
						}
					}

					if (!found && (energy == v->f(i, j) + ptr_to_penalty(i , j, ct, data))) 
					{
						stack->push(i, j, 0, v->f(i, j), 1);
						found = true;
					}
					else if (!found && (mod[i] || mod[j]) && can_pair(i, j)) 
					{
						if (energy == v->f(i + 1, j - 1) + ptr_to_penalty(i + 1, j - 1, ct, data) + ptr_to_erg1(i, j, i + 1, j - 1, ct, data)) 
						{
							registerbasepair(ct, i, j);
							stack->push(i + 1, j - 1, 0, v->f(i + 1, j - 1), 1);
							found = true;
						}
					}
				}

				#else
				//start by removing nucleotides:
				while (i<j&&energy==we->f(i+1,j)) i++;              // RMW added i<j condition, removed energy assignment.
				while (i<j&&energy==we->f(i,j-1)) j--;

					//now check for a stem:	
				if (energy==v->f(i+1,j)+erg4(j,i+1,i,2,ct,data,lfce[i])+
				    penalty(i+1,j,ct,data)) {
				  
				  i++;	
				  stack->push(i,j,0,v->f(i,j),1);
				  found = true;
				}
				
				else if ((mod[i+1]||mod[j])&&can_pair(i+1,j)) {

				  if (energy==v->f(i+2,j-1)+erg4(j,i+1,i,2,ct,data,lfce[i])+
				      penalty(i+1,j,ct,data)+erg1(i+1,j,i+2,j-1,ct,data)) {

				    registerbasepair(ct,i+1,j);
				    stack->push(i+2,j-1,0,v->f(i+2,j-1),1);
				    found = true;
				  }
				}

				if (!found&&(energy==v->f(i,j-1)+erg4(j-1,i,j,1,ct,data,lfce[j])+
					     penalty(i,j-1,ct,data))) {
				  
				  j--;
				  stack->push(i,j,0,v->f(i,j),1);
				  found = true;
				}
				    else if (!found&&(mod[i]||mod[j-1])&&can_pair(i,j-1)) {
				      
				      if ((energy==v->f(i+1,j-2)+erg4(j-1,i,j,1,ct,data,lfce[j])+
					   penalty(i,j-1,ct,data)+erg1(i,j-1,i+1,j-2,ct,data))) {
					
					registerbasepair(ct,i,j-1);
					stack->push(i+1,j-2,0,v->f(i+1,j-2),1);
					found = true;
				      }
				    }

				    if (!found&&energy==v->f(i+1,j-1)+
					data->tstack[ct->numseq[j-1]][ct->numseq[i+1]]//tstkm?
					[ct->numseq[j]][ct->numseq[i]] + checknp(lfce[i],lfce[j])+
					penalty(i+1,j-1,ct,data)) {
				      
				      i++;
				      j--;
				      stack->push(i,j,0,v->f(i,j),1);
				      found = true;
				    }
				    else if (!found&&(mod[i+1]||mod[j-1])&&can_pair(i+1,j-1)) {
				      if (energy==v->f(i+2,j-2)+
					  data->tstack[ct->numseq[j-1]][ct->numseq[i+1]]
					  [ct->numseq[j]][ct->numseq[i]] + checknp(lfce[i],lfce[j])+
					  penalty(i+1,j-1,ct,data)
					  +erg1(i+1,j-1,i+2,j-2,ct,data)) {

					registerbasepair(ct,i+1,j-1);
					stack->push(i+2,j-2,0,v->f(i+2,j-2),1);
					found = true;
				      }
				    }

				    if (!found&&energy==v->f(i,j)+penalty(i,j,ct,data)) {
				      stack->push(i,j,0,v->f(i,j),1);
				      found = true;
				    }
				    else if (!found&&(mod[i]||mod[j])&&can_pair(i,j)) {
				      
				      if (energy==v->f(i+1,j-1)+penalty(i,j,ct,data)+erg1(i,j,i+1,j-1,ct,data)) {
					
					registerbasepair(ct,i,j);
					stack->push(i+1,j-1,0,v->f(i+1,j-1),1);
					found = true;
				      }
				    }

				#endif//yinghan1
			}
			else 
			{
				//(open==0), a multibranch loop
				if (ct->intermolecular) 
				{
					//intermolcular folding 

					//start by removing nucleotides:
					while (energy == w2->f(i + 1, j)) 
					{
						i++;
						energy = w2->f(i, j);
						if (i == j) break;
					}
					while (energy == w2->f(i, j - 1)) 
					{
						j--;
						energy = w2->f(i, j);
						if (i == j) break;
					}

					//following commented out by DHM (5/17/08) -- this should be handled by WMB2
					//if (fce->f(i,i)&INTER) {
					//	if (energy==w2->f(i+1,j)-INFINITE_ENERGY+data->init) {
					//		i++;
					//		energy = w2->f(i,j);
					//	}
					//}
					//if (fce->f(j,j)&INTER) {
					//	if (energy==w2->f(i,j-1)-INFINITE_ENERGY+data->init) {
					//		j--;
					//		energy = w2->f(i,j);
					//	}
					//}

					while (energy == w2->f(i + 1, j)) 
					{
						i++;
						energy = w2->f(i, j);
						if (i == j) break;
					}
					while (energy == w2->f(i, j - 1)) 
					{
						j--;
						energy = w2->f(i, j);
						if (i == j) break;
					}

					//now check for a stem:					
					if (energy == 
						v->f(i + 1, j)
						+ ptr_to_erg4(j, i + 1, i, 2, ct, data, lfce[i])
						+ ptr_to_penalty(i + 1, j, ct, data)
						) 
					{
						i++;
						stack->push(i, j, 0, v->f(i, j), 1);
						found = true;
					}

					else if ((mod[i + 1]||mod[j]) && can_pair(i + 1, j)) 
					{
						if (energy == 
							v->f(i + 2, j - 1)
							+ ptr_to_erg4(j, i + 1, i, 2, ct, data, lfce[i])
							+ ptr_to_penalty(i + 1, j, ct, data)
							+ ptr_to_erg1(i + 1, j, i + 2, j - 1, ct, data)
							) 
						{
							registerbasepair(ct, i + 1, j);
							stack->push(i + 2, j - 1, 0, v->f(i + 2, j - 1), 1);
							found = true;
						}
					}
					if (!found && (energy ==
									v->f(i, j - 1)
									+ ptr_to_erg4(j - 1, i, j, 1, ct, data, lfce[j])
									+ ptr_to_penalty(i, j - 1, ct, data))
						) 
					{
						j--;
						stack->push(i, j, 0, v->f(i, j), 1);
						found = true;
					}
					else if (!found && (mod[i] || mod[j - 1]) && can_pair(i, j - 1)) 
					{ 
						if (!found && (energy == 
										v->f(i + 1, j - 2) 
										+ ptr_to_erg4(j - 1, i, j, 1, ct, data, lfce[j])
										+ ptr_to_penalty(i, j - 1, ct, data)
										+ ptr_to_erg1(i, j - 1, i + 1, j - 2, ct, data))
							) 
						{
							registerbasepair(ct, i, j);
							stack->push(i + 1, j - 2, 0, v->f(i + 1, j - 2), 1);
							found = true;
						}
					}

					if (!found && (energy ==
									v->f(i + 1, j - 1)
								//	+ data->tstkm[ct->numseq[j - 1]][ct->numseq[i + 1]][ct->numseq[j]][ct->numseq[i]] 
									+ multi_tstkm(j - 1, i + 1, j, i, ct)
									+ ct->SHAPEss_give_value(i)
									+ ct->SHAPEss_give_value(j) 
									+ checknp(lfce[i], lfce[j])
									+ ptr_to_penalty(i + 1, j - 1, ct, data))
						) 
					{
						i++;
						j--;
						stack->push(i, j, 0, v->f(i, j), 1);
						found = true;
					}
					else if (!found && (mod[i + 1] || mod[j-1]) && can_pair(i + 1, j - 1)) 
					{
						if (!found && (energy == 
										v->f(i + 2, j - 2)
										// + data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]] 
										+ multi_tstkm(j - 1, i + 1, j, i, ct)
										+ ct->SHAPEss_give_value(i) 
										+ ct->SHAPEss_give_value(j) 
										+ checknp(lfce[i], lfce[j])
										+ ptr_to_penalty(i + 1, j - 1, ct, data)
										+ ptr_to_erg1(i + 1, j - 1, i + 2, j - 2, ct, data))
							) 
						{
							registerbasepair(ct, i + 1, j - 1);
							stack->push(i + 2, j - 2, 0, v->f(i + 2, j - 2), 1);
							found = true;
						}
					}

					if (!found && (energy == v->f(i, j) + ptr_to_penalty(i, j, ct, data))) 
					{
						stack->push(i, j, 0, v->f(i, j), 1);
						found = true;
					}
					else if (!found && (mod[i] || mod[j]) && can_pair(i, j)) 
					{
						if (!found && (energy == 
										v->f(i + 1, j - 1) 
										+ ptr_to_penalty(i, j, ct, data)
										+ ptr_to_erg1(i, j, i + 1, j - 1, ct, data))) 
						{
							registerbasepair(ct, i, j);
							stack->push(i + 1, j - 1, 0, v->f(i + 1, j - 1), 1);
							found = true;
						}
					}
				}
				if (!found) 
				{
					//start by removing nucleotides:
					while (energy == w->f(i + 1, j) + multi_eparam(6, ct) + ct->SHAPEss_give_value(i)) 
					{
						i++;
						energy = w->f(i, j);
						if (i == j) break;

					}
					while (energy == w->f(i, j - 1) + multi_eparam(6, ct) + ct->SHAPEss_give_value(j) && j!=ct->GetSequenceLength() + 1) 
					{
						j--;
						energy = w->f(i, j);
						if (i == j) break;
					}

					//now check for a stem:
					if (energy == 
									v->f(i + 1, j)
									+ ptr_to_erg4(j, i + 1, i, 2, ct, data, lfce[i])
									+ ptr_to_penalty(i + 1, j, ct, data)
									+ multi_eparam(10, ct)
									+ multi_eparam(6, ct)) 
					{
						i++;
						stack->push(i, j, 0, v->f(i, j), 1);
						found = true;
					}
					else if ((mod[i + 1] || mod[j]) && can_pair(i + 1, j)) 
					{
						if (energy == 
								v->f(i + 2, j - 1)
								+ ptr_to_erg4(j, i + 1, i, 2, ct, data, lfce[i])
								+ ptr_to_penalty(i + 1, j, ct, data) 
								+ multi_eparam(10, ct)
								+ multi_eparam(6, ct)
								+ ptr_to_erg1(i + 1, j, i + 2, j - 1, ct, data)) 
						{
							registerbasepair(ct, i + 1, j);
							stack->push(i + 2, j - 1, 0, v->f(i + 2, j - 1), 1);
							found = true;
						}
					}

					if (!found && (energy == 
									v->f(i, j - 1)
									+ ptr_to_erg4(j - 1, i, j, 1, ct, data, lfce[j])
									+ ptr_to_penalty(i, j - 1, ct, data) 
									+ multi_eparam(10, ct)
									+ multi_eparam(6, ct))) 
					{
						j--;
						stack->push(i, j, 0, v->f(i, j), 1);
						found = true;
					}
					else if (!found && (mod[i] || mod[j-1]) && can_pair(i, j-1)) 
					{
						if ((energy == 
								v->f(i + 1, j - 2)
								+ ptr_to_erg4(j - 1, i, j, 1, ct, data, lfce[j])
								+ ptr_to_penalty(i, j - 1, ct, data)
								+ multi_eparam(10, ct)
								+ multi_eparam(6, ct)
								+ ptr_to_erg1(i, j - 1, i + 1, j - 2, ct, data))) 
						{
							registerbasepair(ct, i, j - 1);
							stack->push(i + 1, j - 2, 0, v->f(i + 1, j - 2), 1);
							found = true;
						}
					}

					if (!found && energy == 
									v->f(i + 1, j - 1)
									// + data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]] 
									+ multi_tstkm(j - 1, i + 1, j, i, ct)
									+ checknp(lfce[i], lfce[j])
									+ ptr_to_penalty(i + 1, j - 1, ct, data)
									+ multi_eparam(10, ct) 
									+ 2 * multi_eparam(6, ct)
									+ ct->SHAPEss_give_value(i)
									+ ct->SHAPEss_give_value(j)) 
					{
						i++;
						j--;
						stack->push(i, j, 0, v->f(i, j), 1);
						found = true;
					}
					else if (!found && (mod[i+1] || mod[j-1]) && can_pair(i + 1, j - 1)) 
					{
						if (energy == 
								v->f(i+2, j-2)
							//	+ data->tstkm[ct->numseq[j - 1]][ct->numseq[i + 1]][ct->numseq[j]][ct->numseq[i]] 
								+ multi_tstkm(j - 1, i + 1, j, i, ct)
								+ checknp(lfce[i], lfce[j])
								+ ptr_to_penalty(i + 1, j - 1, ct, data)
								+ multi_eparam(10, ct)
								+ 2 * multi_eparam(6, ct)
								+ ptr_to_erg1(i + 1, j - 1, i + 2, j - 2, ct, data)
								+ ct->SHAPEss_give_value(i)
								+ ct->SHAPEss_give_value(j)) 
						{
							registerbasepair(ct, i + 1, j - 1);
							stack->push(i + 2, j - 2, 0, v->f(i + 2, j - 2), 1);
							found = true;
						}
					}

					if (!found && energy == v->f(i, j) + ptr_to_penalty(i, j, ct, data) + multi_eparam(10, ct)) 
					{
						stack->push(i, j, 0, v->f(i, j), 1);
						found = true;
					}
					else if (!found && (mod[i]||mod[j]) && can_pair(i, j)) 
					{
						if (energy == 
							v->f(i + 1, j - 1)
							+ ptr_to_penalty(i, j, ct, data)
							+ multi_eparam(10, ct)
							+ ptr_to_erg1(i, j, i + 1, j - 1, ct, data)) 
						{
							registerbasepair(ct, i, j);
							stack->push(i + 1, j - 1, 0, v->f(i + 1, j - 1), 1);
							found = true;
						}
					}
				}
			}
		}
		else if (pair == 2) 
		{
			found = false;
			while (energy == wmb->f(i + 1, j) + multi_eparam(6, ct) + ct->SHAPEss_give_value(i)) 
			{
				i++;
				energy = wmb->f(i, j);
				if (i == j) break;
			}
			while (energy == wmb->f(i, j - 1) + multi_eparam(6, ct) + ct->SHAPEss_give_value(j)) 
			{
				j--;
				energy = wmb->f(i, j);
				if (i == j) break;
			}
			if (ct->intermolecular) 
			{
				while (energy == w2->f(i + 1, j)) 
				{
					i++;
					energy = w2->f(i, j);
					if (i == j) break;
				}
				while (energy == w2->f(i, j - 1)) 
				{
					j--;
					energy = w2->f(i, j);
					if (i == j) break;
				}

				//Commented out by DHM (5/17/08) - Should be handled by wmb2
				//if (fce->f(i,i)&INTER) {
				//	if (energy==w2->f(i+1,j)-INFINITE_ENERGY+data->init) {
				//		i++;
				//		energy = w2->f(i,j);
				//	}
				//}
				//if (fce->f(j,j)&INTER) {
				//	if (energy==w2->f(i,j-1)-INFINITE_ENERGY+data->init) {
				//		j--;
				//		energy = w2->f(i,j);
				//	}
				//}

				while (energy == w2->f(i + 1, j)) 
				{
					i++;
					energy = w2->f(i, j);
					if (i == j) break;
				}
				while (energy == w2->f(i, j - 1)) 
				{
					j--;
					energy = w2->f(i, j);
					if (i == j) break;
				}
			}
		}
		if(pair == 2 || !found) 
		{
			//the structure must bifurcate
			//found = false;
			if (open) 
			{
				#ifndef DYNALIGN_II
				//exterior loop:
				if (i == 1) 
				{
					k = 1;
					while(k <= j - minloop && !found) 
					{
						if (energy == w5[k] + v->f(k + 1, j) + ptr_to_penalty(k + 1, j, ct, data)) 
						{
							stack->push(1, k, 1, w5[k], 0);
							stack->push(k + 1, j, 0, v->f(k + 1, j), 1);
							found = true;
						}
						else if (mod[k + 1] || mod[j] && can_pair(k + 1, j)) 
						{
							if (energy == 
									w5[k] 
									+ v->f(k + 2, j - 1)
									+ ptr_to_penalty(k + 1, j, ct, data)
									+ ptr_to_erg1(k + 1, j, k + 2, j - 1, ct, data)) 
							{
								registerbasepair(ct, k + 1, j);
								stack->push(1, k, 1, w5[k], 0);
								stack->push(k + 2, j - 1, 0, v->f(k + 2, j-1), 1);
								found = true;
							}
						}

						if (!found && (energy == 
										w5[k] 
										+ v->f(k + 2, j) 
										+ ptr_to_penalty(k + 2, j, ct, data)
										+ ptr_to_erg4(j, k + 2, k + 1, 2, ct, data, lfce[k + 1]))) 
						{
							stack->push(1, k, 1, w5[k], 0);
							stack->push(k + 2, j, 0, v->f(k + 2, j), 1);
							found = true;
						}
						else if (!found && (mod[k + 2] || mod[j]) && can_pair(k + 2, j)) 
						{
							if (energy  == 
									w5[k]
									+ v->f(k + 3, j - 1)
									+ ptr_to_penalty(k + 2, j, ct, data)
									+ ptr_to_erg4(j, k + 2, k + 1, 2, ct, data, lfce[k + 1])
									+ ptr_to_erg1(k + 2, j, k + 3, j - 2, ct, data)) 
							{
								registerbasepair(ct, k + 2, j);
								stack->push(1, k, 1, w5[k], 0);
								stack->push(k + 3, j - 1, 0, v->f(k + 3, j - 1), 1);
								found = true;
							}
						}

						if (!found && (energy == 
										w5[k] 
										+ v->f(k + 1, j - 1)
										+ ptr_to_penalty(k + 1, j - 1, ct, data)
										+ ptr_to_erg4(j - 1, k + 1, j, 1, ct, data, lfce[j]))) 
						{

							stack->push(1, k, 1, w5[k], 0);
							stack->push(k + 1, j - 1, 0, v->f(k + 1, j - 1), 1);
							found = true;
						}
						else if (!found && (mod[k + 1] || mod[j - 1]) && can_pair(k + 1, j - 1)) 
						{
							if (!found && (energy == 
											w5[k]
											+ v->f(k + 2, j - 2)
											+ ptr_to_penalty(k + 1, j - 1, ct, data)
											+ ptr_to_erg4(j - 1, k + 1, j, 1, ct, data, lfce[j])
											+ ptr_to_erg1(k + 1, j - 1, k + 2, j - 2, ct, data))) 
							{
								registerbasepair(ct, k + 1, j - 1);
								stack->push(1, k, 1, w5[k], 0);
								stack->push(k + 2, j - 2, 0, v->f(k + 2, j - 2), 1);
								found = true;
							}
						}

						if (!found && (energy == 
										w5[k] 
										+ v->f(k + 2, j - 1)
										+ ptr_to_penalty(k + 2, j - 1, ct, data)
									//	+ data->tstack[ct->numseq[j - 1]][ct->numseq[k + 2]][ct->numseq[j]][ct->numseq[k + 1]]
										+ multi_tstack(j - 1, k + 2, j, k + 1, ct)
										+ checknp(lfce[j], lfce[k + 1])
										+ ct->SHAPEss_give_value(j)
										+ ct->SHAPEss_give_value(k + 1))) 
						{
							stack->push(1, k, 1, w5[k], 0);
							stack->push(k + 2, j-1, 0, v->f(k + 2, j-1), 1);
							found = true;
						}
						else if (!found && (mod[k + 2] || mod[j - 1]) && can_pair(k + 2 ,j - 1)) 
						{
							if ((energy == 
										w5[k]
										+ v->f(k + 3, j - 2)
										+ ptr_to_penalty(k + 2, j - 1, ct, data)
								//		+ data->tstack[ct->numseq[j - 1]][ct->numseq[k + 2]][ct->numseq[j]][ct->numseq[k + 1]]
										+ multi_tstack(j - 1, k + 2, j, k + 1, ct)
										+ checknp(lfce[j], lfce[k + 1])
										+ ptr_to_erg1(k + 2, j - 1, k + 3, j - 2, ct, data)
										+ ct->SHAPEss_give_value(j)
										+ ct->SHAPEss_give_value(k + 1))) 
							{
								registerbasepair(ct, k + 2, j - 1);
								stack->push(1, k, 1, w5[k], 0);
								stack->push(k + 3, j - 2, 0, v->f(k + 3, j - 2), 1);
								found = true;
							}
						}

						if (!found) 
						{
							//check for coaxial stacking:
							a = k + minloop + 1;
							while (a < j - minloop && !found) 
							{
								if (energy == 
									w5[k - 1]
									+ v->f(k, a)
									+ v->f(a + 1, j)
									+ ptr_to_ergcoaxflushbases(k, a, a + 1, j, ct, data)
									+ ptr_to_penalty(k, a, ct, data)
									+ ptr_to_penalty(a + 1, j, ct, data)) 
								{
									if (k > 1) stack->push(1, k - 1, 1, w5[k - 1], 0);
									stack->push(k, a, 0, v->f(k, a), 1);
									stack->push(a + 1, j, 0, v->f(a + 1, j), 1);
									found = true;
								}

								else if (mod[k] || mod[a] || mod[a + 1] || mod[j]) 
								{
									if ((mod[k] || mod[a]) && (mod[a + 1] || mod[j]) && can_pair(k, a) && can_pair(a + 1, j)) 
									{

										if (energy == 
											w5[k - 1]
											+ v->f(k + 1, a - 1)
											+ v->f(a + 2, j - 1)
											+ ptr_to_ergcoaxflushbases(k, a, a + 1, j, ct, data)
											+ ptr_to_penalty(k, a, ct, data)
											+ ptr_to_penalty(a + 1, j, ct, data)
											+ ptr_to_erg1(k, a, k + 1, a - 1, ct, data)
											+ ptr_to_erg1(a + 1, j, a + 2, j - 1, ct, data)) 
										{
											registerbasepair(ct, k, a);
											registerbasepair(ct, a + 1, j);
											if (k>1) stack->push(1, k - 1, 1, w5[k - 1], 0);
											stack->push(k + 1, a - 1, 0, v->f(k+1, a-1), 1);
											stack->push(a + 2, j - 1, 0, v->f(a + 2, j - 1), 1);
											found = true;
										}
									}
									if (!found && (mod[k] || mod[a]) && can_pair(k, a)) 
									{
										if (energy == 
											w5[k - 1] 
											+ v->f(k + 1, a - 1)
											+ v->f(a + 1, j)
											+ ptr_to_ergcoaxflushbases(k, a, a + 1, j, ct, data)
											+ ptr_to_penalty(k, a, ct, data)
											+ ptr_to_penalty(a + 1, j, ct, data)
											+ ptr_to_erg1(k, a, k + 1, a-1, ct, data)) 
										{
											registerbasepair(ct, k, a);
											if (k > 1) stack->push(1, k - 1, 1, w5[k-1], 0);
											stack->push(k + 1, a - 1, 0, v->f(k + 1, a-1), 1);
											stack->push(a + 1, j, 0, v->f(a + 1, j), 1);
											found = true;
										}
									}
									if (!found && (mod[a + 1] || mod[j]) && can_pair(a + 1, j)) 
									{
										if (energy == 
											w5[k - 1]
											+ v->f(k, a)
											+ v->f(a + 2, j - 1)
											+ ptr_to_ergcoaxflushbases(k, a, a + 1, j, ct, data)
											+ ptr_to_penalty(k, a, ct, data)
											+ ptr_to_penalty(a + 1, j, ct, data)
											+ ptr_to_erg1(a + 1, j, a + 2, j - 1, ct, data)) 
										{
											registerbasepair(ct, a + 1, j);
											if (k>1) stack->push(1, k - 1, 1, w5[k - 1], 0);
											stack->push(k, a, 0, v->f(k, a), 1);
											stack->push(a + 2, j - 1, 0, v->f(a + 2, j - 1), 1);
											found = true;
										}
									}
								}

								if (!found && (energy == 
												w5[k - 1] 
												+ v->f(k, a)
												+ v->f(a + 2, j - 1)
												+ ptr_to_ergcoaxinterbases2(k, a, a + 2, j - 1, ct, data)
												+ ptr_to_penalty(k, a, ct, data)
												+ ptr_to_penalty(a + 2, j - 1, ct, data))) 
								{
									if (k > 1) stack->push(1, k - 1, 1, w5[k - 1], 0);
									stack->push(k, a, 0, v->f(k, a), 1);
									stack->push(a + 2, j - 1, 0, v->f(a + 2, j - 1), 1);
									found = true;
								}
								else if (!found && (mod[k] || mod[a] || mod[a + 2] || mod[j - 1])) 
								{
									if ((mod[k] || mod[a]) && (mod[a + 2] || mod[j - 1]) && can_pair(k, a) && can_pair(a + 2, j - 1)) 
									{
										if ((energy == 
												w5[k - 1]
												+ v->f(k + 1, a - 1)
												+ v->f(a + 3, j - 2)
												+ ptr_to_ergcoaxinterbases2(k, a, a + 2, j - 1, ct, data)
												+ ptr_to_penalty(k, a, ct, data)
												+ ptr_to_penalty(a + 2, j - 1, ct, data)
												+ ptr_to_erg1(k, a, k + 1, a - 1, ct, data)
												+ ptr_to_erg1(a + 2, j - 1, a + 3, j - 2, ct, data))) 
										{
											registerbasepair(ct, k, a);
											registerbasepair(ct, a+2, j-1);
											if (k>1) stack->push(1, k - 1, 1, w5[k - 1], 0);
											stack->push(k + 1, a - 1, 0, v->f(k + 1, a - 1), 1);
											stack->push(a + 3, j - 2, 0, v->f(a + 3, j - 2), 1);
											found = true;
										}
									}
									if (!found && (mod[k] || mod[a]) && can_pair(k, a)) 
									{
										if ((energy == 
													w5[k-1]
													+ v->f(k + 1, a - 1)
													+ v->f(a + 2, j - 1)
													+ ptr_to_ergcoaxinterbases2(k, a, a + 2, j - 1, ct, data)
													+ ptr_to_penalty(k, a, ct, data)
													+ ptr_to_penalty(a + 2, j - 1, ct, data)
													+ ptr_to_erg1(k, a, k + 1, a - 1, ct, data))) 
										{
											registerbasepair(ct, k, a);
											if (k > 1) stack->push(1, k - 1, 1, w5[k - 1], 0);
											stack->push(k + 1, a - 1, 0, v->f(k + 1, a - 1), 1);
											stack->push(a + 2, j - 1, 0, v->f(a + 2, j - 1), 1);
											found = true;
										}
									}
									if (!found && (mod[a + 2] || mod[j - 1]) && can_pair(a + 2, j - 1)) 
									{
										if ((energy == 
												w5[k - 1]
												+ v->f(k, a)
												+ v->f(a + 3, j - 2)
												+ ptr_to_ergcoaxinterbases2(k, a, a + 2, j - 1, ct, data)
												+ ptr_to_penalty(k, a, ct, data)
												+ ptr_to_penalty(a + 2, j - 1, ct, data)
												+ ptr_to_erg1(a + 2, j - 1, a + 3, j - 2, ct, data))) 
										{
											registerbasepair(ct, a + 2, j - 1);
											if (k > 1) stack->push(1, k - 1, 1, w5[k - 1], 0);
											stack->push(k, a, 0, v->f(k, a), 1);
											stack->push(a + 3, j - 2, 0, v->f(a + 3, j - 2), 1);
											found = true;
										}
									}
								}

								if (!found && (energy == 
												w5[k - 1]
												+ v->f(k + 1, a)
												+ v->f(a + 2, j)
												+ ptr_to_ergcoaxinterbases1(k + 1, a, a + 2, j, ct, data)
												+ ptr_to_penalty(k + 1, a, ct, data)
												+ ptr_to_penalty(a + 2, j, ct, data))) 
								{
									if (k > 1) stack->push(1, k - 1, 1, w5[k - 1], 0);
									stack->push(k + 1, a, 0, v->f(k + 1, a), 1);
									stack->push(a + 2, j, 0, v->f(a + 2, j), 1);
									found = true;
								}
								else if (!found && (mod[k + 1] || mod[a] || mod[a + 2] || mod[j])) 
								{
									if ((mod[k + 1] || mod[a]) && (mod[a + 2]||mod[j]) && can_pair(k + 1, a) && can_pair(a + 2, j)) 
									{
										if ((energy == 
												w5[k-1]
												+ v->f(k + 2, a - 1)
												+ v->f(a + 3, j - 1)
												+ ptr_to_ergcoaxinterbases1(k + 1, a, a + 2, j, ct, data)
												+ ptr_to_penalty(k + 1, a, ct, data)
												+ ptr_to_penalty(a + 2, j, ct, data)
												+ ptr_to_erg1(k + 1, a, k + 2, a - 1, ct, data)
												+ ptr_to_erg1(a + 2, j, a + 3, j - 1, ct, data))) 
										{
											registerbasepair(ct, k + 1, a);
											registerbasepair(ct, a + 2, j);
											if (k > 1) stack->push(1, k - 1, 1, w5[k - 1], 0);
											stack->push(k + 2, a - 1, 0, v->f(k + 2, a - 1), 1);
											stack->push(a + 3, j - 1, 0, v->f(a + 3, j - 1), 1);
											found = true;
										}
									}

									if (!found && (mod[k+1]||mod[a]) && can_pair(k+1, a)) 
									{
										if ((energy == 
												w5[k-1]
												+ v->f(k + 2, a-1)
												+ v->f(a+2, j)
												+ ptr_to_ergcoaxinterbases1(k+1, a, a+2, j, ct, data)
												+ ptr_to_penalty(k+1, a, ct, data)
												+ ptr_to_penalty(a+2, j, ct, data)
												+ ptr_to_erg1(k+1, a, k+2, a-1, ct, data))) 
										{

											registerbasepair(ct, k+1, a);
											if (k>1) stack->push(1, k-1, 1, w5[k-1], 0);
											stack->push(k+2, a-1, 0, v->f(k+2, a-1), 1);
											stack->push(a+2, j, 0, v->f(a+2, j), 1);
											found = true;
										}
									}

									if (!found && (mod[a+2] || mod[j]) && can_pair(a+2, j)) 
									{
										if ((energy == 
												w5[k-1]
												+ v->f(k+1, a)
												+ v->f(a+3, j-1)
												+ ptr_to_ergcoaxinterbases1(k+1, a, a+2, j, ct, data)
												+ ptr_to_penalty(k+1, a, ct, data)
												+ ptr_to_penalty(a+2, j, ct, data)
												+ ptr_to_erg1(a+2, j, a+3, j-1, ct, data))) 
										{
											registerbasepair(ct, a+2, j);
											if (k>1) stack->push(1, k-1, 1, w5[k-1], 0);
											stack->push(k+1, a, 0, v->f(k+1, a), 1);
											stack->push(a+3, j-1, 0, v->f(a+3, j-1), 1);
											found = true;
										}
									}
								}
								a++;
							}
						}
						k++;
					}
				}
				else 
				{
					//j==n
					k = i + minloop;
					while(k <= j && !found) 
					{
						if (energy == w3[k] + v->f(i, k-1) + ptr_to_penalty(i, k-1, ct, data)) 
						{
							stack->push(k, number, 1, w3[k], 0);
							stack->push(i, k-1, 0, v->f(i, k-1), 1);
							found=true;
						}
						else if (mod[i] || mod[k-1] && can_pair(i, k-1)) 
						{
							if (energy == 
									w3[k] 
									+ v->f(i+1, k-2)
									+ ptr_to_penalty(i, k-1, ct, data)
									+ ptr_to_erg1(i, k-1, i+1, k-2, ct, data)) 
							{
								registerbasepair(ct, i, k-1);
								stack->push(k, number, 1, w3[k], 0);
								stack->push(i+1, k-2, 0, v->f(i+1, k-2), 1);
								found=true;
							}
						}
						if (!found && (energy == 
									w3[k]
									+ v->f(i, k-2)
									+ ptr_to_penalty(i, k-2, ct, data)
									+ ptr_to_erg4(k-2, i, k-1, 1, ct, data, lfce[k-1]))) 
						{
							stack->push(k, number, 1, w3[k], 0);
							stack->push(i, k-2, 0, v->f(i, k-2), 1);
							found = true;
						}
						else if (!found && (mod[i] || mod[k-2]) && can_pair(i, k-2)) 
						{
							if (energy == 
									w3[k]
									+ v->f(i+1, k-3)
									+ ptr_to_penalty(i, k-2, ct, data)
									+ ptr_to_erg4(k-2, i, k-1, 1, ct, data, lfce[k-1])
									+ ptr_to_erg1(i, k-2, i+1, k-3, ct, data)) 
							{
								registerbasepair(ct, i, k-2);
								stack->push(k, number, 1, w3[k], 0);
								stack->push(i+1, k-3, 0, v->f(i+1, k-3), 1);
								found = true;
							}
						}

						if (!found && (energy == 
									w3[k]
									+ v->f(i+1, k-1)
									+ ptr_to_penalty(i+1, k-1, ct, data)
									+ ptr_to_erg4(k-1, i+1, i, 2, ct, data, lfce[i]))) 
						{
							stack->push(k, number, 1, w3[k], 0);
							stack->push(i+1, k-1, 0, v->f(i+1, k-1), 1);
							found = true;
						}
						else if (!found && (mod[i+1] || mod[k-1]) && can_pair(i+1, k-1)) 
						{
							if (energy == 
										w3[k]
										+ v->f(i+2, k-2)
										+ ptr_to_penalty(i+1, k-1, ct, data)
										+ ptr_to_erg4(k-1, i+1, i, 2, ct, data, lfce[i])
										+ ptr_to_erg1(i+1, k-1, i+2, k-2, ct, data)) 
							{
								registerbasepair(ct, i+1, k-1);
								stack->push(k, number, 1, w3[k], 0);
								stack->push(i+2, k-2, 0, v->f(i+2, k-2), 1);
								found=true;
							}
						}
						if (!found && (energy == 
										w3[k]
										+ v->f(i+1, k-2)
										+ ptr_to_penalty(i+1, k-2, ct, data)
									//	+ data->tstack[ct->numseq[k-2]][ct->numseq[i+1]][ct->numseq[k-1]][ct->numseq[i]]
										+ multi_tstack(k - 2, i + 1, k - 1, i, ct)
										+ checknp(lfce[k-1], lfce[i])
										+ ct->SHAPEss_give_value(k-1)
										+ ct->SHAPEss_give_value(i))) 
						{
							stack->push(k, number, 1, w3[k], 0);
							stack->push(i+1, k-2, 0, v->f(i+1, k-2), 1);
							found = true;
						}
						else if (!found && (mod[i+1]||mod[k-2]) && can_pair(i+1, k-2)) 
						{
							if (energy == 
										w3[k]
										+ v->f(i+2, k-3)
										+ ptr_to_penalty(i+1, k-2, ct, data)
									//	+ data->tstack[ct->numseq[k-2]][ct->numseq[i+1]][ct->numseq[k-1]][ct->numseq[i]]
										+ multi_tstack(k - 2, i + 1, k - 1, i, ct)
										+ checknp(lfce[k-1], lfce[i])
										+ ptr_to_erg1(i+1, k-2, i+2, k-3, ct, data)
										+ ct->SHAPEss_give_value(k-1)
										+ ct->SHAPEss_give_value(i)) 
							{
								registerbasepair(ct, i+1, k-2);
								stack->push(k, number, 1, w3[k], 0);
								stack->push(i+2, k-3, 0, v->f(i+2, k-3), 1);
								found = true;
							}
						}
						if (!found) 
						{
							//check for coaxial stacking:
							a = i + minloop + 1;
							while (a < k && !found) 
							{
								if (energy == 
										w3[k+1] 
										+ v->f(i, a)
										+ v->f(a+1, k)
										+ ptr_to_ergcoaxflushbases(i, a, a+1, k, ct, data)
										+ ptr_to_penalty(i, a, ct, data)
										+ ptr_to_penalty(a+1, k, ct, data)) 
								{
									if (k<number) stack->push(k+1, number, 1, w3[k+1], 0);
									stack->push(i, a, 0, v->f(i, a), 1);
									stack->push(a+1, k, 0, v->f(a+1, k), 1);
									found =true;
								}
								else if (mod[i] || mod[a] || mod[a+1] || mod[k]) 
								{
									if ((mod[i] || mod[a]) 
										&& (mod[a+1] || mod[k]) 
										&& can_pair(i, a)
										&& can_pair(a + 1, k)) 
									{
										if (energy == 
												w3[k+1]
												+ v->f(i+1, a-1)
												+ v->f(a+2, k-1)
												+ ptr_to_ergcoaxflushbases(i, a, a+1, k, ct, data)
												+ ptr_to_penalty(i, a, ct, data)
												+ ptr_to_penalty(a+1, k, ct, data)
												+ ptr_to_erg1(i, a, i+1, a-1, ct, data)
												+ ptr_to_erg1(a+1, k, a+2, k-1, ct, data)) 
										{
											registerbasepair(ct, i, a);
											registerbasepair(ct, a+1, k);
											if (k<number) stack->push(k+1, number, 1, w3[k+1], 0);
											stack->push(i+1, a-1, 0, v->f(i+1, a-1), 1);
											stack->push(a+2, k-1, 0, v->f(a+2, k-1), 1);
											found = true;
										}	
									}

									if ((mod[i] || mod[a]) && !found && can_pair(i, a)) 
									{
										if (energy == 
												w3[k+1]
												+ v->f(i+1, a-1)
												+ v->f(a+1, k)
												+ ptr_to_ergcoaxflushbases(i, a, a+1, k, ct, data)
												+ ptr_to_penalty(i, a, ct, data)
												+ ptr_to_penalty(a+1, k, ct, data)
												+ ptr_to_erg1(i, a, i+1, a-1, ct, data)) 
										{
											registerbasepair(ct, i, a);
											if (k<number) stack->push(k+1, number, 1, w3[k+1], 0);
											stack->push(i+1, a-1, 0, v->f(i+1, a-1), 1);
											stack->push(a+1, k, 0, v->f(a+1, k), 1);
											found = true;
										}	
									}

									if (!found && (mod[a+1] || mod[k]) && can_pair(a+1, k)) 
									{
										if (energy == 
												w3[k+1]
												+ v->f(i, a)
												+ v->f(a+2, k-1)
												+ ptr_to_ergcoaxflushbases(i, a, a+1, k, ct, data)
												+ ptr_to_penalty(i, a, ct, data)
												+ ptr_to_penalty(a+1, k, ct, data)
												+ ptr_to_erg1(a+1, k, a+2, k-1, ct, data)) 
										{
											registerbasepair(ct, a+1, k);
											if (k < number) stack->push(k+1, number, 1, w3[k+1], 0);
											stack->push(i, a, 0, v->f(i, a), 1);
											stack->push(a+2, k-1, 0, v->f(a+2, k-1), 1);
											found = true;
										}	
									}
								}
								if (!found && (energy == 
											w3[k+1]
											+ v->f(i, a)											
											+ v->f(a+2, k-1)
											+ ptr_to_ergcoaxinterbases2(i, a, a+2, k-1, ct, data)
											+ ptr_to_penalty(i, a, ct, data)
											+ ptr_to_penalty(a+2, k-1, ct, data))) 
								{

									if (k < number) stack->push(k+1, number, 1, w3[k+1], 0);
									stack->push(i, a, 0, v->f(i, a), 1);
									stack->push(a+2, k-1, 0, v->f(a+2, k-1), 1);
									found = true;
								}
								else if (!found&&(mod[i]||mod[a]||mod[a+2]||mod[k-1])) 
								{
									if ((mod[i] || mod[a]) 
										&& (mod[a+2] || mod[k - 1])										 
										&& can_pair(i, a)
										&& can_pair(a+2, k-1))
									{
										if ((energy == 
													w3[k+1]
													+ v->f(i+1, a-1)
													+ v->f(a+3, k-2)
													+ ptr_to_ergcoaxinterbases2(i, a, a+2, k-1, ct, data)
													+ ptr_to_penalty(i, a, ct, data)
													+ ptr_to_penalty(a+2, k-1, ct, data)
													+ ptr_to_erg1(i, a, i+1, a-1, ct, data)
													+ ptr_to_erg1(a+2, k-1, a+3, k-2, ct, data))) 
										{
											registerbasepair(ct, i, a);
											registerbasepair(ct, a+2, k-1);
											if (k < number) stack->push(k+1, number, 1, w3[k+1], 0);
											stack->push(i+1, a-1, 0, v->f(i+1, a-1), 1);
											stack->push(a+3, k-2, 0, v->f(a+3, k-2), 1);
											found = true;
										}
									}

									if ((mod[i] || mod[a]) && !found && can_pair(i, a)) 
									{
										if ((energy == 
												w3[k+1]
												+ v->f(i+1, a-1)
												+ v->f(a+2, k-1)
												+ ptr_to_ergcoaxinterbases2(i, a, a+2, k-1, ct, data)
												+ ptr_to_penalty(i, a, ct, data)
												+ ptr_to_penalty(a+2, k-1, ct, data)
												+ ptr_to_erg1(i, a, i+1, a-1, ct, data)
										    )) 
										{
											registerbasepair(ct, i, a);
											if (k<number) stack->push(k+1, number, 1, w3[k+1], 0);
											stack->push(i+1, a-1, 0, v->f(i+1, a-1), 1);
											stack->push(a+2, k-1, 0, v->f(a+2, k-1), 1);
											found = true;
										}
									}

									if ((!found) && (mod[a+2] || mod[k-1]) && can_pair(a+2, k-1)) 
									{
										if ((energy == 
												w3[k+1]
												+ v->f(i, a)
												+ v->f(a+3, k-2)
												+ ptr_to_ergcoaxinterbases2(i, a, a+2, k-1, ct, data)
												+ ptr_to_penalty(i, a, ct, data)
												+ ptr_to_penalty(a+2, k-1, ct, data)
												+ ptr_to_erg1(a+2, k-1, a+3, k-2, ct, data))) 
										{
											registerbasepair(ct, a+2, k-1);
											if (k<number) stack->push(k+1, number, 1, w3[k+1], 0);
											stack->push(i, a, 0, v->f(i, a), 1);
											stack->push(a+3, k-2, 0, v->f(a+3, k-2), 1);
											found = true;
										}
									}
								}

								if (!found && (energy == 
											w3[k+1]
											+ v->f(i+1, a)
											+ v->f(a+2, k)
											+ ptr_to_ergcoaxinterbases1(i+1, a, a+2, k, ct, data)
											+ ptr_to_penalty(i+1, a, ct, data)
											+ ptr_to_penalty(a+2, k, ct, data)))									    
								{
									if (k < number) stack->push(k+1, number, 1, w3[k+1], 0);
									stack->push(i+1, a, 0, v->f(i+1, a), 1);
									stack->push(a+2, k, 0, v->f(a+2, k), 1);
									found = true;
								}
								else if (!found&&(mod[i+1]||mod[a]||mod[a+2]||mod[k])) 
								{
									if ((mod[i+1]||mod[a])&&(mod[a+2]||mod[k])
											&&can_pair(i+1, a)
											&&can_pair(a+2, k)) 
									{
										if ((energy == 
													w3[k+1]
													+ v->f(i+2, a-1)
													+ v->f(a+3, k-1)
													+ ptr_to_ergcoaxinterbases1(i+1, a, a+2, k, ct, data)
													+ ptr_to_penalty(i+1, a, ct, data)
													+ ptr_to_penalty(a+2, k, ct, data)
													+ ptr_to_erg1(i+1, a, i+2, a-1, ct, data)
													+ ptr_to_erg1(a+2, k, a+3, k-1, ct, data))) 
										{
											registerbasepair(ct, i+1, a);
											registerbasepair(ct, a+2, k);
											if (k < number) stack->push(k+1, number, 1, w3[k+1], 0);
											stack->push(i+2, a-1, 0, v->f(i+2, a-1), 1);
											stack->push(a+3, k-1, 0, v->f(a+3, k-1), 1);
											found = true;
										}	
									}

									if ((mod[i+1] || mod[a]) && (!found) && can_pair(i+1, a)) 
									{
										if ((energy == 
													w3[k+1]
													+ v->f(i+2, a-1)
													+ v->f(a+2, k)
													+ ptr_to_ergcoaxinterbases1(i+1, a, a+2, k, ct, data)
													+ ptr_to_penalty(i+1, a, ct, data)
													+ ptr_to_penalty(a+2, k, ct, data)
													+ ptr_to_erg1(i+1, a, i+2, a-1, ct, data))) 
										{
											registerbasepair(ct, i+1, a);
											if (k<number) stack->push(k+1, number, 1, w3[k+1], 0);
											stack->push(i+2, a-1, 0, v->f(i+2, a-1), 1);
											stack->push(a+2, k, 0, v->f(a+2, k), 1);
											found = true;
										}	
									}

									if ((!found) && (mod[a+2]||mod[k]) && can_pair(a+2, k)) 
									{
										if ((energy ==
													w3[k+1]
													+ v->f(i+1, a)
													+ v->f(a+3, k-1)
													+ ptr_to_ergcoaxinterbases1(i+1, a, a+2, k, ct, data)
													+ ptr_to_penalty(i+1, a, ct, data)
													+ ptr_to_penalty(a+2, k, ct, data)
													+ ptr_to_erg1(a+2, k, a+3, k-1, ct, data))) 
										{
											registerbasepair(ct, a+2, k);
											if (k<number) stack->push(k+1, number, 1, w3[k+1], 0);
											stack->push(i+1, a, 0, v->f(i+1, a), 1);
											stack->push(a+3, k-1, 0, v->f(a+3, k-1), 1);
											found = true;
										}	
									}
								}
								a++;
							}
						}
						k++;
					}
				}
				
				#else
              			  k=i+1;
			  while(k<j&&!found) {
					if (energy==we->f(i,k)+we->f(k+1,j)) {
						stack->push(i,k,1,we->f(i,k),0);
						stack->push(k+1,j,1,we->f(k+1,j),0);
						found = true;
					}

					//also check for coaxial stacking:
					else if (energy==v->f(i,k)+v->f(k+1,j)+penalty(i,k,ct,data)+
					         penalty(k+1,j,ct,data)
					         +ergcoaxflushbases(i,k,k+1,j,ct,data)) {

					  stack->push(i,k,0,v->f(i,k),1);
					  stack->push(k+1,j,0,v->f(k+1,j),1);
					  found = true;
					}

					else if (mod[i]||mod[k]||mod[k+1]||mod[j]) {

					  if ((mod[i]||mod[k])&&(mod[k+1]||mod[j])&&can_pair(i,k)
					      &&can_pair(k+1,j)) {

					    if (energy==v->f(i+1,k-1)+v->f(k+2,j-1)+penalty(i,k,ct,data)+
						penalty(k+1,j,ct,data)
						+ergcoaxflushbases(i,k,k+1,j,ct,data)
						+ erg1(i,k,i+1,k-1,ct,data)
						+ erg1(k+1,j,k+2,j-1,ct,data)) {

					      registerbasepair(ct,i,k);
					      registerbasepair(ct,k+1,j);
					      stack->push(i+1,k-1,0,v->f(i+1,k-1),1);
					      stack->push(k+2,j-1,0,v->f(k+2,j-1),1);
					      found = true;
					    }		
					  }

					  if ((mod[i]||mod[k])&&(!found)&&can_pair(i,k)) {
					    
					    if (energy==v->f(i+1,k-1)+v->f(k+1,j)+penalty(i,k,ct,data)+
						penalty(k+1,j,ct,data)
						+ergcoaxflushbases(i,k,k+1,j,ct,data)
						+erg1(i,k,i+1,k-1,ct,data)) {

					      registerbasepair(ct,i,k);
					      
					      stack->push(i+1,k-1,0,v->f(i+1,k-1),1);
					      stack->push(k+1,j,0,v->f(k+1,j),1);
					      found = true;
					    }		
					  }
					  if ((!found)&&(mod[k+1]||mod[j])&&can_pair(k+1,j)) {

					    if (energy==v->f(i,k)+v->f(k+2,j-1)+penalty(i,k,ct,data)+
						penalty(k+1,j,ct,data)
						+ergcoaxflushbases(i,k,k+1,j,ct,data)
						+erg1(k+1,j,k+2,j-1,ct,data)) {

							
					      registerbasepair(ct,k+1,j);
					      stack->push(i,k,0,v->f(i,k),1);
					      stack->push(k+2,j-1,0,v->f(k+2,j-1),1);
					      found = true;
					    }		
					  }
					}
					if (!found&&(energy==v->f(i,k)+v->f(k+2,j-1)+penalty(i,k,ct,data)+
					             penalty(k+2,j-1,ct,data)
					             +ergcoaxinterbases2(i,k,k+2,j-1,ct,data))) {

					  stack->push(i,k,0,v->f(i,k),1);
					  stack->push(k+2,j-1,0,v->f(k+2,j-1),1);
					  found = true;
					}
					else if (!found&&(mod[i]||mod[k]||mod[k+2]||mod[j-1])) {

					  if ((mod[i]||mod[k])&&(mod[k+2]||mod[j-1])&&can_pair(i,k)
					      &&can_pair(k+2,j-1)) {

					    if ((energy==v->f(i+1,k-1)+v->f(k+3,j-2)+penalty(i,k,ct,data)+
						 penalty(k+2,j-1,ct,data)
						 +ergcoaxinterbases2(i,k,k+2,j-1,ct,data)
						 +erg1(i,k,i+1,k-1,ct,data)
						 +erg1(k+2,j-1,k+3,j-2,ct,data))) {

					      registerbasepair(ct,i,k);
					      registerbasepair(ct,k+2,j-1);
					      stack->push(i+1,k-1,0,v->f(i+1,k-1),1);
					      stack->push(k+3,j-2,0,v->f(k+3,j-2),1);
					      found = true;
					    }	
					  }

					  if ((mod[i]||mod[k])&&(!found)&&can_pair(i,k)) {
					    
					    if ((energy==v->f(i+1,k-1)+v->f(k+2,j-1)+penalty(i,k,ct,data)+
						 penalty(k+2,j-1,ct,data)
						 +ergcoaxinterbases2(i,k,k+2,j-1,ct,data)
						 +erg1(i,k,i+1,k-1,ct,data))) {
					      
					      registerbasepair(ct,i,k);
								
					      stack->push(i+1,k-1,0,v->f(i+1,k-1),1);
					      stack->push(k+2,j-1,0,v->f(k+2,j-1),1);
					      found = true;
					    }	
					  }

					  if ((!found)&&(mod[k+2]||mod[j-1])&&can_pair(k+2,j-1)) {
					    
					    if ((energy==v->f(i,k)+v->f(k+3,j-2)+penalty(i,k,ct,data)+
						 penalty(k+2,j-1,ct,data)+
						 ergcoaxinterbases2(i,k,k+2,j-1,ct,data)
						 +erg1(k+2,j-1,k+3,j-2,ct,data))) {

								
					      registerbasepair(ct,k+2,j-1);
					      stack->push(i,k,0,v->f(i,k),1);
					      stack->push(k+3,j-2,0,v->f(k+3,j-2),1);
					      found = true;
					    }	
					  }
					}
					if (!found&&(energy==v->f(i+1,k)+v->f(k+2,j)+penalty(i+1,k,ct,data)+
					             penalty(k+2,j,ct,data)+
					             ergcoaxinterbases1(i+1,k,k+2,j,ct,data))) {

					  stack->push(i+1,k,0,v->f(i+1,k),1);
					  stack->push(k+2,j,0,v->f(k+2,j),1);
					  found = true;
					}

					else if (!found&&(mod[i+1]||mod[k]||mod[k+2]||mod[j])) {

					  if ((mod[i+1]||mod[k])&&(mod[k+2]||mod[j])&&can_pair(i+1,k)
					      &&can_pair(k+2,j)) {
					    if ((energy==v->f(i+2,k-1)+v->f(k+3,j-1)+penalty(i+1,k,ct,data)+
						 penalty(k+2,j,ct,data)+
						 ergcoaxinterbases1(i+1,k,k+2,j,ct,data)
						 +erg1(i+1,k,i+2,k-1,ct,data)
						 +erg1(k+2,j,k+3,j-1,ct,data))) {

					      registerbasepair(ct,i+1,k);
					      registerbasepair(ct,k+2,j);
					      stack->push(i+2,k-1,0,v->f(i+2,k-1),1);
					      stack->push(k+3,j-1,0,v->f(k+3,j-1),1);
					      found = true;
					    }
					  }

					  if ((mod[i+1]||mod[k])&&(!found)&&can_pair(i+1,k)) {
					    if ((energy==v->f(i+2,k-1)+v->f(k+2,j)+penalty(i+1,k,ct,data)+
						 penalty(k+2,j,ct,data)+
						 ergcoaxinterbases1(i+1,k,k+2,j,ct,data)
						 +erg1(i+1,k,i+2,k-1,ct,data))) {

					      registerbasepair(ct,i+1,k);
								
					      stack->push(i+2,k-1,0,v->f(i+2,k-1),1);
					      stack->push(k+2,j,0,v->f(k+2,j),1);
					      found = true;
					    }
					  }

					  if ((!found)&&(mod[k+2]||mod[j])&&can_pair(k+2,j)) {
					    if ((energy==v->f(i+1,k)+v->f(k+3,j-1)+penalty(i+1,k,ct,data)+
						 penalty(k+2,j,ct,data)+
						 ergcoaxinterbases1(i+1,k,k+2,j,ct,data)
								
						 +erg1(k+2,j,k+3,j-1,ct,data))) {

					      
					      registerbasepair(ct,k+2,j);
					      stack->push(i+1,k,0,v->f(i+1,k),1);
					      stack->push(k+3,j-1,0,v->f(k+3,j-1),1);
					      found = true;
					    }
					  }
					}
					k++;
			  }

				#endif                  
			}
			else 
			{
				//open == 0, multiloop
				k = i + 1;
				while(k < j && !found) 
				{
					if (energy == w->f(i, k) + w->f(k + 1, j)) 
					{
						stack->push(i, k, 0, w->f(i, k), 0);
						stack->push(k + 1, j, 0, w->f(k + 1, j), 0);
						found = true;
					}
					//also check for coaxial stacking:
					else if (energy == 
								v->f(i, k) 
								+ v->f(k + 1, j)
								+ ptr_to_penalty(i, k, ct, data)
								+ ptr_to_penalty(k + 1, j, ct, data)
								+ 2 * multi_eparam(10, ct)
								+ ptr_to_ergcoaxflushbases(i, k, k + 1, j, ct, data)) 
					{
						stack->push(i, k, 0, v->f(i, k), 1);
						stack->push(k + 1, j, 0, v->f(k + 1, j), 1);
						found = true;
					}
					else if (mod[i] || mod[k] || mod[k + 1] || mod[j]) 
					{
						if ((mod[i] || mod[k]) && (mod[k + 1] || mod[j]) && can_pair(i, k) && can_pair(k + 1, j)) 
						{
							if (energy == 
									v->f(i + 1, k - 1)
									+ v->f(k + 2, j-1)
									+ ptr_to_penalty(i, k, ct, data)
									+ ptr_to_penalty(k + 1, j, ct, data)
									+ 2 * multi_eparam(10, ct)
									+ ptr_to_ergcoaxflushbases(i, k, k + 1, j, ct, data)
									+ ptr_to_erg1(i, k, i + 1, k - 1, ct, data)
									+ ptr_to_erg1(k + 1, j, k + 2, j - 1, ct, data)) 
							{
								registerbasepair(ct, i, k);
								registerbasepair(ct, k+1, j);
								stack->push(i+1, k-1, 0, v->f(i+1, k-1), 1);
								stack->push(k+2, j-1, 0, v->f(k+2, j-1), 1);
								found = true;
							}		
						}
						if ((mod[i] || mod[k]) && (!found) && can_pair(i, k)) 
						{
							if (energy == 
									v->f(i+1, k-1)
									+ v->f(k+1, j)
									+ ptr_to_penalty(i, k, ct, data)
									+ ptr_to_penalty(k+1, j, ct, data)
									+ 2 * multi_eparam(10, ct)
									+ ptr_to_ergcoaxflushbases(i, k, k+1, j, ct, data)
									+ ptr_to_erg1(i, k, i+1, k-1, ct, data)) 
							{
								registerbasepair(ct, i, k);
								stack->push(i+1, k-1, 0, v->f(i+1, k-1), 1);
								stack->push(k+1, j, 0, v->f(k+1, j), 1);
								found = true;
							}		
						}
						if ((!found) && (mod[k+1] || mod[j]) && can_pair(k + 1, j)) 
						{
							if (energy == 
								v->f(i, k)
								+ v->f(k + 2, j - 1)
								+ ptr_to_penalty(i, k, ct, data)
								+ ptr_to_penalty(k + 1, j, ct, data)
								+ 2 * multi_eparam(10, ct)
								+ ptr_to_ergcoaxflushbases(i, k, k + 1, j, ct, data)
								+ ptr_to_erg1(k + 1, j, k + 2, j - 1, ct, data)) 
							{
								registerbasepair(ct, k+1, j);
								stack->push(i, k, 0, v->f(i, k), 1);
								stack->push(k+2, j-1, 0, v->f(k+2, j-1), 1);
								found = true;
							}		
						}
					}
					if (!found && (energy == 
									v->f(i, k)
									+ v->f(k+2, j-1)
									+ ptr_to_penalty(i, k, ct, data)
									+ ptr_to_penalty(k+2, j-1, ct, data)
									+ 2 * multi_eparam(10, ct)
									+ 2 * multi_eparam(6, ct)
									+ ptr_to_ergcoaxinterbases2(i, k, k+2, j-1, ct, data))) 
					{
						stack->push(i, k, 0, v->f(i, k), 1);
						stack->push(k+2, j-1, 0, v->f(k+2, j-1), 1);
						found = true;
					}
					else if (!found && (mod[i] || mod[k] || mod[k+2] || mod[j-1])) 
					{
						if ((mod[i] || mod[k]) && (mod[k+2] || mod[j-1]) && can_pair(i, k) && can_pair(k+2, j-1)) 
						{
							if ((energy == 
									v->f(i+1,  k-1)
									+ v->f(k+3, j-2)
									+ ptr_to_penalty(i, k, ct, data)
									+ ptr_to_penalty(k+2, j-1, ct, data)
								    + 2 * multi_eparam(10, ct)
								    + 2 * multi_eparam(6, ct)
								    + ptr_to_ergcoaxinterbases2(i, k, k+2, j-1, ct, data)
									+ ptr_to_erg1(i, k, i+1, k-1, ct, data)
									+ ptr_to_erg1(k+2, j-1, k+3, j-2, ct, data))) 
							{
								registerbasepair(ct, i, k);
								registerbasepair(ct, k+2, j-1);
								stack->push(i+1, k-1, 0, v->f(i+1, k-1), 1);
								stack->push(k+3, j-2, 0, v->f(k+3, j-2), 1);
								found = true;
							}	
						}
						if ((mod[i]||mod[k])&&(!found)&&can_pair(i, k)) 
						{
							if ((energy == 
									v->f(i+1, k-1)
									+ v->f(k+2, j-1)
									+ ptr_to_penalty(i, k, ct, data)
									+ ptr_to_penalty(k+2, j-1, ct, data)
									+ 2*multi_eparam(10, ct)
									+ 2*multi_eparam(6, ct)
								    + ptr_to_ergcoaxinterbases2(i, k, k+2, j-1, ct, data)
									+ ptr_to_erg1(i, k, i+1, k-1, ct, data))) 
							{
								registerbasepair(ct, i, k);
								stack->push(i+1, k-1, 0, v->f(i+1, k-1), 1);
								stack->push(k+2, j-1, 0, v->f(k+2, j-1), 1);
								found = true;
							}	
						}
						if ((!found) && (mod[k+2]||mod[j-1]) && can_pair(k+2, j-1)) 
						{
							if ((energy == 
									v->f(i, k)
									+ v->f(k+3, j-2)
									+ ptr_to_penalty(i, k, ct, data)
									+ ptr_to_penalty(k+2, j-1, ct, data)
									+ 2*multi_eparam(10, ct)
									+ 2*multi_eparam(6, ct)
									+ ptr_to_ergcoaxinterbases2(i, k, k+2, j-1, ct, data)
									+ ptr_to_erg1(k+2, j-1, k+3, j-2, ct, data))) 
							{
								registerbasepair(ct, k+2, j-1);
								stack->push(i, k, 0, v->f(i, k), 1);
								stack->push(k+3, j-2, 0, v->f(k+3, j-2), 1);
								found = true;
							}	
						}
					}

					if (!found && (energy == 
									v->f(i+1, k) 
									+ v->f(k+2, j) 
									+ ptr_to_penalty(i+1, k, ct, data)
									+ ptr_to_penalty(k+2, j, ct, data)
									+ 2*multi_eparam(10, ct)
									+ 2*multi_eparam(6, ct)
									+ ptr_to_ergcoaxinterbases1(i+1, k, k+2, j, ct, data))) 
					{
						stack->push(i+1, k, 0, v->f(i+1, k), 1);
						stack->push(k+2, j, 0, v->f(k+2, j), 1);
						found = true;
					}
					else if (!found && (mod[i+1]||mod[k]||mod[k+2]||mod[j])) 
					{
						if ((mod[i+1]||mod[k])&&(mod[k+2]||mod[j])&&can_pair(i+1, k) &&can_pair(k+2, j)) 
						{
							if ((energy == 
									v->f(i+2, k-1)
									+ v->f(k+3, j-1)
									+ ptr_to_penalty(i+1, k, ct, data)
									+ ptr_to_penalty(k+2, j, ct, data)
									+ 2*multi_eparam(10, ct)
									+ 2*multi_eparam(6, ct)
									+ ptr_to_ergcoaxinterbases1(i+1, k, k+2, j, ct, data)
									+  ptr_to_erg1(i+1, k, i+2, k-1, ct, data)
									+  ptr_to_erg1(k+2, j, k+3, j-1, ct, data))) 
							{
								registerbasepair(ct, i+1, k);
								registerbasepair(ct, k+2, j);
								stack->push(i+2, k-1, 0, v->f(i+2, k-1), 1);
								stack->push(k+3, j-1, 0, v->f(k+3, j-1), 1);
								found = true;
							}
						}
						if ((mod[i+1] || mod[k]) && (!found) && can_pair(i+1, k)) 
						{
							if ((energy == 
									v->f(i+2, k-1)
									+ v->f(k+2, j)
									+ ptr_to_penalty(i+1, k, ct, data)
									+ ptr_to_penalty(k+2, j, ct, data)
									+ 2*multi_eparam(10, ct)
									+ 2*multi_eparam(6, ct)
									+ ptr_to_ergcoaxinterbases1(i+1, k, k+2, j, ct, data)
									+ ptr_to_erg1(i+1, k, i+2, k-1, ct, data))) 
							{
								registerbasepair(ct, i+1, k);
								stack->push(i+2, k-1, 0, v->f(i+2, k-1), 1);
								stack->push(k+2, j, 0, v->f(k+2, j), 1);
								found = true;
							}
						}
						if ((!found) && (mod[k+2]||mod[j]) && can_pair(k+2, j)) 
						{
							if ((energy == 
									v->f(i+1, k)
									+ v->f(k+3, j-1)
									+ ptr_to_penalty(i+1, k, ct, data)
									+ ptr_to_penalty(k+2, j, ct, data)
									+ 2*multi_eparam(10, ct)
									+ 2*multi_eparam(6, ct)
									+ ptr_to_ergcoaxinterbases1(i+1, k, k+2, j, ct, data)
									+ ptr_to_erg1(k+2, j, k+3, j-1, ct, data))) 
							{
								registerbasepair(ct, k+2, j);
								stack->push(i+1, k, 0, v->f(i+1,  k), 1);
								stack->push(k+3, j-1, 0, v->f(k+3,  j-1), 1);
								found = true;
							}
						}
					}
					k++;
				}
				if (ct->intermolecular) 
				{
					k = i + 1;
					while(k<j && !found) 
					{
						if (energy == w2->f(i, k) + w2->f(k+1, j)) 
						{
							stack->push(i, k, 0, w2->f(i, k), 0);
							stack->push(k+1, j, 0, w2->f(k+1, j), 0);
							found = true;
						}
						k++;
					}
					if (i!=number) 
					{
						if (!(fce->f(i, i) & INTER)) 
						{
							if (energy == wmb2->f(i+1, j) ) 
							{
								stack->push(i+1, j, 0, w2->f(i+1, j), 0); 
								found = true;
							}
						}
						else  if (energy == wmb2->f(i+1, j) + data->init - INFINITE_ENERGY) 
						{
							stack->push(i+1, j, 0, w2->f(i+1, j), 0);
							found = true;
						}
					}
					if ((j != number + 1) && (!found)) 
					{
						if (!(fce->f(j, j) & INTER)) 
						{
							if (energy == wmb2->f(i, j-1)) 
							{
								stack->push(i, j-1, 0, wmb2->f(i, j-1), 0);
								found = true;
							}
						}
						else if (energy == wmb2->f(i, j-1) + data->init-INFINITE_ENERGY) 
						{
							stack->push(i, j-1, 0, wmb2->f(i, j-1), 0);
							found = true;
						}
					}
				}
			}
			if (i==j) found = true;
			if (!found) 
			{
				errmsg (100, 2);
				tracebackerror = 1;
			}
		}
		else if (pair == 1) 
		{
			//energy = v->f(i, j) So we have a basepair
			while (energy == v->f(i, j)) 
			{
				//record the found base pair
				registerbasepair(ct, i, j);
			 	energy = energy + covar(i, j, ct);
				//check for stacked pair:
				if (energy == ptr_to_erg1(i, j, i + 1, j - 1, ct, data) + v->f(i + 1, j - 1) && i != number && j != number + 1) 
				{
					i++;
					j--;
					energy = v->f(i, j);
				}
				else break;
			}

			//now past the helical region:

			//check for hairpin loop:
			//if it is a hairpin loop, stop and return to the stack, otherwise
			//continue to define the loop closed by i and j
			if (energy != ptr_to_erg3(i, j, ct, data, fce->f(i, j))) 
			{
				found = false; //use the flag found to keep track as to whether a 
				//multiloop or exterior loop was found
				if (i == number) 
				{
					//place the 5' fragment on the stack if big enough to have structure:
					if (energy == ptr_to_penalty(i, j, ct, data) + w5[j-number-1]) 
					{
						if (j-number-1>minloop+1) 
						{
							stack->push(1, j-number-1, 1, w5[j-number-1], 0);
							found = true;
						}
						else found=true;
					}
					if (!found && energy == ptr_to_penalty(i, j, ct, data) + w5[j - number - 2] + 
						ptr_to_erg4(i, j, j - 1, 2, ct, data, lfce[j-1])) 
					{
						if (j - number - 2 > minloop + 1) 
						{
							stack->push(1, j-number-2, 1, w5[j-number-2], 0);
							found = true;
						}
						else found = true;
					}

					//now consider stacking to the 3' fragment:
					//first consider coaxial stacking to the 5' fragment:
					a = j - number - minloop - 2;
					
					while (a > 0 && !found) 
					{
						if (energy == w5[a-1] 
									+ ptr_to_penalty(i, j, ct, data) 
									+ ptr_to_penalty(a, j-number-1, ct, data)
									+ v->f(a, j-number-1)
									+ ptr_to_ergcoaxflushbases(a, j-number-1, j-number, i, ct, data)) 
						{
							if (a-1>minloop+1) stack->push(1, a-1, 1, w5[a-1], 0);
							if (i+1<number-minloop-1) stack->push(i+1, number, 1, w3[i+1], 0);
							stack->push(a, j-number-1, 0, v->f(a, j-number-1), 1);
							found = true;
						}
						else if (mod[a]||mod[j-number-1]&&can_pair(a, j-number-1)) 
						{

							if (energy == w5[a-1]
										+ ptr_to_penalty(i, j, ct, data)
										+ ptr_to_penalty(a, j-number-1, ct, data)
										+ v->f(a+1, j-number-2)
										+ ptr_to_ergcoaxflushbases(a, j-number-1, j-number, i, ct, data)
										+ ptr_to_erg1(a, j-number-1, a+1, j-number-2, ct, data)) 
							{
								registerbasepair(ct, a, j-number-1);
								if (a-1>minloop+1) stack->push(1, a-1, 1, w5[a-1], 0);
								if (i+1<number-minloop-1) stack->push(i+1, number, 1, w3[i+1], 0);
								stack->push(a+1, j-number-2, 0, v->f(a+1, j-number-2), 1);
								found = true;
							}
						}

						if (!found && (energy == w5[a-1]
												+ ptr_to_penalty(i, j, ct, data)
												+ ptr_to_penalty(a+1, j-number-2, ct, data)
												+ v->f(a+1, j-number-2)
												+ ptr_to_ergcoaxinterbases1(a+1, j-number-2, j-number, i, ct, data))) 
						{
							if (a-1>minloop+1) stack->push(1, a-1, 1, w5[a-1], 0);
							if (i+1<number-minloop-1) stack->push(i+1, number, 1, w3[i+1], 0);
							stack->push(a+1, j-number-2, 0, v->f(a+1, j-number-2), 1);
							found = true;
						}
						else if (!found && (mod[a+1]||mod[j-number-2]) && can_pair(i+1, j-number-2)) 
						{
							if (energy == w5[a-1]
										+ ptr_to_penalty(i, j, ct, data)
										+ ptr_to_penalty(a+1, j-number-2, ct, data)
										+ v->f(a+2, j-number-3)
										+ ptr_to_ergcoaxinterbases1(a+1, j-number-2, j-number, i, ct, data)
										+ ptr_to_erg1(a+1, j-number-2, a+2, j-number-3, ct, data))
							{
								registerbasepair(ct, a+1, j-number-2);
								if (a-1 > minloop+1) stack->push(1, a-1, 1, w5[a-1], 0);
								if (i+1<number-minloop-1) stack->push(i+1, number, 1, w3[i+1], 0);
								stack->push(a+2, j-number-3, 0, v->f(a+2, j-number-3), 1);
								found = true;
							}
						}
						if (!found && i<number) 
						{
							if (energy == w5[a - 1]									
										+ ptr_to_penalty(i, j, ct, data)
										+ ptr_to_penalty(a, j-number-2, ct, data)
										+ v->f(a, j-number-2)
										+ ptr_to_ergcoaxinterbases2(a, j-number-2, j-number, i, ct, data)) 
							{
								if (a-1>minloop+1) stack->push(1, a-1, 1, w5[a-1], 0);
								if (i+2<number-minloop-1) stack->push(i+2, number, 1, w3[i+2], 0);
								stack->push(a, j-number-2, 0, v->f(a, j-number-2), 1);
								found = true;
							}
							else if (mod[a] || mod[j-number-2] && can_pair(a, j-number-2)) 
							{
								if (energy == w5[a-1]
											+ ptr_to_penalty(i, j, ct, data)
											+ ptr_to_penalty(a, j-number-2, ct, data)
											+ v->f(a+1, j-number-3)
											+ ptr_to_ergcoaxinterbases2(a, j-number-2, j-number, i, ct, data)
											+ ptr_to_erg1(a, j-number-2, a+1, j-number-3, ct, data)) 
								{
									registerbasepair(ct, a, j-number-2);
									if (a-1>minloop+1) stack->push(1, a-1, 1, w5[a-1], 0);
									if (i+2<number-minloop-1) stack->push(i+2, number, 1, w3[i+2], 0);
									stack->push(a+1, j-number-3, 0, v->f(a+1, j-number-3), 1);
									found = true;
								}
							}
						}
						a--;
					}
				}
				else 
				{
					//this means that (j>i+2*minloop+2) 

					//check for multiloop closed by i and j

					//check for the case where the bifuraction is in a multibranch loop:

					if (energy == 
								wmb->f(i+1, j-1) 
								+ multi_eparam(10, ct) 
								+ multi_eparam(5, ct) 
								+ ptr_to_penalty(i, j, ct, data)) 
					{
						stack->push(i+1, j-1, 0, wmb->f(i+1, j-1), 2);
						found = true;
					}
					else if (energy == 
								wmb->f(i+2, j-1)
								+ ptr_to_penalty(i, j, ct, data)
								+ multi_eparam(10, ct)
								+ multi_eparam(5, ct)
								+ multi_eparam(6, ct)
								+ ptr_to_erg4(i, j, i+1, 1, ct, data, lfce[i+1])) 
					{
						stack->push(i+2, j-1, 0, wmb->f(i+2, j-1), 2);
						found = true;
					}
					else if (energy == 
								wmb->f(i+1, j-2)
								+ ptr_to_penalty(i, j, ct, data)
								+ multi_eparam(10, ct)
								+ multi_eparam(5, ct)
								+ multi_eparam(6, ct)
								+ ptr_to_erg4(i, j, j-1, 2, ct, data, lfce[j-1])) 
					{
						stack->push(i+1, j-2, 0, wmb->f(i+1, j-2), 2);
						found = true;
					}
					else if (energy == 
								wmb->f(i+2, j-2)
								+ ptr_to_penalty(i, j, ct, data)
								+ multi_eparam(10, ct)
								+ multi_eparam(5, ct)
								+ 2 * multi_eparam(6, ct)
							//	+ data->tstkm[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]]
								+ multi_tstkm(i, j, i + 1, j - 1, ct)
								+ checknp(lfce[i+1], lfce[j-1])
								+ ct->SHAPEss_give_value(i+1)
								+ ct->SHAPEss_give_value(j-1)) 
					{
						stack->push(i+2, j-2, 0, wmb->f(i+2, j-2), 2);
						found = true;
					}
					k = i + 1;
					while (k<j && !found) 
					{
						
						if (k+1<j-1 && 
								energy == ptr_to_penalty(i, j, ct, data) 
								+ ptr_to_penalty(i+1, k, ct, data)
								+ v->f(i+1, k)
								+ ptr_to_ergcoaxflushbases(j, i, i+1, k, ct, data)
								+ w->f(k+1, j-1)
								+ multi_eparam(5, ct)
								+ 2 * multi_eparam(10, ct)) 
						{
							stack->push(i+1, k, 0, v->f(i+1, k), 1);
							stack->push(k+1, j-1, 0, w->f(k+1, j-1), 0);
							found = true;
						}
						else if (k+1 < j-1&&mod[i+1] || mod[k] && can_pair(i+1, k)) 
						{
							if (energy == 
									ptr_to_penalty(i, j, ct, data)
									+ ptr_to_penalty(i+1, k, ct, data)+v->f(i+2, k-1)
									+ ptr_to_ergcoaxflushbases(j, i, i+1, k, ct, data)
									+ w->f(k+1, j-1)
									+ multi_eparam(5, ct)
									+ 2*multi_eparam(10, ct)
									+ ptr_to_erg1(i+1, k, i+2, k-1, ct, data)) 
							{
								registerbasepair(ct, i+1, k);
								stack->push(i+2, k-1, 0, v->f(i+2, k-1), 1);
								stack->push(k+1, j-1, 0, w->f(k+1, j-1), 0);
								found = true;
							}
						}

						if (!found && (k+1 < j-2 && energy == ptr_to_penalty(i, j, ct, data)+
									ptr_to_penalty(i+2, k, ct, data) 
									+ v->f(i+2, k) 
									+ ptr_to_ergcoaxinterbases1(j, i, i+2, k, ct, data)
									+ w->f(k+1, j-2)
									+ multi_eparam(5, ct)
									+ 2*multi_eparam(10, ct)
									+ 2*multi_eparam(6, ct))) 
						{
							stack->push(i+2, k, 0, v->f(i+2, k), 1);
							stack->push(k+1, j-2, 0, w->f(k+1, j-2), 0);
							found = true;
						}
						else if (!found && (k+1<j-2 && (mod[i+2]||mod[k])) && can_pair(i+2, k)) 
						{
							if (energy == 
									ptr_to_penalty(i, j, ct, data)
									+ ptr_to_penalty(i+2, k, ct, data)
									+ v->f(i+3, k-1)
									+ ptr_to_ergcoaxinterbases1(j, i, i+2, k, ct, data)
									+ w->f(k+1, j-2)
									+ multi_eparam(5, ct)
									+ 2*multi_eparam(10, ct)
									+ 2*multi_eparam(6, ct)
									+ ptr_to_erg1(i+2, k, i+3, k-1, ct, data)) 
							{
								registerbasepair(ct, i+2, k);
								stack->push(i+3, k-1, 0, v->f(i+3, k-1), 1);
								stack->push(k+1, j-2, 0, w->f(k+1, j-2), 0);
								found = true;
							}
						}
						if (!found && k+1 < j-1 && energy == 
								ptr_to_penalty(i, j, ct, data)
								+ ptr_to_penalty(i+2, k-1, ct, data)
								+ v->f(i+2, k-1)
								+ ptr_to_ergcoaxinterbases2(j, i, i+2, k-1, ct, data)
								+ w->f(k+1, j-1)
								+ multi_eparam(5, ct)
								+ 2*multi_eparam(10, ct)
								+ 2*multi_eparam(6, ct)) 
						{
							stack->push(i+2, k-1, 0, v->f(i+2, k-1), 1);
							stack->push(k+1, j-1, 0, w->f(k+1, j-1), 0);
							found = true;
						}
						else if (!found && (k+1<j-1) && (mod[i+2]||mod[k-1]) && can_pair(i+2, k-1)) 
						{
							if (energy == 
									ptr_to_penalty(i, j, ct, data)
									+ ptr_to_penalty(i+2, k-1, ct, data)
									+v->f(i+3, k-2)
									+ ptr_to_ergcoaxinterbases2(j, i, i+2, k-1, ct, data)
									+ w->f(k+1, j-1)
									+ multi_eparam(5, ct)
									+ 2*multi_eparam(10, ct)
									+ 2*multi_eparam(6, ct)
									+ ptr_to_erg1(i+2, k-1, i+3, k-2, ct, data)) 
							{
								registerbasepair(ct, i+2, k-1);
								stack->push(i+3, k-2, 0, v->f(i+3, k-2), 1);
								stack->push(k+1, j-1, 0, w->f(k+1, j-1), 0);
								found = true;
							}
						}

						if (!found && i+1<k-1 && energy == 
								ptr_to_penalty(i, j, ct, data)
								+ ptr_to_penalty(k, j-1, ct, data)
								+ v->f(k, j-1)
								+ ptr_to_ergcoaxflushbases(k, j-1, j, i, ct, data)
								+ w->f(i+1, k-1)
								+ multi_eparam(5, ct)
								+ 2*multi_eparam(10, ct)) 
						{
							stack->push(k, j-1, 0, v->f(k, j-1), 1);
							stack->push(i+1, k-1, 0, w->f(i+1, k-1), 0);
							found = true;
						}
						else if (!found && i+1<k-1 && (mod[k]||mod[j-1]) && can_pair(k,  j-1)) 
						{
							if (energy == 
									ptr_to_penalty(i, j, ct, data)
									+ ptr_to_penalty(k, j-1, ct, data)
									+ v->f(k+1, j-2)
									+ ptr_to_ergcoaxflushbases(k, j-1, j, i, ct, data)
									+ w->f(i+1, k-1)
									+ multi_eparam(5, ct)
									+ 2*multi_eparam(10, ct)
									+ ptr_to_erg1(k, j-1, k+1, j-2, ct, data)) 
							{
								registerbasepair(ct, k, j-1);
								stack->push(k+1, j-2, 0, v->f(k+1, j-2), 1);
								stack->push(i+1, k-1, 0, w->f(i+1, k-1), 0);
								found = true;
							}
						}

						if (!found && i+2<k-1 && energy == 
								ptr_to_penalty(i, j, ct, data)
								+ ptr_to_penalty(k, j-2, ct, data)
								+ v->f(k, j-2)
								+ ptr_to_ergcoaxinterbases2(k, j-2, j, i, ct, data)
								+ w->f(i+2, k-1)
								+ multi_eparam(5, ct)
								+ 2*multi_eparam(10, ct)
								+ 2*multi_eparam(6, ct))
						{
							stack->push(k, j-2, 0, v->f(k, j-2), 1);
							stack->push(i+2, k-1, 0, w->f(i+2, k-1), 0);
							found = true;
						}

						else if (!found && i+2<k-1 && (mod[k]||mod[j-2]) && can_pair(k, j-2)) 
						{
							if (energy == 
									ptr_to_penalty(i, j, ct, data)
									+ ptr_to_penalty(k, j-2, ct, data)
									+ v->f(k+1, j-3)
									+ ptr_to_ergcoaxinterbases2(k, j-2, j, i, ct, data)
									+ w->f(i+2, k-1)
									+ multi_eparam(5, ct)
									+ 2*multi_eparam(10, ct)
									+ 2*multi_eparam(6, ct)
									+ ptr_to_erg1(k, j-2, k+1, j-3, ct, data))
							{
								registerbasepair(ct, k, j-2);
								stack->push(k+1, j-3, 0, v->f(k+1, j-3), 1);
								stack->push(i+2, k-1, 0, w->f(i+2, k-1), 0);
								found = true;
							}
						}

						if (!found && i+1<k-1 && energy == 
								ptr_to_penalty(i, j, ct, data)
								+ ptr_to_penalty(k+1, j-2, ct, data)
								+ v->f(k+1, j-2)
								+ ptr_to_ergcoaxinterbases1(k+1, j-2, j, i, ct, data)
								+ w->f(i+1, k-1)
								+ multi_eparam(5, ct)
								+ 2*multi_eparam(10, ct)
								+ 2*multi_eparam(6, ct)) 
						{
							stack->push(k+1, j-2, 0, v->f(k+1, j-2), 1);
							stack->push(i+1, k-1, 0, w->f(i+1, k-1), 0);
							found = true;
						}
						else if (!found && i+1<k-1 && (mod[k+1]||mod[j-2]) && can_pair(k+1, j-2)) 
						{
							if (energy == 
									ptr_to_penalty(i, j, ct, data)
									+ ptr_to_penalty(k+1, j-2, ct, data)
									+ v->f(k+2, j-3)
									+ ptr_to_ergcoaxinterbases1(k+1, j-2, j, i, ct, data)
									+ w->f(i+1, k-1)
									+ multi_eparam(5, ct)
									+ 2*multi_eparam(10, ct)
									+ 2*multi_eparam(6, ct)
									+ ptr_to_erg1(k+1, j-2, k+2, j-3, ct, data)) 
							{
								registerbasepair(ct, k+1, j-2);
								stack->push(k+2, j-3, 0, v->f(k+2, j-3), 1);
								stack->push(i+1, k-1, 0, w->f(i+1, k-1), 0);
								found = true;
							}
						}

						//}
						if (k == number && !found) 
						{
							if (energy == w3[i+1] + w5[j-number-1] + ptr_to_penalty(i, j, ct, data)) 
							{
								if (i<number-minloop-2) stack->push(i+1, number, 1, w3[i+1], 0);
								if (j-number>minloop+3) stack->push(1, j-number-1, 1, w5[j-number-1], 0);
								found = true;
							}
							else if (energy == 
										w3[i+2] 
										+ w5[j-number-1] 
										+ ptr_to_penalty(i, j, ct, data) 
										+ ptr_to_erg4(i, j, i+1, 1, ct, data, lfce[i+1])) 
							{
								if (i<number-minloop-3) stack->push(i+2, number, 1, w3[i+2], 0);
								if (j-number>minloop+3) stack->push(1, j-number-1, 1, w5[j-number-1], 0);
								found = true;
							}
							else if (j-number-2>-1) 
							{
								if (energy == 
										w3[i+1]
										+ w5[j-number-2]
										+ ptr_to_penalty(i, j, ct, data) 
										+ ptr_to_erg4(i, j, j-1, 2, ct, data, lfce[j-1])) 
								{
									if (i<number-minloop-2) stack->push(i+1, number, 1, w3[i+1], 0);
									if (j-number>minloop+4) stack->push(1, j-number-2, 1, w5[j-number-2], 0);
									found = true;
								}
								else if (energy == 
										w3[i+2]
										+ w5[j-number-2]
										+ ptr_to_penalty(i, j, ct, data) 
									//	+ data->tstack[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]]
										+ multi_tstack(i, j, i + 1, j - 1, ct)
										+ checknp(lfce[i+1], lfce[j-1])
										+ ct->SHAPEss_give_value(i+1)
										+ ct->SHAPEss_give_value(j-1)) 
								{
									if (i<number-minloop-3) stack->push(i+2, number, 1, w3[i+2], 0);
									if (j-number>minloop+4) stack->push(1, j-number-2, 1, w5[j-number-2], 0);
									found = true;
								}
							}

							//also consider coaxial stacking:

							//first consider coaxial stacking to the 5' fragment:
							a = j-number-minloop-2;
							while (a>0&&!found) 
							{
								if (energy == 
										w5[a-1]
										+ w3[i+1]
										+ ptr_to_penalty(i, j, ct, data)
										+ ptr_to_penalty(a, j-number-1, ct, data)
										+ v->f(a, j-number-1)
										+ ptr_to_ergcoaxflushbases(a, j-number-1, j-number, i, ct, data)) 
								{
									if (a-1>minloop+1) stack->push(1, a-1, 1, w5[a-1], 0);
									if (i+1<number-minloop-1) stack->push(i+1, number, 1, w3[i+1], 0);
									stack->push(a, j-number-1, 0, v->f(a, j-number-1), 1);
									found = true;
								}
								else if (mod[a] || mod[j-number-1] && can_pair(a,  j-number-1)) 
								{
									if (energy == 
											w5[a-1]
											+ w3[i+1]
											+ ptr_to_penalty(i, j, ct, data)
											+ ptr_to_penalty(a, j-number-1, ct, data)
											+ v->f(a+1, j-number-2)
											+ ptr_to_ergcoaxflushbases(a, j-number-1, j-number, i, ct, data)
											+ ptr_to_erg1(a, j-number-1, a+1, j-number-2, ct, data)) 
									{
										registerbasepair(ct, a, j-number-1);
										if (a-1>minloop+1) stack->push(1, a-1, 1, w5[a-1], 0);
										if (i+1<number-minloop-1) stack->push(i+1, number, 1, w3[i+1], 0);
										stack->push(a+1, j-number-2, 0, v->f(a+1, j-number-2), 1);
										found = true;
									}
								}

								if (!found && energy == 
										w5[a-1]
										+ w3[i+1]
										+ ptr_to_penalty(i, j, ct, data)
										+ ptr_to_penalty(a+1, j-number-2, ct, data)
										+ v->f(a+1, j-number-2)
										+ ptr_to_ergcoaxinterbases1(a+1, j-number-2, j-number, i, ct, data)) 
								{
									if (a-1>minloop+1) stack->push(1, a-1, 1, w5[a-1], 0);
									if (i+1<number-minloop-1) stack->push(i+1, number, 1, w3[i+1], 0);
									stack->push(a+1, j-number-2, 0, v->f(a+1, j-number-2), 1);
									found = true;
								}
								else if (!found && (mod[a+1]||mod[j-number-2]) && can_pair(a+1, j-number-2)) 
								{
									if (energy == 
											w5[a-1]
											+ w3[i+1]
											+ ptr_to_penalty(i, j, ct, data)
											+ ptr_to_penalty(a+1, j-number-2, ct, data)
											+ v->f(a+2, j-number-3)
											+ ptr_to_ergcoaxinterbases1(a+1, j-number-2, j-number, i, ct, data)
											+ ptr_to_erg1(a+1, j-number-2, a+2, j-number-3, ct, data)) 
									{
										registerbasepair(ct, a+1, j-number-2);
										if (a-1>minloop+1) stack->push(1, a-1, 1, w5[a-1], 0);
										if (i+1<number-minloop-1) stack->push(i+1, number, 1, w3[i+1], 0);
										stack->push(a+2, j-number-3, 0, v->f(a+2, j-number-3), 1);
										found = true;
									}
								}

								if (!found && energy == 
												w5[a-1]
												+ w3[i+2]
												+ ptr_to_penalty(i, j, ct, data)
												+ ptr_to_penalty(a, j-number-2, ct, data)
												+ v->f(a, j-number-2)
												+ ptr_to_ergcoaxinterbases2(a, j-number-2, j-number, i, ct, data)) 
								{
									if (a-1>minloop+1) stack->push(1, a-1, 1, w5[a-1], 0);
									if (i+2<number-minloop-1) stack->push(i+2, number, 1, w3[i+2], 0);
									stack->push(a, j-number-2, 0, v->f(a, j-number-2), 1);
									found = true;
								}
								else if (!found && (mod[a]||mod[j-number-2]) && can_pair(a, j-number-2)) 
								{
									if (energy == 
											w5[a-1]
											+ w3[i+2]
											+ ptr_to_penalty(i, j, ct, data)
											+ ptr_to_penalty(a, j-number-2, ct, data)
											+ v->f(a+1, j-number-3)
											+ ptr_to_ergcoaxinterbases2(a, j-number-2, j-number, i, ct, data)
											+ ptr_to_erg1(a, j-number-2, a+1, j-number-3, ct, data)) 
									{
										registerbasepair(ct, a, j-number-2);
										if (a-1>minloop+1) stack->push(1, a-1, 1, w5[a-1], 0);
										if (i+2<number-minloop-1) stack->push(i+2, number, 1, w3[i+2], 0);
										stack->push(a+1, j-number-3, 0, v->f(a+1, j-number-3), 1);
										found = true;
									}
								}
								a--;
							}

							//now consider stacking to the 3' fragment:
							a = i + minloop + 2;
							while (a <= number && !found) 
							{
								if (energy == 
										w5[j-number-1]
										+ w3[a+1]
										+ ptr_to_penalty(i, j, ct, data)
										+ ptr_to_penalty(i+1, a, ct, data) 
										+ v->f(i+1, a)
										+ ptr_to_ergcoaxflushbases(j-number, i, i+1, a, ct, data)) 
								{

									if (j-number-1>minloop+1) stack->push(1, j-number-1, 1, w5[j-number-1], 0);
									if (a+1<number-minloop-1) stack->push(a+1, number, 1, w3[a+1], 0);

									stack->push(i+1, a, 0, v->f(i+1, a), 1);
									found = true;
								}
								else if (mod[i+1] || mod[a] && can_pair(i+1, a)) 
								{
									if (energy == 
											w5[j-number-1]
											+ w3[a+1]
											+ ptr_to_penalty(i, j, ct, data)
											+ ptr_to_penalty(i+1, a, ct, data) 
											+ v->f(i+2, a-1)
											+ ptr_to_ergcoaxflushbases(j-number, i, i+1, a, ct, data)
											+ ptr_to_erg1(i+1, a, i+2, a-1, ct, data)) 
									{
										registerbasepair(ct, i+1, a);
										if (j-number-1>minloop+1) stack->push(1, j-number-1, 1, w5[j-number-1], 0);
										if (a+1<number-minloop-1) stack->push(a+1, number, 1, w3[a+1], 0);
										stack->push(i+2, a-1, 0, v->f(i+2, a-1), 1);
										found = true;
									}
								}

								if (j-number-2>-1) 
								{
									if (!found && energy == 
													w5[j-number-2]
													+ w3[a+1]
													+ ptr_to_penalty(i, j, ct, data)
													+ ptr_to_penalty(i+2, a, ct, data)
													+ v->f(i+2, a)
													+ ptr_to_ergcoaxinterbases1(j-number, i, i+2, a, ct, data)) 
									{
										if (j-number-2>minloop+1) stack->push(1, j-number-2, 1, w5[j-number-2], 0);
										if (a+1<number-minloop-1) stack->push(a+1, number, 1, w3[a+1], 0);
										stack->push(i+2, a, 0, v->f(i+2, a), 1);
										found = true;
									}

									else if (!found && (mod[i+2]||mod[a]) && can_pair(i+2, a)) 
									{
										if (!found && energy == 
														w5[j-number-2]
														+ w3[a+1]
														+ ptr_to_penalty(i, j, ct, data)
														+ ptr_to_penalty(i+2, a, ct, data)
														+ v->f(i+3, a-1)
														+ ptr_to_ergcoaxinterbases1(j-number, i, i+2, a, ct, data)
														+ ptr_to_erg1(i+2, a, i+3, a-1, ct, data)) 
										{
											registerbasepair(ct, i+2, a);
											if (j-number-2>minloop+1) stack->push(1, j-number-2, 1, w5[j-number-2], 0);
											if (a+1<number-minloop-1) stack->push(a+1, number, 1, w3[a+1], 0);
											stack->push(i+3, a-1, 0, v->f(i+3, a-1), 1);
											found = true;
										}
									}
								}

								if (!found && energy == 
												w5[j-number-1]
												+ w3[a+1]+ptr_to_penalty(i, j, ct, data)
												+ ptr_to_penalty(i+2, a-1, ct, data)
												+ v->f(i+2, a-1)
												+ ptr_to_ergcoaxinterbases2(j-number, i, i+2, a-1, ct, data)) 
								{
									if (j-number-1>minloop+1) stack->push(1, j-number-1, 1, w5[j-number-1], 0);
									if (a+1<number-minloop-1) stack->push(a+1, number, 1, w3[a+1], 0);

									stack->push(i+2, a-1, 0, v->f(i+2, a-1), 1);
									found = true;
								}

								else if (!found && (mod[i+2]||mod[a-1]) && can_pair(i+2,  a-1)) 
								{
									if (!found && energy == 
													w5[j-number-1]
													+ w3[a+1]
													+ ptr_to_penalty(i, j, ct, data)
													+ ptr_to_penalty(i+2, a-1, ct, data)
													+ v->f(i+3, a-2)
													+ ptr_to_ergcoaxinterbases2(j-number, i, i+2, a-1, ct, data)
													+ ptr_to_erg1(i+2, a-1, i+3, a-2, ct, data)) 
									{
										registerbasepair(ct, i+2, a-1);
										if (j-number-1>minloop+1) stack->push(1, j-number-1, 1, w5[j-number-1], 0);
										if (a+1<number-minloop-1) stack->push(a+1, number, 1, w3[a+1], 0);
										stack->push(i+3, a-2, 0, v->f(i+3, a-2), 1);
										found = true;
									}
								}
								a++;
							}
						}
						if (found) break;
						k++;
						}
					}
				if (ct->intermolecular) 
				{
						if (energy == 
									wmb2->f(i+1, j-1)
									+ ptr_to_penalty(i, j, ct, data)
									+ INFINITE_ENERGY) 
						{
							stack->push(i+1, j-1, 0, wmb2->f(i+1, j-1), 2);
							found = true;
						}
						else if (energy == 
									wmb2->f(i+2, j-1)
									+ ptr_to_penalty(i, j, ct, data)
									+ ptr_to_erg4(i, j, i+1, 1, ct, data, lfce[i+1])
									+ INFINITE_ENERGY) 
						{
							stack->push(i+2, j-1, 0, wmb2->f(i+2, j-1), 2);
							found = true;
						}
						else if (energy == 
									wmb2->f(i+1, j-2)
									+ ptr_to_penalty(i, j, ct, data)
									+ ptr_to_erg4(i, j, j-1, 2, ct, data, lfce[j-1])
									+ INFINITE_ENERGY) 
						{
							stack->push(i+1, j-2, 0, wmb2->f(i+1, j-2), 2);
							found = true;
						}
						else if (energy == 
									wmb2->f(i+2, j-2)
									+ ptr_to_penalty(i, j, ct, data)
									+ INFINITE_ENERGY
								//	+ data->tstkm[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]]
									+ multi_tstkm(i, j, i + 1, j - 1, ct)
									+ checknp(lfce[i+1], lfce[j-1])) 
						{
							stack->push(i+2, j-2, 0, wmb2->f(i+2, j-2), 2);
							found = true;
						}
				}
				if (!found) 
				{
					//neither a multiloop nor an exterior loop were found, 
					//there must be a bulge/interior loop
					for (d = j-i-minloop;(d>=1)&&!found;d--) 
					{
						for (a = i+1;(a<=j-1-d) && !found; a++) 
						{
							b = d + a;
							if (energy == (ptr_to_erg2(i, j, a, b, ct, data, fce->f(i, a), fce->f(b, j))+ v->f(a,  b)))											
							{
								i = a;
								j = b;
								stack->push(i, j, 0, v->f(i, j), 1);
								found = true;
							}
							if (mod[a] || mod[b] && can_pair(a, b)) 
							{
								if (energy == 
											(ptr_to_erg2(i, j, a, b, ct, data, fce->f(i, a), fce->f(b, j))
											+ v->f(a+1, b-1)
											+ ptr_to_erg1(a, b, a+1, b-1, ct, data))) 
								{
									i = a + 1;
									j = b - 1;
									registerbasepair(ct, a, b);
									stack->push(i, j, 0, v->f(i, j), 1);
									found = true;
								}
							}
						}
					}
				}
				if (!found) 
				{	//something went wrong					
					errmsg (100, 2);
					tracebackerror=1;
				}
			}
		}
	}

	delete stack;
	return tracebackerror;
	#undef can_pair
}

// Starting point for the maximum number of basepairs within %cntrl8 of the minimum free energy.
#define maxsort 90000 

//  Traceback predicts a set of low free energy structures using the mfold heuristic with the fill step information.
//  cntrl8 is the maximum % difference in free energy of suboptimal structures if > 0
//	otherwise,  cntrl8 is a maximum energy difference in kcal/mol*factor
//  This returns an error code,  where zero is no error and non-zero indicates a traceback error.
int traceback(	structure *ct,  
				datatable *data, 
				DynProgArray<integersize> *v, 
				DynProgArray<integersize> *w, 
				DynProgArray<integersize> *wmb, 
				DynProgArray<integersize> *w2,
				DynProgArray<integersize> *wmb2, 
				integersize *w3, integersize *w5, 
				forceclass *fce,
				bool *lfce,
				integersize vmin, 
				int cntrl6, 
				int cntrl8, 
				int cntrl9, 
				bool *mod	) 
{	
	
	integersize(*ptr_to_erg1)(int i, int j, int ip, int jp, structure * ct, datatable * data);
	integersize(*ptr_to_erg2)(int i, int j, int ip, int jp, structure * ct, datatable * data, char a, char b);
	integersize(*ptr_to_erg3)(int i, int j, structure * ct, datatable * data, char dbl);
	integersize(*ptr_to_erg4)(int i, int j, int ip, int jp, structure * ct, datatable * data, bool lfce);
	integersize(*ptr_to_ergcoaxflushbases) (int i, int j, int ip, int jp, structure * ct, datatable * data);
	integersize(*ptr_to_ergcoaxinterbases1) (int i, int j, int ip, int jp, structure * ct, datatable * data);
	integersize(*ptr_to_ergcoaxinterbases2) (int i, int j, int ip, int jp, structure * ct, datatable * data);
	integersize(*ptr_to_penalty) (int i, int j, structure * ct, datatable * data);

	if (ct->GetNumberofSequences() == 1) {
		ptr_to_erg1 = erg1;
		ptr_to_erg2 = erg2;
		ptr_to_erg3 = erg3;
		ptr_to_erg4 = erg4;
		ptr_to_ergcoaxflushbases = ergcoaxflushbases;
		ptr_to_ergcoaxinterbases1 = ergcoaxinterbases1;
		ptr_to_ergcoaxinterbases2 = ergcoaxinterbases2;
		ptr_to_penalty = penalty;
	}
	else if (ct->GetNumberofSequences() > 1) {
		ptr_to_erg1 = multi_erg1;
		ptr_to_erg2 = multi_erg2;
		ptr_to_erg3 = multi_erg3;
		ptr_to_erg4 = multi_erg4;
		ptr_to_ergcoaxflushbases = multi_ergcoaxflushbases;
		ptr_to_ergcoaxinterbases1 = multi_ergcoaxinterbases1;
		ptr_to_ergcoaxinterbases2 = multi_ergcoaxinterbases2;
		ptr_to_penalty = multi_penalty;
	}

	int ii;
	int iret, jret;
	int i, j, num, count, count2; //local
	int numbp, k1, k2, crit;
	int cur, c, k, cntr;
	
	int ji; //(added during debugging)

	int tracebackerror = 0;
	int sort = maxsort;
	int number = ct->GetSequenceLength();
	bool flag = true;
	
	#if defined(debugmode)
	char filename[maxfil];
	char temp[20];
	#endif

	// This code previously prevented tracebacks for structures with positive free energy change.
	// Now, just trap and forbid tracebacks for very high folding free energy changes.  
	// These, for example, 	are structures that have no viable pairs.
	// Return the mfe as the unpaired structure, however.

	if (vmin >= 0) 
	{ 
		//Add an empty structure to ct:
		ct->AddStructure();
		//Set the energy to zero, a random coil:
		ct->SetEnergy(1, 0);
		if (vmin >= INFINITE_ENERGY/2) return 0;//no viable structure found
	}

	// mark keeps track of which pairs have been formed by the suboptimal routine
	// declare and initialize mark
	bool **mark = new bool *[number + 1];
	for (i = 0; i <= number; i++) 
	{
		mark[i] = new bool[number + 1];
	}
	for (count = 0; count <= number; count++) 
	{
		for (count2 = 0; count2 <= number; count2++) 
		{
			mark[count][count2] = false;
		}
	}

	// if cntr18 is 20 and vmin 100, crit will be 120
	if (cntrl8 > 0) crit = (short) ( abs(vmin) * ( float(cntrl8) / 100.0 ) );
	else crit = -cntrl8;
	crit = crit + vmin;

	//  Construct heapi and heapj composed of pairs ij such that a structure
	//	containing ij has energy less than a given percent (ie:cntrl8) of 
	//  the minimum folding energy
	integersize *energy = new integersize[sort + 1];
	int *heapi = new int [sort + 1];
	int *heapj = new int [sort + 1];

	num = 0;
	i = 1;
	j = 2;

	while (i < number) 
	{
		if (num == sort) 
		{
			//allocate more space for the heap
			delete[] heapi;
			delete[] heapj;
			delete[] energy;
			sort = 10 * sort;
			heapi = new int [sort + 1];
			heapj = new int [sort + 1];
			energy = new integersize [sort + 1];
			i = 1;
			j = 2;
			num = 0;
		}
		if ((v->f(i, j) + v->f(j, i + number) + covar(i, j, ct)) <= crit) 
		{
			num++;
			heapi[num] = i;
			heapj[num] = j;
			energy[num] = v->f(i, j) + v->f(j, i + number) + covar(i, j ,ct);
			j = j + 1;
			if (j > number) 
			{
				i++;
				j = i + 1;
			}
		}
		else if (mod[i] || mod[j]) 
		{
			//add i-j pair to heap if it is adjacent to unpaired nucs in one direction
			if (v->f(i, j) < INFINITE_ENERGY) 
			{
				cur = v->f(i, j) + v->f(j + 1, i + number - 1) + ptr_to_erg1(j, i + number, j + 1, i + number - 1, ct, data);
				if (cur <= crit) 
				{
					num++;
					heapi[num] = i;
					heapj[num] = j;
					energy[num] = crit;
					j++;	
					if (j > number) 
					{
						i++;
						j = i + 1;
					}
				}
				else 
				{
					j++;
					if (j > number) 
					{
						i++;
						j = i + 1;
					}
				}
			}
			else if (v->f(j, i + number) == INFINITE_ENERGY) 
			{
				cur = v->f(i + 1, j - 1) + v->f(j, i + number) + ptr_to_erg1(i, j, i + 1, j - 1, ct, data);
				if (cur <= crit) 
				{
					num++;
					heapi[num] = i;
					heapj[num] = j;
					energy[num] = crit;
					j = j++;
					if (j > number) 
					{
						i++;
						j = i + 1;
					}
				}		
				else 
				{
					j++;
					if (j > number) 
					{
						i++;
						j = i + 1;
					}
				}
			}	
			else 
			{
				j++;
				if (j > number) 
				{
					i++;
					j = i + 1;
				}
			}
		}
		else 
		{
			j++;
			if (j > number) 
			{
				i++;
				j = i + 1;
			}
		}
	}

	//sort the base pair list:

	//energy[0] is never used, but gets tested in the case that up=0;
	energy[0] = 0; 

	//make a heap:
	int q, up;
	for (q = 2; q <= num; q++) 
	{
		cur = q;
		up = cur/2;
		while ( (energy[cur] < energy[up]) && up >= 1) 
		{
			swap(heapi[cur], heapi[up]);
			swap(heapj[cur], heapj[up]);
			swap(energy[cur], energy[up]);
			cur = cur/2;
			up = up/2;
		}
	}

	//sort the heap:
	for (int ir = num - 1; ir >= 1; ir--) 
	{
		swap(heapi[ir + 1], heapi[1]);
		swap(heapj[ir + 1], heapj[1]);
		swap(energy[ir + 1], energy[1]);

		up = 1;
		c = 2;
		while (c <= ir) 
		{
			if (c != ir) 
			{
				if (energy[c + 1] < energy[c]) 
					c++;
			}
			if (energy[c] < energy[up]) 
			{
				swap(heapi[c], heapi[up]);
				swap(heapj[c], heapj[up]);
				swap(energy[c], energy[up]);
				up = c;
				c = 2 * c;
			}
			else c = ir + 1;
		}
	}

	cntr = num;

	// keep track if a structure was allocated, but has been rejected with failedprevious.
	bool failedprevious = false;
	while (flag) 
	{
		// This is routine to select the region of the structure to be
		// folded, it allows for sub-optimal structure predictions
		// err=0;
		// Select the next valid unmarked basepair
		while (mark[heapi[cntr]][heapj[cntr]]) 
		{
			if (cntr == 1) 
			{
				flag = false;
				goto sub900;
			}
			cntr--;
		}

		iret = heapi[cntr];
		jret = heapj[cntr];

		// No prior structure was allocated, so add a structure now:
		if (!failedprevious) ct->AddStructure();
		//The prior failed structure needs to be cleaned, i.e. all base pairs removed:
		else ct->CleanStructure(ct->GetNumberofStructures());
		
		// Reset failedprevious as a new traceback is started:
		failedprevious = false;

		// Traceback to find best structure on included fragment (ie: iret to jret)
		if (flag) 
		{
			for (count = 1;count <= 2; count++) 
			{
				if (count == 1) 
				{
					ii = iret;
					ji = jret;
				}
				if (count == 2) 
				{
					ii = jret;
					ji = iret + number;
				}
				if (mod[iret] || mod[jret]) 
				{
					if (v->f(ii, ji) == INFINITE_ENERGY) 
					{
						ii += 1;
						ji -= 1;
					}
				}

				#ifndef DYNALIGN_II
				if (trace(ct, data, ii, ji, v, w, wmb, w2, wmb2, lfce, fce, w3, w5, mod) != 0) tracebackerror = 1;
				#else
				if(trace(ct, data, ii, ji, v, w, wmb, w2, wmb2, lfce, fce, w3, w5, mod, NULL)!=0) tracebackerror=1;
				#endif

				if (count == 2) 
				{
					ct->SetEnergy(ct->GetNumberofStructures(), energy[cntr]);
					
					//count the number of new base pairs not within window of existing base pairs
					numbp = 0;
					for (k = 1; k <= number; k++) 
					{
						if (k<(ct->GetPair(k, ct->GetNumberofStructures()))) 
						{
							if (!(mark[k][ct->GetPair(k, ct->GetNumberofStructures())])) numbp++;
						}
					}
					for (k = 1; k <= number; k++) 
					{
						if (k < ct->GetPair(k, ct->GetNumberofStructures())) 
						{
							// Mark "traced back" base pairs and also base pairs
							//	which are within a window of cntrl9
							mark[k][ct->GetPair(k, ct->GetNumberofStructures())] = true;
							if (cntrl9 > 0) 
							{
								for (k1 = max(1, k - cntrl9); k1 <= min(number, k + cntrl9); k1++) 
								{
									for (k2 = max(k1, ct->GetPair(k, ct->GetNumberofStructures()) - cntrl9); k2 <= min(number, ct->GetPair(k, ct->GetNumberofStructures()) + cntrl9); k2++) 
									{    	
										mark[k1][k2] = true;
									}
								}
							}
						}
					}
					if (numbp <= cntrl9 && ct->GetNumberofStructures() > 1) 
					{
						failedprevious = true;
						//goto sub900;
					}					
					else 
					{
						//place the structure name (from ctlabel[1]) into each structure
						ct->SetCtLabel(ct->GetSequenceLabel(), ct->GetNumberofStructures());
				
						#if defined(debugmode)
					
						strcpy(filename,"energydump");
						itoa(ct->GetNumberofStructures(),temp,10);
						strcat(filename,temp);
						strcat(filename,".out");
						energydump (ct, v, data, ct->GetNumberofStructures(),filename,iret,jret);

						#endif
					}
					if (ct->GetNumberofStructures() == cntrl6 && !failedprevious) 
					{
						flag = false;
					}
				}
			}
			sub900:
				continue;
		}
	}

	if (failedprevious) ct->RemoveLastStructure();

	de_allocate(mark, number + 1);
	delete[] energy;
	delete[] heapi;
	delete[] heapj;
	return tracebackerror;
} // traceback
#else
#endif

#ifndef DYNALIGN_II
	//This is the dynamic algorithm of Zuker:
	//cntrl6 = #tracebacks
	//cntrl8 = percent sort
	//cntrl9 = window

	//ProgressHandler is an interface for returning the progress of the calculation (TProgressDialog is a sub-class of ProgressHandler suitable for text interfaces)
	//Savfile is for creating a file with arrays and parameters for refolding with different 
	//suboptimal tracebacks
	//quickenergy indicates whether to find the lowest free energy for the sequence without a structure
#ifndef INSTRUMENTED
	int dynamic(structure* ct,datatable* data,int cntrl6, int cntrl8,int cntrl9,
			ProgressHandler* update, bool quickenergy, char* save, int maxinter, bool quickstructure, bool simple_iloops, bool disablecoax, bool allow_isolated)


#else //INSTRUMENTED IS DEFINED
		void dynamic(structure* ct,datatable* data,int cntrl6, int cntrl8,int cntrl9,
				DynProgArray<integersize> *v, DynProgArray<integersize> *vmb/*tracks MB loops*/, DynProgArray<integersize> *vext/*tracks exterior loops*/,
				ProgressHandler* update, bool quickenergy, char* save, int maxinter, bool quickstructure, bool simple_iloops, bool disablecoax, bool allow_isolated)
#endif //END of INSTRUMENTED iS DEFINED
		{		
			int number;		
			int vmin;
			//int i,j,k,l;
			int i, j;
			//			int m,n,o,p;

			bool *lfce,*mod;//[maxbases+1][maxbases+1];

			//integersize **work,**work2;
			integersize *w5,*w3;//**wca,**curE,**prevE,**tempE;

			if (ct->GetThermodynamicDataTable()!=data) {
				cerr << "In dynamic ("<<__FILE__<<"): The structure's datatable ("<<ct->GetThermodynamicDataTable()<<") does not match the passed-in datatable ("<<data<<"). This can cause problems with IsNuc etc." << endl;
				ct->SetThermodynamicDataTable(data);
			}

	// #define can_pair(i, j) (multi_can_pair(i, j, ct))
	#define can_pair(i, j) (ct->can_pair(i, j, ct))
	
	int tracebackerror = 0;

	#ifdef SMP
		if(!simple_iloops) return 27;		
	#endif

	// inc is an array that saves time by showing which bases can pair before erg is called
		number = (ct->GetSequenceLength());//place the number of bases in a registered integer
	// time the algorithm execution
	#ifdef timer 
		#include <time.h>
		ofstream timeout;
		int seconds;
		char timerstring[100];
		char timelength[10];
		strcpy(timerstring,"time_pred_");
		sprintf(timelength,"%i",ct->GetSequenceLength());
		strcat(timerstring,timelength);
		strcat(timerstring,".out");
		timeout.open(timerstring);
		timeout<<time(NULL)<<"\n";
		seconds = time(NULL);
	#endif

	//  declare v, w and wmb arrays:
	//	v[i][j] is the best energy for subsequence i to j when i and j are paired
	//	w[i][j] is the best energy for subsequence i to j
	//	for i<j<n, it is the interior framgment between nucleotides i and j
	//	for i<n<j, it is the exterior structure fragmnet that is 5' to j-n and 3' to i
	DynProgArray<integersize> w(number);
	DynProgArray<integersize> wmb(number);
	#ifndef INSTRUMENTED // if INSTRUMENTED compiler flag is defined
	DynProgArray<integersize> v(number);
	#endif
	
	forceclass fce(number);

	// declare and initialize w2 and wmb2 arrays for internmolecular folding:
	DynProgArray<integersize> *w2, *wmb2;
	if (ct->intermolecular) 
	{
		w2 = new DynProgArray<integersize>(number);
		wmb2 = new DynProgArray<integersize>(number);
	}
	else 
	{
		wmb2 = NULL;
		w2 = NULL;
	}

	// declare and initialize mod and lfce
	lfce = new bool [2 * number + 1];
	mod = new bool [2 * number + 1];
	for (i = 0; i <= 2 * number; i++) 
	{
		lfce[i] = false;
		mod[i] = false;
	}
	for (i = 0; i < ct->GetNumberofModified(); i++) 
	{
		if (ct->GetModified(i) > 1 && ct->GetModified(i) < ct->GetSequenceLength()) 
		{
			mod[ct->GetModified(i)] = true;
			mod[ct->GetModified(i) + ct->GetSequenceLength()] = true;
		}
	}
	
	// declare and initialize w5 and w3
			w5 = new integersize [number+1];
			w3 = new integersize [number+2];

			for (i=0;i<=number;i++) {
				w5[i] = 0;
				w3[i] = 0;
			}

			w3[number+1] = 0;
			force(ct,&fce,lfce);
			vmin=INFINITE_ENERGY;

#ifndef INSTRUMENTED//If pre-compiler flag INSTRUMENTED is not defined, compile the following code
#ifndef DYNALIGN_II
			//perform the fill steps:(i.e. fill arrays v and w.)
#ifndef _CUDA_CALC_DYNALIGN_
			fill(ct, v, w, wmb, fce, vmin,lfce, mod,w5, w3, quickenergy, data, w2, wmb2, update, maxinter, quickstructure, simple_iloops, disablecoax, allow_isolated);
#endif
#else//ifdef DYNALIGN_II
#ifndef _CUDA_CALC_DYNALIGN_
                        fill(ct, v, w, wmb, fce, vmin,lfce, mod,w5, w3, quickenergy, data, w2, wmb2, NULL, update, maxinter);
#endif
#endif

#ifdef _CUDA_CALC_DYNALIGN_
        //read data
        struct fparam par;
        const char *path = getDataPath();
        if (!path) die("%s: need to set environment variable $DATAPATH", "no nearest neighbor parameters found");
        fparam_read_from_text(path, 0, &par);
        //copy the sequence, correcting for 1-indexing
        const int len = ct->GetSequenceLength();
        char* buf = new char[len+1];
        for(int i=0;i<len;i++)
		{
                buf[i] = ct->nucs[i+1];
        }
        buf[len] = '\0'; //make it null-terminated

        // do fill routine
        frna_t gpu = frna_new(buf, &par);

        //copy arrays
        vmin = INFINITE_ENERGY;
        for(int i=1;i<=len;i++)
		{
			const int ip = i-1;
			for(int j=i+1;j<=len;j++)
			{
                        const int jp = j-1;
                        v.f(i,j) = (integersize) (gpu->v[ip*len+jp]);
                        w.f(i,j) = (integersize) (gpu->w[ip*len+jp]);
                        wmb.f(i,j) = (integersize) (gpu->wm[ip*len+jp]);
                        v.f(j,i+len) = (integersize) (gpu->v[jp*len+ip]);
                        w.f(j,i+len) = (integersize) (gpu->w[jp*len+ip]);
                        wmb.f(j,i+len) = (integersize) (gpu->wm[jp*len+ip]);
						const int temp = gpu->v[ip*len+jp]+gpu->v[jp*len+ip];
						vmin = min(temp,vmin);
			}
			w5[i] = gpu->w5[ip];
			w3[i] = gpu->w3[ip];
		}

        for(int i=1;i<=len;i++)
		{
			for(int j=i+1;j<=len;j++)
			{
				int tmp = v.f(i,j)+v.f(j,i+len);
				vmin = min(vmin,tmp);
			}
		}
        //delete data
        frna_delete(gpu);
        delete[] buf;

	#endif	//_CUDA_CALC_DYNALIGN_
	
	 
	
	if (!(update && update->canceled())) 
	{	// if the user cancelled the operation, skip to the end to do delete[] cleanup
		if (save != 0) 
		{
			ofstream sav(save, ios::binary);
	
			// write the save file information so that the sequence can be re-folded,
			// include thermodynamic data to prevent traceback errors
			short vers = safiversion;
			
			// save a version of the save file, start with structure information
			write(&sav, &vers); 
			
			int sequencelength = ct->GetSequenceLength();
			
			write(&sav, &(sequencelength));

			write(&sav, &(ct->intermolecular));

			int pairs = ct->GetNumberofPairs();
			
			write(&sav, &(pairs));

			for (i = 0; i < ct->GetNumberofPairs(); i++) 
			{
				pairs=ct->GetPair5(i);
				write(&sav, &(pairs));
				pairs=ct->GetPair3(i);
				write(&sav, &(pairs));
			}

			pairs = ct->GetNumberofForbiddenPairs();
			
			write(&sav, &(pairs));

			for (i = 0; i < ct->GetNumberofForbiddenPairs(); i++) 
			{
				pairs = ct->GetForbiddenPair5(i);
				write(&sav, &(pairs));
				pairs = ct->GetForbiddenPair3(i);
				write(&sav, &(pairs));
			}

			for (i = 0; i <= ct->GetSequenceLength(); i++) 
			{
				write(&sav, &(ct->hnumber[i]));
				sav.write(&(ct->nucs[i]), 1);
			}

			for (i = 0; i <= 2 * ct->GetSequenceLength(); i++) write(&sav, &(ct->numseq[i]));
	
			int doubles = ct->GetNumberofDoubles();

			write(&sav, &(doubles));

			for (i = 0; i < ct->GetNumberofDoubles(); i++) 
			{
				doubles = ct->GetDouble(i);
				write(&sav, &(doubles));
			}

			if (ct->intermolecular) 
			{
				for (i = 0; i < 3; i++) write(&sav, &(ct->inter[i]));
			}

			int singles = ct->GetNumberofSingles();
			
			write(&sav, &singles);
			
			for (i = 0; i < ct->GetNumberofSingles(); i++) 
			{
				singles = ct->GetSingle(i);	
				write(&sav, &(singles));
			}

			int modified;
			
			modified = ct->GetNumberofModified();

			write(&sav, &modified);
			
			for (i = 0; i < ct->GetNumberofModified(); i++) 
			{
				modified = ct->GetModified(i);
				write(&sav, &(modified));
			}
	
			modified = ct->GetNumberofGU();
			
			write(&sav, &(modified));
			
			for (i = 0; i < ct->GetNumberofGU(); i++) 
			{
				modified = ct->GetGUpair(i);
				write(&sav, &(modified));
			}
	
			string label = ct->GetSequenceLabel();
			
			write(&sav, &(label));

			write(&sav, &(ct->templated));
			
			if (ct->templated) 
			{
				for (i = 0; i <= ct->GetSequenceLength(); i++) 
				{
					for (j = 0; j <= i; j++) write(&sav, &(ct->tem[i][j]));	
				}
			}

			// write the SHAPE data (for pseudo-free energy constraints)
			write(&sav, &(ct->shaped));
			
			if (ct->shaped) 
			{
				for (i = 0; i <= 2 * ct->GetSequenceLength(); i++) write(&sav, &(ct->SHAPE[i]));
			}

			//now write the array class data for v, w, and wmb:
			for (i = 0; i<= ct->GetSequenceLength(); i++) 
			{
				write(&sav, &(w3[i]));
				write(&sav, &(w5[i]));
				
				for (j = 0; j <= ct->GetSequenceLength(); j++) 
				{
					write(&sav, &(v.dg[i][j+i]));
					write(&sav, &(w.dg[i][j+i]));
					write(&sav, &(wmb.dg[i][j+i]));
					writesinglechar(&sav, &(fce.dg[i][j]));

					if (ct->intermolecular) 
					{
						write(&sav, &(w2->dg[i][j+i]));
						write(&sav, &(wmb2->dg[i][j+i]));
					}
				}	
			}

			write(&sav, &(w3[ct->GetSequenceLength() + 1]));
			
			for (i = 0; i <= 2 * ct->GetSequenceLength(); i++) 
			{
				write(&sav, &(lfce[i]));
				write(&sav, &(mod[i]));
			}

			write(&sav, &vmin);

			// write the alphabet information
			// now write the thermodynamic data:
			write(&sav, data);
			sav.close();
		}
		if (quickenergy) 
		{
			// Don't do traceback, just return energy
			ct->AddStructure();
			ct->SetEnergy(1, w5[number]);
			tracebackerror = 0;
		}
		else if (quickstructure)
		{
			// Calculate only the lowest free energy structure
			tracebackerror = trace(ct, data, 1, ct->GetSequenceLength(), &v, &w, &wmb, w2, wmb2, lfce, &fce, w3, w5, mod, true);
		}
		else tracebackerror = traceback(ct, data, &v, &w, &wmb, w2, wmb2, w3, w5, &fce, lfce, vmin, cntrl6, cntrl8, cntrl9, mod);
	}
	
	//arraydump(v, w, wmb, w5, w3, ct->GetSequenceLength(), "post_traceback.out");

	#else//INSTRUMENTED IS DEFINED

	fill(ct, *v, *vmb, *vext, w, wmb, fce, vmin,lfce, mod,w5, w3, quickenergy, data, w2, wmb2, update, maxinter);
	
	/*
	//write a dot plot file that tracks MB loops and exterior loops
	
	ofstream sav(save);
	double pairs;
	double mbpairs;
	double extpairs;
		
	for (i=1;i<=ct->GetSequenceLength();i++) {
		for (j=(i+1);j<=ct->GetSequenceLength();j++) {
				
			if ((v.f(i,j)+v.f(j,i+ct->GetSequenceLength()))<0) {
					
				pairs= ((double)(v.f(i,j)+v.f(j,i+ct->GetSequenceLength()))/10.0);
				mbpairs=min( ((double)(vmb->f(i,j)+v.f(j,i+ct->GetSequenceLength()))/10.0), ((double)(v.f(i,j)+vmb->f(j,i+ct->GetSequenceLength()))/10.0));
				extpairs= ((double)(v.f(i,j)+vext->f(j,i+ct->GetSequenceLength()))/10.0);
				sav << i << "\t" << j << "\t" << pairs<<"\t"<< mbpairs<< "\t"<<extpairs<<"\n"; 
					
					
			}
			}
			}

			sav.close();
	*/


	#endif	//END INSTRUMENTED 

	delete[] lfce;
	delete[] mod;
	delete[] w5;
	delete[] w3;

	if (ct->intermolecular) 
	{
		delete w2;
		delete wmb2;
	}

	#undef can_pair
	#ifdef timer
		timeout << time(NULL)<<"\n";
		timeout << time(NULL) - seconds;
		timeout.close();
	#endif

	#ifndef INSTRUMENTED // if INSTRUMENTED is not defined
	return tracebackerror;
	#endif // if INSTRUMENTED is defined
}

#else  // #ifndef DYNALIGN_II
#endif

//	The fill routine is encapsulated in function fill.
//	This was separated from dynamic on 3/12/06 by DHM.  This provides greater flexibility
//	for use of the arrays for other tasks than secondary structure prediction, e.g. dot plots.
#if defined DYNALIGN_II
void fill(structure *ct, DynProgArray<integersize> &v, DynProgArray<integersize> &w, DynProgArray<integersize> &wmb, forceclass &fce, int &vmin,bool *lfce, bool *mod,
          integersize *w5, integersize *w3, bool quickenergy,
          datatable *data, DynProgArray<integersize> *w2, DynProgArray<integersize> *wmb2, DynProgArray<integersize> *we,ProgressHandler* update, int maxinter, bool quickstructure, bool simple_iloops, bool disablecoax, bool allow_isolated)

#elif !defined INSTRUMENTED//If pre-compiler flag INSTRUMENTED is not defined, compile the following code
	void fill(structure *ct, DynProgArray<integersize> &v, DynProgArray<integersize> &w, DynProgArray<integersize> &wmb, forceclass &fce, int &vmin,bool *lfce, bool *mod,
			integersize *w5, integersize *w3, bool quickenergy,
			datatable *data, DynProgArray<integersize> *w2, DynProgArray<integersize> *wmb2, ProgressHandler* update, int maxinter,bool quickstructure, bool simple_iloops, bool disablecoax, bool allow_isolated)

#else //IF DEFINED INSTRUMENTED
		void fill(structure *ct, DynProgArray<integersize> &v, DynProgArray<integersize> &vmb, DynProgArray<integersize> &vext, DynProgArray<integersize> &w, DynProgArray<integersize> &wmb, forceclass &fce, int &vmin,bool *lfce, bool *mod,
				integersize *w5, integersize *w3, bool quickenergy,
				datatable *data, DynProgArray<integersize> *w2, DynProgArray<integersize> *wmb2, ProgressHandler* update, int maxinter,bool quickstructure, bool simple_iloops, bool disablecoax, bool allow_isolated)

#endif //end !INTRUMENTED

{
	integersize(*ptr_to_erg1)(int i, int j, int ip, int jp, structure * ct, datatable * data);
	integersize(*ptr_to_erg2)(int i, int j, int ip, int jp, structure * ct, datatable * data, char a, char b);
	integersize(*ptr_to_erg3)(int i, int j, structure * ct, datatable * data, char dbl);
	integersize(*ptr_to_erg4)(int i, int j, int ip, int jp, structure * ct, datatable * data, bool lfce);
	integersize(*ptr_to_ergcoaxflushbases) (int i, int j, int ip, int jp, structure * ct, datatable * data);
	integersize(*ptr_to_ergcoaxinterbases1) (int i, int j, int ip, int jp, structure * ct, datatable * data);
	integersize(*ptr_to_ergcoaxinterbases2) (int i, int j, int ip, int jp, structure * ct, datatable * data);
	integersize(*ptr_to_penalty) (int i, int j, structure * ct, datatable * data);

	if (ct->GetNumberofSequences() == 1) {
		ptr_to_erg1 = erg1;
		ptr_to_erg2 = erg2;
		ptr_to_erg3 = erg3;
		ptr_to_erg4 = erg4;
		ptr_to_ergcoaxflushbases = ergcoaxflushbases;
		ptr_to_ergcoaxinterbases1 = ergcoaxinterbases1;
		ptr_to_ergcoaxinterbases2 = ergcoaxinterbases2;
		ptr_to_penalty = penalty;
	}
	else if (ct->GetNumberofSequences() > 1) {
		ptr_to_erg1 = multi_erg1;
		ptr_to_erg2 = multi_erg2;
		ptr_to_erg3 = multi_erg3;
		ptr_to_erg4 = multi_erg4;
		ptr_to_ergcoaxflushbases = multi_ergcoaxflushbases;
		ptr_to_ergcoaxinterbases1 = multi_ergcoaxinterbases1;
		ptr_to_ergcoaxinterbases2 = multi_ergcoaxinterbases2;
		ptr_to_penalty = multi_penalty;
	}

	int h, maximum, d, maxj;
	integersize **wca, **curE, **prevE, **tempE;
	int number = ct->GetSequenceLength();

	// initialize the multibranch loop parameters, multiply by the number of sequences in CT for alignmentfold calculations:
	
	// the free energy per branch in a multibranch loop:
	integersize mbl_branch = multi_eparam(10, ct); 

	// the free energy cost applied to each multibranch loop:
	integersize mbl_closure = multi_eparam(5, ct); 

	// the free energycost per nucleotide in a a
	integersize mbl_nucleotide = multi_eparam(6, ct); 

	if (ct->GetThermodynamicDataTable() != data) 
	{
		cerr << "In fill (algorithm.cpp): The structure's datatable does not match the passed-in datatable. This can cause problems with IsNuc etc." << endl;
		ct->SetThermodynamicDataTable(data);
	}

	// transposed v and w arrays
	// allows fast traversal of columns
	// vT is used for coaxial stacking search along columns
	// wT is used for coaxial stacking and multibranch loops
	DynProgArrayT<integersize> vT = DynProgArrayT<integersize>(number);
	DynProgArrayT<integersize> wT = DynProgArrayT<integersize>(number);

	// can_pair is used to information about base pairing
	vector<vector<bool>> can_pair_array(2*number+1, vector<bool>(2*number+1, false));
	
	for (int r = 0; r < 2 * number; r++) {
		for (int s = 0; s < 2 * number; s++) {
			if (ct->GetNumberofSequences() == 1) {
				can_pair_array[r + 1][s + 1] = data->can_pair(r + 1, s + 1, ct->numseq);
			}
			else if (ct->GetNumberofSequences() > 1) {
				can_pair_array[r + 1][s + 1] = ct->pairing_matrix[r+1][s+1];
			}
		}
	}
	
	

	if (!ct->intermolecular)
	{
		//This code is needed for O(N^3) prediction of internal loops
		wca = new integersize * [number + 1];
		curE = new integersize * [number + 1];
		prevE = new integersize * [number + 1];

		for (int locali = 0; locali <= number; locali++)
		{
			wca[locali] = new integersize[number + 1];
			curE[locali] = new integersize[number + 1];
			prevE[locali] = new integersize[number + 1];

			for (int localj = 0; localj <= number; localj++)
			{
				wca[locali][localj] = INFINITE_ENERGY;
				curE[locali][localj] = INFINITE_ENERGY;
				prevE[locali][localj] = INFINITE_ENERGY;
			}
		}
	}
	else // intermolecular folding
	{	
		// initialize wca only
		wca = new integersize * [number + 1];
		for (int locali = 0; locali <= number; locali++)
		{
			wca[locali] = new integersize[number + 1];
			for (int localj = 0; localj <= number; localj++)
			{
				wca[locali][localj] = INFINITE_ENERGY;
			}
		}
	}

	if (quickenergy || quickstructure) maximum = number;
	else maximum = (2 * number - 1);

	vmin = INFINITE_ENERGY;
	
	// h or (h - number + 1) is the distance between i and j
	for (h = 0; h < maximum; h++) 
	{	
		// d = j - i;
		// d = (h <= (number - 1)) ? h : (h - number + 1);
		d = (h <= number-1) ? h : (h - number + 1);
		
		if (((h % 10) == 0) && update) update->update((100 * h) / (maximum + 1));
		if (update && update->canceled()) break; // If user cancelled calculation, then break out of for loop and do cleanup at end.

		if (h == number && !ct->intermolecular)
		{
			for (int locali = 0; locali <= number; locali++)
			{
				for (int localj = 0; localj <= number; localj++)
				{
					curE[locali][localj] = INFINITE_ENERGY;
					prevE[locali][localj] = INFINITE_ENERGY;
				}
			}
		}

		// These variables for start and end (as opposed the fancy syntax that was here before) are needed because of openmp.
		int startme, endme;//start and end for loop over i
		// if (h <= (number - 1))
		if (h <= number-1)
		{
			startme = 1;
			endme = number - h;
		}
		else
		{
			startme = 2 * number - h;
			endme = number;
		}

		#ifdef SMP
		#pragma omp parallel for
		#endif

		for (int i = startme; i <= endme; i++)
		{	
		//	cout << h << " " << i << endl;
			int j = i + d;
			register int rarray;
			int dp, ll, jpf, jf, bl, maxasym;
			int before, after;
			int e[6];
			int k, p;
			int ip, jp, ii, jj, di;

		//  #define can_pair(i, j) (multi_can_pair(i, j, ct))
		//	#define can_pair(i, j) (ct->can_pair(i, j, ct))
			#define can_pair(i, j) can_pair_array[i][j]
		
			#ifndef disablecoax
			int castack; //variable for coaxial stacking, not needed if coaxial stacking is disabled
			#endif
			// int maxi = ((h <= (number - 1)) ? (number - h) : number);
			int maxi = (h < number) ? (number - h) : number;
			if (i == maxi)	maxj = j;
			maxasym = maxinter;

			if (ct->templated)
			{
				if (i > ct->GetSequenceLength()) ii = i - ct->GetSequenceLength();
				else ii = i;
				if (j > ct->GetSequenceLength()) jj = j - ct->GetSequenceLength();
				else jj = j;
				if (jj < ii)
				{
					p = jj;
					jj = ii;
					ii = p;
				}
				if (!ct->tem[jj][ii]) goto sub2;
			}

			// Compute v[i][j], the minimum energy of the substructure from i to j,
			// inclusive, where i and j are base paired

			if (fce.f(i, j) & SINGLE)
			{
				//i or j is forced single-stranded
				v.f(i, j) = vT.f(i, j) = INFINITE_ENERGY + 50;
				goto sub2;
			}
			if (fce.f(i, j) & NOPAIR)
			{
				//i or j is forced into a pair elsewhere
				v.f(i, j) = vT.f(i, j) = INFINITE_ENERGY + 50;
				goto sub2;
			}
			if (j <= number)
			{
				if ((j - i) <= minloop) goto sub3;
			}
			
			v.f(i, j) = vT.f(i, j) = INFINITE_ENERGY;

			if (!can_pair(i, j)) goto sub2;

			//force u's into gu pairs
			for (ip = 0; ip < ct->GetNumberofGU(); ip++)
			{
				if (ct->GetGUpair(ip) == i)
				{
					if (ct->numseq[j] != 3)
					{
						v.f(i, j) = vT.f(i, j) = INFINITE_ENERGY;
						goto sub2;
					}
				}
				else if (ct->GetGUpair(ip) == j)
				{
					if (ct->numseq[i] != 3)
					{
						v.f(i, j) = vT.f(i, j) = INFINITE_ENERGY;
						goto sub2;
					}
				}
				else if ((ct->GetGUpair(ip) + number) == j)
				{
					if (ct->numseq[i] != 3)
					{
						v.f(i, j) = vT.f(i, j) = INFINITE_ENERGY;
						goto sub2;
					}
				}
			}

			//now check to make sure that this isn't an isolated pair:
			//	(consider a pair separated by a bulge as not! stacked)

			//before = 0 if a stacked pair cannot form 5' to i
			if (allow_isolated) before = 1;
			else if ((i > 1 && j < (2 * number) && j != number)) {
				if ((j > number && ((i - j + number) > minloop + 2)) || j < number) {
					before = can_pair(i - 1, j + 1);
				}
				else before = 0;
			}
			else before = 0;

			//after = 0 if a stacked pair cannot form 3' to i
			if (allow_isolated) after = 1;
			else if ((((j - i) > minloop + 2) && (j <= number) || (j > number + 1)) && (i != number)) {
				after = can_pair(i + 1, j - 1);
			}
			else after = 0;
			
			if (ct->GetNumberofSequences() > 1) {
				if (before == 0 && after == 0) {
					if (j > number)
					{
						before = ct->can_pair_isolated(j - number, i, ct);
					}
					else
						before = ct->can_pair_isolated(i, j, ct);
				}
			}						
			//if there are no stackable pairs to i.j then don't allow a pair i,j
			if ((before == 0) && (after == 0))
			{
				//v.f(i,j)= 0;
				goto sub2;
			}


				rarray = INFINITE_ENERGY;

			if (i == (number) || j == ((number) + 1)) goto sub1;

			//Perhaps i and j close a hairpin:
			rarray = min(rarray, ptr_to_erg3(i, j, ct, data, fce.f(i, j)));

			if ((j - i - 1) >= (minloop + 2) || j > (number))
			{
				// Perhaps i,j stacks over i+1,j-1
				if (!mod[i] && !mod[j])
				{	
					int xx = ptr_to_erg1(i, j, i + 1, j - 1, ct, data);
					int yy = v.f(i + 1, j - 1);
					//make sure this is not a site of chemical modification
					rarray = min(rarray, (ptr_to_erg1(i, j, i + 1, j - 1, ct, data) + v.f(i + 1, j - 1)));
				}
				else
				{
					//allow G-U to be modified or a pair next to a G-U to be modified
					if (((ct->IsNuc(i, 'G') || ct->IsNuc(i, 'g')) && (ct->IsNuc(j, 'U') || ct->IsNuc(j, 'u'))) ||
						((ct->IsNuc(i, 'U') || ct->IsNuc(i, 'u')) && (ct->IsNuc(j, 'G') || ct->IsNuc(j, 'g'))))
					{
						rarray = min(rarray, (ptr_to_erg1(i, j, i + 1, j - 1, ct, data) + v.f(i + 1, j - 1)));
					}
					else if (((ct->IsNuc(i + 1, 'G') || ct->IsNuc(i + 1, 'g')) && (ct->IsNuc(j - 1, 'U') || ct->IsNuc(j - 1, 'u'))) ||
						((ct->IsNuc(i + 1, 'U') || ct->IsNuc(i + 1, 'u')) && (ct->IsNuc(j - 1, 'G') || ct->IsNuc(j - 1, 'g'))))
					{
						rarray = min(rarray, (ptr_to_erg1(i, j, i + 1, j - 1, ct, data) + v.f(i + 1, j - 1)));
					}
					else if (i - 1 > 0 && j + 1 < 2 * number)
					{
						if (((ct->IsNuc(i - 1, 'G') || ct->IsNuc(i - 1, 'g')) && (ct->IsNuc(j + 1, 'U') || ct->IsNuc(j + 1, 'u'))) ||
							((ct->IsNuc(i - 1, 'U') || ct->IsNuc(i - 1, 'u')) && (ct->IsNuc(j + 1, 'G') || ct->IsNuc(j + 1, 'g'))))
						{
							rarray = min(rarray, (ptr_to_erg1(i, j, i + 1, j - 1, ct, data) + v.f(i + 1, j - 1)));
						}
					}
				}
			}

			// Perhaps i,j closes an interior or bulge loop, search for the best possibility
			// If this is intermolecular folding, or if running in parallel with openMP, revert to the old O(N^4) algorithm
			// For now, also use O(N^4) internal loops for applying SHAPE to single-stranded regions
			if (simple_iloops || ct->intermolecular || (ct->shaped && (ct->SHAPEintercept_ss != 0 || ct->SHAPEslope_ss != 0)) || ct->ssoffset)
			{
				if (((j - i - 1) >= (minloop + 3)) || (j > (number)))
				{
					for (di = (j - i - 3); di >= 1; di--)
					{
						for (ip = (i + 1); ip <= (j - 1 - di); ip++)
						{
							jp = di + ip;
							if (!can_pair(ip, jp)) continue;
							if ((j - i - 2 - di) > (multi_eparam(7, ct))) goto sub1;
							if (abs(ip - i + j - jp) <= (maxinter))
							{
								if (ip > number)
								{
									rarray = min( 
												rarray, 
												(ptr_to_erg2(i, j, ip, jp, ct, data, fce.f(i, ip - number), fce.f(jp - number, j - number)) 
												+ v.f(ip - (number), jp - (number))));

									if (mod[ip] || mod[jp])
									{
										//ip or jp is modified

										rarray = min(
													rarray, 
													ptr_to_erg2(i, j, ip, jp, ct, data, fce.f(i, ip - number), fce.f(jp - number, j - number)) 
													+ v.f(ip - (number)+1, jp - (number)-1) 
													+ ptr_to_erg1(ip - number, jp - number, ip + 1 - number, jp - 1 - number, ct, data));
									}
								}
								else
								{
									if (jp <= number)
									{
										rarray = min(
													rarray, 
													ptr_to_erg2(i, j, ip, jp, ct, data, fce.f(i, ip), fce.f(jp, j)) 
													+ v.f(ip, jp));

										if (mod[ip] || mod[jp])
										{
											//i or j is modified
											rarray = min(
														rarray, 
														ptr_to_erg2(i, j, ip, jp, ct, data, fce.f(i, ip), fce.f(jp, j)) 
														+ v.f(ip + 1, jp - 1) 
														+ ptr_to_erg1(ip, jp, ip + 1, jp - 1, ct, data));
										}
									}
									else
									{
										rarray = min(
													rarray, 
													ptr_to_erg2(i, j, ip, jp, ct, data, fce.f(i, ip), fce.f(jp - number, j - number)) 
													+ v.f(ip, jp));

										if (mod[ip] || mod[jp])
										{
											//i or j is modified
											rarray = min(
												rarray,
														ptr_to_erg2(i, j, ip, jp, ct, data, fce.f(i, ip), fce.f(jp - number, j - number)) 
														+ v.f(ip + 1, jp - 1) 
														+ ptr_to_erg1(ip, jp, ip + 1, jp - 1, ct, data));
										}
									}
								}
							}
						}
					}
				}
			}

			/*

				Perhaps i,j closes an interior or bulge loop, search for the best possibility
				fill the interior loops' energy rarray first
				calculate the small loop (size<=5) first
				the larger loop is prefilled with curE[dp][i] (after sub2)

				d = j - i, dp = jp - ip (interior loop)
				i<ip<number<jp<j or i<ip<jp<j<number
			*/

			//Below is the O(N^3) algorithm
			//this is no longer used by default because time benchmarks show it is slower
			//when the maximum internal loop size is 30
			else
			{
				if ((d - 1) >= (minloop + 3) || j > number)
				{
					for (dp = d - 3; dp >= ((j > number) ? 1 : minloop + 1); dp--)
					{
						ll = d - dp - 2;
						//ll is the loop length in terms of number of unpaired nucs
						//calculate every ip,jp when ll <=5: this includes loops with special rules:
						//0x1,0x2,0x3,0x4,0x5,1x1,1x2,1x3,1x4,2x2,2x3
						if (ll >= 1 && ll <= 5)
						{
							for (ip = i + 1; ip <= j - 1 - dp; ip++)
							{
								jp = ip + dp;
								if (can_pair(ip, jp))
								{
									if ((ip <= number && jp > number) || (ip <= number && j <= number))
									{
										//using jpf and jf and  instead of jp and j when j,jp>number
										jpf = ((jp <= number) ? jp : jp - number);
										jf = ((j <= number) ? j : j - number);
										rarray = min(
													rarray,
													ptr_to_erg2(i, j, ip, jp, ct, data, fce.f(i, ip), fce.f(jpf, jf))
													+ v.f(ip, jp));
										
										//i or j is modified
										if ((mod[ip] || mod[jp]) && can_pair(ip, jp) && !(fce.f(ip, jp) & SINGLE))
											rarray = min(
														rarray,
														ptr_to_erg2(i, j, ip, jp, ct, data, fce.f(i, ip), fce.f(jpf, jf))
														+ v.f(ip + 1, jp - 1)
														+ ptr_to_erg1(ip, jp, ip + 1, jp - 1, ct, data));
									}
								}
							}
						}
						// when size >= 6 and <= 30;
						// else if ( ll>=6 && ll <= (data->eparam[7]) )
						else if (ll >= 6 && ll <= maxinter)
						{

							//interior energy prefilled (after sub2:) + exterior engergy (the function for erg2in and erg2ex is found in rna_library_inter.cpp
							//The interior energy, curE, is the lowest free energy for a interior fragment with a pair closing dp nucs
							//the interior fragment and asymetry contriubutions are included

							rarray = min(rarray, (curE[dp][i] + erg2ex(i, j, ll, ct, data)));

							//considering loop 1x(n-1) (bl=1)and 0 x n(bulge)(bl=0) as stacking bonus on 1x(n-1) and bulge is not allowed
							for (bl = 0; bl <= 1; bl++)
							{
								ip = i + 1 + bl;
								jp = ip + dp;
								jpf = ((jp <= number) ? jp : jp - number);
								jf = ((j <= number) ? j : j - number);
								if ((ip <= number && jp > number) || (ip <= number && j <= number))
								{
									if (can_pair(ip, jp))
									{
										if (abs(ip - i + jp - j) <= maxasym)
										{
											rarray = min(
														rarray,
														ptr_to_erg2(i, j, ip, jp, ct, data, fce.f(i, ip), fce.f(jpf, jf))
														+ v.f(ip, jp));
											//i or j is modified
											if ((mod[ip] || mod[jp]) && can_pair(ip, jp) && !(fce.f(ip, jp) & SINGLE))
											{
												rarray = min(
															rarray,
															ptr_to_erg2(i, j, ip, jp, ct, data, fce.f(i, ip), fce.f(jpf, jf))
															+ v.f(ip + 1, jp - 1)
															+ ptr_to_erg1(ip, jp, ip + 1, jp - 1, ct, data));
											}
										}
									}
								}

								jp = j - 1 - bl;
								ip = jp - dp;
								jpf = ((jp <= number) ? jp : jp - number);
								jf = ((j <= number) ? j : j - number);
								if ((ip <= number && jp > number) || (ip <= number && j <= number))
								{
									if (can_pair(ip, jp))
									{
										if (abs(ip - i + jp - j) <= maxasym)
										{
											rarray = min(
														rarray, 
														ptr_to_erg2(i, j, ip, jp, ct, data, fce.f(i, ip), fce.f(jpf, jf)) 
														+ v.f(ip, jp));
											if ((mod[ip] || mod[jp]) && can_pair(ip, jp) && !(fce.f(ip, jp) & SINGLE))
											{
												rarray = min(
															rarray,
															ptr_to_erg2(i, j, ip, jp, ct, data, fce.f(i, ip), fce.f(jpf, jf))
															+ v.f(ip + 1, jp - 1)
															+ ptr_to_erg1(ip, jp, ip + 1, jp - 1, ct, data));
											}
										}
									}
								}
							}
						}
					}
				}
			}


			

			//   ----- Indentation ----- //
			// Perhaps i,j closes a multibranch or exterior loop, search for the best possibility
			sub1:

			#ifndef INSTRUMENTED//

			#else//If INSTRUMENTED IS DEFINED
			v.f(i, j) = vT.f(i, j) = rarray;
			rarray = INFINITE_ENERGY;
			#endif//END INSTRUMENTED


			if (((j - i - 1) >= (2 * minloop + 4)) || (j > (number))) 
			{
				for (ii = 1; ii <= 4; ii++) e[ii] = INFINITE_ENERGY;
				//consider the exterior loop closed by i,j
				if (j > number)
				{
					rarray = min(
							rarray,
							w3[i + 1]
							+ w5[j - number - 1]
							+ ptr_to_penalty(i, j, ct, data));
					if (i != number)
					{
						rarray = min(	
									rarray, 
									ptr_to_erg4(i, j, i + 1, 1, ct, data, lfce[i + 1]) 
									+ ptr_to_penalty(i, j, ct, data) 
									+ w3[i + 2] 
									+ w5[j - number - 1]);							
					}
					if (j != number + 1)
					{
						rarray = min(
									rarray,
									ptr_to_erg4(i, j, j - 1, 2, ct, data, lfce[j - 1])
									+ ptr_to_penalty(i, j, ct, data)
									+ w3[i + 1]
									+ w5[j - number - 2]);
					}
					if (i != number && j != number + 1) 
					{
						rarray = min(
									rarray, 
								//	data->tstack[ct->numseq[i]][ct->numseq[j]][ct->numseq[i + 1]][ct->numseq[j - 1]]
									multi_tstack(i, j, i + 1, j - 1, ct)
									+ checknp(lfce[i + 1], lfce[j - 1]) 
									+ w3[i + 2] 
									+ w5[j - number - 2]
									+ ptr_to_penalty(i, j, ct, data) 
									+ ct->SHAPEss_give_value(i + 1) 
									+ ct->SHAPEss_give_value(j - 1));				
					}

					//consider the coaxial stacking of a helix from i to ip onto helix ip+1 or ip+2 to j:
					//first consider a helix stacking from the 5' sequence fragment:
					#ifndef disablecoax
					if (!disablecoax) 
					{
						for (ip = j - number - minloop - 1; ip > 0; ip--) 
						{
							//first consider flush stacking
							rarray = min(								
										rarray, 
										w3[i + 1] 
										+ w5[ip - 1] 
										+ ptr_to_penalty(i, j, ct, data) 
										+ ptr_to_penalty(j - number - 1, ip, ct, data) 
										+ ptr_to_ergcoaxflushbases(ip, j - number - 1, j - number, i, ct, data) 
										+ vT.f(ip, j - number - 1));

							if ((mod[ip] || mod[j - number - 1]) && j - number - 2 > 0 && !(fce.f(ip, j - number - 1) & SINGLE)) 
							{
								if (can_pair(ip + 1, j - number - 2)) 
								{
									rarray = min(	
												rarray, 
												w3[i + 1] 
												+ w5[ip - 1] 
												+ ptr_to_penalty(i, j, ct, data) 
												+ ptr_to_penalty(j - number - 1, ip, ct, data) 
												+ ptr_to_ergcoaxflushbases(ip, j - number - 1, j - number, i, ct, data) + vT.f(ip + 1, j - number - 2) 
												+ ptr_to_erg1(ip, j - number - 1, ip + 1, j - number - 2, ct, data));											
								}
							}
							if (j - number - 2 > 0) 
							{
								//now consider an intervening nuc
								if (i < number) 
								{
									rarray = min(
												rarray, 
												w3[i + 2] 
												+ w5[ip - 1] 
												+ ptr_to_penalty(i, j, ct, data) 
												+ ptr_to_penalty(ip, j - number - 2, ct, data) 
												+ ptr_to_ergcoaxinterbases2(ip, j - number - 2, j - number, i, ct, data) + vT.f(ip, j - number - 2) 
												+ checknp(lfce[i + 1], lfce[j - number - 1]));
											
									if ((mod[ip] || mod[j - number - 2]) && can_pair(ip + 1, j - number - 3) != 0 && !(fce.f(ip, j - number - 2) & SINGLE)) 
									{
										rarray = min(
													rarray, 
													w3[i + 2] 
													+ w5[ip - 1] 
													+ ptr_to_penalty(i, j, ct, data) 
													+ ptr_to_penalty(ip, j - number - 2, ct, data) 
													+ ptr_to_ergcoaxinterbases2(ip, j - number - 2, j - number, i, ct, data) 
													+ vT.f(ip + 1, j - number - 3) 
													+ ptr_to_erg1(ip, j - number - 2, ip + 1, j - number - 3, ct, data) 
													+ checknp(lfce[i + 1], lfce[j - number - 1]));											
									}
								}

								//consider the other possibility for an intervening nuc
								rarray = min(
											rarray,
											w3[i + 1] 
											+ w5[ip - 1] 
											+ ptr_to_penalty(i, j, ct, data) 
											+ ptr_to_penalty(ip + 1, j - number - 2, ct, data) 
											+ ptr_to_ergcoaxinterbases1(ip + 1, j - number - 2, j - number, i, ct, data) 
											+ vT.f(ip + 1, j - number - 2) 
											+ checknp(lfce[ip], lfce[j - number - 1]));
									
								if ((mod[ip + 1] || mod[j - number - 2]) && can_pair(ip + 2, j - number - 3) && !(fce.f(ip + 1, j - number - 2) & SINGLE)) 
								{
									rarray = min(										
												rarray,
												w3[i + 1] 
												+ w5[ip - 1] 
												+ ptr_to_penalty(i, j, ct, data) 
												+ ptr_to_penalty(ip + 1, j - number - 2, ct, data) 
												+ ptr_to_ergcoaxinterbases1(ip + 1, j - number - 2, j - number, i, ct, data) 
												+ vT.f(ip + 2, j - number - 3)
												+ ptr_to_erg1(ip + 1, j - number - 2, ip + 2, j - number - 3, ct, data) 
												+ checknp(lfce[ip], lfce[j - number - 1]));											
								}
							}
						}
						//now consider a helix stacking from the 3' sequence fragment:
						for (ip = i + minloop + 1; ip <= number; ip++) 
						{
							//first consider flush stacking
							rarray = min(								
										rarray,
										w3[ip + 1] 
										+ w5[j - number - 1] 
										+ ptr_to_penalty(i, j, ct, data) 
										+ ptr_to_penalty(ip, i + 1, ct, data) 
										+ ptr_to_ergcoaxflushbases(j - number, i, i + 1, ip, ct, data) 
										+ v.f(i + 1, ip));
								
							if (mod[i + 1] || mod[ip] && can_pair(i + 2, ip - 1) && !(fce.f(i + 1, ip) & SINGLE)) 
							{
								rarray = min(									
											rarray,
											w3[ip + 1] 
											+ w5[j - number - 1] 
											+ ptr_to_penalty(i, j, ct, data) 
											+ ptr_to_penalty(ip, i + 1, ct, data) 
											+ ptr_to_ergcoaxflushbases(j - number, i, i + 1, ip, ct, data) 
											+ v.f(i + 2, ip - 1)
											+ ptr_to_erg1(i + 1, ip, i + 2, ip - 1, ct, data));										
							}
							//now consider an intervening nuc
							if (j - number > 1) 
							{
								rarray = min(
											rarray,
											w3[ip + 1] 
											+ w5[j - number - 2] 
											+ ptr_to_penalty(i, j, ct, data) 
											+ ptr_to_penalty(ip, i + 2, ct, data) 
											+ ptr_to_ergcoaxinterbases1(j - number, i, i + 2, ip, ct, data) 
											+ v.f(i + 2, ip) 
											+ checknp(lfce[i + 1], lfce[j - number - 1]));
		
								if ((mod[i + 2] || mod[ip]) && can_pair(i + 3, ip - 1) && !(fce.f(i + 2, ip) & SINGLE)) 
								{
									rarray = min(										
												rarray,
												w3[ip + 1] 
												+ w5[j - number - 2] 
												+ ptr_to_penalty(i, j, ct, data) 
												+ ptr_to_penalty(ip, i + 2, ct, data) 
												+ ptr_to_ergcoaxinterbases1(j - number, i, i + 2, ip, ct, data) + v.f(i + 3, ip - 1)
												+ ptr_to_erg1(i + 2, ip, i + 3, ip - 1, ct, data) 
												+ checknp(lfce[i + 1], lfce[j - number - 1]));											
								}
							}
							//consider the other possibility for an intervening nuc
							rarray = min(								
										rarray,
										w3[ip + 1] 
										+ w5[j - number - 1] 
										+ ptr_to_penalty(i, j, ct, data) 
										+ ptr_to_penalty(ip - 1, i + 2, ct, data) 
										+ ptr_to_ergcoaxinterbases2(j - number, i, i + 2, ip - 1, ct, data) 
										+ v.f(i + 2, ip - 1) 
										+ checknp(lfce[i + 1], lfce[ip]));
						
							if ((mod[i + 2] || mod[ip - 1]) && can_pair(i + 3, ip - 2) && !(fce.f(i + 2, ip - 1) & SINGLE)) 
							{
								rarray = min(									
											rarray,
											w3[ip + 1] 
											+ w5[j - number - 1] 
											+ ptr_to_penalty(i, j, ct, data) 
											+ ptr_to_penalty(ip - 1, i + 2, ct, data) 
											+ ptr_to_ergcoaxinterbases2(j - number, i, i + 2, ip - 1, ct, data) 
											+ v.f(i + 3, ip - 2)
											+ ptr_to_erg1(i + 2, ip - 1, i + 3, ip - 2, ct, data) 
											+ checknp(lfce[i + 1], lfce[ip]));									
							}
						}
					}
					#endif //ifndef disablecoax
				}

				#ifndef INSTRUMENTED
				#else//IF INSTRUMENTED IS DEFINED
				v.f(i, j) = vT.f(i, j) = min(rarray, v.f(i, j));//QUESTION: DO WE NEED THIS? WAS IN algorithm.napss.cpp
				vext.f(i, j) = rarray;
				rarray = INFINITE_ENERGY;//QUESTION: DO WE NEED THIS? WAS IN algorithm.napss.cpp
				#endif//END INSTRUMENTED

				// consider the multiloop closed by i,j
				if ((j - i) > (2 * minloop + 4) && i != number) 
				{
					//no dangling ends on i-j pair:
					if (j - 1 != number) 
					{
						rarray = min(							
									rarray, 
									wmb.f(i + 1, j - 1) 
									+ mbl_closure + mbl_branch
									+ ptr_to_penalty(i, j, ct, data));							

						//i+1 dangles on i-j pair:
						if (i + 1 != number)
						{
							rarray = min(
										rarray, 
										ptr_to_erg4(i, j, i + 1, 1, ct, data, lfce[i + 1]) 
										+ ptr_to_penalty(i, j, ct, data) 
										+ wmb.f(i + 2, j - 1) 
										+ mbl_closure 
										+ mbl_nucleotide 
										+ mbl_branch);								
						}
					}
					if (j - 2 != number) 
					{
						//j-1 dangles
						if (j != (number + 1)) 
						{
							rarray = min(
										rarray, 
										ptr_to_erg4(i, j, j - 1, 2, ct, data, lfce[j - 1]) 
										+ ptr_to_penalty(i, j, ct, data) 
										+ wmb.f(i + 1, j - 2) 
										+ mbl_closure 
										+ mbl_nucleotide 
										+ mbl_branch);								
						}
						//both i+1 and j-1 dangle
						if ((i + 1 != number) && (j != (number + 1))) 
						{
							rarray = min(
										rarray,
									//	data->tstkm[ct->numseq[i]][ct->numseq[j]][ct->numseq[i + 1]][ct->numseq[j - 1]] 
										+ multi_tstkm(i, j, i + 1, j - 1, ct)
										+ checknp(lfce[i + 1], lfce[j - 1]) 
										+ wmb.f(i + 2, j - 2) 
										+ mbl_closure 
										+ 2 * mbl_nucleotide 
										+ mbl_branch
										+ ptr_to_penalty(i, j, ct, data) 
										+ ct->SHAPEss_give_value(i + 1) 
										+ ct->SHAPEss_give_value(j - 1));						
						}
					}

					//consider the coaxial stacking of a helix from i to j onto helix i+1 or i+2 to ip:
					#ifndef disablecoax
					if (!disablecoax) 
					{
						for (ip = i + 1; (ip < j); ip++) 
						{
							//first consider flush stacking
							//conditions guarantee that the coaxial stacking isn't considering an exterior loop 
							//if ((i!=number)/*&&(i+1!=number)*//*&&((j>number)||(ip!=number)&&(ip+1!=number))&&(j-1!=number)*/) 
							if (i != number && ip != number && j - 1 != number) 
							{
								if (can_pair(i + 1, ip)) 
								{
									//only proceed if i+1 can pair with ip
									rarray = min(										
												rarray, 
												ptr_to_penalty(i, j, ct, data) 
												+ v.f(i + 1, ip) 
												+ ptr_to_penalty(i + 1, ip, ct, data) 
												+ mbl_closure
												+ 2 * mbl_branch 
												+ wT.f(ip + 1, j - 1) 
												+ ptr_to_ergcoaxflushbases(j, i, i + 1, ip, ct, data));
										
									if ((mod[i + 1] || mod[ip]) && can_pair(i + 2, ip - 1) && !(fce.f(i + 1, ip) & SINGLE)) 
									{
										rarray = min(
													rarray, 
													ptr_to_penalty(i, j, ct, data) 
													+ v.f(i + 2, ip - 1) 
													+ ptr_to_penalty(i + 1, ip, ct, data) 
													+ mbl_closure
													+ 2 * mbl_branch 
													+ wT.f(ip + 1, j - 1) 
													+ ptr_to_ergcoaxflushbases(j, i, i + 1, ip, ct, data)
													+ ptr_to_erg1(i + 1, ip, i + 2, ip - 1, ct, data));
									}
								}
								//if ((ip<j-1)&&(i+2!=number)) 
								if (can_pair(i + 2, ip) && ip + 2 < j - 1 && i + 1 != number && ip + 1 != number) 
								{
									//now consider an intervening nuc
									if ((ip + 2 < j - 1)/*&&(j>number||ip+2!=number)*/)
									{
										rarray = min(	
													rarray,
													ptr_to_penalty(i, j, ct, data)
													+ v.f(i + 2, ip)
													+ ptr_to_penalty(i + 2, ip, ct, data)
													+ mbl_closure
													+ 2 * mbl_nucleotide
													+ 2 * mbl_branch
													+ wT.f(ip + 2, j - 1)
													+ ptr_to_ergcoaxinterbases2(j, i, i + 2, ip, ct, data)
													+ checknp(lfce[i + 1], lfce[ip + 1]));
											
									}
									if ((mod[i + 2] || mod[ip]) && can_pair(i + 3, ip - 1) && !(fce.f(i + 2, ip) & SINGLE)) 
									{
										rarray = min(											
													rarray, 
													ptr_to_penalty(i, j, ct, data) 
													+ v.f(i + 3, ip - 1) 
													+ ptr_to_penalty(i + 2, ip, ct, data) 
													+ mbl_closure
													+ 2 * mbl_nucleotide 
													+ 2 * mbl_branch 
													+ wT.f(ip + 2, j - 1)
													+ ptr_to_ergcoaxinterbases2(j, i, i + 2, ip, ct, data)
													+ ptr_to_erg1(i + 2, ip, i + 3, ip - 1, ct, data) 
													+ checknp(lfce[i + 1], lfce[ip + 1]));
											
									}
									if (ip + 1 < j - 2 && j - 2 != number)
									{
										rarray = min(																						
													rarray,
													ptr_to_penalty(i, j, ct, data)
													+ v.f(i + 2, ip)
													+ ptr_to_penalty(i + 2, ip, ct, data)
													+ mbl_closure
													+ 2 * mbl_nucleotide
													+ 2 * mbl_branch
													+ wT.f(ip + 1, j - 2)
													+ ptr_to_ergcoaxinterbases1(j, i, i + 2, ip, ct, data)
													+ checknp(lfce[i + 1], lfce[j - 1]));
											
									}
									if ((mod[i + 2] || mod[ip]) && can_pair(i + 3, ip - 1) && !(fce.f(i + 2, ip) & SINGLE)) 
									{
										rarray = min(																						
													rarray, 
													ptr_to_penalty(i, j, ct, data) 
													+ v.f(i + 3, ip - 1) 
													+ ptr_to_penalty(i + 2, ip, ct, data) 
													+ mbl_closure
													+ 2 * mbl_nucleotide 
													+ 2 * mbl_branch 
													+ wT.f(ip + 1, j - 2)
													+ ptr_to_ergcoaxinterbases1(j, i, i + 2, ip, ct, data)
													+ ptr_to_erg1(i + 2, ip, i + 3, ip - 1, ct, data) 
													+ checknp(lfce[i + 1], lfce[j - 1]));											
									}
								} //end (ip+2<j-1&&i+1!=number&&ip+1!=number)
							}
						}

						//consider the coaxial stacking of a helix from i to j onto helix ip to j-2 or j-1:
						for (ip = j - 1; ip > i; ip--) 
						{
						//	cout << ip << "ip	i" << i << endl;
							//conditions guarantee that the coaxial stacking isn't considering an exterior loop
							//if ((i!=number)&&(i+1!=number)&&((j>number)||(ip!=number)&&(ip-1!=number))&&(j-1!=number)) {}
							if (j - 1 != number && ip - 1 != number && i != number) 
							{
								//first consider flush stacking
								if (can_pair(j - 1, ip)) 
								{
									//only proceed if j-1 can pair to ip
									rarray = min(										
												rarray, 
												ptr_to_penalty(i, j, ct, data) 
												+ vT.f(ip, j - 1) 
												+ ptr_to_penalty(j - 1, ip, ct, data) 
												+ mbl_closure
												+ 2 * mbl_branch 
												+ w.f(i + 1, ip - 1) 
												+ ptr_to_ergcoaxflushbases(ip, j - 1, j, i, ct, data));
										
									if ((mod[ip] || mod[j - 1]) && can_pair(ip + 1, j - 2) && !(fce.f(ip, j - 1) & SINGLE)) 
									{
										rarray = min(											
													rarray, 
													ptr_to_penalty(i, j, ct, data) 
													+ vT.f(ip + 1, j - 2) 
													+ ptr_to_penalty(j - 1, ip, ct, data) 
													+ mbl_closure
													+ 2 * mbl_branch 
													+ w.f(i + 1, ip - 1) 
													+ ptr_to_ergcoaxflushbases(ip, j - 1, j, i, ct, data)
													+ ptr_to_erg1(ip, j - 1, ip + 1, j - 2, ct, data));												
									}
								}
								if (can_pair(j - 2, ip) && j - 2 != number) 
								{
									//now consider an intervening nuc
									//if ((ip>i+1)&&(j>number||ip-2!=number))
									if ((ip - 2 > i + 1) && ip - 2 != number) 
									{
										rarray = min(
													rarray, 
													ptr_to_penalty(i, j, ct, data) 
													+ vT.f(ip, j - 2) 
													+ ptr_to_penalty(j - 2, ip, ct, data) 
													+ mbl_closure
													+ 2 * mbl_nucleotide 
													+ 2 * mbl_branch 
													+ w.f(i + 1, ip - 2)
													+ ptr_to_ergcoaxinterbases1(ip, j - 2, j, i, ct, data) 
													+ checknp(lfce[j - 1], lfce[ip - 1]));
											
										if ((mod[ip] || mod[j - 2]) && can_pair(ip + 1, j - 3) && !(fce.f(ip, j - 2) & SINGLE)) 
										{
											rarray = min(												
														rarray, 
														ptr_to_penalty(i, j, ct, data) 
														+ vT.f(ip + 1, j - 3) 
														+ ptr_to_penalty(j - 2, ip, ct, data) 
														+ mbl_closure
														+ 2 * mbl_nucleotide 
														+ 2 * mbl_branch 
														+ w.f(i + 1, ip - 2)
														+ ptr_to_ergcoaxinterbases1(ip, j - 2, j, i, ct, data)
														+ ptr_to_erg1(ip, j - 2, ip + 1, j - 3, ct, data) 
														+ checknp(lfce[j - 1], lfce[ip - 1]));												
										}
									}
									if ((ip - 1 > i + 2) && i + 1 != number) 
									{
										rarray = min(
													rarray, 
													ptr_to_penalty(i, j, ct, data) 
													+ vT.f(ip, j - 2) 
													+ ptr_to_penalty(j - 2, ip, ct, data) 
													+ mbl_closure
													+ 2 * mbl_nucleotide 
													+ 2 * mbl_branch 
													+ w.f(i + 2, ip - 1)
													+ ptr_to_ergcoaxinterbases2(ip, j - 2, j, i, ct, data) 
													+ checknp(lfce[j - 1], lfce[i + 1]));
											
										if ((mod[ip] || mod[j - 2]) && can_pair(ip + 1, j - 3) && !(fce.f(ip, j - 2) & SINGLE)) 
										{
											rarray = min(																								
														rarray, 
														ptr_to_penalty(i, j, ct, data) 
														+ vT.f(ip + 1, j - 3) 
														+ ptr_to_penalty(j - 2, ip, ct, data) 
														+ mbl_closure
														+ 2 * mbl_nucleotide 
														+ 2 * mbl_branch 
														+ w.f(i + 2, ip - 1)
														+ ptr_to_ergcoaxinterbases2(ip, j - 2, j, i, ct, data)
														+ ptr_to_erg1(ip, j - 2, ip + 1, j - 3, ct, data) 
														+ checknp(lfce[j - 1], lfce[i + 1]));												
										}
									}
								}
							}
						}
					}
					#endif //ifndef disablecoax
					
					if (ct->intermolecular) 
					{
						//intermolecular, so consider wmb2,
						//don't add the multiloop penalties because this is a exterior loop
						rarray = min(
									rarray, 
									wmb2->f(i + 1, j - 1) 
									+ ptr_to_penalty(i, j, ct, data) 
									+ INFINITE_ENERGY);
							
						//i+1 dangles on i-j pair:
						if (i != number) 
						{
							rarray = min(
										rarray,
										ptr_to_erg4(i, j, i + 1, 1, ct, data, lfce[i + 1])
										+ ptr_to_penalty(i, j, ct, data)
										+ wmb2->f(i + 2, j - 1)
										+ INFINITE_ENERGY);							
						}
						//j-1 dangles
						if (j != (number + 1))
						{
							rarray = min(
										rarray, 
										ptr_to_erg4(i, j, j - 1, 2, ct, data, lfce[j - 1]) 
										+ ptr_to_penalty(i, j, ct, data) 
										+ wmb2->f(i + 1, j - 2) 
										+ INFINITE_ENERGY);							
						}
						//both i+1 and j-1 dangle
						if ((i != number) && (j != (number + 1))) 
						{
							rarray = min(
										rarray,
									//	data->tstkm[ct->numseq[i]][ct->numseq[j]][ct->numseq[i + 1]][ct->numseq[j - 1]] 
										multi_tstkm(i, j, i + 1, j - 1, ct)
										+ checknp(lfce[i + 1], lfce[j - 1]) 
										+ wmb2->f(i + 2, j - 2) 
										+ ptr_to_penalty(i, j, ct, data) 
										+ INFINITE_ENERGY);							
						}
					}
				}
			}

			#ifndef INSTRUMENTED//IF INSTRUMENTED IS NOT DEFINED
			v.f(i, j) = vT.f(i, j) = rarray - covar(i, j, ct);
			#else//IF INSTRUMENTED IS DEFINED
			vmb.f(i, j) = rarray;
			v.f(i, j) = vT.f(i, j) = min(rarray, v.f(i, j));//QUESTION: THIS IS DIFFERENT FROM THE MAIN algorithm.cpp (see #ifndef above)
			#endif//END INSTRUMENTED
			
			sub2:

			/*prefill curE[i] and prev[i] for the first two diagonals as ll =4 and 5
			  As d =10, only fill curE (ll=4, d=10)
			  As d =11, only fill prevE (ll=4||5, d=11)
			  As d>11, fill curE(ll=4||5, d>11)
			  exchange curE and prevE after d >11  as curE[h][i][dp]=curE=curE[h-2][i+1][dp]
			  (curE[i][j][ll]=curE[i+1][j-1][ll-2])

				note: dp = jp-ip, where ip and jp are the interior pairs in interior loop
				d = j - i;
			 */
					 //this is part of the ON^3 internal loop calculation(disabled by default)
			if (!simple_iloops) 
			{
				if (((d - 1) >= (minloop + 3) || j > number) && !ct->intermolecular) 
				{
					for (dp = d - 3; dp >= ((j > number) ? 1 : minloop + 1); dp--)
					{
						ll = d - dp - 2;
						//calculate every ip>ip+1,jp<jp-1 when ll ==5 ||4
						if (ll == 4 || ll == 5)
						{
							for (ip = i + 2; ip <= j - 2 - dp; ip++)
							{
								jp = ip + dp;
								if (can_pair(ip, jp))
								{
									if ((ip <= number && jp > number) || (ip <= number && j <= number))
									{
										if (d == ((j > number) ? 7 : 10) || d > ((j > number) ? 8 : 11))
											//fill the first diagonal of d and first two of ll for every larger d 
										{
											curE[dp][i] = min(curE[dp][i], (erg2in(i, j, ip, jp, ct, data, fce.f(i, ip), fce.f(jp, j)) + v.f(ip, jp)));
											//i or j is modified
											if ((mod[ip] || mod[jp]) && can_pair(ip, jp) && !(fce.f(ip, jp) & SINGLE)) 
											{
												curE[dp][i] = min(													
																curE[dp][i],
																erg2in(i, j, ip, jp, ct, data, fce.f(i, ip), fce.f(jp, j))
																+ v.f(ip + 1, jp - 1)
																+ ptr_to_erg1(ip, jp, ip + 1, jp - 1, ct, data));												
											}
										}
										else if (d == ((j > number) ? 8 : 11))  //fill the second diagonal of d
										{
											prevE[dp][i] = min(prevE[dp][i], (erg2in(i, j, ip, jp, ct, data, fce.f(i, ip), fce.f(jp, j)) + v.f(ip, jp)));
											if ((mod[ip] || mod[jp]) && can_pair(ip, jp) && !(fce.f(ip, jp) & SINGLE))
											{
												prevE[dp][i] = min(													
																	prevE[dp][i], 
																	erg2in(i, j, ip, jp, ct, data, fce.f(i, ip), fce.f(jp, j)) 
																	+ v.f(ip + 1, jp - 1) 
																	+ ptr_to_erg1(ip, jp, ip + 1, jp - 1, ct, data));													
											}
										}
									}
								}
							}
						}
						//when size >=6 and <=30;
						//   else if (ll>=6&&ll<=(data->eparam[7]))
						else if (ll >= 6 && ll <= maxinter) 
						{
							//calculate minimum curE[dp][i] of 1 x (n-1) for next step's 2 x (n-2)
							ip = i + 2;
							jp = ip + dp;
							if ((ip <= number && jp > number) || (ip <= number && j <= number)) 
							{
								if (abs(ip - i + jp - j) <= maxasym) 
								{
									if (can_pair(ip, jp)) 
									{
										curE[dp][i] = min(
														curE[dp][i], 
														erg2in(i, j, ip, jp, ct, data, fce.f(i, ip), fce.f(jp, j)) 
														+ v.f(ip, jp));											

										if ((mod[ip] || mod[jp]) && can_pair(ip, jp) && !(fce.f(ip, jp) & SINGLE)) 
										{
											curE[dp][i] = min(
															curE[dp][i], 
															erg2in(i, j, ip, jp, ct, data, fce.f(i, ip), fce.f(jp, j)) 
															+ v.f(ip + 1, jp - 1) 
															+ ptr_to_erg1(ip, jp, ip + 1, jp - 1, ct, data));												
										}
									}
								}
							}

							jp = j - 2;
							ip = jp - dp;
							if ((ip <= number && jp > number) || (ip <= number && j <= number))
							{
								if (abs(ip - i + jp - j) <= maxasym)
								{
									if (can_pair(ip, jp))
									{
										curE[dp][i] = min(	
														curE[dp][i], 
														erg2in(i, j, ip, jp, ct, data, fce.f(i, ip), fce.f(jp, j)) 
														+ v.f(ip, jp));
											
										if ((mod[ip] || mod[jp]) && can_pair(ip, jp) && !(fce.f(ip, jp) & SINGLE))
										{
											curE[dp][i] = min(
															curE[dp][i], 
															erg2in(i, j, ip, jp, ct, data, fce.f(i, ip), fce.f(jp, j)) 
															+ v.f(ip + 1, jp - 1) 
															+ ptr_to_erg1(ip, jp, ip + 1, jp - 1, ct, data));												
										}
									}
								}
							}
						}
					}
				}

				//also block propagation of interior loops that contain nucleotides that need to be double-stranded:
				if ((lfce[i] || lfce[j]) && !ct->intermolecular) for (dp = 1; dp <= d; dp++) curE[dp][i] = INFINITE_ENERGY;//QUESTION: THIS WASN'T IN THE algirithm.napss.cpp
			}


			//Compute w[i][j]: best energy between i and j where i,j does not have
			//	to be a base pair
			//(an exterior loop when it contains n and 1 (ie:n+1)   )


			w.f(i, j) = INFINITE_ENERGY;
		//	cout << i << '\t' << j << '\t' << v.f(i, j) << endl;

			//force a pair between i and j
			if (fce.f(i, j) & PAIR) 
			{
				w.f(i, j) = v.f(i, j) + mbl_branch + ptr_to_penalty(i, j, ct, data);
				goto sub3;
			}

			for (ii = 1; ii <= 5; ii++) e[ii] = INFINITE_ENERGY;

			if (i != number) 
			{
				//calculate the energy of i stacked onto the pair of i+1,j

				e[1] = 
					v.f(i + 1, j) 
					+ mbl_branch 
					+ mbl_nucleotide 
					+ ptr_to_erg4(j, i + 1, i, 2, ct, data, lfce[i]) 
					+ ptr_to_penalty(i + 1, j, ct, data);

				if ((mod[i + 1] || mod[j]) && can_pair(i + 2, j - 1) && !(fce.f(i + 1, j) & SINGLE)) 
				{
					e[1] = min(
							e[1], 
							v.f(i + 2, j - 1) 
							+ mbl_branch 
							+ mbl_nucleotide 
							+ ptr_to_erg4(j, i + 1, i, 2, ct, data, lfce[i]) 
							+ ptr_to_penalty(i + 1, j, ct, data)
							+ ptr_to_erg1(i + 1, j, i + 2, j - 1, ct, data));
				}
				if (!lfce[i]) 
				{
					if (!(fce.f(i, i) & INTER)) 
					{
						//add a nuc to an existing loop:
						e[4] = w.f(i + 1, j) + mbl_nucleotide + ct->SHAPEss_give_value(i);
					}
					//this is for when i represents the center of an intermolecular linker:
					else 
					{ 
						e[4] = w.f(i + 1, j) + mbl_nucleotide + INFINITE_ENERGY; 
					}
				}
			}
			if (j != (number + 1)) 
			{
				//calculate the energy of j stacked onto the pair of i,j-1
				if (j != 1) 
				{
					e[2] = 
						v.f(i, j - 1) 
						+ mbl_branch 
						+ mbl_nucleotide 
						+ ptr_to_erg4(j - 1, i, j, 1, ct, data, lfce[j]) 
						+ ptr_to_penalty(i, j - 1, ct, data);

					if ((mod[i] || mod[j - 1]) && can_pair(i + 1, j - 2) && !(fce.f(i, j - 1) & SINGLE)) 
					{

						e[2] = min(
								e[2], 
								v.f(i + 1, j - 2) 
								+ mbl_branch 
								+ mbl_nucleotide 
								+ ptr_to_erg4(j - 1, i, j, 1, ct, data, lfce[j])
								+ ptr_to_penalty(i, j - 1, ct, data)
								+ ptr_to_erg1(i, j - 1, i + 1, j - 2, ct, data));
					}

					if (!lfce[j]) 
					{
						if (!(fce.f(j, j) & INTER)) 
						{	//add a nuc to an existing loop:
							e[5] = w.f(i, j - 1) + mbl_nucleotide + ct->SHAPEss_give_value(j);
						}
						else e[5] = w.f(i, j - 1) + mbl_nucleotide + INFINITE_ENERGY;
					}
				}
			}
			if ((i != (number)) && (j != ((number)+1))) 
			{
				//calculate i and j stacked onto the pair of i+1,j-1
				if (j != 1 && !lfce[i] && !lfce[j]) 
				{
					e[3] =	v.f(i + 1, j - 1) 
							+ mbl_branch 
							+ 2 * (mbl_nucleotide)
						//	+ data->tstkm[ct->numseq[j - 1]][ct->numseq[i + 1]][ct->numseq[j]][ct->numseq[i]]
							+ multi_tstkm(j - 1, i + 1, j, i, ct)
							+ ptr_to_penalty(j - 1, i + 1, ct, data) 
							+ ct->SHAPEss_give_value(i) 
							+ ct->SHAPEss_give_value(j);

					if ((mod[i + 1] || mod[j - 1]) && (j - 2 > 0) && !(fce.f(i + 1, j - 1) & SINGLE)) 
					{
						if (can_pair(i + 2, j - 2)) 
						{
							e[3] = min(
									e[3], 
									v.f(i + 2, j - 2) 
									+ mbl_branch 
									+ 2 * (mbl_nucleotide)
								//	+ data->tstkm[ct->numseq[j - 1]][ct->numseq[i + 1]][ct->numseq[j]][ct->numseq[i]]
									+ multi_tstkm(j - 1, i + 1, j, i, ct)
									+ ptr_to_penalty(j - 1, i + 1, ct, data) 
									+ ptr_to_erg1(i + 1, j - 1, i + 2, j - 2, ct, data) 
									+ ct->SHAPEss_give_value(i) 
									+ ct->SHAPEss_give_value(j));
						}
					}
				}
			}

			//fragment with i paired to j
			e[1] = min((mbl_branch + v.f(i, j) + ptr_to_penalty(j, i, ct, data)), e[1]);

			if ((mod[i] || mod[j]) && can_pair(i + 1, j - 1) && !(fce.f(i, j) & SINGLE)) 
			{
				e[1] = min
				(
					mbl_branch
					+ v.f(i + 1, j - 1) 
					+ ptr_to_penalty(j, i, ct, data) 
					+ ptr_to_erg1(i, j, i + 1, j - 1, ct, data), 
					e[1]
				);
			}
			w.f(i, j) = min(e[1], e[2]);
			w.f(i, j) = min(w.f(i, j), e[3]);
			w.f(i, j) = min(w.f(i, j), e[4]);
			w.f(i, j) = min(w.f(i, j), e[5]);

			if (ct->intermolecular) 
			{

				// wmb2[i][j%3] = INFINITE_ENERGY;
				// keep track of w2:
				for (ii = 1; ii <= 5; ii++) e[ii] = 2 * INFINITE_ENERGY;
				if (i != number) 
				{
					//calculate the energy of i stacked onto the pair of i + 1, j
					e[1] = 
						v.f(i + 1, j) 
						+ ptr_to_erg4(j, i + 1, i, 2, ct, data, lfce[i]) 
						+ ptr_to_penalty(i + 1, j, ct, data);

					if ((mod[i + 1] || mod[j]) && can_pair(i + 2, j - 1) && !(fce.f(i + 1, j) & SINGLE)) 
					{
						e[1] = min
						(
							e[1], 
							v.f(i + 2, j - 1) 
							+ ptr_to_erg4(j, i + 1, i, 2, ct, data, lfce[i]) 
							+ ptr_to_penalty(i + 1, j, ct, data)
							+ ptr_to_erg1(i + 1, j, i + 2, j - 1, ct, data)
						);
					}

					//if (!lfce[i]) {
					//if (!(fce.f(i,i)&DUBLE))
					e[4] = w2->f(i + 1, j);
					//this is for when i represents the center of an intermolecular linker:
					//else e[4] = w2->f(i+1,j) - INFINITE_ENERGY + data->init;
					//}
				}
				if (j != (number + 1)) 
				{
					//calculate the energy of j stacked onto the pair of i, j - 1
					if (j != 1) 
					{
						e[2] = 
							v.f(i, j - 1) 
							+ ptr_to_erg4(j - 1, i, j, 1, ct, data, lfce[j]) 
							+ ptr_to_penalty(i, j - 1, ct, data);

						if ((mod[i] || mod[j - 1]) && can_pair(i + 1, j - 2) && !(fce.f(i, j - 1) & SINGLE)) 
						{
							e[2] = 
								min(
								e[2], 
								v.f(i + 1, j - 2) 
								+ ptr_to_erg4(j - 1, i, j, 1, ct, data, lfce[j]) 
								+ ptr_to_penalty(i, j - 1, ct, data)
								+ ptr_to_erg1(i, j - 1, i + 1, j - 2, ct, data));
						}
						e[5] = w2->f(i, j - 1);						
					}
				}
				if ((i != number) && (j != (number + 1))) 
				{
					//calculate i and j stacked onto the pair of i+1,j-1
					if (j != 1) 
					{
						e[3] = 
							v.f(i + 1, j - 1) 
						//	+ data->tstkm[ct->numseq[j - 1]][ct->numseq[i + 1]][ct->numseq[j]][ct->numseq[i]]
							+ multi_tstkm(j - 1, i + 1, j, i, ct)
							+ checknp(lfce[i], lfce[j])
							+ ptr_to_penalty(j - 1, i + 1, ct, data);

						if ((mod[i + 1] || mod[j - 1]) && can_pair(i + 2, j - 2) && !(fce.f(i + 1, j - 1) & SINGLE)) 
						{
							e[3] = 
								min
								(
									e[3], 
									v.f(i + 2, j - 2) 
						//			+ data->tstkm[ct->numseq[j - 1]][ct->numseq[i + 1]][ct->numseq[j]][ct->numseq[i]]
									+ multi_tstkm(j - 1, i + 1, j, i, ct)
									+ checknp(lfce[i], lfce[j])
									+ ptr_to_penalty(j - 1, i + 1, ct, data) 
									+ ptr_to_erg1(i + 1, j - 1, i + 2, j - 2, ct, data)
								);
						}
					}
				}

				e[1] = min(e[1], (v.f(i, j) + ptr_to_penalty(j, i, ct, data)));

				if (mod[i] || mod[j] && can_pair(i + 1, j - 1) && !(fce.f(i, j) & SINGLE)) 
				{
					e[1] = 
						min
						(
							e[1],
							v.f(i + 1, j - 1) 
							+ ptr_to_penalty(j, i, ct, data) 
							+ ptr_to_erg1(i, j, i + 1, j - 1, ct, data)	
						);
				}

				w2->f(i, j) = min(e[1], e[2]);
				w2->f(i, j) = min(w2->f(i, j), e[3]);
				w2->f(i, j) = min(w2->f(i, j), e[4]);
				w2->f(i, j) = min(w2->f(i, j), e[5]);

			}

			/*if (((j-i-1)>(2*minloop+2))||(j>(number))) {
			// search for an open bifuraction:
			for (k=i;k<=j-1;k++) {
			if (k==(number)) w.f(i,j)=min(w.f(i,j),
			w3[i]+w5[j-(number)]);
			else w.f(i,j) = min(w.f(i,j),
			w.f(i,k)+work[k+1][j%3]);
			}
			}  */
			////fill wmb:

			if (((j - i - 1) > (2 * minloop + 2)) || j > number) 
			{

				//search for an open bifurcation:
				//the underlying w arrays are accessed directly
				//instead of through the f function
				//to avoid the overhead of the conditions in f
				int end = min(number, j);
				int tmp = wmb.f(i, j);
				for (k = i; k < end; ++k) 
				{
					tmp = min(tmp, w.dg[i][k] + wT.dg[j][k + 1]);
				}
				for (k = number; k < j; ++k) 
				{
					tmp = min(tmp, w.dg[i][k] + wT.dg[j - number][k + 1 - number]);
				}
			
				wmb.f(i, j) = tmp;

				//for the sake of coaxial stacking, also consider the addition of nucs
				//to a previously calculated wmb
				if (i != number)
					if (!lfce[i]) wmb.f(i, j) = min(wmb.f(i, j), wmb.f(i + 1, j) + mbl_nucleotide + ct->SHAPEss_give_value(i));
				if (j != number + 1)
					if (!lfce[j]) wmb.f(i, j) = min(wmb.f(i, j), wmb.f(i, j - 1) + mbl_nucleotide + ct->SHAPEss_give_value(j));

				e[1] = INFINITE_ENERGY;
				e[2] = INFINITE_ENERGY;
				//also consider the coaxial stacking of two helixes

				#ifndef disablecoax
				if (!disablecoax) 
				{
					for (ip = i + minloop + 1; ip < j - minloop - 1; ip++) 
					{
						//first consider flush stacking
						if (ip != number) 
						{
							if (can_pair(i, ip) && can_pair(j, ip + 1)) 
							{
								//ony proceed if i and ip and ip+1 and j can pair
								e[1] = 
									min
									(
										e[1], 
										v.f(i, ip) 
										+ vT.f(ip + 1, j) 
										+ ptr_to_penalty(i, ip, ct, data)
										+ ptr_to_penalty(ip + 1, j, ct, data) 
										+ ptr_to_ergcoaxflushbases(i, ip, ip + 1, j, ct, data)
									);

								if (mod[i] || mod[ip] || mod[ip + 1] || mod[j]) 
								{
									if ((mod[i] || mod[ip]) && (mod[ip + 1] || mod[j]) && can_pair(i + 1, ip - 1)
										&& can_pair(ip + 2, j - 1) && !(fce.f(i, ip) & SINGLE) && !(fce.f(ip + 1, j) & SINGLE)) 
									{
										e[1] = min
										(
											e[1], 
											v.f(i + 1, ip - 1) 
											+ vT.f(ip + 2, j - 1) 
											+ ptr_to_penalty(i, ip, ct, data)
											+ ptr_to_penalty(ip + 1, j, ct, data) 
											+ ptr_to_ergcoaxflushbases(i, ip, ip + 1, j, ct, data)
											+ ptr_to_erg1(i, ip, i + 1, ip - 1, ct, data) 
											+ ptr_to_erg1(ip + 1, j, ip + 2, j - 1, ct, data)
										);
									}
									if ((mod[i] || mod[ip]) && can_pair(i + 1, ip - 1) && !(fce.f(i, ip) & SINGLE)) 
									{
										e[1] = min
										(
											e[1], 
											v.f(i + 1, ip - 1) + vT.f(ip + 1, j) 
											+ ptr_to_penalty(i, ip, ct, data)
											+ ptr_to_penalty(ip + 1, j, ct, data) 
											+ ptr_to_ergcoaxflushbases(i, ip, ip + 1, j, ct, data)
											+ ptr_to_erg1(i, ip, i + 1, ip - 1, ct, data)
										);
									}
									if ((mod[ip + 1] || mod[j]) && can_pair(ip + 2, j - 1) && !(fce.f(ip + 1, j) & SINGLE)) 
									{
										e[1] = min
										(
											e[1], 
											v.f(i, ip) 
											+ vT.f(ip + 2, j - 1) 
											+ ptr_to_penalty(i, ip, ct, data)
											+ ptr_to_penalty(ip + 1, j, ct, data) 
											+ ptr_to_ergcoaxflushbases(i, ip, ip + 1, j, ct, data)
											+ ptr_to_erg1(ip + 1, j, ip + 2, j - 1, ct, data)
										);
									}
								}
							}//end (can_pair(i,ip)&&can_pair(j,ip+1)) 
							if (ip + 1 != number) 
							{
								if (can_pair(i, ip) && can_pair(j - 1, ip + 2)) 
								{
									//only proceed if i and ip and j-1 and ip+2 can pair
									if (!lfce[ip + 1] && !lfce[j]) 
									{
										//now consider an intervening mismatch

										e[2] = 
											min
											(
												e[2], 
												v.f(i, ip) 
												+ vT.f(ip + 2, j - 1) 
												+ ptr_to_penalty(i, ip, ct, data)
												+ ptr_to_penalty(ip + 2, j - 1, ct, data) 
												+ ptr_to_ergcoaxinterbases2(i, ip, ip + 2, j - 1, ct, data)
											);

										if (mod[i] || mod[ip] || mod[ip + 2] || mod[j - 1]) 
										{
											if ((mod[i] || mod[ip]) && (mod[ip + 2] || mod[j - 1]) && can_pair(i + 1, ip - 1)
												&& can_pair(ip + 3, j - 2) && !(fce.f(i, ip) & SINGLE) && !(fce.f(ip + 2, j - 1) & SINGLE)) 
											{
												if (!lfce[ip + 1] && !lfce[j]) 
												{
													e[2] =
														min
														(
															e[2],
															v.f(i + 1, ip - 1)
															+ vT.f(ip + 3, j - 2)
															+ ptr_to_penalty(i, ip, ct, data)
															+ ptr_to_penalty(ip + 2, j - 1, ct, data)
															+ ptr_to_ergcoaxinterbases2(i, ip, ip + 2, j - 1, ct, data)
															+ ptr_to_erg1(i, ip, i + 1, ip - 1, ct, data)
															+ ptr_to_erg1(ip + 2, j - 1, ip + 3, j - 2, ct, data)
														);
												}
											}
											if ((mod[i] || mod[ip]) && can_pair(i + 1, ip - 1) && !(fce.f(i, ip) & SINGLE)) 
											{
												if (!lfce[ip + 1] && !lfce[j]) 
												{
													e[2] = 
														min
														(
															e[2], 
															v.f(i + 1, ip - 1) 
															+ vT.f(ip + 2, j - 1) 
															+ ptr_to_penalty(i, ip, ct, data)
															+ ptr_to_penalty(ip + 2, j - 1, ct, data) 
															+ ptr_to_ergcoaxinterbases2(i, ip, ip + 2, j - 1, ct, data)
															+ ptr_to_erg1(i, ip, i + 1, ip - 1, ct, data)
														);
												}
											}
											if ((mod[ip + 2] || mod[j - 1]) && can_pair(ip + 3, j - 2) && !(fce.f(ip + 2, j - 1) & SINGLE))
											{
												if (!lfce[ip + 1] && !lfce[j]) 
												{
													e[2] =
														min
														(
															e[2],
															v.f(i, ip)
															+ vT.f(ip + 3, j - 2)
															+ ptr_to_penalty(i, ip, ct, data)
															+ ptr_to_penalty(ip + 2, j - 1, ct, data)
															+ ptr_to_ergcoaxinterbases2(i, ip, ip + 2, j - 1, ct, data)
															+ ptr_to_erg1(ip + 2, j - 1, ip + 3, j - 2, ct, data)
														);
												}
											}
										}
									} // end (!lfce[ip+1]&&!lfce[j])
								} // end (can_pair(i,ip)&&can_pair(j-1,ip+2))

								if (can_pair(i + 1, ip) && can_pair(j, ip + 2)) 
								{
									//only proceed if i and ip j and ip+2 can pair
									if (!lfce[i] && !lfce[ip + 1] && i != number) 
									{
										e[2] = 
											min
											(
												e[2], 
												v.f(i + 1, ip) 
												+ vT.f(ip + 2, j) 
												+ ptr_to_penalty(i + 1, ip, ct, data)
												+ ptr_to_penalty(ip + 2, j, ct, data) 
												+ ptr_to_ergcoaxinterbases1(i + 1, ip, ip + 2, j, ct, data)
											);

										if (mod[i + 1] || mod[ip] || mod[ip + 2] || mod[j]) 
										{
											if ((mod[i + 1] || mod[ip]) && (mod[ip + 2] || mod[j]) && can_pair(i + 2, ip - 1)
												&& can_pair(ip + 3, j - 1) && !(fce.f(i + 1, ip) & SINGLE) && !(fce.f(ip + 2, j) & SINGLE)) 
											{
												e[2] = 
													min
													(
														e[2], 
														v.f(i + 2, ip - 1) 
														+ vT.f(ip + 3, j - 1) 
														+ ptr_to_penalty(i + 1, ip, ct, data)
														+ ptr_to_penalty(ip + 2, j, ct, data)
														+ ptr_to_ergcoaxinterbases1(i + 1, ip, ip + 2, j, ct, data)
														+ ptr_to_erg1(i + 1, ip, i + 2, ip - 1, ct, data) 
														+ ptr_to_erg1(ip + 2, j, ip + 3, j - 1, ct, data)
													);
											}
											if ((mod[i + 1] || mod[ip]) && can_pair(i + 2, ip - 1) && !(fce.f(i + 1, ip) & SINGLE)) 
											{
												e[2] = 
													min
													(
														e[2], 
														v.f(i + 2, ip - 1) + vT.f(ip + 2, j) 
														+ ptr_to_penalty(i + 1, ip, ct, data)
														+ ptr_to_penalty(ip + 2, j, ct, data) 
														+ ptr_to_ergcoaxinterbases1(i + 1, ip, ip + 2, j, ct, data)
														+ ptr_to_erg1(i + 1, ip, i + 2, ip - 1, ct, data)
													);
											}
											if ((mod[ip + 2] || mod[j]) && can_pair(ip + 3, j - 1) && !(fce.f(ip + 2, j) & SINGLE)) 
											{
												e[2] = 
													min
													(
														e[2], 
														v.f(i + 1, ip) 
														+ vT.f(ip + 3, j - 1) 
														+ ptr_to_penalty(i + 1, ip, ct, data)
														+ ptr_to_penalty(ip + 2, j, ct, data) 
														+ ptr_to_ergcoaxinterbases1(i + 1, ip, ip + 2, j, ct, data)
														+ ptr_to_erg1(ip + 2, j, ip + 3, j - 1, ct, data)
													);
											}
										}
									}
								}//end (can_pair(i+1,ip)&&can_pair(j,ip+2)) 
							}
						}
					}
				}
				#endif //ifndef disablecoax

				wmb.f(i, j) = min(wmb.f(i, j), e[1] + 2 * mbl_branch);
				wmb.f(i, j) = min(wmb.f(i, j), e[2] + 2 * mbl_branch + 2 * mbl_nucleotide);

				if (j <= number) wca[i][j] = min(e[1], e[2]);

				w.f(i, j) = min(w.f(i, j), wmb.f(i, j));

				if (ct->intermolecular) 
				{
					//intermolecular folding:
					//search for an open bifurcation:
					for (k = i; k <= j; k++) 
					{
						if (k != number) wmb2->f(i, j) = min(wmb2->f(i, j), w2->f(i, k) + w2->f(k + 1, j));
					}
					if (i != number) 
					{
						if (!(fce.f(i, i) & INTER)) wmb2->f(i, j) = min(wmb2->f(i, j), wmb2->f(i + 1, j));
						else  wmb2->f(i, j) = min(wmb2->f(i, j), wmb2->f(i + 1, j) + data->init - INFINITE_ENERGY);
					}
					if (j != number + 1) 
					{
						if (!(fce.f(j, j) & INTER)) wmb2->f(i, j) = min(wmb2->f(i, j), wmb2->f(i, j - 1));
						else wmb2->f(i, j) = min(wmb2->f(i, j), wmb2->f(i, j - 1) + data->init - INFINITE_ENERGY);
					}
					w2->f(i, j) = min(w2->f(i, j), wmb2->f(i, j));
				}
			}
			sub3:
				wT.f(i, j) = w.f(i, j);

			// Calculate vmin, the best energy for the entire sequence
			if (j > (number)) 
			{
				vmin = min(vmin, v.f(i, j) + v.f(j - (number), i));
			}
	
			// compute w5[i], the energy of the best folding from 1->i, and
			// w3[i], the energy of the best folding from i-->GetSequenceLength()
			if (i == 1 && j <= number)
			{
				if (j <= minloop + 1) 
				{
					if (lfce[j]) w5[j] = INFINITE_ENERGY;
					else  w5[j] = w5[j - 1] + ct->SHAPEss_give_value(j);
				}

				else 
				{
					if (lfce[j]) w5[j] = INFINITE_ENERGY;
					else w5[j] = w5[j - 1] + ct->SHAPEss_give_value(j);

					for (k = 1; k <= 5; k++) e[k] = INFINITE_ENERGY;//e[k]=0;
					
					#ifndef disablecoax
					castack = INFINITE_ENERGY;
					#endif //ifndef disablecoax

					for (k = 0; k <= (j - 4); k++) 
					{
						e[1] = min(e[1], (w5[k] + v.f(k + 1, j) + ptr_to_penalty(j, k + 1, ct, data)));

						if (mod[k + 1] || mod[j] && can_pair(k + 2, j - 1) && !(fce.f(k + 1, j) & SINGLE)) 
						{
							e[1] = min(
									e[1], 
									(w5[k] 
									+ v.f(k + 2, j - 1) 
									+ ptr_to_penalty(j, k + 1, ct, data)
									+ ptr_to_erg1(k + 1, j, k + 2, j - 1, ct, data)));								
						}

						e[2] = min(
								e[2], 
								(w5[k] 
								+ ptr_to_erg4(j, k + 2, k + 1, 2, ct, data, lfce[k + 1]) 
								+ v.f(k + 2, j) 
								+ ptr_to_penalty(j, k + 2, ct, data)));
							
						if (mod[k + 2] || mod[j] && can_pair(k + 3, j - 1) && !(fce.f(k + 2, j) & SINGLE)) 
						{
							e[2] = min(
									e[2], 
									(w5[k] 
									+ ptr_to_erg4(j, k + 2, k + 1, 2, ct, data, lfce[k + 1]) 
									+ v.f(k + 3, j - 1)
									+ ptr_to_penalty(j, k + 2, ct, data) 
									+ ptr_to_erg1(k + 2, j, k + 3, j - 1, ct, data)));								
						}

						e[3] = min(														
								e[3], 
								w5[k] 
								+ ptr_to_erg4(j - 1, k + 1, j, 1, ct, data, lfce[j]) 
								+ v.f(k + 1, j - 1) 
								+ ptr_to_penalty(j - 1, k + 1, ct, data));
							
						if (mod[k + 1] || mod[j - 1] && can_pair(k + 2, j - 2) && !(fce.f(k + 1, j - 1) & SINGLE)) 
						{
							e[3] = min(
									e[3], 
									w5[k] 
									+ ptr_to_erg4(j - 1, k + 1, j, 1, ct, data, lfce[j]) 
									+ v.f(k + 2, j - 2)
									+ ptr_to_penalty(j - 1, k + 1, ct, data) 
									+ ptr_to_erg1(k + 1, j - 1, k + 2, j - 2, ct, data));								
						}

						e[4] = min(														
								e[4], 
								w5[k] 
						//		+ data->tstack[ct->numseq[j - 1]][ct->numseq[k + 2]][ct->numseq[j]][ct->numseq[k + 1]]
								+ multi_tstack(j - 1, k + 2, j, k + 1, ct)
								+ checknp(lfce[j], lfce[k + 1]) 
								+ v.f(k + 2, j - 1) 
								+ ptr_to_penalty(j - 1, k + 2, ct, data) 
								+ ct->SHAPEss_give_value(j) 
								+ ct->SHAPEss_give_value(k + 1));
							

						if (mod[k + 2] || mod[j - 1] && can_pair(k + 3, j - 2) && !(fce.f(k + 2, j - 1) & SINGLE)) 
						{
							e[4] = min(
									e[4], 
									w5[k] 
								//	+ data->tstack[ct->numseq[j - 1]][ct->numseq[k + 2]][ct->numseq[j]][ct->numseq[k + 1]]
									+ multi_tstack(j - 1, k + 2, j, k + 1, ct)
									+ checknp(lfce[j], lfce[k + 1]) 
									+ v.f(k + 3, j - 2) 
									+ ptr_to_penalty(j - 1, k + 2, ct, data) 
									+ ptr_to_erg1(k + 2, j - 1, k + 3, j - 2, ct, data)
									+ ct->SHAPEss_give_value(j) 
									+ ct->SHAPEss_give_value(k + 1));								
						}

						#ifndef disablecoax
						castack = min(castack, w5[k] + wca[k + 1][j]);
						#endif //ifndef disablecoax
					}

					w5[j] = min(w5[j], e[1]);
					w5[j] = min(w5[j], e[2]);
					w5[j] = min(w5[j], e[3]);
					w5[j] = min(w5[j], e[4]);
					
					#ifndef disablecoax
					w5[j] = min(w5[j], castack);
					#endif //ifndef disablecoax
				}
			}
			if (j == number)
			{
				w3[0] = 0;
				w3[number + 1] = 0;
				if (i >= (number - minloop)) 
				{    //number+1 ... number-minloop
					if (lfce[i]) w3[i] = INFINITE_ENERGY;
					else w3[i] = w3[i + 1] + ct->SHAPEss_give_value(i);
				}
				//w3[i]=0;
				if (i >= 1 && i <= (number - minloop - 1)) 
				{

					if (lfce[i]) w3[i] = INFINITE_ENERGY;
					else w3[i] = w3[i + 1] + ct->SHAPEss_give_value(i);

					for (k = 1; k <= 5; k++) e[k] = INFINITE_ENERGY;
					#ifndef disablecoax
					castack = INFINITE_ENERGY;
					#endif //ifndef disablecoax
					for (k = (number + 1); k >= (i + 4); k--) 
					{
						e[1] = min(e[1], (v.f(i, k - 1) + w3[k] + ptr_to_penalty(k - 1, i, ct, data)));

						if ((mod[i] || mod[k - 1]) && can_pair(i + 1, k - 2) && !(fce.f(i, k - 1) & SINGLE)) 
						{
							e[1] = min(
									e[1], 
									(v.f(i + 1, k - 2) 
									+ w3[k] 
									+ ptr_to_penalty(k - 1, i, ct, data) 
									+ ptr_to_erg1(i, k - 1, i + 1, k - 2, ct, data)));
						}

						e[2] = min(
								e[2], 
								(v.f(i + 1, k - 1) 
								+ ptr_to_erg4(k - 1, i + 1, i, 2, ct, data, lfce[i]) 
								+ ptr_to_penalty(k - 1, i + 1, ct, data) 
								+ w3[k]));

						if ((mod[i + 1] || mod[k - 1]) && can_pair(i + 2, k - 2) && !(fce.f(i + 1, k - 1) & SINGLE)) 
						{
							e[2] = min(
									e[2], 
									(v.f(i + 2, k - 2) 
									+ ptr_to_erg4(k - 1, i + 1, i, 2, ct, data, lfce[i]) 
									+ ptr_to_penalty(k - 1, i + 1, ct, data) 
									+ w3[k] 
									+ ptr_to_erg1(i + 1, k - 1, i + 2, k - 2, ct, data)));
						}

						e[3] = min(
								e[3], 
								(v.f(i, k - 2) 
								+ ptr_to_erg4(k - 2, i, k - 1, 1, ct, data, lfce[k - 1]) 
								+ ptr_to_penalty(k - 2, i, ct, data) + w3[k]));

						if ((mod[i] || mod[k - 2]) && can_pair(i + 1, k - 3) && !(fce.f(i, k - 2) & SINGLE)) 
						{
							e[3] = min(
									e[3], 
									(v.f(i + 1, k - 3) 
									+ ptr_to_erg4(k - 2, i, k - 1, 1, ct, data, lfce[k - 1]) 
									+ ptr_to_penalty(k - 2, i, ct, data) 
									+ w3[k] 
									+ ptr_to_erg1(i, k - 2, i + 1, k - 3, ct, data)));
						}

						if (!lfce[i] && !lfce[k - 1]) 
						{
							e[4] = min(
									e[4], 
									(v.f(i + 1, k - 2) 
								//	+ data->tstack[ct->numseq[k - 2]][ct->numseq[i + 1]][ct->numseq[k - 1]][ct->numseq[i]]
									+ multi_tstack(k - 2, i + 1, k - 1, i, ct)
									+ checknp(lfce[k - 1], lfce[i]) 
									+ w3[k] 
									+ ptr_to_penalty(k - 2, i + 1, ct, data)) 
									+ ct->SHAPEss_give_value(i) 
									+ ct->SHAPEss_give_value(k - 1));

							if ((mod[i + 1] || mod[k - 2]) && can_pair(i + 2, k - 3) && !(fce.f(i + 1, k - 2) & SINGLE)) 
							{
								e[4] = min(
										e[4], 
										(v.f(i + 2, k - 3) 
									//	+ data->tstack[ct->numseq[k - 2]][ct->numseq[i + 1]][ct->numseq[k - 1]][ct->numseq[i]]
										+ multi_tstack(k - 2, i + 1, k - 1, i, ct)
										+ checknp(lfce[k - 1], lfce[i]) 
										+ w3[k] 
										+ ptr_to_penalty(k - 2, i + 1, ct, data) 
										+ ptr_to_erg1(i + 1, k - 2, i + 2, k - 3, ct, data))
										+ ct->SHAPEss_give_value(i) 
										+ ct->SHAPEss_give_value(k - 1));
							}
						}

						//also consider coaxial stacking:
						#ifndef disablecoax
						castack = min(wca[i][k - 1] + w3[k], castack);
						#endif //ifndef disablecoax
					}

					w3[i] = min(w3[i], e[1]);
					w3[i] = min(w3[i], e[2]);
					w3[i] = min(w3[i], e[3]);
					w3[i] = min(w3[i], e[4]);
					
					#ifndef disablecoax
					w3[i] = min(w3[i], castack);
					#endif //ifndef disablecoax
				}
			}
			#ifdef DYNALIGN_II

			if (j <= number) {
				for (size_t it = 1; it <= 5; it++)e[it] = 0;

				castack = 0;

				//	we.f(i,j)=INFINITE_ENERGY;
				if (j - i + 1 <= minloop + 1) {
					if (lfce[j] || lfce[i])we->f(i, j) = INFINITE_ENERGY;
					else we->f(i, j) = min(we->f(i, j - 1), we->f(i + 1, j));

				}

				else {
					if (lfce[j] || lfce[i])we->f(i, j) = INFINITE_ENERGY;
					else we->f(i, j) = min(we->f(i, j - 1), we->f(i + 1, j));

					for (k = 1; k <= 5; k++)e[k] = INFINITE_ENERGY;

					castack = INFINITE_ENERGY;

					e[1] = min(e[1], v.f(i, j) + penalty(j, i, ct, data));

					if (mod[i] || mod[j] && can_pair(i + 1, j - 1) && !(fce.f(i, j) & SINGLE)) {
						e[1] = min(e[1], (v.f(i + 1, j - 1) + penalty(j, i, ct, data)
							+ erg1(i, j, i + 1, j - 1, ct, data)));
					}

					e[2] = min(e[2], (v.f(i + 1, j) + penalty(j, i + 1, ct, data) + erg4(j, i + 1, i, 2, ct, data, lfce[i])));//yinghan:erg4, how to set up the fourth parameter? sep.27

					if (mod[i + 1] || mod[j] && can_pair(i + 2, j - 1) && !(fce.f(i + 1, j) & SINGLE)) {
						e[2] = min(e[2], (v.f(i + 2, j - 1) + penalty(j, i + 1, ct, data) + erg4(j, i + 1, i, 2, ct, data, lfce[i])
							+ erg1(i + 1, j, i + 2, j - 1, ct, data)));
					}//yinghan:shoudln't they check mod first? sep.27

					e[3] = min(e[3], (v.f(i, j - 1) + penalty(j - 1, i, ct, data) + erg4(j - 1, i, j, 1, ct, data, lfce[j])));

					if (mod[i] || mod[j - 1] && can_pair(i + 1, j - 2) && !(fce.f(i, j - 1) & SINGLE)) {
						e[3] = min(e[3], (v.f(i + 1, j - 2) + penalty(j - 1, i, ct, data) + erg4(j - 1, i, j, 1, ct, data, lfce[j])
							+ erg1(i, j - 1, i + 1, j - 2, ct, data)));
					}

					e[4] = min(e[4], (data->tstack[ct->numseq[j - 1]][ct->numseq[i + 1]][ct->numseq[j]][ct->numseq[i]]
						+ checknp(lfce[j], lfce[i]) + v.f(i + 1, j - 1) + penalty(j - 1, i + 1, ct, data)));

					if (mod[i + 1] || mod[j - 1] && can_pair(i + 2, j - 2) && !(fce.f(i + 1, j - 1) & SINGLE)) {
						e[4] = min(e[4], (data->tstack[ct->numseq[j - 1]][ct->numseq[i + 1]][ct->numseq[j]][ct->numseq[i]]
							+ checknp(lfce[j], lfce[i]) + v.f(i + 2, j - 2) + erg1(i + 1, j - 1, i + 2, j - 2, ct, data)
							+ penalty(j - 1, i + 1, ct, data)));
					}


					for (k = i; k <= j; k++)e[5] = min(e[5], we->f(i, k) + we->f(k + 1, j));


					castack = min(castack, wca[i][j]);
				}

				we->f(i, j) = min(we->f(i, j), e[1]);
				we->f(i, j) = min(we->f(i, j), e[2]);
				we->f(i, j) = min(we->f(i, j), e[3]);
				we->f(i, j) = min(we->f(i, j), e[4]);
				we->f(i, j) = min(we->f(i, j), e[5]);
				we->f(i, j) = min(we->f(i, j), castack);

			}
			#else
			#endif

			#undef can_pair
		}

		if (!simple_iloops) 
		{
			if (!ct->intermolecular)
			{
				if (d > (maxj > number ? 8 : 11))
				{
					tempE = curE;
					curE = prevE;
					prevE = tempE;
				}
				if (d > (maxj > number ? 7 : 10)) 
				{
					for (int localdp = 1; localdp <= d - 1; localdp++)
					{
						for (int locali = ((h <= (number - 2)) ? 1 : (2 * number - h - 1)); locali <= ((h <= (number - 2)) ? (number - h - 1) : number); locali++)
						{
							if (locali < number) curE[localdp][locali] = curE[localdp][locali + 1];
						}
					}
				}
			}
		}
	}

	// SOMETHING WRONG
	// arraydump(v, w, wmb, w5, w3, ct->GetSequenceLength(), "3.out");

	// output V, W, WMB, and W2V:
	#if defined (debugmode)
	#ifndef INSTRUMENTED
	arraydump(v, w, wmb, w5, w3, ct->GetSequenceLength(), "arrays.out");
	#endif
	#endif

	// clean up memory use:

	for (int locali = 0; locali <= number; locali++)
	{
		delete[] wca[locali];
	}
	delete[] wca;

	if (!ct->intermolecular) 
	{
		for (int locali = 0; locali <= number; locali++) 
		{
			delete[] curE[locali];
			delete[] prevE[locali];
		}
		delete[] curE;
		delete[] prevE;
	}

} //end fill

#ifndef INSTRUMENTED//IF INSTRUMENTED IS NOT DEFINED

void filter(	structure* ct, 
				int percent, 
				int max, 
				int window	) 
{
	//structure temp;
	short int i, j, k1, k2, crit;
	bool **mark;
	//bool keep;

	mark = new bool *[ct->GetSequenceLength() + 1];

	for (i = 0; i <= ct->GetSequenceLength(); i++) 
	{
		mark[i] = new bool [ct->GetSequenceLength() + 1];
	}

	for (i = 1; i <= ct->GetSequenceLength(); i++) 
	{
		for (j = i; j <= ct->GetSequenceLength(); j++) 
		{
			mark[i][j] = false;
		}
	}

	crit = (short) (ct->GetEnergy(1) + abs((int)((float)ct->GetEnergy(1) * ((((float)percent)/100.0)))));
	// temp.numofstructures = 0;

	// check each structure:
	for (i = 1; i <= ct->GetNumberofStructures(); ++i) 
	{
		if (ct->GetEnergy(i) > crit) 
		{
			// none of the remaining structures should be kept bcs the free
			for (j = ct->GetNumberofStructures(); j >= i; j--) ct->RemoveLastStructure();
			//clean up memory:
			de_allocate(mark, ct->GetSequenceLength()+1);
			return; 
			// energy is > than % sort
		}
		else if (i > max) 
		{
			//none of the remaining structures should be kept bcs the max #
			//of structures has been reached
			for (j = ct->GetNumberofStructures(); j >= i; j--) ct->RemoveLastStructure();
			de_allocate(mark, ct->GetSequenceLength() + 1);
			return;
		}

		//now map the baspairs within size < window and verify whether this structure
		//should be kept or discarded
		//keep = false;
		int newpairs=0;

		for (j = 1; j <= ct->GetSequenceLength(); j++) 
		{
			if (ct->GetPair(j, i) > j) 
			{
				if (!mark[j][ct->GetPair(j, i)]) 
				{
					//this base has not been marked so keep the structure:
					++newpairs;
				}
			}
		}

		for (j = 1; j <= ct->GetSequenceLength(); j++) 
		{
			if (ct->GetPair(j, i) > j) 
			{
				//now mark the basepairs:
				for (k1 = j - window; k1 <= j + window; k1++) 
				{
					for (k2 = ct->GetPair(j, i) - window; k2 <= ct->GetPair(j, i) + window; k2++) 
					{
						if ((k1 > 0) && (k2 > 0) && (k1 <= ct->GetSequenceLength()) && (k2 <= ct->GetSequenceLength())) 
						{
							mark[k1][k2] = true;
						}
					}
				}
			}
		}

		if (newpairs <= window) 
		{   //structure needs to be kept, copy it over to temp
			//Get rid of this structure:
			ct->RemoveStructure(i);
			--i;//decrement the counter so that the new ith structure is not missed
		}
	}
	//clean up memory use:
	de_allocate(mark, ct->GetSequenceLength() + 1);
}

void cctout(	structure *ct, 
				char *filename	) 
{
	int i, j;
	ofstream out(filename);
       
	out << "-100\n";
	out << ct->GetSequenceLength() << "\n";
	out << ct->GetNumberofStructures() << " ";
	out << ct->GetCtLabel(1).c_str();

	for (i = 1; i <= ct->GetSequenceLength(); i++) 
	{
		out << ct->numseq[i] << "\n";
	}

	for (i = 1; i <= ct->GetNumberofStructures(); i++) 
	{
		out << ct->GetEnergy(1) << "\n";
		for (j = 1; j <= ct->GetSequenceLength(); j++) 
		{
			out << ct->GetPair(j, i) << "\n";
		}
	}
}

void calcpnum(	dotarray *dots, 
				int *pnum, 
				int increment, 
				short int numofbases,
				ProgressHandler *PD	) 
{
	short int i, j;
	for (i = 1; i <= numofbases; i++) 
	{
		pnum[i] = 0;
		//count the dots in the ith column
		for (j = i + 1; j <= numofbases; j++) 
		{
			if (dots->dot(i, j) <= increment) pnum[i]++;
		}
		//count the dots in the ith row
		for (j = 1; j < i; j++) 
		{
			if (dots->dot(j, i) <= increment) pnum[i]++;
		}
	}
}

//This function will give an energy breakdown for a structure, int n,
//stored in ct, it takes the arrays below as input
void energydump(structure *ct, 
				DynProgArray<integersize> *v, 
				datatable *data, 
				int n,
				char *filename,
				int ii, 
				int ji	) 
{
	short int stack[500], stackpos, i, j, k, temp, count, helix, stacke[500], stackpose, start;
	ofstream out;
	char number[15];
	bool exterior;

	out.open(filename);
	stackpos = 0;

	sprintf(number, "%f", (float (ct->GetEnergy(n)))/conversionfactor);
	out << "Structure:  "<< n <<"\n";
	out << "\n# " << n << "  Total Energy = "<< number << "\n\n";
	out << "pair of origin: "<< ii <<" - "<< ji << "\n\n";
   
	stackpos = 0;

	//Analyze the exterior loop
	i = 0;
	temp = 0;
	while (i < ct->GetSequenceLength()) 
	{
		i++;
		if (ct->GetPair(i, n) > 0) 
		{
			if (i <= ii && ct->GetPair(i, n) >= ji) 
			{
				temp = temp + v->f(ct->GetPair(i, n), i + ct->GetSequenceLength());  
			}
			else temp = temp-v->f(i, ct->GetPair(i, n));
			i = ct->GetPair(i, n);
		}
	}

	sprintf(number, "%f", (float (temp))/conversionfactor);
	out << "Exterior loop energy = "<<number<<"\n";

	//Place the forced pair on the stack
	stackpose = 1;
	stacke[1] = ii;

	//Trace out the exterior fragment
	while (stackpose > 0) 
	{
		i = stacke[stackpose];
		stackpose--;
		helix = 0;

		//follow the helix:
		while (ct->GetPair(i - 1, n) == ct->GetPair(i, n) + 1) 
		{
			//helix continues:
			temp = erg1(ct->GetPair(i, n), i, ct->GetPair(i, n) + 1, i - 1, ct, data);
			//note: v->f(i,ct->basepr[n][i])-v->f(i+1,ct->basepr[n][i+1]);
			//does not work here for cases when a loop is open in min free energy structure
			//but closed in the sub optimal structure being traced
			sprintf(number, "%f", (float (temp))/conversionfactor);
			helix = helix + temp;
			out << "\tStack energy = " << number << "  for " << ct->GetPair(i, n) <<"-"
			    << i << " onto " << ct->GetPair(i, n) + 1 << "-" << i-1 << "\n";
			i--;
		}
		sprintf(number, "%f", (float (helix))/conversionfactor);
		out << "Helix energy = " << number << "\n";

		//now we've come to a loop, what type?
		j = ct->GetPair(i, n);
		start = i;
		exterior = false;
		//check for exterior loop
      
		temp = v->f(j, i+ct->GetSequenceLength());
		
		while (i != j+1 && i != j) 
		{
			i--;
			if (i==0) 
			{
				exterior = true;
				i = ct->GetSequenceLength() + 1; 
			}
			else if (ct->GetPair(i, n) > 0) 
			{
				//we've found another helix
				k = ct->GetPair(i, n);      

				if (i <= start && k < i) 
				{
					stackpos++;
					stack[stackpos] = ct->GetPair(i, n);
					temp = temp - v->f(k, i);
				}
				else if (i > start) 
				{
					stackpos++;
					stack[stackpos] = ct->GetPair(i, n);
					temp = temp - v->f(k, i);
				} 
				else 
				{
					stackpose++;
					stacke[stackpose] = i;
					temp = temp - v->f(k, i + ct->GetSequenceLength());
				}

				i = k;
			}
		}

		if (!exterior) 
		{
			//multi loop
			sprintf(number, "%f", (float (temp))/conversionfactor);
			out << "Multibranch loop energy = " << number << "  for closure by " <<
				ct->GetPair(j, n) << "-" << j << "\n";
		}
	}
   
	//Place the forced pair on the stack
	stackpos++;
	stack[stackpos] = ii;

	//Trace out the interior fragment
	while (stackpos > 0) 
	{
		i = stack[stackpos];
		stackpos--;

		helix = 0;
		//follow the helix:
		while (ct->GetPair(i + 1, n) == ct->GetPair(i, n) - 1) 
		{
			//helix continues:
			temp = erg1(i, ct->GetPair(i, n), i + 1, ct->GetPair(i + 1, n), ct, data);
			//note: v->f(i,ct->basepr[n][i])-v->f(i+1,ct->basepr[n][i+1]);
			//does not work here for cases when a loop is open in min free energy structure
			//but closed in the sub optimal structure being traced
			sprintf(number, "%f", (float (temp))/conversionfactor);
			helix = helix + temp;
			out << "\tStack energy = "<< number<<"  for "<< (i + 1)<< "-"
			    << ct->GetPair(i + 1, n) <<" onto "<< i << "-" << ct->GetPair(i, n) << "\n";
			i++;
		}

		sprintf(number, "%f", (float (helix))/conversionfactor);
		out << "Helix energy = "<< number<<"\n";

		//now we've come to a loop, what type?
		j = ct->GetPair(i, n);

		temp = v->f(i, j);
		count = 0;
		
		while (i < j-1) 
		{
			i++;
			if (ct->GetPair(i, n)>0) 
			{
				//we've found another helix
				count++;
				k = ct->GetPair(i, n);
				temp = temp - v->f(i, k);
				stackpos++;
				stack[stackpos] = i;
				i = k;
			}
		}

		if (count == 0) 
		{
			//hairpin:
			sprintf(number, "%f", (float (temp))/conversionfactor);
			out << "Hairpin energy = " << number << "  for closure by " <<
				ct->GetPair(j, n) << "-" << j << "\n";
		}
		else if (count == 1) 
		{
			//bulge or internal loop:
			sprintf(number, "%f", (float (temp))/conversionfactor);
			out << "Bulge/Internal loop energy = " << number << "  for closure by " <<
				ct->GetPair(j, n) << "-" << j << "\n";
		}
		else 
		{
			//multi loop
			sprintf(number, "%f", (float (temp))/conversionfactor);
			out << "Multibranch loop energy = " << number << "  for closure by " <<
				ct->GetPair(j, n) << "-" << j << "\n";
		}
	}
	out.close();
}

//This function will give an energy breakdown for a structure, int n,
//stored in ct, it takes the arrays below as input
//note:this is not fully debugged
void energydump(structure *ct, 
				datatable *data,
				DynProgArray<integersize> *v, 
				int n,
				char *filename	) 
{
	int stack[500], stackpos, i, j, k, temp, count, helix, auaddition;
	ofstream out;
	char number[6], auend[6];
	bool sbulge;

	sbulge = false;
	out.open(filename);
	stackpos = 0;

  
	sprintf(auend, "%f", (float (data->auend))/conversionfactor);

   
	sprintf(number, "%f", (float (ct->GetEnergy(n)))/conversionfactor);
	out << "Structure:  "<<n<<"\n";
	out <<"\n# "<<n<<"  Total Energy = "<<number<<"\n\n";

	temp = ct->GetEnergy(n);

	//Analyze the exterior loop
	i = 0;
	while (i<ct->GetSequenceLength()) 
	{
		i++;
		if (ct->GetPair(i, n)>0) 
		{
			stackpos++;
			stack[stackpos] = i;
			if (data->pairing[ct->numseq[i]][ct->numseq[ct->GetPair(i, n)]]) 
			{
				temp = temp - data->auend;
			}
			temp = temp - v->f(i, ct->GetPair(i, n));
			i =ct->GetPair(i, n);
		}
	}
	sprintf(number, "%f", (float (temp))/conversionfactor);
	out << "Exterior loop energy = "<<number<<"\n";

	while (stackpos > 0) 
	{
		i = stack[stackpos];
		stackpos--;
		helix = 0;

		if (data->pairing[ct->numseq[i]][ct->numseq[ct->GetPair(i, n)]]&&!sbulge) 
		{
			out << "Non-GC end = "<<auend<<"\n";
			helix = helix + data->auend;
		}
		sbulge = false;

		//follow the helix:
		while (ct->GetPair(i+1, n) == ct->GetPair(i, n)-1) 
		{
			//helix continues:
			temp = erg1(i, ct->GetPair(i, n), i+1, ct->GetPair(i+1, n), ct, data);
			//note: v->f(i, ct->basepr[n][i])-v->f(i+1, ct->basepr[n][i+1]);
			//does not work here for cases when a loop is open in min free energy structure
			//but closed in the sub optimal structure being traced
			sprintf(number, "%f", (float (temp))/conversionfactor);
			helix = helix + temp;
			out << "Stack energy = "<<number<<"  for "<<(i+1)<<"-"
			    <<ct->GetPair(i+1, n)<<" onto "<<i<<"-"<<ct->GetPair(i, n)<<"\n";
			i++;
		}

		//now we've come to a loop, what type?
		auaddition = 0;
		j = ct->GetPair(i, n);
		temp = v->f(i, j);
		count = 0;
		while (i<j-1) 
		{
			i++;
			if (ct->GetPair(i, n)>0) 
			{
				//we've found another helix
				if (data->pairing[ct->numseq[i]][ct->numseq[ct->GetPair(i, n)]])
					auaddition = auaddition - data->auend;
				count++;
				k = ct->GetPair(i, n);
				temp = temp - v->f(i, k);
				stackpos++;
				stack[stackpos] = i;
				i = k+1;
			}
		}
		//check for single bulge:
		
		if (count == 1) 
		{
			if (ct->GetPair(j-1, n)>0) 
			{
				if (ct->GetPair(ct->GetPair(j, n) + 1, n)==0) 
				{
					if (ct->GetPair(ct->GetPair(j, n)+2, n)>0) 
					{
						sbulge = true;
						auaddition = 0;
					}
				}
			}
			else if (ct->GetPair(ct->GetPair(j, n) + 1, n)>0) 
			{
				if (ct->GetPair(j-1, n) == 0) 
				{
					if (ct->GetPair(j-2, n) > 0) 
					{
						sbulge = true;
						auaddition = 0;
					}
				}
			}
		}
		if (data->pairing[ct->numseq[j]][ct->numseq[ct->GetPair(j, n)]] && !sbulge) 
		{
			out << "Non-GC end = "<<auend<<"\n";
			helix = helix + data->auend;
		}

		sprintf(number, "%f", (float (helix))/conversionfactor);
		out << "\tHelix energy = "<<number<<"\n";

		if (count == 0) 
		{
			//hairpin:
			temp = temp + auaddition;
			sprintf(number, "%f", (float (temp))/conversionfactor);
			out << "Hairpin energy = "<<number<<"  for closure by "<<
				ct->GetPair(j, n)<<"-"<<j<<"\n";
		}
		else if (count == 1) 
		{
			//bulge or internal loop:
			temp = temp + auaddition;
			sprintf(number, "%f", (float (temp))/conversionfactor);
			out << "Bulge/Internal loop energy = "<<number<<"  for closure by "<<
				ct->GetPair(j, n)<<"-"<<j<<"\n";
		}
		else 
		{
			//multi loop
			temp = temp + auaddition;
			sprintf(number, "%f", (float (temp))/conversionfactor);
			out << "Multibranch loop energy = "<<number<<"  for closure by "<<
				ct->GetPair(j, n)<<"-"<<j<<"\n";
		}
	}
	out.close();
}

//cntrl6 = #tracebacks
//cntrl8 = percent sort
//cntrl9 = window
void readsav(	const char *filename, 
				structure *ct, 
				DynProgArray<integersize> *w2, 
				DynProgArray<integersize> *wmb2, 
				integersize *w5, integersize *w3, 
				bool *lfce, 
				bool *mod, 
				datatable *data,
				DynProgArray<integersize> *v, 
				DynProgArray<integersize> *w, 
				DynProgArray<integersize> *wmb, 
				forceclass *fce, int *vmin	) 
{
	//	int i,j,k,l,m,n,o,p;
	int i, j;

	short vers;//save file version

	ifstream sav(filename, ios::binary);
	
	//read the save file information so that the sequence can be re-folded,
	//data->pairinglude thermodynamic data to prevent traceback errors

	//read the file version first
	read(&sav, &vers);

	//start with structure information
	int sequencelength;

	read(&sav, &(sequencelength));

	read(&sav, &(ct->intermolecular));

	int pairnumber;
	
	read(&sav, &(pairnumber));

	int pair5, pair3;

	for (i = 0; i < pairnumber; i++) 
	{
		read(&sav, &(pair5));	
		read(&sav, &(pair3));
		ct->AddPair(pair5, pair3);
	}
	
	read(&sav, &(pairnumber));

	for (i = 0; i < pairnumber; i++) 
	{
		read(&sav, &(pair5));
		read(&sav, &(pair3));
		ct->AddForbiddenPair(pair5, pair3);
	}

	for (i = 0; i <= ct->GetSequenceLength(); i++) 
	{	
		read(&sav, &(ct->hnumber[i]));
		sav.read(&(ct->nucs[i]), 1);
	}

	for (i = 0; i <= 2 * ct->GetSequenceLength(); i++) read(&sav, &(ct->numseq[i]));
	
	//Read the number of nucleotides forced double stranded
	read(&sav, &(pairnumber));

	//Read and set double stranded nucleotides
	for (i = 0; i < pairnumber; i++) 
	{
		int doubles;
		read(&sav, &(doubles));
		ct->AddDouble(doubles);
	}

	if (ct->intermolecular) 
	{
		w2 = new DynProgArray<integersize>(ct->GetSequenceLength());
		wmb2 = new DynProgArray<integersize>(ct->GetSequenceLength());		
		for (i = 0; i < 3; i++) read(&sav, &(ct->inter[i]));
	}

	//Read the number of nucleotides that are not allowed to pair
	read(&sav, &(pairnumber));

	//Read the nucleotides that cannot pair and store the information.
	for (i = 0; i < pairnumber; i++) 
	{
		int singles;
		read(&sav, &(singles));
		ct->AddSingle(singles);
	}

	//Read the number of chemically modified nucleotides
	read(&sav, &(pairnumber));

	//Read and store the chemically modified nucleotides
	for (i = 0; i < pairnumber; i++) 
	{
		int modified;
		read(&sav, &(modified));
		ct->AddModified(modified);
	}
	
	//Read the number of Us in GU pairs
	read(&sav, &(pairnumber));

	//Read the Us in GU pairs, and store information
	for (i = 0; i < pairnumber; i++) 
	{	
		int GU;
		read(&sav, &(GU));  
		ct->AddGUPair(GU);
	}
	
	string label;
	read(&sav, &(label));
	ct->SetSequenceLabel(label);

	read(&sav, &(ct->templated));
	if (ct->templated) 
	{
		ct->allocatetem();
		for (i = 0;i <= ct->GetSequenceLength(); i++) 
		{
			for (j = 0; j <= i; j++) read(&sav, &(ct->tem[i][j]));	
		}
	}

	read(&sav, &ct->shaped);
	if (ct->shaped) 
	{
		ct->SHAPE = new double [2*ct->GetSequenceLength() + 1];
		for (i = 0; i <= 2*ct->GetSequenceLength(); i++) read(&sav, &(ct->SHAPE[i]));
	}
	
	//now read the array class data for v, w, and wmb:
	//now write the array class data for v, w, and wmb:
	for (i = 0; i <= ct->GetSequenceLength(); i++) 
	{
		read(&sav, &(w3[i]));
		read(&sav, &(w5[i]));
		for (j = 0; j <= ct->GetSequenceLength(); j++) 
		{
			read(&sav, &(v->dg[i][j+i]));
			read(&sav, &(w->dg[i][j+i]));
			read(&sav, &(wmb->dg[i][j+i]));
			readsinglechar(&sav, &(fce->dg[i][j]));
			if (ct->intermolecular) 
			{
				read(&sav, &(w2->dg[i][j+i]));
				read(&sav, &(wmb2->dg[i][j+i]));
			}
		}
	}

	read(&sav, &(w3[ct->GetSequenceLength() + 1]));
	
	for (i = 0; i<=2 * ct->GetSequenceLength(); i++) 
	{
		read(&sav, &(lfce[i]));
		read(&sav, &(mod[i]));
	}

	read(&sav, vmin);
	
	read(&sav, data);

	//set the location of the datatables inside the structure class:
	ct->SetThermodynamicDataTable(data);
	
	sav.close();
}

//opensav will open a save file created by the fill algorithm (function dynamic)
//	it then runs the traceback routine with the parameters provided
void opensav(	char* filename, 
				structure* ct, 
				int cntrl6, 
				int cntrl8,
				int cntrl9	) 
{
	int i;
	DynProgArray<integersize> *w2, *wmb2;
	integersize *w5, *w3;
	int vmin;
	bool *lfce, *mod;
	datatable *data;

	data = new datatable;
	short vers;

	//peek at the length of the sequence and whether the folding is intermolecular to allocate arrays:
	ifstream sav(filename, ios::binary);
	
	read(&sav, &vers);//could add functionality here to make sure that the correct version is being read
	int sequencelength;
	read(&sav, &(sequencelength));
	read(&sav, &(ct->intermolecular));
	sav.close();

	//allocate everything
	ct->allocate(sequencelength);
	
	DynProgArray<integersize> w(ct->GetSequenceLength());
	DynProgArray<integersize> v(ct->GetSequenceLength());
	DynProgArray<integersize> wmb(ct->GetSequenceLength());	

	forceclass fce(ct->GetSequenceLength());

	lfce = new bool [2*ct->GetSequenceLength()+1];
	mod = new bool [2*ct->GetSequenceLength()+1];
	
	w5 = new integersize [ct->GetSequenceLength()+1];
	w3 = new integersize [ct->GetSequenceLength()+2];

	if (ct->intermolecular) 
	{
		w2 = new DynProgArray<integersize>(ct->GetSequenceLength());
		wmb2 = new DynProgArray<integersize>(ct->GetSequenceLength());	
		for (i=0;i<3;i++) read(&sav, &(ct->inter[i]));
	}
	else 
	{
		w2 = NULL;
		wmb2 = NULL;
	}
	
	readsav(filename, ct, w2, wmb2, w5, w3, lfce, mod, data, &v, &w, &wmb, &fce, &vmin);
	traceback(ct, data, &v, &w, &wmb, w2, wmb2, w3, w5, &fce, lfce, vmin, cntrl6, cntrl8, cntrl9, mod);

	delete[] lfce;
	delete[] mod;	
	delete[] w5;
	delete[] w3;

	if (ct->intermolecular) 
	{
		delete w2;
		delete wmb2;
	}

	delete data;
}

#endif //END INSTRUMENTED

#ifndef INSTRUMENTED
void arraydump(	DynProgArray<integersize>& v,
				DynProgArray<integersize>& w,
				DynProgArray<integersize>& wmb, 
				const integersize* w5, 
				const integersize* w3, 
				const int n, 
				const char* filename)
{
	ofstream foo;
	foo.open(filename);
	foo << "i" << "\t"<<"j"<<"\t"<<"v.f(i, j)"<<"\t"<<"w.f(i, j)"<<"\t"<<"wmb.f(i, j)"<<"\t"<<"v.f(j, i+number)"<<"\t"<<"w.f(j, i+number)"<<"\t"<<"wmb.f(j, i+number)"<<"\n";
	for (int j = 1; j <= n; j++) 
	{
		for (int i = 1; i <= j; i++) 
		{
			foo << i << "\t"<<j<<"\t"<<v.f(i, j)<<"\t"<<w.f(i, j)<<"\t"<<wmb.f(i, j)<<"\t"<<v.f(j, i+n)<<"\t"<<w.f(j, i+n)<<"\t"<<wmb.f(j, i+n)<<"\n";
		}	
	}
	foo <<"\n\n\n";
	foo << "i" << "\t" << "w5[i]" << "\t" << "w3[i]" << "\n";
	for (int i = 1;i <= n; i++) foo << i << "\t" << w5[i] << "\t" << w3[i] << "\n";

	foo.close();
}
#endif

#ifndef INSTRUMENTED
/* bool multi_can_pair(int i, int j, structure *ct) {
	return ct->pairing_matrix[i][j];
}*/

integersize multi_eparam(int i, structure *ct) 
{
	return ct->GetThermodynamicDataTable()->eparam[i] * ct->number_of_sequences;
}

integersize multi_tstack(int a, int b, int c, int d, structure *ct)
{		
	int N = ct->number_of_sequences;
	if (N == 1) {
		return ct->GetThermodynamicDataTable()->tstack[ct->numseq[a]][ct->numseq[b]][ct->numseq[c]][ct->numseq[d]];
	}
	else if (N > 1) {
		integersize energy = 0;
		for (int n = 0; n < N; n++)
		{
			structure* ct_n = ct->get_individual_sequence(n);
			datatable* dt_n = ct_n->GetThermodynamicDataTable();
			energy = energy + dt_n->tstack[ct_n->numseq[a]][ct_n->numseq[b]][ct_n->numseq[c]][ct_n->numseq[d]];
		}
		return energy;
	}	
}

integersize multi_tstkm(int a, int b, int c, int d, structure *ct)
{
	int N = ct->number_of_sequences;
	if (N == 1) {
		return ct->GetThermodynamicDataTable()->tstkm[ct->numseq[a]][ct->numseq[b]][ct->numseq[c]][ct->numseq[d]];
	}
	else if (N > 1)	{
		integersize energy = 0;
		for (int n = 0; n < N; n++)
		{
			structure* ct_n = ct->get_individual_sequence(n);
			datatable* dt_n = ct_n->GetThermodynamicDataTable();
			energy = energy + dt_n->tstkm[ct_n->numseq[a]][ct_n->numseq[b]][ct_n->numseq[c]][ct_n->numseq[d]];
		}
		return energy;
	}
}

integersize covar(int i, int j, structure* ct) {
	return 0;
}
#endif


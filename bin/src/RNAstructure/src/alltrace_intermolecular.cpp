//This version of alltrace is used in intermolecular.* for OligoWalk calculations.
//This is largely redundant to the code in alltrace.cpp and should be removed in the future.

#include "alltrace.h"
#include "alltrace_intermolecular.h"

#ifndef _WINDOWS_GUI
	#include "platform.h"
#endif

#define startingsize 500  //maximum number of structure fragments to start
//#define startingrefinementstacksize 25

#undef pfdebugmode  //flag to indicate debugging
//#define pfdebugmode

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>

using namespace std;


//perecentdelta is the maximum percent energy diff in suboptimal structures from optimal
//absolutedelta is the maximum energy diff in suboptimal structures from optimal
//	the minimum difference from the above is used as the cutoff stored in crit
//Note, if ct->stacking is set to true, the stacking of nucleotides in multibranch and exterior loops will be tracked
//	On the stack, a flush coaxial stack (eg i-j on j+1 k) is stacked as i on k and k on i.
//	Also on the stack, an intervening stack (eg i-j with j+2 k and k+1 intervening) is stacked as i on k+1 and k+1 on k.
//If ctname is set (to other than null), the structures will be written to a ct file as they are produced.  This is helpful for long
//seqeunces, where even a small energy increment (delta) can lead to too many structures to store in memory.
void alltracetraceback(bool writect,int *ctenergy, structure *ct, atDynProgArray *v, atDynProgArray *w, atDynProgArray *wmb, atDynProgArray *wl, atDynProgArray *wmbl, 
	atDynProgArray *wcoax, forceclass *fce, short *w5, bool *lfce, bool *mod, datatable *data, short percentdelta, short absolutedelta, 
	char *ctname=NULL) {
	
	alltracestructurestack stack(ct->GetSequenceLength());
	
	short i,j,k,ip,d,jp,kp;
	integersize crit;
	integersize energy,branch;
	short open,pair;
	long int current;
	int count,pos,pos2;
	bool passed,found;
	ofstream *out;
	
	int energy_size=MaxStructure; //these two will the be used to increase the size of ctenergy array
	int *tempenergy;

	vector< vector<bool> > inc = data->pairing;

	current = 1; //#of structure to be traced


	if (ctname!=NULL) {

		out = new ofstream(ctname);
		*out << "-100\n";
		*out << ct->GetSequenceLength()<<"\n";
		*out << "0"<<" ";//this line of a cct file would normally indicate the number of structures
						//in this case, the number of structures is unknown, so a 0 will be used

		*out << ct->GetSequenceLabel();

		for (i=1;i<=ct->GetSequenceLength();i++) {
			*out << ct->numseq[i]<<"\n";
		}

	}


	stack.pushtorefinement(1,ct->GetSequenceLength(),1,w5[ct->GetSequenceLength()],0);
	stack.energy[0]=w5[ct->GetSequenceLength()];

	crit = w5[ct->GetSequenceLength()]-max(-absolutedelta,(short)((double)w5[ct->GetSequenceLength()]*(double)percentdelta/100.0));

		
	while (stack.current>-1) {
		while (stack.pullfromrefinement(&i,&j,&open,&energy,&pair)) {
			//continue with refinements

			if (open==1) {
				//this is an exterior fragment
				
				
				
				//check all possibilities:
				if (energy!=0) {
					//if energy==0, this is a fragment with no interactions, throw it away from consideration...
					if ((stack.peekatenergy()+w5[j-1]-energy)<=crit) {

						if (j-i>minloop) stack.push(stack.peekatenergy()+w5[j-1]-energy,false,0,0,true,i,j-1,1,w5[j-1],0,false,0,0,false,0,0,0,0,0);
					
					}
					//check for split:

					for (k=0;k<=j-4;k++) {
						
						if ((stack.peekatenergy()+w5[k]+v->f(k+1,j)+penalty(j,k+1,ct,data)-energy)<=crit) {
							stack.push(stack.peekatenergy()+w5[k]+v->f(k+1,j)+penalty(j,k+1,ct,data)-energy,true,k+1,j,true,k+1,j,0,v->f(k+1,j),1,false,0,0,true,1,k,1,w5[k],0);
						}

						if (mod[k+1]||mod[j]&&inc[ct->numseq[k+2]][ct->numseq[j-1]]&&notgu(k+1,j,ct)&&!(fce->f(k+1,j)&SINGLE)) {
							if ((stack.peekatenergy()+w5[k]+v->f(k+2,j-1)+erg1(k+1,j,k+2,j-1,ct,data)+penalty(j,k+1,ct,data)-energy)<=crit) {
								stack.push(stack.peekatenergy()+w5[k]+v->f(k+1,j)+erg1(k+1,j,k+2,j-1,ct,data)+penalty(j,k+1,ct,data)-energy,true,k+1,j,true,k+1,j,0,v->f(k+2,j-1)+erg1(k+1,j,k+2,j-1,ct,data),1,false,0,0,true,1,k,1,w5[k],0);
						
							}

						}

						//now look at a dangling end:
						if ((stack.peekatenergy()+w5[k]+erg4(j,k+2,k+1,2,ct,data,lfce[k+1])+v->f(k+2,j)+penalty(j,k+2,ct,data)-energy)<=crit) {
							stack.push(stack.peekatenergy()+w5[k]+erg4(j,k+2,k+1,2,ct,data,lfce[k+1])+v->f(k+2,j)+penalty(j,k+2,ct,data)-energy,true,k+2,j,true,k+2,j,0,v->f(k+2,j),1,false,0,0,true,1,k,1,w5[k],0);
							stack.nstack(k+1,k+2);
						}

						if(mod[k+2]||mod[j]&&inc[ct->numseq[k+3]][ct->numseq[j-1]]&&notgu(k+2,j,ct)
							&&!(fce->f(k+2,j)&SINGLE)) {

							if ((stack.peekatenergy()+w5[k]+erg4(j,k+2,k+1,2,ct,data,lfce[k+1])+v->f(k+3,j-1)+penalty(j,k+2,ct,data)+erg1(k+2,j,k+3,j-1,ct,data)-energy)<=crit) {
								stack.push(stack.peekatenergy()+w5[k]+erg4(j,k+2,k+1,2,ct,data,lfce[k+1])+v->f(k+3,j-1)+penalty(j,k+2,ct,data)+erg1(k+2,j,k+3,j-1,ct,data)-energy,true,k+2,j,true,k+2,j,0,v->f(k+3,j-1)+erg1(k+2,j,k+3,j-1,ct,data),1,false,0,0,true,1,k,1,w5[k],0);
								stack.nstack(k+1,k+2);
							}

						}

						//other dangling end:
						if ((stack.peekatenergy()+w5[k]+erg4(j-1,k+1,j,1,ct,data,lfce[j])+v->f(k+1,j-1)+penalty(j-1,k+1,ct,data)-energy)<=crit) {
							stack.push(stack.peekatenergy()+w5[k]+erg4(j-1,k+1,j,1,ct,data,lfce[j])+v->f(k+1,j-1)+penalty(j-1,k+1,ct,data)-energy,true,k+1,j-1,true,k+1,j-1,0,v->f(k+1,j-1),1,false,0,0,true,1,k,1,w5[k],0);
							stack.nstack(j,j-1);
						}
						
						if (mod[k+1]||mod[j-1]&&inc[ct->numseq[k+2]][ct->numseq[j-2]]&&notgu(k+1,j-1,ct)&&!(fce->f(k+1,j-1)&SINGLE)) {
							if ((stack.peekatenergy()+w5[k]+erg4(j-1,k+1,j,1,ct,data,lfce[j])+v->f(k+2,j-2)+penalty(j-1,k+1,ct,data)+erg1(k+1,j-1,k+2,j-2,ct,data)-energy)<=crit) {
								stack.push(stack.peekatenergy()+w5[k]+erg4(j-1,k+1,j,1,ct,data,lfce[j])+v->f(k+2,j-2)+penalty(j-1,k+1,ct,data)+erg1(k+1,j-1,k+2,j-2,ct,data)-energy,true,k+1,j-1,true,k+1,j-1,0,v->f(k+2,j-2)+erg1(k+1,j-1,k+2,j-2,ct,data),1,false,0,0,true,1,k,1,w5[k],0);
								stack.nstack(j,j-1);
							}

						}

						//terminal stack possibility:
						if ((stack.peekatenergy()+w5[k]+data->tstack[ct->numseq[j-1]][ct->numseq[k+2]][ct->numseq[j]][ct->numseq[k+1]] 
							+checknp(lfce[j],lfce[k+1]) + v->f(k+2,j-1)+
							penalty(j-1,k+2,ct,data)-energy)<=crit) {

							stack.push(stack.peekatenergy()+w5[k]+data->tstack[ct->numseq[j-1]][ct->numseq[k+2]][ct->numseq[j]][ct->numseq[k+1]] 
									+checknp(lfce[j],lfce[k+1]) + v->f(k+2,j-1)+
									penalty(j-1,k+2,ct,data)-energy,true,k+2,j-1,true,k+2,j-1,0,v->f(k+2,j-1),1,false,0,0,true,1,k,1,w5[k],0);
							stack.nstack(j,j-1,k+1,k+2);
						}

						if (mod[k+2]||mod[j-1]&&inc[ct->numseq[k+3]][ct->numseq[j-2]]&&notgu(k+2,j-1,ct)&&!(fce->f(k+2,j-1)&SINGLE)) {

							if ((stack.peekatenergy()+w5[k]+data->tstack[ct->numseq[j-1]][ct->numseq[k+2]][ct->numseq[j]][ct->numseq[k+1]] 
									+checknp(lfce[j],lfce[k+1]) + v->f(k+3,j-2)+
									penalty(j-1,k+2,ct,data)+erg1(k+2,j-1,k+3,j-2,ct,data)-energy)<=crit) {

								stack.push(stack.peekatenergy()+w5[k]+data->tstack[ct->numseq[j-1]][ct->numseq[k+2]][ct->numseq[j]][ct->numseq[k+1]] 
									+checknp(lfce[j],lfce[k+1]) + v->f(k+3,j-2)+
									penalty(j-1,k+2,ct,data)+erg1(k+2,j-1,k+3,j-2,ct,data)-energy,true,k+2,j-1,true,k+2,j-1,0,v->f(k+3,j-2)+erg1(k+2,j-1,k+3,j-2,ct,data),1,false,0,0,true,1,k,1,w5[k],0);
								stack.nstack(j,j-1,k+1,k+2);
							}	

						}

						//now check coaxial stacking options:
						for (ip=k+minloop+1;ip<j-minloop-1;ip++) {
							//first consider flush stacking
							
							if ((stack.peekatenergy()+w5[k]+v->f(k+1,ip)+v->f(ip+1,j)+penalty(k+1,ip,ct,data)
								+penalty(ip+1,j,ct,data)+ergcoax(k+1,ip,ip+1,j,0,ct,data)-energy)<=crit) {

								stack.push(stack.peekatenergy()+w5[k]+v->f(k+1,ip)+v->f(ip+1,j)+penalty(k+1,ip,ct,data)
								+penalty(ip+1,j,ct,data)+ergcoax(k+1,ip,ip+1,j,0,ct,data)-energy,true,k+1,ip,true,k+1,ip,0,v->f(k+1,ip),1,true,ip+1,j,true,ip+1,j,0,v->f(ip+1,j),1,true,1,k,1,w5[k],0);
								stack.nstack(k+1,j,j,k+1);
							}
							
							
							if ((mod[k+1]||mod[ip]||mod[ip+1]||mod[j])) {

								if ((mod[k+1]||mod[ip])&&(mod[ip+1]||mod[j])&&inc[ct->numseq[k+1+1]][ct->numseq[ip-1]]
									&&inc[ct->numseq[ip+2]][ct->numseq[j-1]]&&notgu(k+1,ip,ct)&&notgu(ip+1,j,ct)
										&&!(fce->f(ip+1,j)&SINGLE)&&!(fce->f(i,ip)&SINGLE)) {


									if ((stack.peekatenergy()+w5[k]+v->f(k+2,ip-1)+v->f(ip+2,j-1)+penalty(k+1,ip,ct,data)
										+penalty(ip+1,j,ct,data)+ergcoax(k+1,ip,ip+1,j,0,ct,data)
										+erg1(k+1,ip,k+2,ip-1,ct,data)+erg1(ip+1,j,ip+2,j-1,ct,data)-energy)<=crit) {

										stack.push(stack.peekatenergy()+w5[k]+v->f(k+2,ip-1)+v->f(ip+2,j-1)+penalty(k+1,ip,ct,data)
											+penalty(ip+1,j,ct,data)+ergcoax(k+1,ip,ip+1,j,0,ct,data)
											+erg1(k+1,ip,k+2,ip-1,ct,data)+erg1(ip+1,j,ip+2,j-1,ct,data)-energy,true,k+1,ip,true,k+1,ip,0,v->f(k+2,ip-1)+erg1(k+1,ip,k+2,ip-1,ct,data),1,true,ip+1,j,true,ip+1,j,0,v->f(ip+2,j-1)+erg1(ip+1,j,ip+2,j-1,ct,data),1,true,1,k,1,w5[k],0);
										stack.nstack(k+1,j,j,k+1);
									}
									

								}

								if ((mod[k+1]||mod[ip])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]&&notgu(k+1,ip,ct)&&!(fce->f(k+1,ip)&SINGLE)) {
							

									if ((stack.peekatenergy()+w5[k]+v->f(k+2,ip-1)+v->f(ip+1,j)+penalty(k+1,ip,ct,data)
										+penalty(ip+1,j,ct,data)+ergcoax(k+1,ip,ip+1,j,0,ct,data)
										+erg1(k+1,ip,k+2,ip-1,ct,data)-energy)<=crit) {

										stack.push(stack.peekatenergy()+w5[k]+v->f(k+2,ip-1)+v->f(ip+1,j)+penalty(k+1,ip,ct,data)
										+penalty(ip+1,j,ct,data)+ergcoax(k+1,ip,ip+1,j,0,ct,data)
										+erg1(k+1,ip,k+2,ip-1,ct,data)-energy,true,k+1,ip,true,k+1,ip,0,v->f(k+2,ip-1)+erg1(k+1,ip,k+2,ip-1,ct,data),1,true,ip+1,j,true,ip+1,j,0,v->f(ip+1,j),1,true,1,k,1,w5[k],0);
										stack.nstack(k+1,j,j,k+1);
									}
									
										
								}

								if ((mod[ip+1]||mod[j])&&inc[ct->numseq[ip+2]][ct->numseq[j-1]]&&notgu(ip+1,j,ct)&&!(fce->f(ip+1,j)&SINGLE)) {


									if ((stack.peekatenergy()+w5[k]+v->f(k+1,ip)+v->f(ip+2,j-1)+penalty(k+1,ip,ct,data)
										+penalty(ip+1,j,ct,data)+ergcoax(k+1,ip,ip+1,j,0,ct,data)
										+erg1(ip+1,j,ip+2,j-1,ct,data)-energy)<=crit) {

										stack.push(stack.peekatenergy()+w5[k]+v->f(k+1,ip)+v->f(ip+2,j-1)+penalty(k+1,ip,ct,data)
										+penalty(ip+1,j,ct,data)+ergcoax(k+1,ip,ip+1,j,0,ct,data)
										+erg1(ip+1,j,ip+2,j-1,ct,data)-energy,true,k+1,ip,true,k+1,ip,0,v->f(k+1,ip),1,true,ip+1,j,true,ip+1,j,0,v->f(ip+2,j-1)+erg1(ip+1,j,ip+2,j-1,ct,data),1,true,1,k,1,w5[k],0);
										stack.nstack(k+1,j,j,k+1);
									}
									

								}


							}
						
							if (!lfce[ip+1]&&!lfce[j]) {
								//now consider an intervening mismatch
								if(!lfce[ip+1]&&!lfce[j]) {

									if ((stack.peekatenergy()+w5[k]+v->f(k+1,ip)+v->f(ip+2,j-1)+penalty(k+1,ip,ct,data)
										+penalty(ip+2,j-1,ct,data)+ergcoax(k+1,ip,ip+2,j-1,j,ct,data)-energy)<=crit) {

										stack.push(stack.peekatenergy()+w5[k]+v->f(k+1,ip)+v->f(ip+2,j-1)+penalty(k+1,ip,ct,data)
											+penalty(ip+2,j-1,ct,data)+ergcoax(k+1,ip,ip+2,j-1,j,ct,data)-energy,true,k+1,ip,true,k+1,ip,0,v->f(k+1,ip),1,true,ip+2,j-1,true,ip+2,j-1,0,v->f(ip+2,j-1),1,true,1,k,1,w5[k],0);
										stack.nstack(k+1,j,j,j-1);
									}
								

								}
								if (mod[k+1]||mod[ip]||mod[ip+2]||mod[j-1]) {
									if ((mod[k+1]||mod[ip])&&(mod[ip+2]||mod[j-1])&&inc[ct->numseq[k+1+1]][ct->numseq[ip-1]]
										&&inc[ct->numseq[ip+3]][ct->numseq[j-2]]&&notgu(k+1,ip,ct)&&notgu(ip+2,j-1,ct)
											&&!(fce->f(k+1,ip)&SINGLE)&&!(fce->f(ip+2,j-1)&SINGLE)) {

										if(!lfce[ip+1]&&!lfce[j]) {

											if ((stack.peekatenergy()+w5[k]+v->f(k+1+1,ip-1)+v->f(ip+3,j-2)+penalty(k+1,ip,ct,data)
												+penalty(ip+2,j-1,ct,data)+ergcoax(k+1,ip,ip+2,j-1,j,ct,data)
												+erg1(k+1,ip,k+1+1,ip-1,ct,data)+erg1(ip+2,j-1,ip+3,j-2,ct,data)-energy)<=crit) {

												stack.push(stack.peekatenergy()+w5[k]+v->f(k+1+1,ip-1)+v->f(ip+3,j-2)+penalty(k+1,ip,ct,data)
												+penalty(ip+2,j-1,ct,data)+ergcoax(k+1,ip,ip+2,j-1,j,ct,data)
												+erg1(k+1,ip,k+1+1,ip-1,ct,data)+erg1(ip+2,j-1,ip+3,j-2,ct,data)-energy,true,k+1,ip,true,k+1,ip,0,v->f(k+2,ip-1)+erg1(k+1,ip,k+1+1,ip-1,ct,data),1,true,ip+2,j-1,true,ip+2,j-1,0,v->f(ip+3,j-2)+erg1(ip+2,j-1,ip+3,j-2,ct,data),1,true,1,k,1,w5[k],0);
												stack.nstack(k+1,j,j,j-1);
											}

										}
											


									}

									if ((mod[k+1]||mod[ip])&&inc[ct->numseq[k+1+1]][ct->numseq[ip-1]]&&notgu(k+1,ip,ct)&&!(fce->f(k+1,ip)&SINGLE)) {
							
										if(!lfce[ip+1]&&!lfce[j]) {
											
											if ((stack.peekatenergy()+w5[k]+v->f(k+1+1,ip-1)+v->f(ip+2,j-1)+penalty(k+1,ip,ct,data)
												+penalty(ip+2,j-1,ct,data)+ergcoax(k+1,ip,ip+2,j-1,j,ct,data)
												+erg1(k+1,ip,k+1+1,ip-1,ct,data)-energy)<=crit) {

												stack.push(stack.peekatenergy()+w5[k]+v->f(k+1+1,ip-1)+v->f(ip+2,j-1)+penalty(k+1,ip,ct,data)
													+penalty(ip+2,j-1,ct,data)+ergcoax(k+1,ip,ip+2,j-1,j,ct,data)
													+erg1(k+1,ip,k+1+1,ip-1,ct,data)-energy,true,k+1,ip,true,k+1,ip,0,v->f(k+2,ip-1)+erg1(k+1,ip,k+1+1,ip-1,ct,data),1,true,ip+2,j-1,true,ip+2,j-1,0,v->f(ip+2,j-1),1,true,1,k,1,w5[k],0);
												stack.nstack(k+1,j,j,j-1);
											}
										}
										
									}

									if ((mod[ip+2]||mod[j-1])&&inc[ct->numseq[ip+3]][ct->numseq[j-2]]&&notgu(ip+2,j-1,ct)&&!(fce->f(ip+2,j-1)&SINGLE)) {


										if(!lfce[ip+1]&&!lfce[j]) {
											
											if ((stack.peekatenergy()+w5[k]+v->f(k+1,ip)+v->f(ip+3,j-2)+penalty(k+1,ip,ct,data)
											+penalty(ip+2,j-1,ct,data)+ergcoax(k+1,ip,ip+2,j-1,j,ct,data)
											+erg1(ip+2,j-1,ip+3,j-2,ct,data)-energy)<=crit) {

												stack.push(stack.peekatenergy()+w5[k]+v->f(k+1,ip)+v->f(ip+3,j-2)+penalty(k+1,ip,ct,data)
													+penalty(ip+2,j-1,ct,data)+ergcoax(k+1,ip,ip+2,j-1,j,ct,data)
													+erg1(ip+2,j-1,ip+3,j-2,ct,data)-energy,true,k+1,ip,true,k+1,ip,0,v->f(k+1,ip),1,true,ip+2,j-1,true,ip+2,j-1,0,v->f(ip+3,j-2)+erg1(ip+2,j-1,ip+3,j-2,ct,data),1,true,1,k,1,w5[k],0);
												stack.nstack(k+1,j,j,j-1);	
											}
											
											
										}
										
										
									}
								}
							}
		
							if(!lfce[k+1]&&!lfce[ip+1]) {
								//check the other coaxial stack
								if ((stack.peekatenergy()+w5[k]+v->f(k+1+1,ip)+v->f(ip+2,j)+penalty(k+1+1,ip,ct,data)
									+penalty(ip+2,j,ct,data)+ergcoax(k+1+1,ip,ip+2,j,k+1,ct,data)-energy)<=crit) {

									stack.push(stack.peekatenergy()+w5[k]+v->f(k+1+1,ip)+v->f(ip+2,j)+penalty(k+1+1,ip,ct,data)
									+penalty(ip+2,j,ct,data)+ergcoax(k+1+1,ip,ip+2,j,k+1,ct,data)-energy,true,k+2,ip,true,k+2,ip,0,v->f(k+2,ip),1,true,ip+2,j,true,ip+2,j,0,v->f(ip+2,j),1,true,1,k,1,w5[k],0);
									stack.nstack(k+2,k+1,k+1,j);
								}
								

								if (mod[k+1+1]||mod[ip]||mod[ip+2]||mod[j]) {
									if ((mod[k+1+1]||mod[ip])&&(mod[ip+2]||mod[j])&&inc[ct->numseq[k+1+2]][ct->numseq[ip-1]]
										&&inc[ct->numseq[ip+3]][ct->numseq[j-1]]&&notgu(k+1+1,ip,ct)&&notgu(ip+2,j,ct)
										&&!(fce->f(k+1+1,ip)&SINGLE)&&!(fce->f(ip+2,j)&SINGLE)	) {


										if ((stack.peekatenergy()+w5[k]+v->f(k+1+2,ip-1)+v->f(ip+3,j-1)+penalty(k+1+1,ip,ct,data)
												+penalty(ip+2,j,ct,data)+ergcoax(k+1+1,ip,ip+2,j,k+1,ct,data)
												+erg1(k+1+1,ip,k+1+2,ip-1,ct,data)+erg1(ip+2,j,ip+3,j-1,ct,data)-energy)<=crit) {

											stack.push(stack.peekatenergy()+v->f(k+1+2,ip-1)+v->f(ip+3,j-1)+penalty(k+1+1,ip,ct,data)
												+penalty(ip+2,j,ct,data)+ergcoax(k+1+1,ip,ip+2,j,k+1,ct,data)
												+erg1(k+1+1,ip,k+1+2,ip-1,ct,data)+erg1(ip+2,j,ip+3,j-1,ct,data)-energy,true,k+2,ip,true,k+2,ip,0,v->f(k+3,ip-1)+erg1(k+1+1,ip,k+1+2,ip-1,ct,data),1,true,ip+2,j,true,ip+2,j,0,v->f(ip+3,j-1)+erg1(ip+2,j,ip+3,j-1,ct,data),1,true,1,k,1,w5[k],0);
											stack.nstack(k+2,k+1,k+1,j);
										}
										
												

									}
									if ((mod[k+1+1]||mod[ip])&&inc[ct->numseq[k+1+2]][ct->numseq[ip-1]]&&notgu(k+1+1,ip,ct)&&!(fce->f(k+1+1,ip)&SINGLE)) {
							
										if ((stack.peekatenergy()+w5[k]+v->f(k+1+2,ip-1)+v->f(ip+2,j)+penalty(k+1+1,ip,ct,data)
											+penalty(ip+2,j,ct,data)+ergcoax(k+1+1,ip,ip+2,j,k+1,ct,data)
											+erg1(k+1+1,ip,k+1+2,ip-1,ct,data)-energy)<=crit) {

											stack.push(stack.peekatenergy()+v->f(k+1+2,ip-1)+v->f(ip+2,j)+penalty(k+1+1,ip,ct,data)
											+penalty(ip+2,j,ct,data)+ergcoax(k+1+1,ip,ip+2,j,k+1,ct,data)
											+erg1(k+1+1,ip,k+1+2,ip-1,ct,data)-energy,true,k+2,ip,true,k+2,ip,0,v->f(k+3,ip-1)+erg1(k+1+1,ip,k+1+2,ip-1,ct,data),1,true,ip+2,j,true,ip+2,j,0,v->f(ip+2,j),1,true,1,k,1,w5[k],0);
											stack.nstack(k+2,k+1,k+1,j);
										}
										
									}

									if ((mod[ip+2]||mod[j])&&inc[ct->numseq[ip+3]][ct->numseq[j-1]]&&notgu(ip+2,j,ct)&&!(fce->f(ip+2,j)&SINGLE)) {

										if ((stack.peekatenergy()+w5[k]+v->f(k+1+1,ip)+v->f(ip+3,j-1)+penalty(k+1+1,ip,ct,data)
											+penalty(ip+2,j,ct,data)+ergcoax(k+1+1,ip,ip+2,j,k+1,ct,data)
											+erg1(ip+2,j,ip+3,j-1,ct,data)-energy)<=crit) {

											stack.push(stack.peekatenergy()+v->f(k+1+1,ip)+v->f(ip+3,j-1)+penalty(k+1+1,ip,ct,data)
											+penalty(ip+2,j,ct,data)+ergcoax(k+1+1,ip,ip+2,j,k+1,ct,data)
											+erg1(ip+2,j,ip+3,j-1,ct,data)-energy,true,k+2,ip,true,k+2,ip,0,v->f(k+2,ip),1,true,ip+2,j,true,ip+2,j,0,v->f(ip+3,j-1)+erg1(ip+2,j,ip+3,j-1,ct,data),1,true,1,k,1,w5[k],0);
											stack.nstack(k+2,k+1,k+1,j);
										}


									}
								}
							}
							
							

						}


					}

					stack.flushbullpen();
				}
				//after checking all possibilities, record the last one found into the stack:
				

			

			}
			else {
				//open == 0 -- so this is not an exterior loop fragment
				if (pair==1) {
					//i-j is a pair
					if (mod[i]||mod[j]) {
						//this is a modified pair that was chosen because it is stacked on the next pair and closed a loop
						if (energy==v->f(i+1,j-1)+erg1(i,j,i+1,j-1,ct,data)) {

							stack.push(stack.peekatenergy(),true,i+1,j-1,true,i+1,j-1,0,v->f(i+1,j-1),1,false,0,0,false,0,0,0,0,0);
							

						}


					}
					else {
						//we have a pair, decide what options can branch from here:

						//check for hairpin
						if ((stack.peekatenergy()+erg3(i,j,ct,data,fce->f(i,j))-energy)<=crit) {
							//hairpin is acceptable
							stack.push((stack.peekatenergy()+erg3(i,j,ct,data,fce->f(i,j)))-energy,false,0,0,false,0,0,0,0,0,false,0,0,false,0,0,0,0,0);
							
						}

						//check for stack:
						if ((stack.peekatenergy()+erg1(i,j,i+1,j-1,ct,data)+v->f(i+1,j-1)-energy)<=crit) {
							//stack is acceptable:
							stack.push((stack.peekatenergy()+erg1(i,j,i+1,j-1,ct,data)+v->f(i+1,j-1))-energy,true,i+1,j-1,true,i+1,j-1,0,v->f(i+1,j-1),1,false,0,0,false,0,0,0,0,0);
				
						}
						//check all possible internal/bulge loops:
						      
						if ((j-i-1)>=(minloop+3)) {
      						for (d=(j-i-3);d>=1;d--) {
         						for (ip=(i+1);ip<=(j-1-d);ip++) {
            						jp = d+ip;
									if ((j-i-2-d)>(data->eparam[7])) goto skip;
									if (abs(ip-i+jp-j)<=(data->eparam[8])) {
               						

										if ((stack.peekatenergy()+erg2(i,j,ip,jp,ct,data,fce->f(i,ip),fce->f(jp,j))+v->f(ip,jp)-energy)<=crit) {
											//loop is acceptable:
											stack.push((stack.peekatenergy()+erg2(i,j,ip,jp,ct,data,fce->f(i,ip),fce->f(jp,j))+v->f(ip,jp)-energy),true,ip,jp,true,ip,jp,0,v->f(ip,jp),1,false,0,0,false,0,0,0,0,0);
				
										}

                  							

										if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]&&notgu(ip,jp,ct)&&!(fce->f(ip,jp)&SINGLE)) {
											//i or j is modified

											if ((stack.peekatenergy()+erg2(i,j,ip,jp,ct,data,fce->f(i,ip),fce->f(jp,j))+
                  								v->f(ip+1,jp-1)+erg1(ip,jp,ip+1,jp-1,ct,data)-energy)<=crit) {
											//loop is acceptable:
												stack.push((stack.peekatenergy()+erg2(i,j,ip,jp,ct,data,fce->f(i,ip),fce->f(jp,j))+v->f(ip,jp)-energy),true,ip,jp,true,ip,jp,0,v->f(ip+1,jp-1)+erg1(ip,jp,ip+1,jp-1,ct,data),1,false,0,0,false,0,0,0,0,0);
				
											}
											

										}
						

									
									}


			   
								}
							 }
						}

						skip:

						

						//consider the multiloop closed by i,j
						if ((j-i)>(2*minloop+4)) {
          					//no dangling ends on i-j pair:
							if ((stack.peekatenergy()+wmb->f(i+1,j-1)+data->eparam[5]+data->eparam[10]+penalty(i,j,ct,data)-energy)<=crit) {
								//loop is acceptable:
								stack.push((stack.peekatenergy()+wmb->f(i+1,j-1)+data->eparam[5]+data->eparam[10]+penalty(i,j,ct,data)-energy),false,0,0,true,i+1,j-1,0,wmb->f(i+1,j-1),2,false,0,0,false,0,0,0,0,0);
				
							}

							//i+1 dangles on i-j pair:
							if ((stack.peekatenergy()+erg4(i,j,i+1,1,ct,data,lfce[i+1])+penalty(i,j,ct,data)+
            					wmb->f(i+2,j-1)+ data->eparam[5] + data->eparam[6] + data->eparam[10]-energy)<=crit) {
								//loop is acceptable:
								stack.push((stack.peekatenergy()+erg4(i,j,i+1,1,ct,data,lfce[i+1])+penalty(i,j,ct,data)+
            						wmb->f(i+2,j-1)+ data->eparam[5] + data->eparam[6] + data->eparam[10]-energy),false,0,0,true,i+2,j-1,0,wmb->f(i+2,j-1),2,false,0,0,false,0,0,0,0,0);
								stack.nstack(i+1,i);
							}
							
							if ((stack.peekatenergy()+erg4(i,j,j-1,2,ct,data,lfce[j-1]) + penalty(i,j,ct,data) +
            					wmb->f(i+1,j-2) + data->eparam[5] + data->eparam[6] + data->eparam[10]-energy)<=crit) {
								//loop is acceptable:
								stack.push((stack.peekatenergy()+erg4(i,j,j-1,2,ct,data,lfce[j-1]) + penalty(i,j,ct,data) +
            						wmb->f(i+1,j-2) + data->eparam[5] + data->eparam[6] + data->eparam[10]-energy),false,0,0,true,i+1,j-2,0,wmb->f(i+1,j-2),2,false,0,0,false,0,0,0,0,0);
								stack.nstack(j-1,j);
							}
							
							
							//both i+1 and j-1 dangle (a terminal mismatch)
							if ((stack.peekatenergy()+data->tstkm[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]]+
									checknp(lfce[i+1],lfce[j-1])+
									wmb->f(i+2,j-2) + data->eparam[5] + data->eparam[6] + data->eparam[6]+ data->eparam[10]
									+penalty(i,j,ct,data)-energy)<=crit) {
								//loop is acceptable:
								stack.push((stack.peekatenergy()+data->tstkm[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]]+
									checknp(lfce[i+1],lfce[j-1])+
									wmb->f(i+2,j-2) + data->eparam[5] + data->eparam[6] + data->eparam[6]+ data->eparam[10]
									+penalty(i,j,ct,data)-energy),false,0,0,true,i+2,j-2,0,wmb->f(i+2,j-2),2,false,0,0,false,0,0,0,0,0);
								stack.nstack(i+1,i,j-1,j);
							}
							
		
							//consider the coaxial stacking of a helix from i to j onto helix i+1 or i+2 to ip:
			
							for (ip=i+1;(ip<j);ip++) {
								//first consider flush stacking
												
								if ((stack.peekatenergy()+penalty(i,j,ct,data)+v->f(i+1,ip)+
									penalty(i+1,ip,ct,data)+data->eparam[5]
									+data->eparam[10]+data->eparam[10]+(min(w->f(ip+1,j-1),wmb->f(ip+1,j-1)))+ergcoax(j,i,i+1,ip,0,ct,data)-energy)<=crit) {

									//loop is acceptable:
									stack.push((stack.peekatenergy()+penalty(i,j,ct,data)+v->f(i+1,ip)+
										penalty(i+1,ip,ct,data)+data->eparam[5]
										+data->eparam[10]+data->eparam[10]+(min(w->f(ip+1,j-1),wmb->f(ip+1,j-1)))+ergcoax(j,i,i+1,ip,0,ct,data)-energy),true,i+1,ip,true,i+1,ip,0,v->f(i+1,ip),1,false,0,0,true,ip+1,j-1,0,min(w->f(ip+1,j-1),wmb->f(ip+1,j-1)),3);
									stack.nstack(j,ip,ip,j);
								}

	
								if((mod[i+1]||mod[ip])&&inc[ct->numseq[i+2]][ct->numseq[ip-1]]&&notgu(i+1,ip,ct)&&!(fce->f(i+1,ip)&SINGLE)) {

									if ((stack.peekatenergy()+penalty(i,j,ct,data)+v->f(i+2,ip-1)+
										penalty(i+1,ip,ct,data)+data->eparam[5]
										+data->eparam[10]+data->eparam[10]+(min(w->f(ip+1,j-1),wmb->f(ip+1,j-1)))+ergcoax(j,i,i+1,ip,0,ct,data)
										+erg1(i+1,ip,i+2,ip-1,ct,data)-energy)<=crit) {

										//loop is acceptable:
										stack.push((stack.peekatenergy()+penalty(i,j,ct,data)+v->f(i+1,ip)+
											penalty(i+1,ip,ct,data)+data->eparam[5]
											+data->eparam[10]+data->eparam[10]+(min(w->f(ip+1,j-1),wmb->f(ip+1,j-1)))+ergcoax(j,i,i+1,ip,0,ct,data)-energy),true,i+1,ip,true,i+1,ip,0,v->f(i+2,ip-1)+erg1(i+1,ip,i+2,ip-1,ct,data),1,false,0,0,true,ip+1,j-1,0,min(w->f(ip+1,j-1),wmb->f(ip+1,j-1)),3);
										stack.nstack(j,ip,ip,j);
									}


								

								}
				


									
								if (ip+2<j-1) {
								//now consider an intervening nuc
								if ((ip+2<j-1))


									if ((stack.peekatenergy()+penalty(i,j,ct,data)+v->f(i+2,ip)+
										penalty(i+2,ip,ct,data)+data->eparam[5]
										+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]+
										(min(w->f(ip+2,j-1),wmb->f(ip+2,j-1)))
										+ergcoax(j,i,i+2,ip,ip+1,ct,data)+checknp(lfce[i+1],lfce[ip+1])-energy)<=crit) {

										//loop is acceptable:
										stack.push((stack.peekatenergy()+penalty(i,j,ct,data)+v->f(i+2,ip)+
											penalty(i+2,ip,ct,data)+data->eparam[5]
											+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]+
											(min(w->f(ip+2,j-1),wmb->f(ip+2,j-1)))
											+ergcoax(j,i,i+2,ip,ip+1,ct,data)+checknp(lfce[i+1],lfce[ip+1])-energy),true,i+2,ip,true,i+2,ip,0,v->f(i+2,ip),1,false,0,0,true,ip+2,j-1,0,min(w->f(ip+2,j-1),wmb->f(ip+2,j-1)),3);
										stack.nstack(j,ip+1,ip+1,ip);
									}




									if((mod[i+2]||mod[ip])&&inc[ct->numseq[i+3]][ct->numseq[ip-1]]&&notgu(i+2,ip,ct)
										&&!(fce->f(i+2,ip)&SINGLE)) {

										if ((stack.peekatenergy()+penalty(i,j,ct,data)+v->f(i+3,ip-1)+
											penalty(i+2,ip,ct,data)+data->eparam[5]
											+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]+
											(min(w->f(ip+2,j-1),wmb->f(ip+2,j-1)))
											+ergcoax(j,i,i+2,ip,ip+1,ct,data)
											+erg1(i+2,ip,i+3,ip-1,ct,data)+checknp(lfce[i+1],lfce[ip+1])-energy)<=crit) {

											//loop is acceptable:
											stack.push((stack.peekatenergy()+penalty(i,j,ct,data)+v->f(i+3,ip-1)+
												penalty(i+2,ip,ct,data)+data->eparam[5]
												+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]+
												(min(w->f(ip+2,j-1),wmb->f(ip+2,j-1)))
												+ergcoax(j,i,i+2,ip,ip+1,ct,data)
												+erg1(i+2,ip,i+3,ip-1,ct,data)+checknp(lfce[i+1],lfce[ip+1])-energy),true,i+2,ip,true,i+2,ip,0,v->f(i+3,ip-1)+erg1(i+2,ip,i+3,ip-1,ct,data),1,false,0,0,true,ip+2,j-1,0,min(w->f(ip+2,j-1),wmb->f(ip+2,j-1)),3);
											stack.nstack(j,ip+1,ip+1,ip);
										}


									}
		


									if (ip+1<j-2) {

										if ((stack.peekatenergy()+penalty(i,j,ct,data)+v->f(i+2,ip)+
											penalty(i+2,ip,ct,data)+data->eparam[5]
											+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]+
											(min(w->f(ip+1,j-2),wmb->f(ip+1,j-2)))
											+ergcoax(j,i,i+2,ip,j-1,ct,data)+checknp(lfce[i+1],lfce[j-1])-energy)<=crit) {

											//loop is acceptable:
											stack.push((stack.peekatenergy()+penalty(i,j,ct,data)+v->f(i+2,ip)+
												penalty(i+2,ip,ct,data)+data->eparam[5]
												+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]+
												(min(w->f(ip+1,j-2),wmb->f(ip+1,j-2)))
												+ergcoax(j,i,i+2,ip,j-1,ct,data)+checknp(lfce[i+1],lfce[j-1])-energy),true,i+2,ip,true,i+2,ip,0,v->f(i+2,ip),1,false,0,0,true,ip+1,j-2,0,min(w->f(ip+1,j-2),wmb->f(ip+1,j-2)),3);
											stack.nstack(j,j-1,j-1,ip);
										}

										if((mod[i+2]||mod[ip])&&inc[ct->numseq[i+3]][ct->numseq[ip-1]]&&notgu(i+2,ip,ct)
											&&!(fce->f(i+2,ip)&SINGLE)) {

											if ((stack.peekatenergy()+penalty(i,j,ct,data)+v->f(i+3,ip-1)+
												penalty(i+2,ip,ct,data)+data->eparam[5]
												+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]+
												(min(w->f(ip+1,j-2),wmb->f(ip+1,j-2)))
												+ergcoax(j,i,i+2,ip,j-1,ct,data)
												+erg1(i+2,ip,i+3,ip-1,ct,data)+checknp(lfce[i+1],lfce[j-1])-energy)<=crit) {

												//loop is acceptable:
												stack.push((stack.peekatenergy()+penalty(i,j,ct,data)+v->f(i+3,ip-1)+
													penalty(i+2,ip,ct,data)+data->eparam[5]
													+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]+
													(min(w->f(ip+1,j-2),wmb->f(ip+1,j-2)))
													+ergcoax(j,i,i+2,ip,j-1,ct,data)
													+erg1(i+2,ip,i+3,ip-1,ct,data)+checknp(lfce[i+1],lfce[j-1])-energy),true,i+2,ip,true,i+2,ip,0,v->f(i+3,ip-1)+erg1(i+2,ip,i+3,ip-1,ct,data),1,false,0,0,true,ip+1,j-2,0,min(w->f(ip+1,j-2),wmb->f(ip+1,j-2)),3);
												stack.nstack(j,j-1,j-1,ip);
											}


										}

									}

				


								}



							}

							//consider the coaxial stacking of a helix from i to j onto helix ip to j-2 or j-1:
							for (ip=j-1;ip>i;ip--) {


								//conditions guarantee that the coaxial stacking isn't considering an exterior loop
								
								
								//first consider flush stacking
								if ((stack.peekatenergy()+penalty(i,j,ct,data)+v->f(ip,j-1)+
									penalty(j-1,ip,ct,data)+data->eparam[5]
									+data->eparam[10]+data->eparam[10]+(min(w->f(i+1,ip-1),wmb->f(i+1,ip-1)))+ergcoax(ip,j-1,j,i,0,ct,data)-energy)<=crit) {

									//loop is acceptable:
									stack.push((stack.peekatenergy()+penalty(i,j,ct,data)+v->f(ip,j-1)+
										penalty(j-1,ip,ct,data)+data->eparam[5]
										+data->eparam[10]+data->eparam[10]+(min(w->f(i+1,ip-1),wmb->f(i+1,ip-1)))+ergcoax(ip,j-1,j,i,0,ct,data)-energy),
										true,ip,j-1,true,ip,j-1,0,v->f(ip,j-1),1,false,0,0,true,i+1,ip-1,0,min(w->f(i+1,ip-1),wmb->f(i+1,ip-1)),3);
									stack.nstack(ip,i,i,ip);
								}


								if((mod[ip]||mod[j-1])&&inc[ct->numseq[ip+1]][ct->numseq[j-2]]&&notgu(ip,j-1,ct)&&!(fce->f(ip,j-1)&SINGLE)) {
							
									if ((stack.peekatenergy()+penalty(i,j,ct,data)+v->f(ip+1,j-2)+
										penalty(j-1,ip,ct,data)+data->eparam[5]
										+data->eparam[10]+data->eparam[10]+(min(w->f(i+1,ip-1),wmb->f(i+1,ip-1)))
										+ergcoax(ip,j-1,j,i,0,ct,data)
										+erg1(ip,j-1,ip+1,j-2,ct,data)-energy)<=crit) {

										//loop is acceptable:
										stack.push((stack.peekatenergy()+penalty(i,j,ct,data)+v->f(ip+1,j-2)+
											penalty(j-1,ip,ct,data)+data->eparam[5]
											+data->eparam[10]+data->eparam[10]+(min(w->f(i+1,ip-1),wmb->f(i+1,ip-1)))
											+ergcoax(ip,j-1,j,i,0,ct,data)
											+erg1(ip,j-1,ip+1,j-2,ct,data)-energy),
											true,ip,j-1,true,ip,j-1,0,v->f(ip+1,j-2)+erg1(ip,j-1,ip+1,j-2,ct,data),1,false,0,0,true,i+1,ip-1,0,min(w->f(i+1,ip-1),wmb->f(i+1,ip-1)),3);
										stack.nstack(ip,i,i,ip);
									}
									

								}

				

		

								
								//now consider an intervening nuc
									
								if (ip-2>i+1) {

									if ((stack.peekatenergy()+penalty(i,j,ct,data)+v->f(ip,j-2)+
										penalty(j-2,ip,ct,data)+data->eparam[5]
										+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]+
										(min(w->f(i+1,ip-2),wmb->f(i+1,ip-2)))
										+ergcoax(ip,j-2,j,i,ip-1,ct,data)+checknp(lfce[j-1],lfce[ip-1])-energy)<=crit) {

										//loop is acceptable:
										stack.push((stack.peekatenergy()+penalty(i,j,ct,data)+v->f(ip,j-2)+
											penalty(j-2,ip,ct,data)+data->eparam[5]
											+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]+
											(min(w->f(i+1,ip-2),wmb->f(i+1,ip-2)))
											+ergcoax(ip,j-2,j,i,ip-1,ct,data)+checknp(lfce[j-1],lfce[ip-1])-energy),
											true,ip,j-2,true,ip,j-2,0,v->f(ip,j-2),1,false,0,0,true,i+1,ip-2,0,min(w->f(i+1,ip-2),wmb->f(i+1,ip-2)),3);
										stack.nstack(ip,ip-1,ip-1,i);
									}
									

									if((mod[ip]||mod[j-2])&&inc[ct->numseq[ip+1]][ct->numseq[j-3]]&&notgu(ip,j-2,ct)&&!(fce->f(ip,j-2)&SINGLE)) {
										
										if ((stack.peekatenergy()+penalty(i,j,ct,data)+v->f(ip+1,j-3)+
											penalty(j-2,ip,ct,data)+data->eparam[5]
											+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]+
											(min(w->f(i+1,ip-2),wmb->f(i+1,ip-2)))
											+ergcoax(ip,j-2,j,i,ip-1,ct,data)
											+erg1(ip,j-2,ip+1,j-3,ct,data)+checknp(lfce[j-1],lfce[ip-1])-energy)<=crit) {

												//loop is acceptable:
												stack.push((stack.peekatenergy()+penalty(i,j,ct,data)+v->f(ip+1,j-3)+
													penalty(j-2,ip,ct,data)+data->eparam[5]
													+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]+
													(min(w->f(i+1,ip-2),wmb->f(i+1,ip-2)))
													+ergcoax(ip,j-2,j,i,ip-1,ct,data)
													+erg1(ip,j-2,ip+1,j-3,ct,data)+checknp(lfce[j-1],lfce[ip-1])-energy),
													true,ip,j-2,true,ip,j-2,0,v->f(ip+1,j-3)+erg1(ip,j-2,ip+1,j-3,ct,data),1,false,0,0,true,i+1,ip-2,0,min(w->f(i+1,ip-2),wmb->f(i+1,ip-2)),3);
												stack.nstack(ip,ip-1,ip-1,i);
										}
										
									}
								}



								if ((ip-1>i+2)) {

									if ((stack.peekatenergy()+penalty(i,j,ct,data)+v->f(ip,j-2)+
										penalty(j-2,ip,ct,data)+data->eparam[5]
										+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]+
										(min(w->f(i+2,ip-1),wmb->f(i+2,ip-1)))
										+ergcoax(ip,j-2,j,i,i+1,ct,data)+checknp(lfce[j-1],lfce[i+1])-energy)<=crit) {

										//loop is acceptable:
										stack.push((stack.peekatenergy()+penalty(i,j,ct,data)+v->f(ip,j-2)+
											penalty(j-2,ip,ct,data)+data->eparam[5]
											+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]+
											(min(w->f(i+2,ip-1),wmb->f(i+2,ip-1)))
											+ergcoax(ip,j-2,j,i,i+1,ct,data)+checknp(lfce[j-1],lfce[i+1])-energy),
											true,ip,j-2,true,ip,j-2,0,v->f(ip,j-2),1,false,0,0,true,i+2,ip-1,0,min(w->f(i+2,ip-1),wmb->f(i+2,ip-1)),3);
										stack.nstack(ip,i+1,i+1,i);
									}

									
									if((mod[ip]||mod[j-2])&&inc[ct->numseq[ip+1]][ct->numseq[j-3]]&&notgu(ip,j-2,ct)
										&&!(fce->f(ip,j-2)&SINGLE)) {

										if ((stack.peekatenergy()+penalty(i,j,ct,data)+v->f(ip+1,j-3)+
											penalty(j-2,ip,ct,data)+data->eparam[5]
											+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]+
											(min(w->f(i+2,ip-1),wmb->f(i+2,ip-1)))
											+ergcoax(ip,j-2,j,i,i+1,ct,data)
											+erg1(ip,j-2,ip+1,j-3,ct,data)+checknp(lfce[j-1],lfce[i+1])-energy)<=crit) {

											//loop is acceptable:
											stack.push((stack.peekatenergy()+penalty(i,j,ct,data)+v->f(ip+1,j-3)+
												penalty(j-2,ip,ct,data)+data->eparam[5]
												+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]+
												(min(w->f(i+2,ip-1),wmb->f(i+2,ip-1)))
												+ergcoax(ip,j-2,j,i,i+1,ct,data)
												+erg1(ip,j-2,ip+1,j-3,ct,data)+checknp(lfce[j-1],lfce[i+1])-energy),
												true,ip,j-2,true,ip,j-2,0,v->f(ip+1,j-3)+erg1(ip,j-2,ip+1,j-3,ct,data),1,false,0,0,true,i+2,ip-1,0,min(w->f(i+2,ip-1),wmb->f(i+2,ip-1)),3);
											stack.nstack(ip,i+1,i+1,i);
										}


									}
								}					
				
							}

						}


      
	 

						
					}

				}
				//pair==5 -> a branch that has no extra added nucs - v
				//pair==4 -> a single branch - w or wl or v
				//pair==3 -> a single or double branch - wmb or wmbl or w or wl or v or wcoax
				//pair==2 -> a double branch - wmb or wmbl or wcoax
				//pair==6 -> a double branch from a coaxial stack - wcoax
				//pair==8 -> a double branch that can only have 5' unpaired nucs - wmbl or wcoax (not wmb)
				//pair==7 -> a single branch that can only have 5' unpaired nucs - v or wl (not w)

				if (pair==5||pair==4||pair==3||pair==7) {
					
					//look for the single branch
					if ((stack.peekatenergy()+data->eparam[10]+v->f(i,j)+penalty(j,i,ct,data)-energy)<=crit) {

						//loop is acceptable:
						stack.push((stack.peekatenergy()+data->eparam[10]+v->f(i,j)+penalty(j,i,ct,data)-energy),
							true,i,j,true,i,j,0,v->f(i,j),1,false,0,0,false,0,0,0,0,0);
				
					}
					
					if ((mod[i]||mod[j])&&inc[ct->numseq[i+1]][ct->numseq[j-1]]&&notgu(i,j,ct)) {

						if ((stack.peekatenergy()+data->eparam[10]+v->f(i+1,j-1)+penalty(j,i,ct,data)+erg1(i,j,i+1,j-1,ct,data)-energy)<=crit) {

							//loop is acceptable:
							stack.push((stack.peekatenergy()+data->eparam[10]+v->f(i,j)+penalty(j,i,ct,data)-energy),
								true,i,j,true,i,j,0,v->f(i+1,j-1)+erg1(i,j,i+1,j-1,ct,data),1,false,0,0,false,0,0,0,0,0);
				
						}
									

					} 

	
      				//calculate the energy of i stacked onto the pair of i+1,j
					if ((stack.peekatenergy()+v->f(i+1,j)+data->eparam[10]+data->eparam[6]+
         				erg4(j,i+1,i,2,ct,data,lfce[i])+penalty(i+1,j,ct,data)-energy)<=crit) {

						//loop is acceptable:
						stack.push((stack.peekatenergy()+v->f(i+1,j)+data->eparam[10]+data->eparam[6]+
         					erg4(j,i+1,i,2,ct,data,lfce[i])+penalty(i+1,j,ct,data)-energy),
							true,i+1,j,true,i+1,j,0,v->f(i+1,j),1,false,0,0,false,0,0,0,0,0);
						stack.nstack(i,i+1);
				
					}

					
					if ((mod[i+1]||mod[j])&&inc[ct->numseq[i+2]][ct->numseq[j-1]]&&notgu(i+1,j,ct)&&!(fce->f(i+1,j)&SINGLE)) {

						if ((stack.peekatenergy()+v->f(i+2,j-1) + data->eparam[10] +data->eparam[6] +
         					erg4(j,i+1,i,2,ct,data,lfce[i])+penalty(i+1,j,ct,data)
							+erg1(i+1,j,i+2,j-1,ct,data)-energy)<=crit) {

							//loop is acceptable:
							stack.push((stack.peekatenergy()+v->f(i+2,j-1) + data->eparam[10] +data->eparam[6] +
         						erg4(j,i+1,i,2,ct,data,lfce[i])+penalty(i+1,j,ct,data)
								+erg1(i+1,j,i+2,j-1,ct,data)-energy),
								true,i+1,j,true,i+1,j,0,v->f(i+2,j-1)+erg1(i+1,j,i+2,j-1,ct,data),1,false,0,0,false,0,0,0,0,0);
							stack.nstack(i,i+1);
						}
						


					}
         
					
					
					//calculate the energy of j stacked onto the pair of i,j-1
					if (j!=1) {

						if ((stack.peekatenergy()+v->f(i,j-1)+ data->eparam[10] + data->eparam[6] +
         					erg4(j-1,i,j,1,ct,data,lfce[j])+penalty(i,j-1,ct,data)-energy)<=crit) {

							//loop is acceptable:
							stack.push((stack.peekatenergy()+v->f(i,j-1)+ data->eparam[10] + data->eparam[6] +
         					erg4(j-1,i,j,1,ct,data,lfce[j])+penalty(i,j-1,ct,data)-energy),
								true,i,j-1,true,i,j-1,0,v->f(i,j-1),1,false,0,0,false,0,0,0,0,0);
							stack.nstack(j,j-1);
				
						}
         				
						if ((mod[i]||mod[j-1])&&inc[ct->numseq[i+1]][ct->numseq[j-2]]&&notgu(i,j-1,ct)&&!(fce->f(i,j-1)&SINGLE)) {

							if ((stack.peekatenergy()+v->f(i+1,j-2) + data->eparam[10] + data->eparam[6] +
         							erg4(j-1,i,j,1,ct,data,lfce[j])+penalty(i,j-1,ct,data)
									+erg1(i,j-1,i+1,j-2,ct,data)-energy)<=crit) {

								//loop is acceptable:
								stack.push((stack.peekatenergy()+v->f(i+1,j-2) + data->eparam[10] + data->eparam[6] +
         							erg4(j-1,i,j,1,ct,data,lfce[j])+penalty(i,j-1,ct,data)
									+erg1(i,j-1,i+1,j-2,ct,data)-energy),
									true,i,j-1,true,i,j-1,0,v->f(i+1,j-2)+erg1(i,j-1,i+1,j-2,ct,data),1,false,0,0,false,0,0,0,0,0);
								stack.nstack(j,j-1);
							}

						}

         	
					}
					
					
      				//calculate i and j stacked onto the pair of i+1,j-1
					if (j!=1&&!lfce[i]&&!lfce[j]) {
						if ((stack.peekatenergy()+v->f(i+1,j-1) +data->eparam[10] + (data->eparam[6]+data->eparam[6]) +
         					data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
							+penalty(j-1,i+1,ct,data)-energy)<=crit) {

							//loop is acceptable:
							stack.push((stack.peekatenergy()+v->f(i+1,j-1) +data->eparam[10] + (data->eparam[6]+data->eparam[6]) +
         						data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
								+penalty(j-1,i+1,ct,data)-energy),
								true,i+1,j-1,true,i+1,j-1,0,v->f(i+1,j-1),1,false,0,0,false,0,0,0,0,0);
							stack.nstack(j,j-1,i,i+1);
						}

         				
						if ((mod[i+1]||mod[j-1])&&(j-2>0)&&!(fce->f(i+1,j-1)&SINGLE)) {
							if(inc[ct->numseq[i+2]][ct->numseq[j-2]]&&notgu(i+1,j-1,ct)) {

								if ((stack.peekatenergy()+v->f(i+2,j-2) + data->eparam[10] + (data->eparam[6]+data->eparam[6]) +
         							data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
									+penalty(j-1,i+1,ct,data)+erg1(i+1,j-1,i+2,j-2,ct,data)-energy)<=crit) {

									//loop is acceptable:
									stack.push((stack.peekatenergy()+v->f(i+2,j-2) + data->eparam[10] + (data->eparam[6]+data->eparam[6]) +
         								data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
										+penalty(j-1,i+1,ct,data)+erg1(i+1,j-1,i+2,j-2,ct,data)-energy),
										true,i+1,j-1,true,i+1,j-1,0,v->f(i+2,j-2)+erg1(i+1,j-1,i+2,j-2,ct,data),1,false,0,0,false,0,0,0,0,0);
									stack.nstack(j,j-1,i,i+1);
								}

							}
						}
					}
				}
				if (pair==4||pair==3||pair==7) {
					//this is a fragment of a multibranch loop requiring one stem or more



					if (!lfce[i]) {
						if ((stack.peekatenergy()+wl->f(i+1,j)+data->eparam[6]-energy)<=crit) {

							//loop is acceptable:
							stack.push((stack.peekatenergy()+wl->f(i+1,j)+data->eparam[6]-energy),
								false,0,0,true,i+1,j,0,wl->f(i+1,j),7,false,0,0,false,0,0,0,0,0);
							
						}
         					
            					
					}

					
					if (!lfce[j]&&pair!=7) {
             					
						if ((stack.peekatenergy()+w->f(i,j-1) + data->eparam[6]-energy)<=crit) {

							//loop is acceptable:
							stack.push((stack.peekatenergy()+w->f(i,j-1) + data->eparam[6]-energy),
								false,0,0,true,i,j-1,0,w->f(i,j-1),4,false,0,0,false,0,0,0,0,0);
				
						}
               					

					}


					
				}
				if (pair==3||pair==2||pair==6||pair==8) {
					for (ip=i+minloop+1;ip<j-minloop-1;ip++) {
					//first consider flush stacking
		
						if ((stack.peekatenergy()+v->f(i,ip)+v->f(ip+1,j)+penalty(i,ip,ct,data)
							+penalty(ip+1,j,ct,data)+ergcoax(i,ip,ip+1,j,0,ct,data)+data->eparam[10]+data->eparam[10]-energy)<=crit) {

							//loop is acceptable:
							stack.push((stack.peekatenergy()+v->f(i,ip)+v->f(ip+1,j)+penalty(i,ip,ct,data)
							+penalty(ip+1,j,ct,data)+ergcoax(i,ip,ip+1,j,0,ct,data)+data->eparam[10]+data->eparam[10]-energy),
								true,i,ip,true,i,ip,0,v->f(i,ip),1,true,ip+1,j,true,ip+1,j,0,v->f(ip+1,j),1);
							stack.nstack(i,j,j,i);
						}
						


						if ((mod[i]||mod[ip]||mod[ip+1]||mod[j])) {

							if ((mod[i]||mod[ip])&&(mod[ip+1]||mod[j])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]
								&&inc[ct->numseq[ip+2]][ct->numseq[j-1]]&&notgu(i,ip,ct)&&notgu(ip+1,j,ct)
								&&!(fce->f(ip+1,j)&SINGLE)&&!(fce->f(i,ip)&SINGLE)) {

								if ((stack.peekatenergy()+v->f(i+1,ip-1)+v->f(ip+2,j-1)+penalty(i,ip,ct,data)
									+penalty(ip+1,j,ct,data)+ergcoax(i,ip,ip+1,j,0,ct,data)
									+erg1(i,ip,i+1,ip-1,ct,data)+erg1(ip+1,j,ip+2,j-1,ct,data)+data->eparam[10]+data->eparam[10]-energy)<=crit) {

									//loop is acceptable:
									stack.push((stack.peekatenergy()+v->f(i,ip)+v->f(ip+1,j)+penalty(i,ip,ct,data)
										+penalty(ip+1,j,ct,data)+ergcoax(i,ip,ip+1,j,0,ct,data)+data->eparam[10]+data->eparam[10]-energy),
										true,i,ip,true,i,ip,0,v->f(i+1,ip-1)+erg1(i,ip,i+1,ip-1,ct,data),1,true,ip+1,j,true,ip+1,j,0,v->f(ip+2,j-1)+erg1(ip+1,j,ip+2,j-1,ct,data),1);
									stack.nstack(i,j,j,i);
								}
								
							}

							if ((mod[i]||mod[ip])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]&&notgu(i,ip,ct)&&!(fce->f(i,ip)&SINGLE)) {
							
								if ((stack.peekatenergy()+v->f(i+1,ip-1)+v->f(ip+1,j)+penalty(i,ip,ct,data)
									+penalty(ip+1,j,ct,data)+ergcoax(i,ip,ip+1,j,0,ct,data)
									+erg1(i,ip,i+1,ip-1,ct,data)+data->eparam[10]+data->eparam[10]-energy)<=crit) {

									//loop is acceptable:
									stack.push((stack.peekatenergy()+v->f(i+1,ip-1)+v->f(ip+1,j)+penalty(i,ip,ct,data)
									+penalty(ip+1,j,ct,data)+ergcoax(i,ip,ip+1,j,0,ct,data)
									+erg1(i,ip,i+1,ip-1,ct,data)+data->eparam[10]+data->eparam[10]-energy),
										true,i,ip,true,i,ip,0,v->f(i+1,ip-1)+erg1(i,ip,i+1,ip-1,ct,data),1,true,ip+1,j,true,ip+1,j,0,v->f(ip+1,j),1);
									stack.nstack(i,j,j,i);
								}
								
							}

							if ((mod[ip+1]||mod[j])&&inc[ct->numseq[ip+2]][ct->numseq[j-1]]&&notgu(ip+1,j,ct)&&!(fce->f(ip+1,j)&SINGLE)) {


								if ((stack.peekatenergy()+v->f(i,ip)+v->f(ip+2,j-1)+penalty(i,ip,ct,data)
									+penalty(ip+1,j,ct,data)+ergcoax(i,ip,ip+1,j,0,ct,data)
									+erg1(ip+1,j,ip+2,j-1,ct,data)+data->eparam[10]+data->eparam[10]-energy)<=crit) {

									//loop is acceptable:
									stack.push((stack.peekatenergy()+v->f(i,ip)+v->f(ip+2,j-1)+penalty(i,ip,ct,data)
										+penalty(ip+1,j,ct,data)+ergcoax(i,ip,ip+1,j,0,ct,data)
										+erg1(ip+1,j,ip+2,j-1,ct,data)+data->eparam[10]+data->eparam[10]-energy),
										true,i,ip,true,i,ip,0,v->f(i,ip),1,true,ip+1,j,true,ip+1,j,0,v->f(ip+2,j-1)+erg1(ip+1,j,ip+2,j-1,ct,data),1);
									stack.nstack(i,j,j,i);
								}
								

							}


						}
				


						
						if (!lfce[ip+1]&&!lfce[j]) {
							//now consider an intervening mismatch
							if(!lfce[ip+1]&&!lfce[j]) 
								
								
								if ((stack.peekatenergy()+v->f(i,ip)+v->f(ip+2,j-1)+penalty(i,ip,ct,data)
									+penalty(ip+2,j-1,ct,data)+ergcoax(i,ip,ip+2,j-1,j,ct,data)+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]-energy)<=crit) {

									//loop is acceptable:
									stack.push((stack.peekatenergy()+v->f(i,ip)+v->f(ip+2,j-1)+penalty(i,ip,ct,data)
										+penalty(ip+2,j-1,ct,data)+ergcoax(i,ip,ip+2,j-1,j,ct,data)+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]-energy),
										true,i,ip,true,i,ip,0,v->f(i,ip),1,true,ip+2,j-1,true,ip+2,j-1,0,v->f(ip+2,j-1),1);
									stack.nstack(i,j,j,j-1);
								}
							
							
							if (mod[i]||mod[ip]||mod[ip+2]||mod[j-1]) {
								if ((mod[i]||mod[ip])&&(mod[ip+2]||mod[j-1])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]
									&&inc[ct->numseq[ip+3]][ct->numseq[j-2]]&&notgu(i,ip,ct)&&notgu(ip+2,j-1,ct)
									&&!(fce->f(i,ip)&SINGLE)&&!(fce->f(ip+2,j-1)&SINGLE)) {

									if(!lfce[ip+1]&&!lfce[j]) 
										if ((stack.peekatenergy()+v->f(i+1,ip-1)+v->f(ip+3,j-2)+penalty(i,ip,ct,data)
										+penalty(ip+2,j-1,ct,data)+ergcoax(i,ip,ip+2,j-1,j,ct,data)
										+erg1(i,ip,i+1,ip-1,ct,data)+erg1(ip+2,j-1,ip+3,j-2,ct,data)+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]-energy)<=crit) {

											//loop is acceptable:
											stack.push((stack.peekatenergy()+v->f(i+1,ip-1)+v->f(ip+3,j-2)+penalty(i,ip,ct,data)
												+penalty(ip+2,j-1,ct,data)+ergcoax(i,ip,ip+2,j-1,j,ct,data)
												+erg1(i,ip,i+1,ip-1,ct,data)+erg1(ip+2,j-1,ip+3,j-2,ct,data)+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]-energy),
												true,i,ip,true,i,ip,0,v->f(i+1,ip-1)+erg1(i,ip,i+1,ip-1,ct,data),1,true,ip+2,j-1,true,ip+2,j-1,0,v->f(ip+3,j-2)+erg1(ip+2,j-1,ip+3,j-2,ct,data),1);
											stack.nstack(i,j,j,j-1);
										}
										
										
								}

								if ((mod[i]||mod[ip])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]&&notgu(i,ip,ct)&&!(fce->f(i,ip)&SINGLE)) {
							
									if(!lfce[ip+1]&&!lfce[j]) 
										
										
										if ((stack.peekatenergy()+v->f(i+1,ip-1)+v->f(ip+2,j-1)+penalty(i,ip,ct,data)
											+penalty(ip+2,j-1,ct,data)+ergcoax(i,ip,ip+2,j-1,j,ct,data)
											+erg1(i,ip,i+1,ip-1,ct,data)+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]-energy)<=crit) {

											//loop is acceptable:
											stack.push((stack.peekatenergy()+v->f(i+1,ip-1)+v->f(ip+2,j-1)+penalty(i,ip,ct,data)
												+penalty(ip+2,j-1,ct,data)+ergcoax(i,ip,ip+2,j-1,j,ct,data)
												+erg1(i,ip,i+1,ip-1,ct,data)+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]-energy),
												true,i,ip,true,i,ip,0,v->f(i+1,ip-1)+erg1(i,ip,i+1,ip-1,ct,data),1,true,ip+2,j-1,true,ip+2,j-1,0,v->f(ip+2,j-1),1);
											stack.nstack(i,j,j,j-1);
										}
									
								}

								if ((mod[ip+2]||mod[j-1])&&inc[ct->numseq[ip+3]][ct->numseq[j-2]]&&notgu(ip+2,j-1,ct)&&!(fce->f(ip+2,j-1)&SINGLE)) {


									if(!lfce[ip+1]&&!lfce[j]) 
										if ((stack.peekatenergy()+v->f(i,ip)+v->f(ip+3,j-2)+penalty(i,ip,ct,data)
										+penalty(ip+2,j-1,ct,data)+ergcoax(i,ip,ip+2,j-1,j,ct,data)
										+erg1(ip+2,j-1,ip+3,j-2,ct,data)+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]-energy)<=crit) {

											//loop is acceptable:
											stack.push((stack.peekatenergy()+v->f(i,ip)+v->f(ip+3,j-2)+penalty(i,ip,ct,data)
												+penalty(ip+2,j-1,ct,data)+ergcoax(i,ip,ip+2,j-1,j,ct,data)
												+erg1(ip+2,j-1,ip+3,j-2,ct,data)+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]-energy),
												true,i,ip,true,i,ip,0,v->f(i,ip),1,true,ip+2,j-1,true,ip+2,j-1,0,v->f(ip+3,j-2)+erg1(ip+2,j-1,ip+3,j-2,ct,data),1);
											stack.nstack(i,j,j,j-1);
										}
										
										
								}
							}
						}
		
						if(!lfce[i]&&!lfce[ip+1]) {

							if ((stack.peekatenergy()+v->f(i+1,ip)+v->f(ip+2,j)+penalty(i+1,ip,ct,data)
								+penalty(ip+2,j,ct,data)+ergcoax(i+1,ip,ip+2,j,i,ct,data)+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]-energy)<=crit) {

									//loop is acceptable:
									stack.push((stack.peekatenergy()+v->f(i+1,ip)+v->f(ip+2,j)+penalty(i+1,ip,ct,data)
										+penalty(ip+2,j,ct,data)+ergcoax(i+1,ip,ip+2,j,i,ct,data)+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]-energy),
										true,i+1,ip,true,i+1,ip,0,v->f(i+1,ip),1,true,ip+2,j,true,ip+2,j,0,v->f(ip+2,j),1);
									stack.nstack(i+1,i,i,j);
							}
							
							if (mod[i+1]||mod[ip]||mod[ip+2]||mod[j]) {
								if ((mod[i+1]||mod[ip])&&(mod[ip+2]||mod[j])&&inc[ct->numseq[i+2]][ct->numseq[ip-1]]
									&&inc[ct->numseq[ip+3]][ct->numseq[j-1]]&&notgu(i+1,ip,ct)&&notgu(ip+2,j,ct)
									&&!(fce->f(i+1,ip)&SINGLE)&&!(fce->f(ip+2,j)&SINGLE)	) {


									if ((stack.peekatenergy()+v->f(i+2,ip-1)+v->f(ip+3,j-1)+penalty(i+1,ip,ct,data)
										+penalty(ip+2,j,ct,data)+ergcoax(i+1,ip,ip+2,j,i,ct,data)
										+erg1(i+1,ip,i+2,ip-1,ct,data)+erg1(ip+2,j,ip+3,j-1,ct,data)+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]-energy)<=crit) {

										//loop is acceptable:
										stack.push((stack.peekatenergy()+v->f(i+2,ip-1)+v->f(ip+3,j-1)+penalty(i+1,ip,ct,data)
										+penalty(ip+2,j,ct,data)+ergcoax(i+1,ip,ip+2,j,i,ct,data)
										+erg1(i+1,ip,i+2,ip-1,ct,data)+erg1(ip+2,j,ip+3,j-1,ct,data)+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]-energy),
											true,i+1,ip,true,i+1,ip,0,v->f(i+2,ip-1)+erg1(i+1,ip,i+2,ip-1,ct,data),1,true,ip+2,j,true,ip+2,j,0,v->f(ip+3,j-1)+erg1(ip+2,j,ip+3,j-1,ct,data),1);
										stack.nstack(i+1,i,i,j);
									}
									
								}
								if ((mod[i+1]||mod[ip])&&inc[ct->numseq[i+2]][ct->numseq[ip-1]]&&notgu(i+1,ip,ct)&&!(fce->f(i+1,ip)&SINGLE)) {
							
									if ((stack.peekatenergy()+v->f(i+2,ip-1)+v->f(ip+2,j)+penalty(i+1,ip,ct,data)
										+penalty(ip+2,j,ct,data)+ergcoax(i+1,ip,ip+2,j,i,ct,data)
										+erg1(i+1,ip,i+2,ip-1,ct,data)+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]-energy)<=crit) {

										//loop is acceptable:
										stack.push((stack.peekatenergy()+v->f(i+2,ip-1)+v->f(ip+2,j)+penalty(i+1,ip,ct,data)
											+penalty(ip+2,j,ct,data)+ergcoax(i+1,ip,ip+2,j,i,ct,data)
											+erg1(i+1,ip,i+2,ip-1,ct,data)+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]-energy),
											true,i+1,ip,true,i+1,ip,0,v->f(i+2,ip-1)+erg1(i+1,ip,i+2,ip-1,ct,data),1,true,ip+2,j,true,ip+2,j,0,v->f(ip+2,j),1);
										stack.nstack(i+1,i,i,j);
									}
									
								}

								if ((mod[ip+2]||mod[j])&&inc[ct->numseq[ip+3]][ct->numseq[j-1]]&&notgu(ip+2,j,ct)&&!(fce->f(ip+2,j)&SINGLE)) {
									if ((stack.peekatenergy()+v->f(i+1,ip)+v->f(ip+3,j-1)+penalty(i+1,ip,ct,data)
										+penalty(ip+2,j,ct,data)+ergcoax(i+1,ip,ip+2,j,i,ct,data)
										+erg1(ip+2,j,ip+3,j-1,ct,data)+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]-energy)<=crit) {

										//loop is acceptable:
										stack.push((stack.peekatenergy()+v->f(i+1,ip)+v->f(ip+3,j-1)+penalty(i+1,ip,ct,data)
											+penalty(ip+2,j,ct,data)+ergcoax(i+1,ip,ip+2,j,i,ct,data)
											+erg1(ip+2,j,ip+3,j-1,ct,data)+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]-energy),
											true,i+1,ip,true,i+1,ip,0,v->f(i+1,ip),1,true,ip+2,j,true,ip+2,j,0,v->f(ip+3,j-1)+erg1(ip+2,j,ip+3,j-1,ct,data),1);
										stack.nstack(i+1,i,i,j);
									}

									
								}
							}
						}
						
					}



				}
				if (pair==3||pair==2||pair==8) {
					//this is a multibranch loop fragment requiring or allowing two stems 
					
					
					//search for a bifurcation:
					for (k=i+1;k<j;k++) {
						

						branch= data->eparam[10]+v->f(i,k)+penalty(k,i,ct,data);

						if ((mod[i]||mod[k])&&inc[ct->numseq[i+1]][ct->numseq[k-1]]&&notgu(i,k,ct)) {

								branch=min(branch, data->eparam[10]+v->f(i+1,k-1)+penalty(k,i,ct,data)+erg1(i,k,i+1,k-1,ct,data));			

						} 



						
      					//calculate the energy of i stacked onto the pair of i+1,k

						branch=min(branch, v->f(i+1,k)+data->eparam[10]+data->eparam[6]+
         					erg4(k,i+1,i,2,ct,data,lfce[i])+penalty(i+1,k,ct,data));

						if ((mod[i+1]||mod[k])&&inc[ct->numseq[i+2]][ct->numseq[k-1]]&&notgu(i+1,k,ct)&&!(fce->f(i+1,k)&SINGLE)) {

								branch=min(branch, v->f(i+2,k-1) + data->eparam[10] +data->eparam[6] +
         							erg4(k,i+1,i,2,ct,data,lfce[i])+penalty(i+1,k,ct,data)
									+erg1(i+1,k,i+2,k-1,ct,data));


						}
         
						
						
      					//calculate the energy of k stacked onto the pair of i,k-1
						if (k!=1) {
         					branch=min(branch, v->f(i,k-1)+ data->eparam[10] + data->eparam[6] +
         						erg4(k-1,i,k,1,ct,data,lfce[k])+penalty(i,k-1,ct,data));

							if ((mod[i]||mod[k-1])&&inc[ct->numseq[i+1]][ct->numseq[k-2]]&&notgu(i,k-1,ct)&&!(fce->f(i,k-1)&SINGLE)) {

									branch=min(branch, v->f(i+1,k-2) + data->eparam[10] + data->eparam[6] +
         								erg4(k-1,i,k,1,ct,data,lfce[k])+penalty(i,k-1,ct,data)
										+erg1(i,k-1,i+1,k-2,ct,data));
				


							}

         	
						}
						
						
      					//calculate i and k stacked onto the pair of i+1,k-1
						if (k!=1&&!lfce[i]&&!lfce[k]) {
         					branch=min(branch, v->f(i+1,k-1) +data->eparam[10] + (data->eparam[6]+data->eparam[6]) +
         						data->tstkm[ct->numseq[k-1]][ct->numseq[i+1]][ct->numseq[k]][ct->numseq[i]]
							+penalty(k-1,i+1,ct,data));



							if ((mod[i+1]||mod[k-1])&&(k-2>0)&&!(fce->f(i+1,k-1)&SINGLE)) {
								if(inc[ct->numseq[i+2]][ct->numseq[k-2]]&&notgu(i+1,k-1,ct)) {

									branch=min(branch, v->f(i+2,k-2) + data->eparam[10] + (data->eparam[6]+data->eparam[6]) +
         								data->tstkm[ct->numseq[k-1]][ct->numseq[i+1]][ct->numseq[k]][ct->numseq[i]]
										+penalty(k-1,i+1,ct,data)+erg1(i+1,k-1,i+2,k-2,ct,data));

								}
							}
						}
						


						

						
						if ((stack.peekatenergy()+branch + wl->f(k+1,j)-energy)<=crit) {

								//loop is acceptable:
								stack.push((stack.peekatenergy()+branch + wl->f(k+1,j)-energy),
									false,0,0,true,i,k,0,branch,5,false,0,0,true,k+1,j,0,wl->f(k+1,j),7);
				
						}
						if ((stack.peekatenergy()+wcoax->f(i,k) + wl->f(k+1,j)-energy)<=crit) {

								//loop is acceptable:
								stack.push((stack.peekatenergy()+wcoax->f(i,k) + wl->f(k+1,j)-energy),
									false,0,0,true,i,k,0,wcoax->f(i,k),6,false,0,0,true,k+1,j,0,wl->f(k+1,j),7);
				
						}
						if ((stack.peekatenergy()+branch + wmbl->f(k+1,j)-energy)<=crit) {

								//loop is acceptable:
								stack.push((stack.peekatenergy()+branch + wmbl->f(k+1,j)-energy),
									false,0,0,true,i,k,0,branch,5,false,0,0,true,k+1,j,0,wmbl->f(k+1,j),8);
				
						}
						if ((stack.peekatenergy()+wcoax->f(i,k) + wmbl->f(k+1,j)-energy)<=crit) {

								//loop is acceptable:
								stack.push((stack.peekatenergy()+wcoax->f(i,k) + wmbl->f(k+1,j)-energy),
									false,0,0,true,i,k,0,wcoax->f(i,k),6,false,0,0,true,k+1,j,0,wmbl->f(k+1,j),8);
				
						}
								

							/*else {
							
							if ((stack.peekatenergy()+wl->f(i,k) + wl->f(k+1,j)-energy)<=crit) {

									//loop is acceptable:
									stack.push((stack.peekatenergy()+wl->f(i,k) + wl->f(k+1,j)-energy),
										false,0,0,true,i,k,0,wl->f(i,k),5,false,0,0,true,k+1,j,0,wl->f(k+1,j),4);
				
							}
							if ((stack.peekatenergy()+wcoax->f(i,k) + wl->f(k+1,j)-energy)<=crit) {

									//loop is acceptable:
									stack.push((stack.peekatenergy()+wcoax->f(i,k) + wl->f(k+1,j)-energy),
										false,0,0,true,i,k,0,wcoax->f(i,k),6,false,0,0,true,k+1,j,0,wl->f(k+1,j),4);
				
							}
							if ((stack.peekatenergy()+wl->f(i,k) + wmbl->f(k+1,j)-energy)<=crit) {

									//loop is acceptable:
									stack.push((stack.peekatenergy()+wl->f(i,k) + wmbl->f(k+1,j)-energy),
										false,0,0,true,i,k,0,wl->f(i,k),5,false,0,0,true,k+1,j,0,wmbl->f(k+1,j),2);
				
							}
							if ((stack.peekatenergy()+wcoax->f(i,k) + wmbl->f(k+1,j)-energy)<=crit) {

									//loop is acceptable:
									stack.push((stack.peekatenergy()+wcoax->f(i,k) + wmbl->f(k+1,j)-energy),
										false,0,0,true,i,k,0,wcoax->f(i,k),6,false,0,0,true,k+1,j,0,wmbl->f(k+1,j),2);
				
							}


						}*/
							
							
							

						
					}

			
			

					
					if (!lfce[i]) 
						if ((stack.peekatenergy()+wmbl->f(i+1,j)+data->eparam[6]-energy)<=crit) {

							//loop is acceptable:
							stack.push((stack.peekatenergy()+wmbl->f(i+1,j)+data->eparam[6]-energy),
								false,0,0,true,i+1,j,0,wmbl->f(i+1,j),8,false,0,0,false,0,0,0,0,0);
				
						}
						
					
					
					if (!lfce[j]&&pair!=8) 
						
						if ((stack.peekatenergy()+wmb->f(i,j-1)+data->eparam[6]-energy)<=crit) {

							//loop is acceptable:
							stack.push((stack.peekatenergy()+wmb->f(i,j-1)+data->eparam[6]-energy),
								false,0,0,true,i,j-1,0,wmb->f(i,j-1),2,false,0,0,false,0,0,0,0,0);
				
						}
						
						

				}
				stack.flushbullpen();

			}
		
		
		}
		//no refinements left to be made with the current structure:

		
	//ct->numofstructures=current;
	ct->AddStructure();

	if (!writect) {//not write ct, just output the energy
		//the size of ctenery cannot be too large to overflow the memory 4x100M
		if (current >= 100000000)	{
			cout<< "Stopped: ctenergy[] is overflowed, as the structure number from alltrace() is larger than Max Structure number: ";
			cout<< current;
			cout <<"\n";
			return;
		}
		//increase the size of ctenergy array if the size was exceeded, add 4x10M each time
		if ( current >= energy_size ) {
			energy_size+=10000000;
			tempenergy=new int[energy_size+1];
			for (ip=0;ip<=energy_size-10000000;ip++)	{	tempenergy[ip]=ctenergy[ip];	}
			delete ctenergy;
			ctenergy=tempenergy;
		}
		//get the energy
		ctenergy[current]= stack.peekatenergy();
		current++;
		stack.pull();
		}
	else {//for now, write to ct:	

		//ct->checknumberofstructures();//make sure there is enough space allocated in ct for this structure
		for (ip=1;ip<=ct->GetSequenceLength();ip++) {
			ct->SetPair(ip,stack.readpair(ip),current);	
		}
		//if (ct->stacking) {
		//	for (ip=1;ip<=ct->GetSequenceLength();ip++) {
		//		ct->basepr[current][ip+ct->GetSequenceLength()]=stack.readstacking(ip);	
		//	}

		//}
		ct->SetEnergy(current,stack.peekatenergy());
		ct->SetCtLabel(ct->GetSequenceLabel(),current);
		
		

		current++;
		//if there are NMR constraints, check them here
		if (ct->min_gu>0||ct->min_g_or_u>0||ct->nneighbors>0||ct->nregion>0) {
			passed = true;//assume the structure will pass the constraints
			if (ct->min_gu>0) {
				//count the GU pairs
				count = 0;
				for (pos=1;pos<=ct->GetSequenceLength();pos++) {
					if (ct->numseq[pos]==3) {
						if (ct->GetPair(pos,current-1)>0) {
							if (ct->numseq[ct->GetPair(pos,current-1)]==4) count++;
						}

					}
				}
				if (ct->min_gu>count) passed=false;

			}
			if (ct->min_g_or_u>0&&passed) {
				//count the Gs and Us paired
				count = 0;
				for (pos=1;pos<=ct->GetSequenceLength();pos++) {
					if (ct->numseq[pos]==3||ct->numseq[pos]==4) {
						if (ct->GetPair(pos,current-1)>0) count++;

					}
				}
				if (ct->min_g_or_u>count) passed=false;


			}

			if (ct->nneighbors>0&&passed) {
				//check the neighbors
				for (ip=0;ip<ct->nneighbors&&passed;ip++) {
					found = false;
					pos = 1;
					
					while (pos<=ct->GetSequenceLength()&&!found) {
						jp = 0;
						pos2=pos;
						while ((ct->neighbors[ip][jp]==ct->numseq[pos2]||ct->neighbors[ip][jp]==ct->numseq[ct->GetPair(pos2,current-1)])
							&&ct->GetPair(pos2,current-1)>0&&!found&&pos2<=ct->GetSequenceLength()) {
							pos2++;
							jp++;
							if (ct->neighbors[ip][jp]==0) {
								found = true;
							}
						}
						pos++;
					}
					if (!found) {
						passed = false;
					}

				}

			}
			for (kp=0;kp<ct->nregion&&passed;kp++) {
				
				
				passed = true;//assume the structure will pass the constraints
				if (ct->rmin_gu[kp]>0) {
					//count the GU pairs
					count = 0;
					for (pos=ct->start[kp];pos<=ct->stop[kp];pos++) {
						if (ct->numseq[pos]==3) {
							if (ct->GetPair(pos,current-1)>=ct->start[kp]&&ct->GetPair(pos,current-1)<=ct->stop[kp]) {
								if (ct->numseq[ct->GetPair(pos,current-1)]==4) count++;
							}

						}
					}
					if (ct->rmin_gu[kp]>count) passed=false;

				}
				if (ct->rmin_g_or_u[kp]>0&&passed) {
					//count the Gs and Us paired
					count = 0;
					for (pos=ct->start[kp];pos<=ct->stop[kp];pos++) {
						if (ct->numseq[pos]==3||ct->numseq[pos]==4) {
							if (ct->GetPair(pos,current-1)>0) count++;

						}
					}
					if (ct->rmin_g_or_u[kp]>count) passed=false;


				}

				if (ct->rnneighbors[kp]>0&&passed) {
					//check the neighbors
					for (ip=0;ip<ct->rnneighbors[kp]&&passed;ip++) {
						found = false;
						pos = ct->start[kp];
					
						while (pos<=ct->stop[kp]&&!found) {
							jp = 0;
							pos2=pos;
							while ((((ct->rneighbors[kp][ip][jp]==ct->numseq[pos2])&&ct->GetPair(pos2,current-1)>0)||(ct->rneighbors[kp][ip][jp]==ct->numseq[ct->GetPair(pos2,current-1)]&&ct->GetPair(pos2,current-1)>=ct->start[kp]&&ct->GetPair(pos2,current-1)<=ct->stop[kp]))&&!found&&pos2<=ct->stop[kp]) {
								pos2++;
								jp++;
								if (ct->rneighbors[kp][ip][jp]==0) {
									found = true;
								}
							}
							pos++;
						}
						if (!found) {
							passed = false;
						}

					}

				}


			}

			if (!passed) {
				current--;
				ct->RemoveLastStructure();
			}

		}
		if (ct->nmicroarray>0) {
			//if there are microarray constraints, check them here
			passed = true;

			for (kp=0;kp<=ct->nmicroarray&&passed;kp++) {
				count = 0;
				for (pos=ct->microstart[kp];pos<=ct->microstop[kp];pos++) {
					if (ct->GetPair(pos,current-1)==0) count++;
					

				}
				if (count<ct->microunpair[kp]) passed = false;

			}

			if (!passed) {
				current--;
				ct->RemoveLastStructure();
			}


		}

		//remove that structure from stack
		stack.pull();

		if (ctname!=NULL) {
			//structures should be written to disk as they are found

			//check to see if current==2, if not, the structure was rejected for some reason.
			if (current==2) {
				current =1;
				*out << ct->GetEnergy(1)<<"\n";
				for (kp=1;kp<=ct->GetSequenceLength();j++) {
					*out << ct->GetPair(kp)<<"\n";
				}
				ct->RemoveLastStructure();

			}

		}
	}
	}
	
	if (ctname!=NULL) {
		out->close();
		delete out;
	}
	else if(writect) {
		ct->sort();  //only sort structures as they are stored in ct
	}
}


//This function calcuates the arrays for tracing back all secondary structures 
void alltrace(bool writect,int *ctenergy,structure* ct,datatable* data, short percentdelta, short absolutedelta, ProgressHandler* update, char* save)
{

		
int d,ip,jp,ii,jj;
register int i,j;
int k,l,m,n,o,p;
bool *lfce,*mod;//[maxbases+1][maxbases+1];
int before,after;
register short e;
short *w5,*wca;
atDynProgArray *w2,*wmb2;
register int number;
register short rarray;
short branch;

vector< vector<bool> > inc = data->pairing;

if (ct->GetThermodynamicDataTable()!=data) {
	cerr << "In alltrace ("<<__FILE__<<"): The structure's datatable does not match the passed-in datatable. This can cause problems with IsNuc etc." << endl;
	ct->SetThermodynamicDataTable(data);
}
//array *v,*w;
//inc is an array that saves time by showing which bases can pair before
//	erg is called

//number is the number of bases being folded
//v[i][j] is the best energy for subsequence i to j when i and j are paired
//	for i<j<n, it is the interior framgment between nucleotides i and j
//	for i<n<j, it is the exterior structure fragmnet that is 5' to j-n and 3' to i
//w[i][j] is the best energy for subsequence i to j

number = (ct->GetSequenceLength());//place the number of bases in a registered integer

//scaling is a per nucleotide scale factor for which W and V are divided
//This is necessary to keep the doubles from overflowing:



//allocate space for the v and w arrays:
atDynProgArray w(number);
atDynProgArray v(number);
atDynProgArray wmb(number);
atDynProgArray wl(number);
atDynProgArray wmbl(number);
atDynProgArray wcoax(number);
forceclass fce(number);




//add a second array for intermolecular folding:

if (ct->intermolecular) {
	w2 = new atDynProgArray(number);
	wmb2 = new atDynProgArray(number);
	
   
	

}

else {

	wmb2=NULL;
	w2=NULL;

}


   	


   lfce = new bool [2*number+1];
   mod = new bool [2*number+1];

   for (i=0;i<=2*number;i++) {
	   lfce[i] = false;
	   mod[i] = false;
   }

   for (i=0;i<ct->GetNumberofModified();++i) {

		if (ct->GetModified(i)!=1&&ct->GetModified(i)!=ct->GetSequenceLength()) {
			mod[ct->GetModified(i)]=true;
			mod[ct->GetModified(i)+ct->GetSequenceLength()]=true;
		}
   }




   w5 = new short [number+1];
   
   wca = new short [number+1];
   
   

   
   		w5[0] = 0;
	
	for (i=0;i<=number;i++) {
		wca[i] = 0;
   }
   



	force(ct,&fce,lfce);


//This is the fill routine:



for (j=1;j<=(number);j++) {

	if (((j%10)==0)&&update) update->update((100*j)/(ct->GetSequenceLength()));
	
	for (i=j;i>=1;i--) {

		
	rarray=INFINITE_ENERGY;


   if (ct->templated) {
   	if (i>ct->GetSequenceLength()) ii = i - ct->GetSequenceLength();
      else ii = i;
      if (j>ct->GetSequenceLength()) jj = j - ct->GetSequenceLength();
      else jj = j;
      if (jj<ii) {
         p = jj;
      	jj = ii;
         ii = p;
      }
   	if (!ct->tem[jj][ii]) goto sub2;
   }

    //Compute v[i][j], the minimum energy of the substructure from i to j,
	//inclusive, where i and j are base paired
	if (fce.f(i,j)&SINGLE) {
		//i or j is forced single-stranded
		//v.f(i,j) = 0;
		goto sub2;
	}
	if (fce.f(i,j)&NOPAIR) {
		//i or j is forced into a pair elsewhere
   		//v.f(i,j)= 0;
		
		goto sub2;
   }

   
	if ((j-i)<=minloop) goto sub3;
   
   
	   
   

   if (inc[ct->numseq[i]][ct->numseq[j]]==0) {
	   //v.f(i,j)= 0;
	   goto sub2;
   }

   	
   //force u's into gu pairs
   for (ip=0;ip<ct->GetNumberofGU();++ip) {
   	if (ct->GetGUpair(ip)==i) {
       	if (ct->numseq[j]!=3) {
         	//v.f(i,j) = 0;
         	goto sub2;
         }
	}
      
      else if (ct->GetGUpair(ip)==j) {
       	if (ct->numseq[i]!=3) {
         	//v.f(i,j) = 0;
         	goto sub2;
         }
	  }
      

   }


	//now check to make sure that this isn't an isolated pair:
	//	(consider a pair separated by a bulge as not! stacked)

	//before = 0 if a stacked pair cannot form 5' to i
   before =0;
	if ((i>1&&j!=number)) {
		
			before = inc[ct->numseq[i-1]][ct->numseq[j+1]];
		
	}
	

	//after = 0 if a stacked pair cannot form 3' to i
	if ((j-i)>minloop+2) {
		after = inc[ct->numseq[i+1]][ct->numseq[j-1]];

	}
	else after = 0;

	//if there are no stackable pairs to i.j then don't allow a pair i,j
	if ((before==0)&&(after==0)) {
		//v.f(i,j)= 0;
		goto sub2;
	}
	
	


   	//Perhaps i and j close a hairpin:
      rarray=erg3(i,j,ct,data,fce.f(i,j));

      if ((j-i-1)>=(minloop+2))
      	//Perhaps i,j stacks over i+1,j-1
		if (!mod[i]&&!mod[j])  //make sure this is not a site of chemical modification
			rarray=min(rarray,erg1(i,j,i+1,j-1,ct,data)+v.f(i+1,j-1));
		else {
			//allow G-U to be modified or a pair next to a G-U to be modified
			if ((ct->numseq[i]==3&&ct->numseq[j]==4)||(ct->numseq[i]==4&&ct->numseq[j]==3)) {
				rarray=min(rarray,erg1(i,j,i+1,j-1,ct,data)+v.f(i+1,j-1));	

			}
			else if ((ct->numseq[i+1]==3&&ct->numseq[j-1]==4)||(ct->numseq[i+1]==4&&ct->numseq[j-1]==3)) {

				rarray=min(rarray,erg1(i,j,i+1,j-1,ct,data)+v.f(i+1,j-1));

			}
			else if (i-1>0&&j+1<2*number) {
				if ((ct->numseq[i-1]==3&&ct->numseq[j+1]==4)||(ct->numseq[i-1]==4&&ct->numseq[j+1]==3)) {

					rarray=min(rarray,erg1(i,j,i+1,j-1,ct,data)+v.f(i+1,j-1));
				
				}

			}

		}

		
      //Perhaps i,j closes an interior or bulge loop, enumerate all possibilities
      	//possibility
		if ((j-i-1)>=(minloop+3)) {
      	for (d=(j-i-3);d>=1;d--) {
         	for (ip=(i+1);ip<=(j-1-d);ip++) {
            	jp = d+ip;
               if ((j-i-2-d)>(data->eparam[7])) goto sub1;
               if (abs(ip-i+jp-j)<=(data->eparam[8])) {
               	
                  		rarray=min(rarray,erg2(i,j,ip,jp,ct,data,fce.f(i,ip),fce.f(jp,j))+v.f(ip,jp));

						if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]&&notgu(ip,jp,ct)&&!(fce.f(ip,jp)&SINGLE)) {
							//i or j is modified
							rarray=min(rarray,erg2(i,j,ip,jp,ct,data,fce.f(i,ip),fce.f(jp,j))+
                  				v.f(ip+1,jp-1)+erg1(ip,jp,ip+1,jp-1,ct,data));

						}
						



			   }
			}


			   
		}
		}
	  

      //Perhaps i,j closes a multibranch or exterior loop, enumerate all possibilities


      sub1:
		
  
      if ((j-i-1)>=(2*minloop+4)) {
     		


		//consider the multiloop closed by i,j
         if ((j-i)>(2*minloop+4)) {
          	//no dangling ends on i-j pair:
             
			rarray=min(rarray,wmb.f(i+1,j-1)+data->eparam[5]+data->eparam[10]
            	+penalty(i,j,ct,data));


				//i+1 dangles on i-j pair:
			
			rarray=min(rarray,erg4(i,j,i+1,1,ct,data,lfce[i+1])+penalty(i,j,ct,data)+
            		wmb.f(i+2,j-1)+ data->eparam[5] + data->eparam[6] + data->eparam[10]);
			
			
				//j-1 dangles
			rarray=min(rarray,erg4(i,j,j-1,2,ct,data,lfce[j-1]) + penalty(i,j,ct,data) +
            	wmb.f(i+1,j-2) + data->eparam[5] + data->eparam[6] + data->eparam[10]);
				
            rarray=min(rarray,data->tstkm[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]]+
					checknp(lfce[i+1],lfce[j-1])+
					wmb.f(i+2,j-2) + data->eparam[5] + data->eparam[6] + data->eparam[6]+ data->eparam[10]
					+penalty(i,j,ct,data));
				
			

		
////////////////////			

		
			//consider the coaxial stacking of a helix from i to j onto helix i+1 or i+2 to ip:
			
			for (ip=i+1;(ip<j);ip++) {
				//first consider flush stacking


				

				//conditions guarantee that the coaxial stacking isn't considering an exterior loop 
				//if ((i!=number)/*&&(i+1!=number)*//*&&((j>number)||(ip!=number)&&(ip+1!=number))&&(j-1!=number)*/) {
				if (i!=number&&ip!=number&&j-1!=number) {
					rarray=min(rarray,penalty(i,j,ct,data)+v.f(i+1,ip)+
						penalty(i+1,ip,ct,data)+data->eparam[5]
						+data->eparam[10]+data->eparam[10]+(min(w.f(ip+1,j-1),wmb.f(ip+1,j-1)))+ergcoax(j,i,i+1,ip,0,ct,data));

			
					if((mod[i+1]||mod[ip])&&inc[ct->numseq[i+2]][ct->numseq[ip-1]]&&notgu(i+1,ip,ct)&&!(fce.f(i+1,ip)&SINGLE)) {

						rarray=min(rarray,penalty(i,j,ct,data)+v.f(i+2,ip-1)+
							penalty(i+1,ip,ct,data)+data->eparam[5]
							+data->eparam[10]+data->eparam[10]+(min(w.f(ip+1,j-1),wmb.f(ip+1,j-1)))+ergcoax(j,i,i+1,ip,0,ct,data)
							+erg1(i+1,ip,i+2,ip-1,ct,data));

					}
				


					//if ((ip<j-1)&&(i+2!=number)) {
					if (ip+2<j-1&&i+1!=number&&ip+1!=number) {
					//now consider an intervening nuc
						if ((ip+2<j-1)/*&&(j>number||ip+2!=number)*/)
						rarray=min(rarray,penalty(i,j,ct,data)+v.f(i+2,ip)+
							penalty(i+2,ip,ct,data)+data->eparam[5]
							+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]+
							(min(w.f(ip+2,j-1),wmb.f(ip+2,j-1)))
							+ergcoax(j,i,i+2,ip,ip+1,ct,data)+checknp(lfce[i+1],lfce[ip+1]));

						if((mod[i+2]||mod[ip])&&inc[ct->numseq[i+3]][ct->numseq[ip-1]]&&notgu(i+2,ip,ct)
							&&!(fce.f(i+2,ip)&SINGLE)) {

							rarray=min(rarray,penalty(i,j,ct,data)+v.f(i+3,ip-1)+
								penalty(i+2,ip,ct,data)+data->eparam[5]
								+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]+
								(min(w.f(ip+2,j-1),wmb.f(ip+2,j-1)))
								+ergcoax(j,i,i+2,ip,ip+1,ct,data)
								+erg1(i+2,ip,i+3,ip-1,ct,data)+checknp(lfce[i+1],lfce[ip+1]));

						}
		


						if (ip+1<j-2&&j-2!=number)
						rarray=min(rarray,penalty(i,j,ct,data)+v.f(i+2,ip)+
							penalty(i+2,ip,ct,data)+data->eparam[5]
							+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]+
							(min(w.f(ip+1,j-2),wmb.f(ip+1,j-2)))
							+ergcoax(j,i,i+2,ip,j-1,ct,data)+checknp(lfce[i+1],lfce[j-1]));

						if((mod[i+2]||mod[ip])&&inc[ct->numseq[i+3]][ct->numseq[ip-1]]&&notgu(i+2,ip,ct)
							&&!(fce.f(i+2,ip)&SINGLE)) {

							rarray=min(rarray,penalty(i,j,ct,data)+v.f(i+3,ip-1)+
								penalty(i+2,ip,ct,data)+data->eparam[5]
								+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]+
								(min(w.f(ip+1,j-2),wmb.f(ip+1,j-2)))
								+ergcoax(j,i,i+2,ip,j-1,ct,data)
								+erg1(i+2,ip,i+3,ip-1,ct,data)+checknp(lfce[i+1],lfce[j-1]));

						}

				


					}




			
				}


			}

			//consider the coaxial stacking of a helix from i to j onto helix ip to j-2 or j-1:
			for (ip=j-1;ip>i;ip--) {


				//conditions guarantee that the coaxial stacking isn't considering an exterior loop
				//if ((i!=number)&&(i+1!=number)&&((j>number)||(ip!=number)&&(ip-1!=number))&&(j-1!=number)) {
				if (j-1!=number&&ip-1!=number&&i!=number) {
					//first consider flush stacking
					rarray=min(rarray,penalty(i,j,ct,data)+v.f(ip,j-1)+
						penalty(j-1,ip,ct,data)+data->eparam[5]
						+data->eparam[10]+data->eparam[10]+(min(w.f(i+1,ip-1),wmb.f(i+1,ip-1)))+ergcoax(ip,j-1,j,i,0,ct,data));


					if((mod[ip]||mod[j-1])&&inc[ct->numseq[ip+1]][ct->numseq[j-2]]&&notgu(ip,j-1,ct)&&!(fce.f(ip,j-1)&SINGLE)) {
						rarray=min(rarray,penalty(i,j,ct,data)+v.f(ip+1,j-2)+
							penalty(j-1,ip,ct,data)+data->eparam[5]
							+data->eparam[10]+data->eparam[10]+(min(w.f(i+1,ip-1),wmb.f(i+1,ip-1)))
							+ergcoax(ip,j-1,j,i,0,ct,data)
							+erg1(ip,j-1,ip+1,j-2,ct,data));

					}

				

		

					if (j-2!=number) {
						//now consider an intervening nuc
						//if ((ip>i+1)&&(j>number||ip-2!=number))
						if (ip-2>i+1&&ip-2!=number) {
							rarray=min(rarray,penalty(i,j,ct,data)+v.f(ip,j-2)+
								penalty(j-2,ip,ct,data)+data->eparam[5]
								+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]+
								(min(w.f(i+1,ip-2),wmb.f(i+1,ip-2)))
								+ergcoax(ip,j-2,j,i,ip-1,ct,data)+checknp(lfce[j-1],lfce[ip-1]));
				


							if((mod[ip]||mod[j-2])&&inc[ct->numseq[ip+1]][ct->numseq[j-3]]&&notgu(ip,j-2,ct)&&!(fce.f(ip,j-2)&SINGLE)) {
								rarray=min(rarray,penalty(i,j,ct,data)+v.f(ip+1,j-3)+
									penalty(j-2,ip,ct,data)+data->eparam[5]
									+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]+
									(min(w.f(i+1,ip-2),wmb.f(i+1,ip-2)))
									+ergcoax(ip,j-2,j,i,ip-1,ct,data)
									+erg1(ip,j-2,ip+1,j-3,ct,data)+checknp(lfce[j-1],lfce[ip-1]));

							}
						}



						if ((ip-1>i+2)&&i+1!=number) {
							rarray=min(rarray,penalty(i,j,ct,data)+v.f(ip,j-2)+
								penalty(j-2,ip,ct,data)+data->eparam[5]
								+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]+
								(min(w.f(i+2,ip-1),wmb.f(i+2,ip-1)))
								+ergcoax(ip,j-2,j,i,i+1,ct,data)+checknp(lfce[j-1],lfce[i+1]));

							if((mod[ip]||mod[j-2])&&inc[ct->numseq[ip+1]][ct->numseq[j-3]]&&notgu(ip,j-2,ct)
								&&!(fce.f(ip,j-2)&SINGLE)) {
								rarray=min(rarray,penalty(i,j,ct,data)+v.f(ip+1,j-3)+
									penalty(j-2,ip,ct,data)+data->eparam[5]
									+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10]+
									(min(w.f(i+2,ip-1),wmb.f(i+2,ip-1)))
									+ergcoax(ip,j-2,j,i,i+1,ct,data)
									+erg1(ip,j-2,ip+1,j-3,ct,data)+checknp(lfce[j-1],lfce[i+1]));

							}
						}

						
					}
				
				}

				
				
				
			}



            /*if (ct->intermolecular) {

            	//intermolecular, so consider wmb2,
               //don't add the multiloop penalties because this is a exterior loop

            	e[1] = min(e[1],wmb2->f(i+1,j-1) + penalty(i,j,ct,data)+infinity);


            	//i+1 dangles on i-j pair:
            	if (i!=number) e[2] = min(e[2],erg4(i,j,i+1,1,ct,data,lfce[i+1]) + penalty(i,j,ct,data) +
            		wmb2->f(i+2,j-1)+infinity);
            	//j-1 dangles
            	if (j!=(number+1)) e[3] = min(e[3],erg4(i,j,j-1,2,ct,data,lfce[j-1]) + penalty(i,j,ct,data) +
            		wmb2->f(i+1,j-2)+infinity);
            	//both i+1 and j-1 dangle
            	if ((i!=number)&&(j!=(number+1))) {
            		e[4] = min(e[4],
            		data->tstkm[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] +
							pfchecknp(lfce[i+1],lfce[j-1]) +
               				wmb2->f(i+2,j-2) + penalty(i,j,ct,data)+infinity);

				}

				


			}*/


         }


      
	  }

      
      sub2:

	  v.f(i,j) = rarray;

	  if (fce.f(i,j)&PAIR)  {//force a pair between i and j
	  		w.f(i,j) = v.f(i,j)+data->eparam[10]+penalty(i,j,ct,data);
			wl.f(i,j) = v.f(i,j)+data->eparam[10]+penalty(i,j,ct,data);
	  		goto sub3;
      }
	  
      ////fill wmb:
	rarray = INFINITE_ENERGY;
	  if (((j-i-1)>(2*minloop+2))||j>number) {
         


         


			//also consider the coaxial stacking of two helixes in wv
			e = INFINITE_ENERGY;
			for (ip=i+minloop+1;ip<j-minloop-1;ip++) {
				//first consider flush stacking
		
				if (ip!=number) {
					rarray=min(rarray,v.f(i,ip)+v.f(ip+1,j)+penalty(i,ip,ct,data)
						+penalty(ip+1,j,ct,data)+ergcoax(i,ip,ip+1,j,0,ct,data));


					if ((mod[i]||mod[ip]||mod[ip+1]||mod[j])) {

						if ((mod[i]||mod[ip])&&(mod[ip+1]||mod[j])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]
							&&inc[ct->numseq[ip+2]][ct->numseq[j-1]]&&notgu(i,ip,ct)&&notgu(ip+1,j,ct)
								&&!(fce.f(ip+1,j)&SINGLE)&&!(fce.f(i,ip)&SINGLE)) {

							rarray=min(rarray,v.f(i+1,ip-1)+v.f(ip+2,j-1)+penalty(i,ip,ct,data)
								+penalty(ip+1,j,ct,data)+ergcoax(i,ip,ip+1,j,0,ct,data)
								+erg1(i,ip,i+1,ip-1,ct,data)+erg1(ip+1,j,ip+2,j-1,ct,data));


						}

						if ((mod[i]||mod[ip])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]&&notgu(i,ip,ct)&&!(fce.f(i,ip)&SINGLE)) {
							
							rarray=min(rarray,v.f(i+1,ip-1)+v.f(ip+1,j)+penalty(i,ip,ct,data)
								+penalty(ip+1,j,ct,data)+ergcoax(i,ip,ip+1,j,0,ct,data)
								+erg1(i,ip,i+1,ip-1,ct,data));


						}

						if ((mod[ip+1]||mod[j])&&inc[ct->numseq[ip+2]][ct->numseq[j-1]]&&notgu(ip+1,j,ct)&&!(fce.f(ip+1,j)&SINGLE)) {


							rarray=min(rarray,v.f(i,ip)+v.f(ip+2,j-1)+penalty(i,ip,ct,data)
								+penalty(ip+1,j,ct,data)+ergcoax(i,ip,ip+1,j,0,ct,data)
								+erg1(ip+1,j,ip+2,j-1,ct,data));

						}


					}
				


					if (ip+1!=number&&j!=number+1) {
						if (!lfce[ip+1]&&!lfce[j]) {
							//now consider an intervening mismatch
							if(!lfce[ip+1]&&!lfce[j]) e=min(e,v.f(i,ip)+v.f(ip+2,j-1)+penalty(i,ip,ct,data)
								+penalty(ip+2,j-1,ct,data)+ergcoax(i,ip,ip+2,j-1,j,ct,data));

							if (mod[i]||mod[ip]||mod[ip+2]||mod[j-1]) {
								if ((mod[i]||mod[ip])&&(mod[ip+2]||mod[j-1])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]
									&&inc[ct->numseq[ip+3]][ct->numseq[j-2]]&&notgu(i,ip,ct)&&notgu(ip+2,j-1,ct)
										&&!(fce.f(i,ip)&SINGLE)&&!(fce.f(ip+2,j-1)&SINGLE)) {

									if(!lfce[ip+1]&&!lfce[j]) e=min(e,v.f(i+1,ip-1)+v.f(ip+3,j-2)+penalty(i,ip,ct,data)
										+penalty(ip+2,j-1,ct,data)+ergcoax(i,ip,ip+2,j-1,j,ct,data)
										+erg1(i,ip,i+1,ip-1,ct,data)+erg1(ip+2,j-1,ip+3,j-2,ct,data));


								}

								if ((mod[i]||mod[ip])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]&&notgu(i,ip,ct)&&!(fce.f(i,ip)&SINGLE)) {
							
									if(!lfce[ip+1]&&!lfce[j]) e=min(e,v.f(i+1,ip-1)+v.f(ip+2,j-1)+penalty(i,ip,ct,data)
										+penalty(ip+2,j-1,ct,data)+ergcoax(i,ip,ip+2,j-1,j,ct,data)
										+erg1(i,ip,i+1,ip-1,ct,data));


								}

								if ((mod[ip+2]||mod[j-1])&&inc[ct->numseq[ip+3]][ct->numseq[j-2]]&&notgu(ip+2,j-1,ct)&&!(fce.f(ip+2,j-1)&SINGLE)) {


									if(!lfce[ip+1]&&!lfce[j]) e=min(e,v.f(i,ip)+v.f(ip+3,j-2)+penalty(i,ip,ct,data)
										+penalty(ip+2,j-1,ct,data)+ergcoax(i,ip,ip+2,j-1,j,ct,data)
										+erg1(ip+2,j-1,ip+3,j-2,ct,data));

								}
							}
						}
		
						if(!lfce[i]&&!lfce[ip+1]&&i!=number) {
							e=min(e,v.f(i+1,ip)+v.f(ip+2,j)+penalty(i+1,ip,ct,data)
								+penalty(ip+2,j,ct,data)+ergcoax(i+1,ip,ip+2,j,i,ct,data));

							if (mod[i+1]||mod[ip]||mod[ip+2]||mod[j]) {
								if ((mod[i+1]||mod[ip])&&(mod[ip+2]||mod[j])&&inc[ct->numseq[i+2]][ct->numseq[ip-1]]
									&&inc[ct->numseq[ip+3]][ct->numseq[j-1]]&&notgu(i+1,ip,ct)&&notgu(ip+2,j,ct)
									&&!(fce.f(i+1,ip)&SINGLE)&&!(fce.f(ip+2,j)&SINGLE)	) {

									e=min(e,v.f(i+2,ip-1)+v.f(ip+3,j-1)+penalty(i+1,ip,ct,data)
										+penalty(ip+2,j,ct,data)+ergcoax(i+1,ip,ip+2,j,i,ct,data)
										+erg1(i+1,ip,i+2,ip-1,ct,data)+erg1(ip+2,j,ip+3,j-1,ct,data));	


							
								}
								if ((mod[i+1]||mod[ip])&&inc[ct->numseq[i+2]][ct->numseq[ip-1]]&&notgu(i+1,ip,ct)&&!(fce.f(i+1,ip)&SINGLE)) {
							
									e=min(e,v.f(i+2,ip-1)+v.f(ip+2,j)+penalty(i+1,ip,ct,data)
										+penalty(ip+2,j,ct,data)+ergcoax(i+1,ip,ip+2,j,i,ct,data)
										+erg1(i+1,ip,i+2,ip-1,ct,data));


								}

								if ((mod[ip+2]||mod[j])&&inc[ct->numseq[ip+3]][ct->numseq[j-1]]&&notgu(ip+2,j,ct)&&!(fce.f(ip+2,j)&SINGLE)) {


									e=min(e,v.f(i+1,ip)+v.f(ip+3,j-1)+penalty(i+1,ip,ct,data)
										+penalty(ip+2,j,ct,data)+ergcoax(i+1,ip,ip+2,j,i,ct,data)
										+erg1(ip+2,j,ip+3,j-1,ct,data));

								}
							}
						}
					}
				}

			



			}

			wca[i] = min(rarray,e);
			rarray =rarray+data->eparam[10]+data->eparam[10];
			e = e+data->eparam[6]+data->eparam[6]+data->eparam[10]+data->eparam[10];
			 rarray = min(rarray,e);
			 wcoax.f(i,j) = rarray;
			

			//search for an open bifurcation:
			 for (k=i+1;k<j;k++) {
				//e = 0;
				
						branch= data->eparam[10]+v.f(i,k)+penalty(k,i,ct,data);

						if ((mod[i]||mod[k])&&inc[ct->numseq[i+1]][ct->numseq[k-1]]&&notgu(i,k,ct)) {

								branch=min(branch, data->eparam[10]+v.f(i+1,k-1)+penalty(k,i,ct,data)+erg1(i,k,i+1,k-1,ct,data));			

						} 



						
      					//calculate the energy of i stacked onto the pair of i+1,k

						branch=min(branch, v.f(i+1,k)+data->eparam[10]+data->eparam[6]+
         					erg4(k,i+1,i,2,ct,data,lfce[i])+penalty(i+1,k,ct,data));

						if ((mod[i+1]||mod[k])&&inc[ct->numseq[i+2]][ct->numseq[k-1]]&&notgu(i+1,k,ct)&&!(fce.f(i+1,k)&SINGLE)) {

								branch=min(branch, v.f(i+2,k-1) + data->eparam[10] +data->eparam[6] +
         							erg4(k,i+1,i,2,ct,data,lfce[i])+penalty(i+1,k,ct,data)
									+erg1(i+1,k,i+2,k-1,ct,data));


						}
         
						
						
      					//calculate the energy of k stacked onto the pair of i,k-1
						if (k!=1) {
         					branch=min(branch, v.f(i,k-1)+ data->eparam[10] + data->eparam[6] +
         						erg4(k-1,i,k,1,ct,data,lfce[k])+penalty(i,k-1,ct,data));

							if ((mod[i]||mod[k-1])&&inc[ct->numseq[i+1]][ct->numseq[k-2]]&&notgu(i,k-1,ct)&&!(fce.f(i,k-1)&SINGLE)) {

									branch=min(branch, v.f(i+1,k-2) + data->eparam[10] + data->eparam[6] +
         								erg4(k-1,i,k,1,ct,data,lfce[k])+penalty(i,k-1,ct,data)
										+erg1(i,k-1,i+1,k-2,ct,data));
				


							}

         	
						}
						
						
      					//calculate i and k stacked onto the pair of i+1,k-1
						if (k!=1&&!lfce[i]&&!lfce[k]) {
         					branch=min(branch, v.f(i+1,k-1) +data->eparam[10] + (data->eparam[6]+data->eparam[6]) +
         						data->tstkm[ct->numseq[k-1]][ct->numseq[i+1]][ct->numseq[k]][ct->numseq[i]]
							+penalty(k-1,i+1,ct,data));



							if ((mod[i+1]||mod[k-1])&&(k-2>0)&&!(fce.f(i+1,k-1)&SINGLE)) {
								if(inc[ct->numseq[i+2]][ct->numseq[k-2]]&&notgu(i+1,k-1,ct)) {

									branch=min(branch, v.f(i+2,k-2) + data->eparam[10] + (data->eparam[6]+data->eparam[6]) +
         								data->tstkm[ct->numseq[k-1]][ct->numseq[i+1]][ct->numseq[k]][ct->numseq[i]]
										+penalty(k-1,i+1,ct,data)+erg1(i+1,k-1,i+2,k-2,ct,data));

								}
							}
						}
						


						rarray=min(rarray,(min(branch,wcoax.f(i,k)))+(min(wl.f(k+1,j),wmbl.f(k+1,j))));

					

			 }
			
			

			if (i!=number)
				if (!lfce[i]) rarray=min(rarray,wmbl.f(i+1,j)+data->eparam[6]);

			wmbl.f(i,j) = rarray;

			wmb.f(i,j) = rarray;
			if (j!=number+1)
				if (!lfce[j]) wmb.f(i,j)=min(wmb.f(i,j),wmb.f(i,j-1)+data->eparam[6]);	
			

		/*if (ct->intermolecular) {
         	//intermolecular folding:


			

         	//search for an open bifurcation:
         	for (k=i;k<=j;k++) {

				if (k!=number) wmb2->f(i,j) = min(wmb2->f(i,j),w2->f(i,k)+work2[k+1][jmt]);

				
         	}

		
			
			if (i!=number)
				if (!(fce.f(i,i)&INTER)) wmb2->f(i,j) = min(wmb2->f(i,j) ,wmb2->f(i+1,j) );
				else  wmb2->f(i,j) = min(wmb2->f(i,j) ,wmb2->f(i+1,j) + data->init - infinity);

			if (j!=number+1)
				if (!(fce.f(j,j)&INTER)) wmb2->f(i,j)  = min(wmb2->f(i,j) ,wmb2->f(i,j-1));
				else wmb2->f(i,j)  = min(wmb2->f(i,j) ,wmb2->f(i,j-1) +data->init-infinity);


         	
			w2->f(i,j) = min(w2->f(i,j),wmb2->f(i,j) );

		}*/


	  } 


	  //Compute w[i][j]

     wl.f(i,j)= data->eparam[10]+v.f(i,j)+penalty(j,i,ct,data);

	 if ((mod[i]||mod[j])&&inc[ct->numseq[i+1]][ct->numseq[j-1]]&&notgu(i,j,ct)) {

			wl.f(i,j)=min(wl.f(i,j), data->eparam[10]+v.f(i+1,j-1)+penalty(j,i,ct,data)+erg1(i,j,i+1,j-1,ct,data));			

	} 



      if (i!=number) {
      	//calculate the energy of i stacked onto the pair of i+1,j

         wl.f(i,j)=min(wl.f(i,j), v.f(i+1,j)+data->eparam[10]+data->eparam[6]+
         	erg4(j,i+1,i,2,ct,data,lfce[i])+penalty(i+1,j,ct,data));

		 if ((mod[i+1]||mod[j])&&inc[ct->numseq[i+2]][ct->numseq[j-1]]&&notgu(i+1,j,ct)&&!(fce.f(i+1,j)&SINGLE)) {

				wl.f(i,j)=min(wl.f(i,j), v.f(i+2,j-1) + data->eparam[10] +data->eparam[6] +
         			erg4(j,i+1,i,2,ct,data,lfce[i])+penalty(i+1,j,ct,data)
					+erg1(i+1,j,i+2,j-1,ct,data));


		 }
         
      }
      if (j!=((number)+1)) {
      	//calculate the energy of j stacked onto the pair of i,j-1
         if (j!=1) {
         	wl.f(i,j)=min(wl.f(i,j), v.f(i,j-1)+ data->eparam[10] + data->eparam[6] +
         		erg4(j-1,i,j,1,ct,data,lfce[j])+penalty(i,j-1,ct,data));

			if ((mod[i]||mod[j-1])&&inc[ct->numseq[i+1]][ct->numseq[j-2]]&&notgu(i,j-1,ct)&&!(fce.f(i,j-1)&SINGLE)) {

					wl.f(i,j)=min(wl.f(i,j), v.f(i+1,j-2) + data->eparam[10] + data->eparam[6] +
         				erg4(j-1,i,j,1,ct,data,lfce[j])+penalty(i,j-1,ct,data)
						+erg1(i,j-1,i+1,j-2,ct,data));
				


			}

         	
         }
      }
      if ((i!=(number))&&(j!=((number)+1))) {
      	//calculate i and j stacked onto the pair of i+1,j-1
         if (j!=1&&!lfce[i]&&!lfce[j]) {
         	wl.f(i,j)=min(wl.f(i,j), v.f(i+1,j-1) +data->eparam[10] + (data->eparam[6]+data->eparam[6]) +
         		data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
               +penalty(j-1,i+1,ct,data));



			if ((mod[i+1]||mod[j-1])&&(j-2>0)&&!(fce.f(i+1,j-1)&SINGLE)) {
				if(inc[ct->numseq[i+2]][ct->numseq[j-2]]&&notgu(i+1,j-1,ct)) {

					wl.f(i,j)=min(wl.f(i,j), v.f(i+2,j-2) + data->eparam[10] + (data->eparam[6]+data->eparam[6]) +
         				data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
						+penalty(j-1,i+1,ct,data)+erg1(i+1,j-1,i+2,j-2,ct,data));

				}
			}
         }
	  }

	  
	  

	  if (i!=number&&!lfce[i]) {
         	//if (!(fce.f(i,i)&INTER))
               //add a nuc to an existing loop:
         		wl.f(i,j)=min(wl.f(i,j),wl.f(i+1,j)+data->eparam[6]);
            	//this is for when i represents the center of an intermolecular linker:
           // else e[4] = w.f(i+1,j) + data->eparam[6] + infinity;
      }

	  w.f(i,j) = wl.f(i,j);
	  if (j!=number+1&&!lfce[j]) {
             	//if (!(fce.f(j,j)&INTER)) {
               	//add a nuc to an existing loop:
               	w.f(i,j)=min(w.f(i,j), w.f(i,j-1) + data->eparam[6]);
               //}
               //else e[5] = w.f(i,j-1) + data->eparam[6] + infinity;

      }
     

	 /* if (ct->intermolecular) {

      		//wmb2[i][j%3] = infinity;
      		//keep track of w2:
			for (ii=1;ii<=5;ii++) e[ii] = 2*infinity;


			if (i!=number) {
      			//calculate the energy of i stacked onto the pair of i+1,j

         		e[1] = v.f(i+1,j) +
         			erg4(j,i+1,i,2,ct,data,lfce[i])+penalty(i+1,j,ct,data);

				if ((mod[i+1]||mod[j])&&inc[ct->numseq[i+2]][ct->numseq[j-1]]) {

					e[1] = min(e[1],v.f(i+2,j-1) +
         				erg4(j,i+1,i,2,ct,data,lfce[i])+penalty(i+1,j,ct,data)
						+erg1(i+1,j,i+2,j-1,ct,data));

				}


         		
         		e[4] = w2->f(i+1,j);
            	
			}
      		if (j!=((number)+1)) {
      		//calculate the energy of j stacked onto the pair of i,j-1
         		if (j!=1) {
         			e[2] = v.f(i,j-1)   +
         				erg4(j-1,i,j,1,ct,data,lfce[j])+penalty(i,j-1,ct,data);

					if ((mod[i]||mod[j-1])&&inc[ct->numseq[i+1]][ct->numseq[j-2]]) {

						e[2] = min(e[2],v.f(i+1,j-2) +
         					erg4(j-1,i,j,1,ct,data,lfce[j])+penalty(i,j-1,ct,data)
							+erg1(i,j-1,i+1,j-2,ct,data));

					}

         			
               		e[5] = w2->f(i,j-1);
                
         		}
      		}
      		if ((i!=(number))&&(j!=((number)+1))) {
      			//calculate i and j stacked onto the pair of i+1,j-1
         		if (j!=1) {
         			e[3] = v.f(i+1,j-1)   +
         				data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]]
									[ct->numseq[j]][ct->numseq[i]]
					+pfchecknp(lfce[i+1],lfce[j-1])
               		+penalty(j-1,i+1,ct,data);



					if ((mod[i+1]||mod[j-1])&&inc[ct->numseq[i+2]][ct->numseq[j-2]]) {

						e[3] = min(e[3],v.f(i+2,j-2) +
         					data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
							+pfchecknp(lfce[i+1],lfce[j-1])
							+penalty(j-1,i+1,ct,data)+erg1(i+1,j-1,i+2,j-2,ct,data));

					}
         		}
      		}

			e[1] = min(e[1],(v.f(i,j)+penalty(j,i,ct,data)));

			if (mod[i]||mod[j]&&inc[ct->numseq[i+1]][ct->numseq[j-1]]) {

				e[1] = min((v.f(i+1,j-1)+penalty(j,i,ct,data)+erg1(i,j,i+1,j-1,ct,data)),e[1]);			

			}


      		w2->f(i,j) = min(e[1],e[2]);
      		w2->f(i,j) = min(w2->f(i,j),e[3]);
      		w2->f(i,j) = min(w2->f(i,j),e[4]);
      		w2->f(i,j) = min(w2->f(i,j),e[5]);




		}*/

        

      




      
      sub3:

      

      
   


      //Compute w5[i], the energy of the best folding from 1->i, and
      	//w3[i], the energy of the best folding from i-->GetSequenceLength()
	if (i==(1)) {

	   if (j<=minloop+1) {
		   if (lfce[j]) w5[j]= INFINITE_ENERGY;
		   else  w5[j] = w5[j-1];
		}
	   
		else {
      		if (lfce[j]) rarray = INFINITE_ENERGY;
			
			else rarray = w5[j-1];
         
      		
			
      		for (k=0;k<=(j-4);k++) {



      			rarray=min(rarray,w5[k]+v.f(k+1,j)+penalty(j,k+1,ct,data));

				if ((mod[k+1]||mod[j])&&inc[ct->numseq[k+2]][ct->numseq[j-1]]&&notgu(k+1,j,ct)&&!(fce.f(k+1,j)&SINGLE)) {

					rarray=min(rarray,w5[k]+v.f(k+2,j-1)+penalty(j,k+1,ct,data)
						+erg1(k+1,j,k+2,j-1,ct,data));
				}

			

				rarray=min(rarray,w5[k]+erg4(j,k+2,k+1,2,ct,data,lfce[k+1])+v.f(k+2,j)+penalty(j,k+2,ct,data));
			
				if((mod[k+2]||mod[j])&&inc[ct->numseq[k+3]][ct->numseq[j-1]]&&notgu(k+2,j,ct)
					&&!(fce.f(k+2,j)&SINGLE)) {
					rarray=min(rarray,w5[k]+erg4(j,k+2,k+1,2,ct,data,lfce[k+1])+v.f(k+3,j-1)+penalty(j,k+2,ct,data)+erg1(k+2,j,k+3,j-1,ct,data));

				}


         		rarray=min(rarray,w5[k]+erg4(j-1,k+1,j,1,ct,data,lfce[j])+v.f(k+1,j-1)+penalty(j-1,k+1,ct,data));

				if ((mod[k+1]||mod[j-1])&&inc[ct->numseq[k+2]][ct->numseq[j-2]]&&notgu(k+1,j-1,ct)&&!(fce.f(k+1,j-1)&SINGLE)) {

					rarray=min(rarray,w5[k]+erg4(j-1,k+1,j,1,ct,data,lfce[j])+v.f(k+2,j-2)
						+penalty(j-1,k+1,ct,data)+erg1(k+1,j-1,k+2,j-2,ct,data));
				}



				rarray=min(rarray,w5[k]+data->tstack[ct->numseq[j-1]][ct->numseq[k+2]][ct->numseq[j]][ct->numseq[k+1]] 
									+checknp(lfce[j],lfce[k+1]) + v.f(k+2,j-1)+
									penalty(j-1,k+2,ct,data));

				if ((mod[k+2]||mod[j-1])&&inc[ct->numseq[k+3]][ct->numseq[j-2]]&&notgu(k+2,j-1,ct)&&!(fce.f(k+2,j-1)&SINGLE)) {

					rarray=min(rarray,w5[k]+data->tstack[ct->numseq[j-1]][ct->numseq[k+2]][ct->numseq[j]][ct->numseq[k+1]] 
									+checknp(lfce[j],lfce[k+1]) + v.f(k+3,j-2)+
									penalty(j-1,k+2,ct,data)+erg1(k+2,j-1,k+3,j-2,ct,data));

				}



				
				
		


				rarray=min(rarray,w5[k]+wca[k+1]);
		
	
			}
			

		w5[j] = rarray;
			
		}

	}
  }


   
}


delete[] wca;

//////////////////////////
//output V, W, WMB, and W2V:
#if defined (pfdebugmode)
ofstream foo;
foo.open("arrays.out");
foo << "i" << "\t"<<"j"<<"\t"<<"v.f(i,j)"<<"\t"<<"w.f(i,j)"<<"\t"<<"wmb.f(i,j)\twmbl.f(i,j)\twcoax.f(i,j)"<<"\t"<<"wl.f(i,j)"<<"\t"<<"v.f(j,i+number)"<<"\t"<<"w.f(j,i+number)"<<"\t"<<"wmb.f(j,i+number)"<<"\t"<<"wl.f(j,i+number)"<<"\t"<<"wmbl.f(j,i+numer)\twcoax.f(j,i+number)"<<"\n";
for (j=1;j<=number;j++) {
	for (i=1;i<=j;i++) {

		foo << i << "\t"<<j<<"\t"<<v.f(i,j)<<"\t"<<w.f(i,j)<<"\t"<<wmb.f(i,j)<<"\t"<<wmbl.f(i,j)<<"\t"<<wcoax.f(i,j)<<"\t"<<wl.f(i,j)<<"\t"<<v.f(j,i+number)<<"\t"<<w.f(j,i+number)<<"\t"<<wmb.f(j,i+number)<<"\t"<<wl.f(j,i+number)<<"\t"<<wmbl.f(j,i+number)<<"\t"<<wcoax.f(j,i+number)<<"\n";

	}
}

foo <<"\n\n\n";
foo << "i" << "\t" << "w5[i]" <<  "\n";
for (i=0;i<=number;i++) foo << i << "\t" << w5[i]  << "\n";

foo.close();

#endif


if (save!=0) {
	ofstream sav(save,ios::binary);
	
	//write the save file information so that the sequence can be re-folded,
		//include thermodynamic data to prevent traceback errors

	//start with structure information
	int localint;
	localint= ct->GetSequenceLength();
	write(&sav,&(localint));
	write(&sav,&(ct->intermolecular));

	writestructuresave(&sav,ct);

	//now write the array class data for v, w, and wmb:
	for (i=0;i<=ct->GetSequenceLength();i++) {
		write(&sav,&(w5[i]));
		for (j=0;j<=ct->GetSequenceLength();j++) {
			write(&sav,&(v.dg[i][j]));
			write(&sav,&(w.dg[i][j]));
			write(&sav,&(wmb.dg[i][j]));
			write(&sav,&(wmbl.dg[i][j]));
			write(&sav,&(wcoax.dg[i][j]));
			write(&sav,&(wl.dg[i][j]));
			writesinglechar(&sav,&(fce.dg[i][j]));
			if (ct->intermolecular) {
				write(&sav,&(w2->dg[i][j]));
				write(&sav,&(wmb2->dg[i][j]));

			}


		}	
	}

	for (i=0;i<=2*ct->GetSequenceLength();i++) {
		write(&sav,&(lfce[i]));
		write(&sav,&(mod[i]));

	}

	//now write the thermodynamic data:
	write(&sav,data);

	sav.close();
}


//do the tracebacks:
alltracetraceback(writect,ctenergy,ct,&v,&w,&wmb,&wl,&wmbl,&wcoax,&fce,w5,lfce,mod,data,percentdelta,absolutedelta);

delete[] lfce;
delete[] mod;



delete[] w5;



if (ct->intermolecular) {
	delete w2;
	delete wmb2;
}



return;
}










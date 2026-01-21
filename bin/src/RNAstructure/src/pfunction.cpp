/*Code for partition function calculation of base pair probabilities.
Copyright 2003,2004,2005,2006 by David H. Mathews

Contributions by Zhi John Lu, 2006, and Michael Sloma, 2013.

*/

//NOTE: Two compiler flags exist to modify the behaviour of MB Loops:
//SIMPLEMBLOOP uses a model in which every helix end has a 5' and 3' dangle (Except at sequence ends)
//disablecoax is a flag that disables coaxial stacking, but leaves the rest of the model intact

// #define ENABLE_DEBUG_LOGS
// #undef DEBUG_LOGGING
#include "debug_logging.h"

#include "pfunction.h"
#include "boltzmann.h" //for boltzman
#include <math.h>
#include <cstdlib>
#include <iomanip>

#ifdef SMP
	#include <omp.h>
#endif

#ifdef PF_LOG_CLASS
const log_double xlog_zero(0); // i.e. xlog(0)
const log_double xlog_one(1); // i.e. xlog(1)
#endif

using namespace std;

#undef pfdebugmode  //flag to indicate debugging
//#define pfdebugmode  //flag to indicate debugging
//#define pfdebugpath "./rescale_alert.out" file where pfdebugmode logging output should go.

#undef equiout
//#define equiout //flag to indicate equilibria should be written to file

#undef timer

//#define timer //flag to indicate the code execution should be timed

// #define maxinter 30   //maximum length of the internal loops
// #define maxasym 30  //maximum asymetry in the internal loops

// int TotalRescaleCount = 0;

void calculatepfunction(structure* ct,pfdatatable* data, ProgressHandler* update, char* save, bool quickQ, PFPRECISION *Q,
	DynProgArray<PFPRECISION> *w, DynProgArray<PFPRECISION> *v, DynProgArray<PFPRECISION> *wmb, DynProgArray<PFPRECISION> *wl, DynProgArray<PFPRECISION> *wlc, DynProgArray<PFPRECISION> *wmbl,
	DynProgArray<PFPRECISION> *wcoax, forceclass *fce,PFPRECISION *w5,PFPRECISION *w3,bool *mod, bool *lfce, bool disablecoax, bool allowisolated, int maxinter) {
int ip,ii,jj;
//int ip, jp, ii, jj, jpf, jf, bl, ll, dp;
int i,j,h,d;
int k,p;
bool **inc;
int number,maxj;
PFPRECISION twoscaling,rarray;
PFPRECISION **wca;
int N = ct->GetNumberofSequences();
// PFPRECISION penaltyij;

#ifdef timer
	#include <time.h>
	ofstream timeout;
	int seconds;
	char timerstring[100];
	char timelength[10];
	strcpy(timerstring,"time_pf_");
	sprintf(timelength,"%i",ct->GetSequenceLength());
	strcat(timerstring,timelength);
	strcat(timerstring,".out");

	timeout.open(timerstring);
	timeout<<time(NULL)<<"\n";
	seconds = time(NULL);
#endif

	PFPRECISION(*ptr_to_erg1)(int i, int j, int ip, int jp, structure * ct, pfdatatable * data);
	PFPRECISION(*ptr_to_erg2)(int i, int j, int ip, int jp, structure * ct, pfdatatable * data, char a, char b);
	PFPRECISION(*ptr_to_erg3)(int i, int j, structure * ct, pfdatatable * data, char dbl);
	PFPRECISION(*ptr_to_erg4)(int i, int j, int ip, int jp, structure * ct, pfdatatable * data, bool lfce);
	PFPRECISION(*ptr_to_ergcoaxflushbases) (int i, int j, int ip, int jp, structure * ct, pfdatatable * data);
	PFPRECISION(*ptr_to_ergcoaxinterbases1) (int i, int j, int ip, int jp, structure * ct, pfdatatable * data);
	PFPRECISION(*ptr_to_ergcoaxinterbases2) (int i, int j, int ip, int jp, structure * ct, pfdatatable * data);
	PFPRECISION(*ptr_to_penalty) (int i, int j, structure * ct, pfdatatable * data);
	bool(*ptr_to_can_pair)(int i, int j, structure * ct, pfdatatable * data);
	PFPRECISION(*ptr_to_tstack)(int a, int b, int c, int d, structure * ct, pfdatatable * data);
	PFPRECISION(*ptr_to_tstkm)(int a, int b, int c, int d, structure * ct, pfdatatable * data);
	//	PFPRECISION(*ptr_to_eparam)(int a, structure * ct, pfdatatable * data);

	if (ct->GetNumberofSequences() == 1) {
		ptr_to_erg1 = erg1;
		ptr_to_erg2 = erg2;
		ptr_to_erg3 = erg3;
		ptr_to_erg4 = erg4;
		ptr_to_ergcoaxflushbases = ergcoaxflushbases;
		ptr_to_ergcoaxinterbases1 = ergcoaxinterbases1;
		ptr_to_ergcoaxinterbases2 = ergcoaxinterbases2;
		ptr_to_penalty = penalty;
		ptr_to_can_pair = can_pair;
		ptr_to_tstack = tstack;
		ptr_to_tstkm = tstkm;
		//		ptr_to_eparam = eparam;
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
		ptr_to_can_pair = multi_can_pair;
		ptr_to_tstack = multi_tstack;
		ptr_to_tstkm = multi_tstkm;
		//		ptr_to_eparam = multi_eparam;
	}


number = (ct->GetSequenceLength());//place the number of bases in a registered integer

DynProgArrayT<PFPRECISION> v_T(number);
DynProgArrayT<PFPRECISION> wl_T(number);
DynProgArrayT<PFPRECISION> wmbl_T(number);

//define and build the inc array:
inc = new bool *[data->alphabet.size()];
for (int letters=0; letters<data->alphabet.size();++letters) {
	inc[letters]=new bool[data->alphabet.size()];
	for (int letterpairs=0;letterpairs<data->alphabet.size();++letterpairs) {
		inc[letters][letterpairs]=data->pairing[letters][letterpairs];

	}
}


wca = new PFPRECISION *[number+1];

for (i=0;i<=number;i++) {
	wca[i]=new PFPRECISION [number+1];
	for(j=0;j<=number;j++) {
		wca[i][j] = (PFPRECISION) ZERO;
	}
}

w5[0] = (PFPRECISION) ONE;  // if PF_LOG_CALC, then ONE becomes log(ONE) = 0.0 initialize the random coil contribution to the partition function
w3[number+1] = (PFPRECISION) ONE;

twoscaling = PROD(data->scaling,data->scaling);

force(ct,fce,lfce);


//This is the fill routine:

if (quickQ) maxj = number;
else maxj = 2*number-1;

for (h=0;h<=( quickQ?(maxj-1):(maxj-1-minloop) );h++){
	
	d=(h<=(number-1))?h:(h-number+1);
	
	if (((h%10)==0)&&update) {
		update->update((100*h)/(2*ct->GetSequenceLength()));
		if (update->canceled()) break; //exit the for-loop because the user canceled the calculation.
	}

	int start;
	int end;
	if (h<=(number-1)) {
		start = 1;
		end = number-h;
	}
	else {
		start = 2*number-h;
		end = number;
	}

	//for (int locali=((h<=(number-1))?1:(2*number-h));locali<=((h<=(number-1))?(number-h):number);locali++){
	#ifdef SMP
		#pragma omp parallel for
	#endif
	for (int locali=start;locali<=end;locali++){
		int localj=locali+d;

		// [Debug_mod]
		// if (locali == 128 && localj == 184){
		// 	SET_DEBUG_LEVEL(DEBUG4)
		// }

		//Test to make sure the fragment is large enough to form pairs:
		if (!(localj<=(number)&&((localj-locali)<=minloop))) {
			//First, calculate V(locali,localj): (The partition function from nucleotides locali to localj, given that locali and localj pair)

			//Start the value of V:
			PFPRECISION localrarray=ZERO;
			PFPRECISION locale;
			PFPRECISION penaltyij;

			//Now, test some conditions as to whether V should be evaluated:

			bool calculatev=true;
			int localii,localjj,localp,localafter;

			if (ct->templated) {
				//Test whether these nucleotides are allowed to pair because there is a pairing mask in ct ("tem").
				if (locali>ct->GetSequenceLength()) localii = locali - ct->GetSequenceLength();
				else localii = locali;
				if (localj>ct->GetSequenceLength()) localjj = localj - ct->GetSequenceLength();
				else localjj = localj;
				if (localjj<localii) {
					localp = localjj;
      				localjj = localii;
					localii = localp;
				}
   				if (!ct->tem[localjj][localii]) calculatev = false;
			}

			if (fce->f(locali,localj)&SINGLE) {
				//locali or localj is forced single-stranded
				calculatev = false;
			}
			if (fce->f(locali,localj)&NOPAIR) {
				//locali or localj is forced into a pair elsewhere
				calculatev = false;
			}

//			if (!inc[ct->numseq[locali]][ct->numseq[localj]]) {
			if (!ptr_to_can_pair(locali, localj, ct, data)) {
				//These are two nucleotides that cannot form a canonical pair
				calculatev = false;
			}

			//force u's into gu pairs, if the user has specified these restraints
			for (int localip=0;localip<ct->GetNumberofGU();localip++) {
				if (ct->GetGUpair(localip)==locali) {
					if (ct->numseq[localj]!=3) {
         				calculatev = false;
					}
				}
				else if (ct->GetGUpair(localip)==localj) {
       				if (ct->numseq[locali]!=3) {
         				calculatev = false;
					}
				}
				else if ((ct->GetGUpair(localip)+number)==localj) {
       				if (ct->numseq[locali]!=3) {
						calculatev = false;
					}
				}
			}

			//now check to make sure that this isn't an isolated pair:
			//	(consider a pair separated by a bulge as not! stacked)

			//before = 0 if a stacked pair cannot form 5' to locali
			int localbefore =0;
			if ((locali>1&&localj<(2*number)&&localj!=number)) {
				if ((localj>number&&((locali-localj+number)>minloop+2))||localj<number) {
	//				localbefore = inc[ct->numseq[locali-1]][ct->numseq[localj+1]];
					localbefore = ptr_to_can_pair(locali - 1, localj + 1, ct, data);
				}
			}

			//after = 0 if a stacked pair cannot form 3' to locali
			if ((((localj-locali)>minloop+2)&&(localj<=number)||(localj>number+1))&&(locali!=number)) {
//				localafter = inc[ct->numseq[locali+1]][ct->numseq[localj-1]];
				localafter = ptr_to_can_pair(locali + 1, localj - 1, ct, data);

			}
			else localafter = 0;

			
			if (ct->GetNumberofSequences() > 1) {
				if (localbefore == 0 && localafter == 0) {
					if (localj > number)
					{
						localbefore = ct->can_pair_isolated(localj - number, locali, ct);
					}
					else
						localbefore = ct->can_pair_isolated(locali, localj, ct);
				}
			}
			
			//if there are no stackable pairs to locali.localj then don't allow a pair locali,localj
			//only implement if !allowisolated, the default behavior to filter these pairs from consideration
			if ((localbefore==0)&&(localafter==0)&&!allowisolated) {
				//v->f(locali,localj)= 0;
				calculatev = false;
			}

			//A large code block for filling V, if all the above conditions pass:
			if (calculatev) {
				penaltyij = ptr_to_penalty(locali,localj,ct,data);

				//Test to make sure this isn't the end of the sequence
				if (!(locali==(number)||localj==((number)+1))) {
   					//Perhaps locali and localj close a hairpin:
					//[X5Y34F]
					localrarray=ptr_to_erg3(locali,localj,ct,data,fce->f(locali,localj));
					LOG_DEBUG4("[X5Y34F]\tlocalrarray_(v(i,j)): " << localrarray);

					if ((localj-locali-1)>=(minloop+2)||localj>(number)) {
      					//Perhaps locali,localj stacks over locali+1,localj-1
						if (!mod[locali]&&!mod[localj]){  //make sure this is not a site of chemical modification
							//[NDYIAA]
							localrarray=SUM(localrarray,PROD(ptr_to_erg1(locali,localj,locali+1,localj-1,ct,data),v->f(locali+1,localj-1)));
							LOG_DEBUG4("[NDYIAA]\tlocalrarray_(v(i,j)): " << localrarray);
						}
						else {
							//allow G-U to be modified or a pair next to a G-U to be modified
							if ((ct->numseq[locali]==3&&ct->numseq[localj]==4)||(ct->numseq[locali]==4&&ct->numseq[localj]==3)) {
								localrarray=SUM(localrarray,PROD(ptr_to_erg1(locali,localj,locali+1,localj-1,ct,data),v->f(locali+1,localj-1)));
								LOG_DEBUG4("[NDYIAA.2]\tlocalrarray_(v(i,j)): " << localrarray);
							}
							else if ((ct->numseq[locali+1]==3&&ct->numseq[localj-1]==4)||(ct->numseq[locali+1]==4&&ct->numseq[localj-1]==3)) {
								localrarray=SUM(localrarray,PROD(ptr_to_erg1(locali,localj,locali+1,localj-1,ct,data),v->f(locali+1,localj-1)));
								LOG_DEBUG4("[NDYIAA.3]\tlocalrarray_(v(i,j)): " << localrarray);
							}
							else if (locali-1>0&&localj+1<2*number) {
								if ((ct->numseq[locali-1]==3&&ct->numseq[localj+1]==4)||(ct->numseq[locali-1]==4&&ct->numseq[localj+1]==3)) {
									localrarray=SUM(localrarray,PROD(ptr_to_erg1(locali,localj,locali+1,localj-1,ct,data),v->f(locali+1,localj-1)));
									LOG_DEBUG4("[NDYIAA.4]\tlocalrarray_(v(i,j)): " << localrarray);
								}
							}
						}
					}

					// Use the O(N^4) recursions for internal loops.  The O(N^3) internal loop code is still available in the log_double branch.
					// Benchmarks showed negligible difference in the performance of partition between the implementations when using a max loop
					// size of 30.
					if (((localj-locali-1)>=(minloop+3))||(localj>(number))) {
						int maxsize;
						if (localj<=number) maxsize = min(maxinter,localj-locali-minloop-3);//interior fragment
						else maxsize = min(localj - locali - 3,maxinter); //exterior fragment
						for (int size=1;size<=maxsize;++size) {
							int localip, localjp;
							for (localip=locali+1,localjp=localj-size-1;localip<=locali+size+1&&localip<=number;++localip,++localjp) {
									if (localjp<=number) {
										localrarray=SUM(localrarray,PROD(ptr_to_erg2(locali,localj,localip,localjp,ct,data,fce->f(locali,localip),
																fce->f(localjp,localj)), v->f(localip,localjp)));
										LOG_TRACE("[FQGH8G.11]\tlocalrarray_(v(i,j)): " << localrarray);

										if ((mod[localip]||mod[localjp])&&ptr_to_can_pair(localip, localjp, ct, data) ) {
											//locali or localj is modified
											localrarray=SUM(localrarray,PROD(ptr_to_erg2(locali,localj,localip,localjp,ct,data,fce->f(locali,localip),
																	fce->f(localjp,localj)),
																v->f(localip+1,localjp-1),ptr_to_erg1(localip,localjp,localip+1,localjp-1,ct,data)));
											LOG_TRACE("[FQGH8G.12]\tlocalrarray_(v(i,j)): " << localrarray);

										}
									}
									else {
										localrarray=SUM(localrarray,PROD(ptr_to_erg2(locali,localj,localip,localjp,ct,data,fce->f(locali,localip),fce->f(localjp-number,localj-number)),
															v->f(localip,localjp)));
										LOG_TRACE("[FQGH8G.13]\tlocalrarray_(v(i,j)): " << localrarray);

										if ((mod[localip]||mod[localjp])&&ptr_to_can_pair(localip, localjp, ct, data)) {
											//locali or localj is modified
											localrarray=SUM(localrarray,PROD(ptr_to_erg2(locali,localj,localip,localjp,ct,data,fce->f(locali,localip),fce->f(localjp-number,localj-number)),
																v->f(localip+1,localjp-1),ptr_to_erg1(localip,localjp,localip+1,localjp-1,ct,data)));
											LOG_TRACE("[FQGH8G.14]\tlocalrarray_(v(i,j)): " << localrarray);
										}
									}
							}//for size
						}//for localip
					}//if ((localj-locali-1)>=(minloop+3))||(localj>(number))
					
					LOG_DEBUG3("[RSYVZ0]\tlocalrarray_(v(i,j)): " << localrarray);
				}//end of condition to make sure this wasn't the end of the sequence

				//Perhaps locali,localj closes a multibranch or exterior loop, enumerate all possibilities
				if (((localj-locali-1)>=(2*minloop+4))||(localj>(number))) {
					//consider the exterior loop closed by locali,localj
					//[9O48X9]
					if (localj>number) {
						#ifdef SIMPLEMBLOOP
						//5' and 3' dangling ends
						if (locali!=number&&localj!=number+1) 
						{
							localrarray=SUM(localrarray, PROD(erg4(locali,localj,locali+1,1,ct,data,lfce[locali+1]),erg4(locali,localj,localj-1,2,ct,data,lfce[localj-1]),
								penaltyij,w3[locali+1],w5[localj-number-1] TWOSCALING));
							LOG_TRACE("[XXXXXX]\tlocalrarray_(v(i,j)): " << localrarray);
						}
						//3' dangling end
						else if (locali!=number) {
							localrarray=SUM(localrarray, PROD(erg4(locali,localj,locali+1,1,ct,data,lfce[locali+1]),
								penaltyij,w3[locali+1],w5[localj-number-1] TWOSCALING));
							LOG_TRACE("[XXXXXX]\tlocalrarray_(v(i,j)): " << localrarray);
						}
						//5' dangling ends
						else if (localj!=number+1) {
							localrarray= SUM(localrarray,PROD(erg4(locali,localj,locali+1,1,ct,data,lfce[locali+1]),
								penaltyij,w3[locali+1],w5[localj-number-1] TWOSCALING));
							LOG_TRACE("[XXXXXX]\tlocalrarray_(v(i,j)): " << localrarray);
						}
						//no dangling ends
						else {
							localrarray= SUM(localrarray,PROD(w3[locali+1],w5[localj-number-1],penaltyij TWOSCALING));
							LOG_TRACE("[XXXXXX]\tlocalrarray_(v(i,j)): " << localrarray);
						}

						#else //not simplembloop
         				localrarray= SUM(localrarray,PROD(w3[locali+1],w5[localj-number-1],penaltyij TWOSCALING));
						LOG_TRACE("[XXXXXX]\tlocalrarray_(v(i,j)): " << localrarray);

						if (locali!=number) {
							localrarray= SUM(localrarray,PROD(ptr_to_erg4(locali,localj,locali+1,1,ct,data,lfce[locali+1]),penaltyij,w3[locali+2],w5[localj-number-1] TWOSCALING));
							LOG_TRACE("[XXXXXX]\tlocalrarray_(v(i,j)): " << localrarray);
						}
						if (localj!=(number+1)) {
							localrarray= SUM(localrarray,PROD(ptr_to_erg4(locali,localj,localj-1,2,ct,data,lfce[localj-1]),penaltyij,w3[locali+1],w5[localj-number-2] TWOSCALING));
							LOG_TRACE("[XXXXXX]\tlocalrarray_(v(i,j)): " << localrarray);
						}
						if ((locali!=number)&&(localj!=(number+1))) {
            				localrarray= SUM(localrarray,PROD(
// original						data->tstack[ct->numseq[locali]][ct->numseq[localj]][ct->numseq[locali+1]][ct->numseq[localj-1]],
								ptr_to_tstack(locali, localj, locali + 1,localj - 1, ct, data),
								pfchecknp(lfce[locali+1],lfce[localj-1]),w3[locali+2],
								w5[localj-number-2],penaltyij TWOSCALING));
							LOG_TRACE("[XXXXXX]\tlocalrarray_(v(i,j)): " << localrarray);
						}

						//consider the coaxial stacking of a helix from locali to localip onto helix localip+1 or localip+2 to localj:
						#ifndef disablecoax //a flag that can turn of coaxial stacking
						if(!disablecoax){
							//first consider a helix stacking from the 5' sequence fragment:
							for (int localip=localj-number-minloop-1;localip>0;localip--) {
								//first consider flush stacking
								localrarray=SUM(localrarray,PROD(
									w3[locali+1],w5[localip-1],penaltyij,ptr_to_penalty(localj-number-1,localip,ct,data),
									ptr_to_ergcoaxflushbases(localip,localj-number-1,localj-number,locali,ct,data),v->f(localip,localj-number-1) TWOSCALING));
								LOG_TRACE("[XXXXXX]\tlocalrarray_(v(i,j)): " << localrarray);

								if ((mod[localip]||mod[localj-number-1])) if (localj-number-2>0&&notgu(localip,localj-number-1,ct)&&!(fce->f(localip,localj-number-1)&SINGLE)) {
									if (ptr_to_can_pair(localip+1, localj-number-2, ct, data)) {
										localrarray=SUM(localrarray,PROD(
											w3[locali+1],w5[localip-1],penaltyij,ptr_to_penalty(localj-number-1,localip,ct,data),
											ptr_to_ergcoaxflushbases(localip,localj-number-1,localj-number,locali,ct,data),v->f(localip+1,localj-number-2),
											ptr_to_erg1(localip,localj-number-1,localip+1,localj-number-2,ct,data) TWOSCALING));
										LOG_TRACE("[XXXXXX]\tlocalrarray_(v(i,j)): " << localrarray);
									}
								}

								if (localj-number-2>0) {
									//now consider an intervening nuc
									if(locali<number) {
										localrarray=SUM(localrarray,PROD(
											w3[locali+2],w5[localip-1],penaltyij,ptr_to_penalty(localip,localj-number-2,ct,data),
											ptr_to_ergcoaxinterbases2(localip,localj-number-2,localj-number,locali,ct,data),v->f(localip,localj-number-2) TWOSCALING,
											pfchecknp(lfce[localj-number-1],lfce[locali+1])));
										LOG_TRACE("[XXXXXX]\tlocalrarray_(v(i,j)): " << localrarray);

										if ((mod[localip]||mod[localj-number-2])) if ( ptr_to_can_pair(localip+1, localj-number-3, ct, data) &&notgu(localip,localj-number-2,ct)
											&&!(fce->f(localip,localj-number-2)&SINGLE)) {
											localrarray=SUM(localrarray,PROD(
												w3[locali+2],w5[localip-1],penaltyij,ptr_to_penalty(localip,localj-number-2,ct,data),
												ptr_to_ergcoaxinterbases2(localip,localj-number-2,localj-number,locali,ct,data),v->f(localip+1,localj-number-3),
												ptr_to_erg1(localip,localj-number-2,localip+1,localj-number-3,ct,data) TWOSCALING,pfchecknp(lfce[localj-number-1],lfce[locali+1])));
											LOG_TRACE("[XXXXXX]\tlocalrarray_(v(i,j)): " << localrarray);
										}
									}

									//consider the other possibility for an intervening nuc
									localrarray=SUM(localrarray,PROD(
										w3[locali+1],w5[localip-1],penaltyij,ptr_to_penalty(localip+1,localj-number-2,ct,data),
										ptr_to_ergcoaxinterbases1(localip+1,localj-number-2,localj-number,locali,ct,data),v->f(localip+1,localj-number-2) TWOSCALING,
										pfchecknp(lfce[localj-number-1],lfce[localip])));
									LOG_TRACE("[XXXXXX]\tlocalrarray_(v(i,j)): " << localrarray);

									if ((mod[localip+1]||mod[localj-number-2])) if (ptr_to_can_pair(localip+2, localj-number-3, ct, data)&&notgu(localip+1,localj-number-2,ct)
										&&!(fce->f(localip+1,localj-number-2)&SINGLE)) {
										localrarray=SUM(localrarray,PROD(
											w3[locali+1],w5[localip-1],penaltyij,ptr_to_penalty(localip+1,localj-number-2,ct,data),
											ptr_to_ergcoaxinterbases1(localip+1,localj-number-2,localj-number,locali,ct,data),v->f(localip+2,localj-number-3),
											ptr_to_erg1(localip+1,localj-number-2,localip+2,localj-number-3,ct,data) TWOSCALING,
											pfchecknp(lfce[localj-number-1],lfce[localip])));
										LOG_TRACE("[XXXXXX]\tlocalrarray_(v(i,j)): " << localrarray);
									}
								}
							}

							//now consider a helix stacking from the 3' sequence fragment:
							for (int localip=locali+minloop+1;localip<=number;localip++) {
								//first consider flush stacking

								localrarray=SUM(localrarray,PROD(
									w3[localip+1],w5[localj-number-1],penaltyij,ptr_to_penalty(localip,locali+1,ct,data),
									ptr_to_ergcoaxflushbases(localj-number,locali,locali+1,localip,ct,data),v->f(locali+1,localip) TWOSCALING));
								LOG_TRACE("[XXXXXX]\tlocalrarray_(v(i,j)): " << localrarray);

								if ((mod[locali+1]||mod[localip])) if (ptr_to_can_pair(locali+2, localip-1, ct, data)&&notgu(locali+1,localip,ct)
									&&!(fce->f(locali+1,localip)&SINGLE)) {
									localrarray=SUM(localrarray,PROD(
										w3[localip+1],w5[localj-number-1],penaltyij,ptr_to_penalty(localip,locali+1,ct,data),
										ptr_to_ergcoaxflushbases(localj-number,locali,locali+1,localip,ct,data),v->f(locali+2,localip-1),
										ptr_to_erg1(locali+1,localip,locali+2,localip-1,ct,data) TWOSCALING));
									LOG_TRACE("[XXXXXX]\tlocalrarray_(v(i,j)): " << localrarray);
								}

								//now consider an intervening nuc
								if (localj-number>1) {
									localrarray=SUM(localrarray,PROD(
										w3[localip+1],w5[localj-number-2],penaltyij,ptr_to_penalty(localip,locali+2,ct,data),
										ptr_to_ergcoaxinterbases1(localj-number,locali,locali+2,localip,ct,data),v->f(locali+2,localip) TWOSCALING,
										pfchecknp(lfce[locali+1],lfce[localj-number-1])));
									LOG_TRACE("[XXXXXX]\tlocalrarray_(v(i,j)): " << localrarray);

									if ((mod[locali+2]||mod[localip])) if (ptr_to_can_pair(locali+3, localip-1, ct, data)&&notgu(locali+2,localip,ct)
										&&!(fce->f(locali+2,localip)&SINGLE)) {
										localrarray=SUM(localrarray,PROD(
											w3[localip+1],w5[localj-number-2],penaltyij,ptr_to_penalty(localip,locali+2,ct,data),
											ptr_to_ergcoaxinterbases1(localj-number,locali,locali+2,localip,ct,data),v->f(locali+3,localip-1),
											ptr_to_erg1(locali+2,localip,locali+3,localip-1,ct,data) TWOSCALING,
											pfchecknp(lfce[locali+1],lfce[localj-number-1])));
										LOG_TRACE("[XXXXXX]\tlocalrarray_(v(i,j)): " << localrarray);
									}
								}

								//consider the other possibility for an intervening nuc
								localrarray=SUM(localrarray,PROD(
									w3[localip+1],w5[localj-number-1],penaltyij,ptr_to_penalty(localip-1,locali+2,ct,data),
									ptr_to_ergcoaxinterbases2(localj-number,locali,locali+2,localip-1,ct,data),v->f(locali+2,localip-1) TWOSCALING,
									pfchecknp(lfce[locali+1],lfce[localip])));
								LOG_TRACE("[XXXXXX]\tlocalrarray_(v(i,j)): " << localrarray);

								if ((mod[locali+2]||mod[localip-1])) if(ptr_to_can_pair(locali+3, localip-2, ct, data)&&notgu(locali+2,localip-1,ct)
									&&!(fce->f(locali+2,localip-1)&SINGLE)) {
									localrarray=SUM(localrarray,PROD(
										w3[localip+1],w5[localj-number-1],penaltyij,ptr_to_penalty(localip-1,locali+2,ct,data),
										ptr_to_ergcoaxinterbases2(localj-number,locali,locali+2,localip-1,ct,data),v->f(locali+3,localip-2),
										ptr_to_erg1(locali+2,localip-1,locali+3,localip-2,ct,data) TWOSCALING,pfchecknp(lfce[locali+1],lfce[localip])));
									LOG_TRACE("[XXXXXX]\tlocalrarray_(v(i,j)): " << localrarray);
								}
							}
						}
						#endif //ifndef disablecoax
						#endif //simplembloop

					}//end of consider exterior loop closed by locali, localj

					//consider the multiloop closed by locali,localj
					//[7TI8Q8]
					if ((localj-locali)>(2*minloop+4)&&locali!=number) {
						#ifdef SIMPLEMBLOOP
						if (localj-1!=number&&ilocal!=number) {
							//5' and 3' dangling ends
							localrarray=SUM(localrarray,PROD(erg4(locali,localj,locali+1,1,ct,data,lfce[locali+1]),
								erg4(locali,localj,localj-1,2,ct,data,lfce[localj-1]),penaltyij,
            					wmb->f(locali+1,localj-1), data->eparam[5] , data->eparam[10] TWOSCALING));
							LOG_DEBUG4("[7TI8Q8]\tlocalrarray_(v(i,j)): " << localrarray);
						}
						else if (localj-1!=number) {//locali==number
							//5' dangling end
							localrarray=SUM(localrarray,PROD(
								erg4(locali,localj,localj-1,2,ct,data,lfce[localj-1]),penaltyij,
            					wmb->f(locali+1,localj-1), data->eparam[5] , data->eparam[10] TWOSCALING));
							LOG_DEBUG4("[7TI8Q8]\tlocalrarray_(v(i,j)): " << localrarray);
						}
						if (locali!=number) {
							//3' dangling ends
							localrarray=SUM(localrarray,PROD(erg4(locali,localj,locali+1,1,ct,data,lfce[locali+1]),
								penaltyij,
            					wmb->f(locali+1,localj-1), data->eparam[5] , data->eparam[10] TWOSCALING));
							LOG_DEBUG4("[7TI8Q8]\tlocalrarray_(v(i,j)): " << localrarray);
						}

						else {
							//no dangling end
							localrarray=SUM(localrarray,PROD(wmb->f(locali+1,localj-1),
								data->eparam[5],data->eparam[10],
            					penaltyij TWOSCALING));
							LOG_DEBUG4("[7TI8Q8]\tlocalrarray_(v(i,j)): " << localrarray);
						}
						#else //!SIMPLEMBLOOP

          				//no dangling ends on locali-localj pair: [M1RT7B]
						if (localj-1!=number) {
							localrarray=SUM(localrarray,PROD(wmb->f(locali+1,localj-1),multi_eparam(5, ct, data),multi_eparam(10, ct, data),
            					penaltyij TWOSCALING));
							LOG_DEBUG4("[M1RT7B]\tlocalrarray_(v(i,j)): " << localrarray);

							//locali+1 dangles on locali-localj pair: [L23JLU]

							if (locali+1!=number) {
								localrarray=SUM(localrarray,PROD(ptr_to_erg4(locali,localj,locali+1,1,ct,data,lfce[locali+1]),penaltyij,
            						wmb->f(locali+2,localj-1), multi_eparam(5, ct, data) , multi_eparam(6, ct, data) , multi_eparam(10, ct, data) TWOSCALING));
								LOG_DEBUG4("[L23JLU]\tlocalrarray_(v(i,j)): " << localrarray);
							}
						}
						if (localj-2!=number) {
							//localj-1 dangles [C4BBSC]
							if (localj!=(number+1)){
								localrarray=SUM(localrarray,PROD(ptr_to_erg4(locali,localj,localj-1,2,ct,data,lfce[localj-1]) , penaltyij ,
            						wmb->f(locali+1,localj-2) , multi_eparam(5, ct, data) , multi_eparam(6, ct, data) , multi_eparam(10, ct, data) TWOSCALING));
								LOG_DEBUG4("[C4BBSC]\tlocalrarray_(v(i,j)): " << localrarray);
							}
							//both locali+1 and localj-1 dangle
							if ((locali+1!=number)&&(localj!=(number+1))) {
            					localrarray=SUM(localrarray,PROD(
// original									data->tstkm[ct->numseq[locali]][ct->numseq[localj]][ct->numseq[locali+1]][ct->numseq[localj-1]],
											ptr_to_tstkm(locali, localj, locali + 1, localj - 1, ct, data),
											pfchecknp(lfce[locali+1],lfce[localj-1]),
											wmb->f(locali+2,localj-2), multi_eparam(5, ct, data), multi_eparam(6, ct, data), multi_eparam(6, ct, data), multi_eparam(10, ct, data),
											penaltyij TWOSCALING));
								LOG_DEBUG4("[LXCS43]\tlocalrarray_(v(i,j)): " << localrarray);
							}
						}

						//consider the coaxial stacking of a helix from locali to localj onto helix locali+1 or locali+2 to localip:
						//
						#ifndef disablecoax //a flag to turn off coaxial stacking
						if(!disablecoax){
							for (int localip=locali+1;(localip<localj);localip++) {
								LOG_DEBUG4("[Coax1.start]\tlocalip: " << localip);
								//first consider flush stacking

								//conditions guarantee that the coaxial stacking isn't considering an exterior loop
								//if ((locali!=number)/*&&(locali+1!=number)*//*&&((localj>number)||(localip!=number)&&(localip+1!=number))&&(localj-1!=number)*/) {
								if (locali!=number&&localip!=number&&localj-1!=number) {
									if (ptr_to_can_pair(locali+1, localip, ct, data)) {
										// #1
										//localrarray=SUM(localrarray,PROD(penaltyij, v->f(locali+1,localip), penalty(locali+1,localip,ct,data), data->eparam[5], data->eparam[10], data->eparam[10], SUM(w->f(localip+1,localj-1),wmb->f(localip+1,localj-1)), ergcoaxflushbases(localj,locali,locali+1,localip,ct,data) TWOSCALING));
										// [SMWF7W]
										// localrarray=SUM(localrarray,PROD(penaltyij,v->f(locali+1,localip), penalty(locali+1,localip,ct,data), data->eparam[5], data->eparam[10], data->eparam[10], SUM(w->f(localip+1,localj-1),wmb->f(localip+1,localj-1)), ergcoaxflushbases(localj,locali,locali+1,localip,ct,data) TWOSCALING));
										localrarray=SUM(localrarray,PROD(penaltyij,v->f(locali+1,localip), ptr_to_penalty(locali+1,localip,ct,data), multi_eparam(5, ct, data), multi_eparam(10, ct, data), multi_eparam(10, ct, data), w->f(localip+1,localj-1), ptr_to_ergcoaxflushbases(localj,locali,locali+1,localip,ct,data) TWOSCALING));
										localrarray=SUM(localrarray,PROD(penaltyij,v->f(locali+1,localip), ptr_to_penalty(locali+1,localip,ct,data), multi_eparam(5, ct, data), multi_eparam(10, ct, data), multi_eparam(10, ct, data), wmb->f(localip+1,localj-1), ptr_to_ergcoaxflushbases(localj,locali,locali+1,localip,ct,data) TWOSCALING));
										LOG_DEBUG4("[SMWF7W.2]\tlocalrarray_(v(i,j)): " << localrarray);

										if((mod[locali+1]||mod[localip])) if (ptr_to_can_pair(locali+2, localip-1, ct, data)&&notgu(locali+1,localip,ct)&&!(fce->f(locali+1,localip)&SINGLE)) {
											// #2
											//localrarray=SUM(localrarray,PROD(penaltyij, v->f(locali+2,localip-1), penalty(locali+1,localip,ct,data), data->eparam[5], data->eparam[10], data->eparam[10], SUM(w->f(localip+1,localj-1),wmb->f(localip+1,localj-1)), ergcoaxflushbases(localj,locali,locali+1,localip,ct,data), erg1(locali+1,localip,locali+2,localip-1,ct,data) TWOSCALING));
											
											localrarray=SUM(localrarray,PROD(penaltyij, v->f(locali+2,localip-1), ptr_to_penalty(locali+1,localip,ct,data), multi_eparam(5, ct, data), multi_eparam(10, ct, data), multi_eparam(10, ct, data), SUM(w->f(localip+1,localj-1), wmb->f(localip+1,localj-1)), ptr_to_ergcoaxflushbases(localj,locali,locali+1,localip,ct,data), ptr_to_erg1(locali+1,localip,locali+2,localip-1,ct,data) TWOSCALING));
											LOG_DEBUG4("[SMWF7W.4]\tlocalrarray_(v(i,j)): " << localrarray);

										}
									}

									//check if locali+2 and localip can pair
									if (ptr_to_can_pair(locali+2, localip, ct, data)) {
										//if ((localip<localj-1)&&(locali+2!=number)) {
										if (localip+2<localj-1&&locali+1!=number&&localip+1!=number) {
										//now consider an intervening nuc
											if ((localip+2<localj-1)/*&&(localj>number||localip+2!=number)*/) {
												// #3 [60M14C]
												// localrarray=SUM(localrarray,PROD(penaltyij,v->f(locali+2,localip), 
												// 				penalty(locali+2,localip,ct,data),data->eparam[5], 
												// 				data->eparam[6], data->eparam[6], data->eparam[10], 
												// 				data->eparam[10], SUM(w->f(localip+2,localj-1), 
												// 				wmb->f(localip+2,localj-1)), 
												// 				ergcoaxinterbases2(localj,locali,locali+2,localip,ct,data) TWOSCALING, 
												// 				pfchecknp(lfce[locali+1],lfce[localip+1])));
												localrarray=SUM(localrarray,PROD(penaltyij,v->f(locali+2,localip), 
																ptr_to_penalty(locali+2,localip,ct,data),multi_eparam(5, ct, data), 
																multi_eparam(6, ct, data), multi_eparam(6, ct, data), multi_eparam(10, ct, data), 
																multi_eparam(10, ct, data), w->f(localip+2,localj-1), 
																ptr_to_ergcoaxinterbases2(localj,locali,locali+2,localip,ct,data) TWOSCALING, 
																pfchecknp(lfce[locali+1],lfce[localip+1])));
												localrarray=SUM(localrarray,PROD(penaltyij,v->f(locali+2,localip), 
																ptr_to_penalty(locali+2,localip,ct,data),multi_eparam(5, ct, data), 
																multi_eparam(6, ct, data), multi_eparam(6, ct, data), multi_eparam(10, ct, data), 
																multi_eparam(10, ct, data),
																wmb->f(localip+2,localj-1), 
																ptr_to_ergcoaxinterbases2(localj,locali,locali+2,localip,ct,data) TWOSCALING, 
																pfchecknp(lfce[locali+1],lfce[localip+1])));
												LOG_DEBUG4("[60M14C.2]\tlocalrarray_(v(i,j)): " << localrarray);

												//localrarray=SUM(localrarray,PROD(penaltyij, v->f(locali+2,localip), penalty(locali+2,localip,ct,data), data->eparam[5], data->eparam[6], data->eparam[6], data->eparam[10], data->eparam[10], SUM(w->f(localip+2,localj-1),wmb->f(localip+2,localj-1)), ergcoaxinterbases2(localj,locali,locali+2,localip,ct,data), twoscaling, pfchecknp(lfce[locali+1],lfce[localip+1]);

												if((mod[locali+2]||mod[localip])) if (ptr_to_can_pair(locali+3, localip-1, ct, data)&&notgu(locali+2,localip,ct)
													&&!(fce->f(locali+2,localip)&SINGLE)) {
														// #4
													localrarray=SUM(localrarray,PROD(penaltyij,v->f(locali+3,localip-1),
														ptr_to_penalty(locali+2,localip,ct,data),multi_eparam(5, ct, data),
														multi_eparam(6, ct, data),multi_eparam(6, ct, data),multi_eparam(10, ct, data),multi_eparam(10, ct, data),
														SUM(w->f(localip+2,localj-1),wmb->f(localip+2,localj-1)),
														ptr_to_ergcoaxinterbases2(localj,locali,locali+2,localip,ct,data),
														ptr_to_erg1(locali+2,localip,locali+3,localip-1,ct,data) TWOSCALING, pfchecknp(lfce[locali+1],lfce[localip+1])));
													LOG_DEBUG4("[60M14C.4]\tlocalrarray_(v(i,j)): " << localrarray);

												}
											}

											if (localip+1<localj-2&&localj-2!=number) {
												// #5 [AWFSPY]
													// localrarray=SUM(localrarray,PROD(penaltyij,v->f(locali+2,localip),
													// 	penalty(locali+2,localip,ct,data),data->eparam[5],
													// 	data->eparam[6],data->eparam[6],data->eparam[10],data->eparam[10],
													// 	SUM(w->f(localip+1,localj-2),wmb->f(localip+1,localj-2)),
													// 	ergcoaxinterbases1(localj,locali,locali+2,localip,ct,data) TWOSCALING,
													// 	pfchecknp(lfce[locali+1],lfce[localj-1])));
													localrarray=SUM(localrarray,PROD(penaltyij,v->f(locali+2,localip),
														ptr_to_penalty(locali+2,localip,ct,data),multi_eparam(5, ct, data),
														multi_eparam(6, ct, data),multi_eparam(6, ct, data),multi_eparam(10, ct, data),multi_eparam(10, ct, data),
														w->f(localip+1,localj-2),
														ptr_to_ergcoaxinterbases1(localj,locali,locali+2,localip,ct,data) TWOSCALING,
														pfchecknp(lfce[locali+1],lfce[localj-1])));
													localrarray=SUM(localrarray,PROD(penaltyij,v->f(locali+2,localip),
														ptr_to_penalty(locali+2,localip,ct,data),multi_eparam(5, ct, data),
														multi_eparam(6, ct, data),multi_eparam(6, ct, data),multi_eparam(10, ct, data),multi_eparam(10, ct, data),
														wmb->f(localip+1,localj-2),
														ptr_to_ergcoaxinterbases1(localj,locali,locali+2,localip,ct,data) TWOSCALING,
														pfchecknp(lfce[locali+1],lfce[localj-1])));
													LOG_DEBUG4("[AWFSPY.2]\tlocalrarray_(v(i,j)): " << localrarray);

													if((mod[locali+2]||mod[localip])) if (ptr_to_can_pair(locali+3, localip-1, ct, data)&&notgu(locali+2,localip,ct)
														&&!(fce->f(locali+2,localip)&SINGLE)) {
															// #6
														localrarray=SUM(localrarray,PROD(penaltyij,v->f(locali+3,localip-1),
															ptr_to_penalty(locali+2,localip,ct,data),multi_eparam(5, ct, data),
															multi_eparam(6, ct, data),multi_eparam(6, ct, data),multi_eparam(10, ct, data),multi_eparam(10, ct, data),
															SUM(w->f(localip+1,localj-2),wmb->f(localip+1,localj-2)),
															ptr_to_ergcoaxinterbases1(localj,locali,locali+2,localip,ct,data),
															ptr_to_erg1(locali+2,localip,locali+3,localip-1,ct,data) TWOSCALING, pfchecknp(lfce[locali+1],lfce[localj-1])));
														LOG_DEBUG4("[AWFSPY.4]\tlocalrarray_(v(i,j)): " << localrarray);

													}
											}
										}
									}//end of check whether locali+2 and loaclip can pair
								}
							}
							LOG_DEBUG4("[Coax1]\tlocalrarray_(v(i,j)): " << localrarray);

							//consider the coaxial stacking of a helix from locali to localj onto helix localip to localj-2 or localj-1:
							for (int localip=localj-1;localip>locali;localip--) {
								LOG_DEBUG4("[Coax2.start]\tlocalip: " << localip);
								//conditions guarantee that the coaxial stacking isn't considering an exterior loop
								//if ((locali!=number)&&(locali+1!=number)&&((localj>number)||(localip!=number)&&(localip-1!=number))&&(localj-1!=number)) {
								if (localj-1!=number&&localip-1!=number&&locali!=number) {
									//first consider flush stacking
									// // #7 [T0FCK6]
									// localrarray=SUM(localrarray,PROD(penaltyij,v->f(localip,localj-1),
									// 	penalty(localj-1,localip,ct,data),data->eparam[5],
									// 	data->eparam[10],data->eparam[10],
									// 	SUM(w->f(locali+1,localip-1),wmb->f(locali+1,localip-1)),
									// 	ergcoaxflushbases(localip,localj-1,localj,locali,ct,data) TWOSCALING));
									localrarray=SUM(localrarray,PROD(penaltyij,v->f(localip,localj-1),
										ptr_to_penalty(localj-1,localip,ct,data),multi_eparam(5, ct, data),
										multi_eparam(10, ct, data),multi_eparam(10, ct, data),
										w->f(locali+1,localip-1),
										ptr_to_ergcoaxflushbases(localip,localj-1,localj,locali,ct,data) TWOSCALING));
									localrarray=SUM(localrarray,PROD(penaltyij,v->f(localip,localj-1),
										ptr_to_penalty(localj-1,localip,ct,data),multi_eparam(5, ct, data),
										multi_eparam(10, ct, data),multi_eparam(10, ct, data),
										wmb->f(locali+1,localip-1),
										ptr_to_ergcoaxflushbases(localip,localj-1,localj,locali,ct,data) TWOSCALING));
									LOG_DEBUG4("[T0FCK6.2]\tlocalrarray_(v(i,j)): " << localrarray);

									if((mod[localip]||mod[localj-1])) if(ptr_to_can_pair(localip+1, localj-2, ct, data)&&notgu(localip,localj-1,ct)&&!(fce->f(localip,localj-1)&SINGLE)) {
										// #8
										localrarray=SUM(localrarray,PROD(penaltyij,v->f(localip+1,localj-2),
											ptr_to_penalty(localj-1,localip,ct,data),multi_eparam(5, ct, data),
											multi_eparam(10, ct, data), multi_eparam(10, ct, data), SUM(w->f(locali+1,localip-1),wmb->f(locali+1,localip-1)),
											ptr_to_ergcoaxflushbases(localip,localj-1,localj,locali,ct,data),
											ptr_to_erg1(localip,localj-1,localip+1,localj-2,ct,data) TWOSCALING));
										LOG_DEBUG4("[T0FCK6.4]\tlocalrarray_(v(i,j)): " << localrarray);

									}

									if (localj-2!=number) {
										//now consider an intervening nuc [332V5Q]
										//if ((localip>locali+1)&&(localj>number||localip-2!=number))
										if (localip-2>locali+1&&localip-2!=number) {
											// #9
											// localrarray=SUM(localrarray,PROD(penaltyij,v->f(localip,localj-2),
											// 	penalty(localj-2,localip,ct,data),data->eparam[5],
											// 	data->eparam[6],data->eparam[6],data->eparam[10],data->eparam[10],
											// 	SUM(w->f(locali+1,localip-2),wmb->f(locali+1,localip-2)),
											// 	ergcoaxinterbases1(localip,localj-2,localj,locali,ct,data) TWOSCALING,pfchecknp(lfce[localj-1],lfce[localip-1])));
											localrarray=SUM(localrarray,PROD(penaltyij,v->f(localip,localj-2),
												ptr_to_penalty(localj-2,localip,ct,data),multi_eparam(5, ct, data),
												multi_eparam(6, ct, data),multi_eparam(6, ct, data),multi_eparam(10, ct, data),multi_eparam(10, ct, data),
												w->f(locali+1,localip-2),
												ptr_to_ergcoaxinterbases1(localip,localj-2,localj,locali,ct,data) TWOSCALING,pfchecknp(lfce[localj-1],lfce[localip-1])));
											localrarray=SUM(localrarray,PROD(penaltyij,v->f(localip,localj-2),
												ptr_to_penalty(localj-2,localip,ct,data),multi_eparam(5, ct, data),
												multi_eparam(6, ct, data),multi_eparam(6, ct, data),multi_eparam(10, ct, data),multi_eparam(10, ct, data),
												wmb->f(locali+1,localip-2),
												ptr_to_ergcoaxinterbases1(localip,localj-2,localj,locali,ct,data) TWOSCALING,pfchecknp(lfce[localj-1],lfce[localip-1])));
											LOG_DEBUG4("[332V5Q.2]\tlocalrarray_(v(i,j)): " << localrarray);
											// LOG_DEBUG4("[332V5Q.2]\tw(i+1, ip-2): " << w->f(locali+1,localip-2))
											// LOG_DEBUG4("[332V5Q.2]\twmb(i+1, ip-2): " << wmb->f(locali+1,localip-2))

											if((mod[localip]||mod[localj-2])) if(ptr_to_can_pair(localip+1, localj-3, ct, data)&&notgu(localip,localj-2,ct)&&!(fce->f(localip,localj-2)&SINGLE)) {
												// #10
												localrarray=SUM(localrarray,PROD(penaltyij,v->f(localip+1,localj-3),
													ptr_to_penalty(localj-2,localip,ct,data),multi_eparam(5, ct, data),
													multi_eparam(6, ct, data),multi_eparam(6, ct, data),multi_eparam(10, ct, data),multi_eparam(10, ct, data),
													SUM(w->f(locali+1,localip-2),wmb->f(locali+1,localip-2)),
													ptr_to_ergcoaxinterbases1(localip,localj-2,localj,locali,ct,data),
													ptr_to_erg1(localip,localj-2,localip+1,localj-3,ct,data) TWOSCALING,pfchecknp(lfce[localj-1],lfce[localip-1])));
												LOG_DEBUG4("[332V5Q.4]\tlocalrarray_(v(i,j)): " << localrarray);

											}
										}

										if ((localip-1>locali+2)&&locali+1!=number) {
											// #11 [Q4LM6F]
											// localrarray=SUM(localrarray,PROD(penaltyij,v->f(localip,localj-2),
											// 	penalty(localj-2,localip,ct,data),data->eparam[5],
											// 	data->eparam[6],data->eparam[6],data->eparam[10],data->eparam[10],
											// 	SUM(w->f(locali+2,localip-1),wmb->f(locali+2,localip-1)),
											// 	ergcoaxinterbases2(localip,localj-2,localj,locali,ct,data) TWOSCALING,pfchecknp(lfce[localj-1],lfce[locali+1])));
											localrarray=SUM(localrarray,PROD(penaltyij,v->f(localip,localj-2),
												ptr_to_penalty(localj-2,localip,ct,data),multi_eparam(5, ct, data),
												multi_eparam(6, ct, data),multi_eparam(6, ct, data),multi_eparam(10, ct, data),multi_eparam(10, ct, data),
												w->f(locali+2,localip-1),
												ptr_to_ergcoaxinterbases2(localip,localj-2,localj,locali,ct,data) TWOSCALING,pfchecknp(lfce[localj-1],lfce[locali+1])));
											localrarray=SUM(localrarray,PROD(penaltyij,v->f(localip,localj-2),
												ptr_to_penalty(localj-2,localip,ct,data),multi_eparam(5, ct, data),
												multi_eparam(6, ct, data),multi_eparam(6, ct, data),multi_eparam(10, ct, data),multi_eparam(10, ct, data),
												wmb->f(locali+2,localip-1),
												ptr_to_ergcoaxinterbases2(localip,localj-2,localj,locali,ct,data) TWOSCALING,pfchecknp(lfce[localj-1],lfce[locali+1])));
											LOG_DEBUG4("[Q4LM6F.2]\tlocalrarray_(v(i,j)): " << localrarray);

											if((mod[localip]||mod[localj-2])) if(ptr_to_can_pair(localip+1, localj-3, ct, data)&&notgu(localip,localj-2,ct)
												&&!(fce->f(localip,localj-2)&SINGLE)) {
													// #12
												localrarray=SUM(localrarray,PROD(penaltyij,v->f(localip+1,localj-3),
													ptr_to_penalty(localj-2,localip,ct,data),multi_eparam(5, ct, data),
													multi_eparam(6, ct, data),multi_eparam(6, ct, data),multi_eparam(10, ct, data),multi_eparam(10, ct, data),
													SUM(w->f(locali+2,localip-1),wmb->f(locali+2,localip-1)),
													ptr_to_ergcoaxinterbases2(localip,localj-2,localj,locali,ct,data),
													ptr_to_erg1(localip,localj-2,localip+1,localj-3,ct,data) TWOSCALING,pfchecknp(lfce[localj-1],lfce[locali+1])));
												LOG_DEBUG4("[Q4LM6F.4]\tlocalrarray_(v(i,j)): " << localrarray);

											}
										}
									}
								}
							}
							LOG_DEBUG4("[Coax2]\tlocalrarray_(v(i,j)): " << localrarray);
						}
						#endif //ifndef disablecoax
						#endif //SIMPLEMBLOOP
					}//end of consider the multibranch loop closed by locali, localj
				}//end of large block for multibranch or exterior closure
			}//end of the large block for calculating V if all conditions are met

//			cout << locali << ",\t" << localj << ":\t" << XLOG_TO_LINEAR(localrarray) << endl;

			v->f(locali,localj) = localrarray;
			DEBUG_LOG(DEBUG2, "[XXXXXX]\tlocali: " << locali << "\tlocalj: " << localj << "\tv->f(locali,localj): " << localrarray);

			//Apply constant (an euilibrium constant for formation of locali-localj pair), if being used:
			if (ct->constant!=NULL) {
				//use ii and jj as the indices for accessing constant array:

				if (locali>ct->GetSequenceLength()) ii = locali - ct->GetSequenceLength();
				else ii = locali;
				if (localj>ct->GetSequenceLength()) jj = localj - ct->GetSequenceLength();
				else jj = localj;
				if (jj<ii) {
					p = jj;
      				jj = ii;
					ii = p;
				}

				v->f(locali,localj) = PROD(v->f(locali,localj),ct->constant[jj][ii]);

			}

			
			if (fce->f(locali,localj)&PAIR)  {//force a pair between locali and localj
	  			w->f(locali,localj) = PROD(v->f(locali,localj),multi_eparam(10, ct, data),penaltyij);
				wl->f(locali,localj) = PROD(v->f(locali,localj),multi_eparam(10, ct, data),penaltyij);
				wl_T.f(locali,localj) = PROD(v->f(locali,localj),multi_eparam(10, ct, data),penaltyij);
				wmb->f(locali,localj) = ZERO;
				wmbl->f(locali,localj) = ZERO;
				wmbl_T.f(locali,localj) = ZERO;
				wcoax->f(locali,localj) = ZERO;
				if (localj<=number) wca[locali][localj] = (PFPRECISION) ZERO;
			}
			else {//not a forced pair

				////fill wmb:
				localrarray = (PFPRECISION) ZERO;
				LOG_DEBUG4("[XXXXXX]\tlocalrarray_(wcoax, wmbl): " << localrarray);
				if (((localj-locali)>(2*minloop+2))||localj>number) {
					#ifdef pfdebugmode
						ofstream dump;
						if (locali==10&&localj==22) {
							dump.open("dump.out");

						}
					#endif

					//also consider the coaxial stacking of two helixes in wv
					locale = (PFPRECISION) ZERO;
					#ifndef SIMPLEMBLOOP
					#ifndef disablecoax //a flag to diable coaxial stacking
					if(!disablecoax){
						for (int localip=locali+minloop+1;localip<localj-minloop-1;localip++) {
							//first consider flush stacking
							if (localip!=number) {
								localrarray=SUM(localrarray,PROD(v->f(locali,localip),v->f(localip+1,localj),ptr_to_penalty(locali,localip,ct,data),
									ptr_to_penalty(localip+1,localj,ct,data), ptr_to_ergcoaxflushbases(locali,localip,localip+1,localj,ct,data)));
								LOG_DEBUG4("[TYG2WQ.1]\tlocalrarray_(wcoax,wmbl): " << localrarray)	;						

								if ((mod[locali]||mod[localip]||mod[localip+1]||mod[localj])) {
									if ((mod[locali]||mod[localip])&&(mod[localip+1]||mod[localj])&&ptr_to_can_pair(locali+1, localip-1, ct, data)
										&&ptr_to_can_pair(localip+2, localj-1, ct, data)&&notgu(locali,localip,ct)&&notgu(localip+1,localj,ct)
											&&!(fce->f(localip+1,localj)&SINGLE)&&!(fce->f(locali,localip)&SINGLE)) {
										localrarray=SUM(localrarray,PROD(v->f(locali+1,localip-1), v->f(localip+2,localj-1), ptr_to_penalty(locali,localip,ct,data),
											ptr_to_penalty(localip+1,localj,ct,data), ptr_to_ergcoaxflushbases(locali,localip,localip+1,localj,ct,data),
											ptr_to_erg1(locali,localip,locali+1,localip-1,ct,data), ptr_to_erg1(localip+1,localj,localip+2,localj-1,ct,data)));
										LOG_DEBUG4("[TYG2WQ.2]\tlocalrarray_(wcoax,wmbl): " << localrarray);
									}

									if ((mod[locali]||mod[localip])&&ptr_to_can_pair(locali+1, localip-1, ct, data)&&notgu(locali,localip,ct)&&!(fce->f(locali,localip)&SINGLE)) {
										localrarray=SUM(localrarray,PROD(v->f(locali+1,localip-1), v->f(localip+1,localj), ptr_to_penalty(locali,localip,ct,data),
											ptr_to_penalty(localip+1,localj,ct,data), ptr_to_ergcoaxflushbases(locali,localip,localip+1,localj,ct,data),
											ptr_to_erg1(locali,localip,locali+1,localip-1,ct,data)));
										LOG_DEBUG4("[TYG2WQ.3]\tlocalrarray_(wcoax,wmbl): " << localrarray);

									}

									if ((mod[localip+1]||mod[localj])&&ptr_to_can_pair(localip+2, localj-1, ct, data)&&notgu(localip+1,localj,ct)&&!(fce->f(localip+1,localj)&SINGLE)) {
										localrarray=SUM(localrarray,PROD(v->f(locali,localip), v->f(localip+2,localj-1), ptr_to_penalty(locali,localip,ct,data),
											ptr_to_penalty(localip+1,localj,ct,data), ptr_to_ergcoaxflushbases(locali,localip,localip+1,localj,ct,data),
											ptr_to_erg1(localip+1,localj,localip+2,localj-1,ct,data)));
										LOG_DEBUG4("[TYG2WQ.4]\tlocalrarray_(wcoax,wmbl): " << localrarray);
									}
								}

								if (localip+1!=number&&localj!=number+1) {
									if (!lfce[localip+1]&&!lfce[localj]) {
										//now consider an intervening mismatch
										locale=SUM(locale,PROD(v->f(locali,localip), v->f(localip+2,localj-1), ptr_to_penalty(locali,localip,ct,data),
											ptr_to_penalty(localip+2,localj-1,ct,data), ptr_to_ergcoaxinterbases2(locali,localip,localip+2,localj-1,ct,data)));

										if (mod[locali]||mod[localip]||mod[localip+2]||mod[localj-1]) {
											if ((mod[locali]||mod[localip])&&(mod[localip+2]||mod[localj-1])&&ptr_to_can_pair(locali+1, localip-1, ct, data)
												&&ptr_to_can_pair(localip+3, localj-2, ct, data)&&notgu(locali,localip,ct)&&notgu(localip+2,localj-1,ct)
													&&!(fce->f(locali,localip)&SINGLE)&&!(fce->f(localip+2,localj-1)&SINGLE)) {
												 locale=SUM(locale,PROD(v->f(locali+1,localip-1), v->f(localip+3,localj-2), ptr_to_penalty(locali,localip,ct,data),
													ptr_to_penalty(localip+2,localj-1,ct,data), ptr_to_ergcoaxinterbases2(locali,localip,localip+2,localj-1,ct,data),
													ptr_to_erg1(locali,localip,locali+1,localip-1,ct,data), ptr_to_erg1(localip+2,localj-1,localip+3,localj-2,ct,data)));

											}

											if ((mod[locali]||mod[localip])&&ptr_to_can_pair(locali+1, localip-1, ct, data)&&notgu(locali,localip,ct)&&!(fce->f(locali,localip)&SINGLE)) {
												locale=SUM(locale,PROD(v->f(locali+1,localip-1), v->f(localip+2,localj-1), ptr_to_penalty(locali,localip,ct,data),
													ptr_to_penalty(localip+2,localj-1,ct,data), ptr_to_ergcoaxinterbases2(locali,localip,localip+2,localj-1,ct,data),
													ptr_to_erg1(locali,localip,locali+1,localip-1,ct,data)));

											}

											if ((mod[localip+2]||mod[localj-1])&&ptr_to_can_pair(localip+3, localj-2, ct, data)&&notgu(localip+2,localj-1,ct)&&!(fce->f(localip+2,localj-1)&SINGLE)) {
												locale=SUM(locale,PROD(v->f(locali,localip), v->f(localip+3,localj-2), ptr_to_penalty(locali,localip,ct,data),
													ptr_to_penalty(localip+2,localj-1,ct,data), ptr_to_ergcoaxinterbases2(locali,localip,localip+2,localj-1,ct,data),
													ptr_to_erg1(localip+2,localj-1,localip+3,localj-2,ct,data)));

											}
										}
									}

									if(!lfce[locali]&&!lfce[localip+1]&&locali!=number) {
										locale=SUM(locale,PROD(v->f(locali+1,localip), v->f(localip+2,localj), ptr_to_penalty(locali+1,localip,ct,data),
											ptr_to_penalty(localip+2,localj,ct,data), ptr_to_ergcoaxinterbases1(locali+1,localip,localip+2,localj,ct,data)));

										if (mod[locali+1]||mod[localip]||mod[localip+2]||mod[localj]) {
											if ((mod[locali+1]||mod[localip])&&(mod[localip+2]||mod[localj])&&ptr_to_can_pair(locali+2, localip-1, ct, data)
												&&ptr_to_can_pair(localip+3, localj-1, ct, data)&&notgu(locali+1,localip,ct)&&notgu(localip+2,localj,ct)
												&&!(fce->f(locali+1,localip)&SINGLE)&&!(fce->f(localip+2,localj)&SINGLE)	) {
												locale=SUM(locale,PROD(v->f(locali+2,localip-1), v->f(localip+3,localj-1), ptr_to_penalty(locali+1,localip,ct,data),
													ptr_to_penalty(localip+2,localj,ct,data), ptr_to_ergcoaxinterbases1(locali+1,localip,localip+2,localj,ct,data),
													ptr_to_erg1(locali+1,localip,locali+2,localip-1,ct,data), ptr_to_erg1(localip+2,localj,localip+3,localj-1,ct,data)));

											}
											if ((mod[locali+1]||mod[localip])&&ptr_to_can_pair(locali+2, localip-1, ct, data)&&notgu(locali+1,localip,ct)&&!(fce->f(locali+1,localip)&SINGLE)) {
												locale=SUM(locale,PROD(v->f(locali+2,localip-1), v->f(localip+2,localj), ptr_to_penalty(locali+1,localip,ct,data),
													ptr_to_penalty(localip+2,localj,ct,data), ptr_to_ergcoaxinterbases1(locali+1,localip,localip+2,localj,ct,data),
													ptr_to_erg1(locali+1,localip,locali+2,localip-1,ct,data)));

											}

											if ((mod[localip+2]||mod[localj])&&ptr_to_can_pair(localip+3, localj-1, ct, data)&&notgu(localip+2,localj,ct)&&!(fce->f(localip+2,localj)&SINGLE)) {
												locale=SUM(locale,PROD(v->f(locali+1,localip), v->f(localip+3,localj-1), ptr_to_penalty(locali+1,localip,ct,data),
													ptr_to_penalty(localip+2,localj,ct,data), ptr_to_ergcoaxinterbases1(locali+1,localip,localip+2,localj,ct,data),
													ptr_to_erg1(localip+2,localj,localip+3,localj-1,ct,data)));

											}
										}
									}
								}
							}
						}
					}
					#endif //ifndef disablecoax
					#endif //ifndef SIMPLEMBLOOP

					if (localj<=number) wca[locali][localj] = SUM(localrarray,locale);

					localrarray =PROD(SUM(localrarray,PROD(locale,multi_eparam(6, ct, data),multi_eparam(6, ct, data))),multi_eparam(10, ct, data),multi_eparam(10, ct, data));
					wcoax->f(locali,localj) = localrarray;
					LOG_DEBUG3("304YK0]\tlocalrarray_(wcoax, wmbl): " << localrarray);

					//search for an open bifurcation:
					int end = min(localj, number);
					if (!lfce[locali]&&locali!=number){
						for (int k=locali+1;k<end;k++) {
							localrarray=SUM(localrarray,PROD(SUM(wlc->dg[locali][k], wcoax->dg[locali][k]), SUM(wl_T.dg[localj][k+1], wmbl_T.dg[localj][k+1])));
							LOG_DEBUG4("[0LNQ9E]\tlocalrarray_(wmbl): " << localrarray);
							// localrarray=SUM(localrarray,PROD(wlc->dg[locali][k], wl_T.dg[localj][k+1]));
							// LOG_DEBUG4("[0LNQ9E.1]\tlocalrarray_(wmbl): " << localrarray)
							// localrarray=SUM(localrarray,PROD(wcoax->dg[locali][k], wl_T.dg[localj][k+1]));
							// LOG_DEBUG4("[0LNQ9E.2]\tlocalrarray_(wmbl): " << localrarray)
							// localrarray=SUM(localrarray,PROD(wlc->dg[locali][k], wmbl_T.dg[localj][k+1]));
							// LOG_DEBUG4("[0LNQ9E.3]\tlocalrarray_(wmbl): " << localrarray)
							// localrarray=SUM(localrarray,PROD(wcoax->dg[locali][k], wmbl_T.dg[localj][k+1]));
							// LOG_DEBUG4("[0LNQ9E.4]\tlocalrarray_(wmbl): " << localrarray)
							//localrarray+=(wl->dg[locali][k] - wl->dg[locali+1][k] * data->eparam[6] * data->scaling + wcoax->dg[locali][k]) * (wl_T.dg[localj][k+1]+wmbl_T.dg[localj][k+1]);
						}
						for (int k=number+1;k<localj;k++) {
							localrarray=SUM(localrarray,PROD(SUM(wlc->dg[locali][k], wcoax->dg[locali][k]),SUM(wl_T.dg[localj-number][k+1-number],wmbl_T.dg[localj-number][k+1-number])));
							LOG_DEBUG4("[XXXXXX]\tlocalrarray_(wmbl): " << localrarray);
							//localrarray+=(wl->dg[locali][k] - wl->dg[locali+1][k] * data->eparam[6] * data->scaling + wcoax->dg[locali][k]) * (wl_T.dg[localj-number][k+1-number]+wmbl_T.dg[localj-number][k+1-number]);
						}
					}
					else {
						for (int k=locali+1;k<end;k++) {
							localrarray=SUM(localrarray,PROD(SUM(wl->dg[locali][k], wcoax->dg[locali][k]), SUM(wl_T.dg[localj][k+1], wmbl_T.dg[localj][k+1])));
							LOG_DEBUG4("[NQ84MW]\tlocalrarray_(wmbl): " << localrarray);
							// localrarray=SUM(localrarray,PROD(wl->dg[locali][k], wl_T.dg[localj][k+1]));
							// LOG_DEBUG4("[NQ84MW.1]\tlocalrarray_(wmbl): " << localrarray)
							// localrarray=SUM(localrarray,PROD(wcoax->dg[locali][k], wl_T.dg[localj][k+1]));
							// LOG_DEBUG4("[NQ84MW.2]\tlocalrarray_(wmbl): " << localrarray)
							// localrarray=SUM(localrarray,PROD(wl->dg[locali][k], wmbl_T.dg[localj][k+1]));
							// LOG_DEBUG4("[NQ84MW.3]\tlocalrarray_(wmbl): " << localrarray)
							// localrarray=SUM(localrarray,PROD(wcoax->dg[locali][k], wmbl_T.dg[localj][k+1]));
							// LOG_DEBUG4("[NQ84MW.4]\tlocalrarray_(wmbl): " << localrarray)
						//	localrarray+=(wl->dg[locali][k] + wcoax->dg[locali][k]) * (wl_T.dg[localj][k+1] + wmbl_T.dg[localj][k+1]);
						}
						for (int k=number+1;k<localj;k++) {
							localrarray=SUM(localrarray,PROD(SUM(wl->dg[locali][k], wcoax->dg[locali][k]), SUM(wl_T.dg[localj-number][k+1-number], wmbl_T.dg[localj-number][k+1-number])));
							LOG_DEBUG4("[XXXXXX]\tlocalrarray_(wmbl): " << localrarray);
						//	localrarray+=(wl->dg[locali][k] + wcoax->dg[locali][k]) * (wl_T.dg[localj-number][k+1-number] + wmbl_T.dg[localj-number][k+1-number]);
						}
					}

					if (locali!=number){
						if (!lfce[locali]) {
							localrarray=SUM(localrarray,PROD(wmbl->f(locali+1,localj),multi_eparam(6, ct, data) SCALING));
							LOG_DEBUG4("[NYT5QE]\tlocalrarray_(wmbl): " << localrarray);
						}
					}
					wmbl->f(locali,localj) = localrarray;
					wmbl_T.f(locali,localj) = localrarray;

					wmb->f(locali,localj) = localrarray;
					if (localj!=number+1)
						if (!lfce[localj]) wmb->f(locali,localj)=SUM(localrarray,PROD(wmb->f(locali,localj-1),multi_eparam(6, ct, data) SCALING));
				}

				//Compute w[locali][localj]
				if (localj>number||(localj-locali>minloop)) {
					#ifdef SIMPLEMBLOOP

  					//calculate the energy of localj stacked onto the pair of locali,localj-1
					if (localj!=N&&locali!=1) {
						//5' and 3' dangling ends
     					wl->f(locali,localj)= PROD(v->f(locali,localj), data->eparam[10],
     						erg4(localj,locali,localj+1,1,ct,data,lfce[localj+1]),
							erg4(localj,locali,locali-1,2,ct,data,lfce[locali-1]),penaltyij);
						LOG_DEBUG4("[XXXXXX]\twl->f(locali,localj): " << wl->f(locali,localj));

						if ((mod[locali]||mod[localj])) if(inc[ct->numseq[locali+1]][ct->numseq[localj-1]]&&notgu(locali,localj,ct)&&!(fce->f(locali,localj)&SINGLE)) {
							wl->f(locali,localj)= SUM(wl->f(locali,localj),PROD(v->f(locali+1,localj-1), data->eparam[10],
     							erg4(localj,locali,localj+1,1,ct,data,lfce[localj+1]),
								erg4(localj,locali,locali-1,2,ct,data,lfce[locali-1]),penalty(locali,localj-1,ct,data),
								erg1(locali,localj,locali+1,localj-1,ct,data)));
							LOG_DEBUG4("[XXXXXX]\twl->f(locali,localj): " << wl->f(locali,localj));

						}
					}
					if (localj!=N) {//locali then ==1
						//3' dangling end
     					wl->f(locali,localj)= PROD(v->f(locali,localj), data->eparam[10],
     						erg4(localj,locali,localj+1,1,ct,data,lfce[localj+1]),
							penaltyij);
						LOG_DEBUG4("[XXXXXX]\twl->f(locali,localj): " << wl->f(locali,localj));

						if ((mod[locali]||mod[localj])) if(inc[ct->numseq[locali+1]][ct->numseq[localj-1]]&&notgu(locali,localj,ct)&&!(fce->f(locali,localj)&SINGLE)) {
							wl->f(locali,localj)= SUM(wl->f(locali,localj),PROD( v->f(locali+1,localj-1), data->eparam[10],
     							erg4(localj,locali,localj+1,1,ct,data,lfce[localj+1]),
								penalty(locali,localj-1,ct,data),
								erg1(locali,localj,locali+1,localj-1,ct,data)));
							LOG_DEBUG4("[XXXXXX]\twl->f(locali,localj): " << wl->f(locali,localj));

						}
					}
					else if (locali!=1) {//then localj==N
						//5' dangling end
     					wl->f(locali,localj)= PROD(v->f(locali,localj), data->eparam[10],
							erg4(localj,locali,locali-1,2,ct,data,lfce[locali-1]),penaltyij);
						LOG_DEBUG4("[XXXXXX]\twl->f(locali,localj): " << wl->f(locali,localj));

						if ((mod[locali]||mod[localj])) if(inc[ct->numseq[locali+1]][ct->numseq[localj-1]]&&notgu(locali,localj,ct)&&!(fce->f(locali,localj)&SINGLE)) {
							wl->f(locali,localj)= SUM(wl->f(locali,localj),PROD( v->f(locali+1,localj-1), data->eparam[10],
								erg4(localj,locali,locali-1,2,ct,data,lfce[locali-1]),penalty(locali,localj-1,ct,data),
								erg1(locali,localj,locali+1,localj-1,ct,data)));
							LOG_DEBUG4("[XXXXXX]\twl->f(locali,localj): " << wl->f(locali,localj));

						}
					}

					else {//locali==1 and localj==N
						wl->f(locali,localj)= PROD(data->eparam[10],v->f(locali,localj),penalty(localj,locali,ct,data));
						LOG_DEBUG4("[XXXXXX]\twl->f(locali,localj): " << wl->f(locali,localj));

						if ((mod[locali]||mod[localj])) if (inc[ct->numseq[locali+1]][ct->numseq[localj-1]]&&notgu(locali,localj,ct)) {
							wl->f(locali,localj)= SUM(wl->f(locali,localj),PROD( data->eparam[10], v->f(locali+1,localj-1), penalty(localj,locali,ct,data), erg1(locali,localj,locali+1,localj-1,ct,data)));
							LOG_DEBUG4("[XXXXXX]\twl->f(locali,localj): " << wl->f(locali,localj));
						}
					}
					#else  //!SIMPLEMBLOOP

					wl->f(locali,localj)= PROD(multi_eparam(10, ct, data),v->f(locali,localj),ptr_to_penalty(localj,locali,ct,data));
					LOG_DEBUG4("[89AGXL]\twl->f(locali,localj): " << wl->f(locali,localj));

					if ((mod[locali]||mod[localj])) if (ptr_to_can_pair(locali+1, localj-1, ct, data)&&notgu(locali,localj,ct)) {
						wl->f(locali,localj)= SUM(wl->f(locali,localj),PROD( multi_eparam(10, ct, data), v->f(locali+1,localj-1), ptr_to_penalty(localj,locali,ct,data), ptr_to_erg1(locali,localj,locali+1,localj-1,ct,data)));
						LOG_DEBUG4("[89AGXL.2]\twl->f(locali,localj): " << wl->f(locali,localj));
					}

					if (locali!=number) {
      					//calculate the energy of locali stacked onto the pair of locali+1,localj
						 // if (v->f(locali+1,localj)*data->eparam[10]*data->eparam[6]*
         				//	erg4(localj,locali+1,locali,2,ct,data,lfce[locali])*penalty(locali+1,localj,ct,data)<0){
							//	 cout << endl;
							//	 cout << v->f(locali+1,localj) << endl;
							//	 cout << wl->f(locali,localj)._value << endl;
							// 	 cout << data->eparam[6] << endl;
							// 	 cout << erg4(localj,locali+1,locali,2,ct,data,lfce[locali]) << endl;
							// 	 cout << penalty(locali+1,localj,ct,data) << endl;
							// 	 cout << v->f(locali+1,localj)*data->eparam[10]*data->eparam[6]*
         					// erg4(localj,locali+1,locali,2,ct,data,lfce[locali])*penalty(locali+1,localj,ct,data) << endl;

						//	 }
						//[J1C6M9]
						wl->f(locali,localj)= SUM(wl->f(locali,localj),PROD( v->f(locali+1,localj),multi_eparam(10, ct, data),multi_eparam(6, ct, data),
         					ptr_to_erg4(localj,locali+1,locali,2,ct,data,lfce[locali]),ptr_to_penalty(locali+1,localj,ct,data)));
						LOG_DEBUG4("[G98UOW]\twl->f(locali,localj): " << wl->f(locali,localj));

						if ((mod[locali+1]||mod[localj])) if(ptr_to_can_pair(locali+2, localj-1,ct,data)&&notgu(locali+1,localj,ct)&&!(fce->f(locali+1,localj)&SINGLE)) {
							wl->f(locali,localj)= SUM(wl->f(locali,localj),PROD( v->f(locali+2,localj-1), multi_eparam(10, ct, data), multi_eparam(6, ct, data),
         						ptr_to_erg4(localj,locali+1,locali,2,ct,data,lfce[locali]), ptr_to_penalty(locali+1,localj,ct,data),
								ptr_to_erg1(locali+1,localj,locali+2,localj-1,ct,data)));
							LOG_DEBUG4("[G98UOW.2]\twl->f(locali,localj): " << wl->f(locali,localj));
						}
					}
					if (localj!=((number)+1)) {
      					//calculate the energy of localj stacked onto the pair of locali,localj-1
						if (localj!=1) {
							// [LX1N9X]
         					wl->f(locali,localj)= SUM(wl->f(locali,localj),PROD(v->f(locali,localj-1), multi_eparam(10, ct, data), multi_eparam(6, ct, data),
         						ptr_to_erg4(localj-1,locali,localj,1,ct,data,lfce[localj]), ptr_to_penalty(locali,localj-1,ct,data)));
							LOG_DEBUG4("[YBN850]\twl->f(locali,localj): " << wl->f(locali,localj));

							if ((mod[locali]||mod[localj-1])) if(ptr_to_can_pair(locali+1, localj-2,ct,data)&&notgu(locali,localj-1,ct)&&!(fce->f(locali,localj-1)&SINGLE)) {
								wl->f(locali,localj)= SUM(wl->f(locali,localj),PROD(v->f(locali+1,localj-2), multi_eparam(10, ct, data), multi_eparam(6, ct, data),
         							ptr_to_erg4(localj-1,locali,localj,1,ct,data,lfce[localj]),ptr_to_penalty(locali,localj-1,ct,data),
									ptr_to_erg1(locali,localj-1,locali+1,localj-2,ct,data)));
								LOG_DEBUG4("[YBN850.2]\twl->f(locali,localj): " << wl->f(locali,localj));
							}
						}
					}
					if ((locali!=(number))&&(localj!=((number)+1))) {
      					//calculate locali and localj stacked onto the pair of locali+1,localj-1
						if (localj!=1&&!lfce[locali]&&!lfce[localj]) {
							//[WERF93]
         					wl->f(locali,localj)= SUM(wl->f(locali,localj),PROD(v->f(locali+1,localj-1), multi_eparam(10, ct, data), multi_eparam(6, ct, data), multi_eparam(6, ct, data),
 // orifinal        			data->tstkm[ct->numseq[localj-1]][ct->numseq[locali+1]][ct->numseq[localj]][ct->numseq[locali]],
								ptr_to_tstkm(localj - 1, locali + 1, localj, locali, ct, data),
								ptr_to_penalty(localj-1,locali+1,ct,data)));
							LOG_DEBUG4("[LX1N9X]\twl->f(locali,localj): " << wl->f(locali,localj));

							if ((mod[locali+1]||mod[localj-1])) if((localj-2>0)&&!(fce->f(locali+1,localj-1)&SINGLE)) {
								if(ptr_to_can_pair(locali+2, localj-2,ct,data)&&notgu(locali+1,localj-1,ct)) {
									wl->f(locali,localj)= SUM(wl->f(locali,localj),PROD(v->f(locali+2,localj-2), multi_eparam(10, ct, data), multi_eparam(6, ct, data), multi_eparam(6, ct, data),
         //	original				data->tstkm[ct->numseq[localj-1]][ct->numseq[locali+1]][ct->numseq[localj]][ct->numseq[locali]],
										ptr_to_tstkm(localj - 1,locali + 1, localj, locali, ct, data),
										ptr_to_penalty(localj-1,locali+1,ct,data), ptr_to_erg1(locali+1,localj-1,locali+2,localj-2,ct,data)));
									LOG_DEBUG4("[LX1N9X.2]\twl->f(locali,localj): " << wl->f(locali,localj));

								}
							}
						}
					}

					#endif  //SIMPLEMBLOOP

					wlc->f(locali, localj) = wl->f(locali,localj);
					
					if (locali!=number&&!lfce[locali]) {
					   //add a nuc to an existing loop:
         				wl->f(locali,localj)= SUM(wl->f(locali,localj),PROD( wl->f(locali+1,localj), multi_eparam(6, ct, data) SCALING));
						LOG_DEBUG4("[YJDEPA]\twl->f(locali,localj): " << wl->f(locali,localj));
            			//this is for when locali represents the center of an intermolecular linker:
						// else e[4] = w->f(locali+1,localj) + data->eparam[6] + infinity;
					}

					w->f(locali,localj) = wl->f(locali,localj);
					if (localj!=number+1&&!lfce[localj]) {
             			//if (!(fce->f(localj,localj)&INTER)) {
               			//add a nuc to an existing loop:
               			w->f(locali,localj)= SUM(w->f(locali,localj),PROD( w->f(locali,localj-1), multi_eparam(6, ct, data) SCALING));
					   //}
					   //else e[5] = w->f(locali,localj-1) + data->eparam[6] + infinity;

					}
					wl_T.f(locali,localj) = wl->f(locali,localj);
				}
			}//end of else "not a forced pair"

		//sub3:
			SET_DEBUG_LEVEL(INFO);
		}

		#ifdef pfdebugmode
		if (twoscaling>PFMAX||twoscaling<PFMIN) {
			//look for positions that will underflow

			ofstream *ufout = new ofstream();
			ufout->open(pfdebugpath,ios::app);
			*ufout<<"locali= "<<locali<<" localj= "<<localj<<" twoscaling= "<<twoscaling<<"\n";
			ufout->close();
			delete ufout;

		}
		#endif//pfdebugmode
	}//end for locali

//*
#ifndef PF_LOG_CALC
	PFPRECISION MAX_PF(PFMAX);
	PFPRECISION MIN_PF(PFMIN);

	int countRecaleUP=0,countRecaleDOWN=0;
	for (i=((h<=(number-1))?1:(2*number-h));i<=((h<=(number-1))?(number-h):number);i++){
		j=i+d;
		bool too_big = (v->f(i,j)>MAX_PF) || (w->f(i,j)>MAX_PF) || (wl->f(i,j)>MAX_PF)
			|| (wcoax->f(i,j)>MAX_PF) || (wmb->f(i,j)>MAX_PF) || (wmbl->f(i,j)>MAX_PF);

		bool too_small = (v->f(i,j)<MIN_PF&&v->f(i,j)>0) || (w->f(i,j)<MIN_PF&&w->f(i,j)>0)
			|| (wl->f(i,j)<MIN_PF&&wl->f(i,j)>0) || (wcoax->f(i,j)<MIN_PF&&wcoax->f(i,j)>0)
			|| (wmb->f(i,j)<MIN_PF&&wmb->f(i,j)>0) || (wmbl->f(i,j)<MIN_PF&&wmbl->f(i,j)>0);

		if(too_big || too_small){
			if (too_big) countRecaleUP++; if (too_small) countRecaleDOWN++; // too_big and too_small are not necessarily mutually exclusive.
			PFPRECISION factor = too_big ? SCALEDOWN : SCALEUP;
			// cout << (too_big?"Too BIG":"Too SMALL") << "\ti=" << i << ", j=" << j << "\tv=" << v->f(i,j) << "\tw=" << w->f(i,j) << "\twl="  << wl->f(i,j) << "\twcoax="  << wcoax->f(i,j) << "\twmb=" << wmb->f(i,j) <<"\twmbl=" <<  wmbl->f(i,j) << endl;
			rescale(h,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,factor);
			twoscaling *= factor*factor;
			copyDPArray(*wl,wl_T);
			copyDPArray(*wmbl,wmbl_T);
		}
	}
	if (countRecaleUP!=0&&countRecaleDOWN!=0) ct->cwarn() << "ERROR in Partition Function: rescale too_big AND too_small!!" << endl;
	if (countRecaleUP>0||countRecaleDOWN>0)   ct->cwarn() << "WARNING in Partition Function: Multiple calls to rescale in for-loop." << endl;

#endif //*/
		//Compute w5[i], the energy of the best folding from 1->i, and

      		//w3[i], the energy of the best folding from i-->numofbases

	if (h<=(number-1)) {
		i = 1;
		j=i+d;

		// if (j == 622)
		// 	SET_DEBUG_LEVEL(DEBUG4)

		if (j<=minloop+1) {
			if (lfce[j]){
				w5[j]= (PFPRECISION) ZERO;				
			}
			else  {
				w5[j] = PROD(w5[j-1] SCALING);
			}
			LOG_DEBUG3("[332V5Q]\tw5[j]: " << w5[j]);
		}

		else {
      		if (lfce[j]){
				rarray = (PFPRECISION) ZERO;
			}

			else {
				rarray = PROD(w5[j-1] SCALING);
			}
			LOG_DEBUG3("[332V5Q]\tw5[j]: " << rarray);

      		for (k=0;k<=(j-4);k++) {
				#ifdef SIMPLEMBLOOP
				//5' and 3' dangling ends
				if (k>0&&j!=N) {
					rarray=SUM(rarray,PROD(w5[k],erg4(j,k+1,j+1,1,ct,data,lfce[j+1]),
						erg4(j,k+1,k,2,ct,data,lfce[k]),
						v->f(k+1,j),penalty(j,k+1,ct,data)));
					LOG_DEBUG4("[XXXXXX]\tw5[j]: " << rarray);

					if ((mod[k+1]||mod[j])) if(inc[ct->numseq[k+2]][ct->numseq[j-1]]&&notgu(k+1,j,ct)&&!(fce->f(k+1,j)&SINGLE)) {
						rarray=SUM(rarray,PROD(w5[k],erg4(j,k+1,j+1,1,ct,data,lfce[j+1]),
							erg4(j,k+1,k,2,ct,data,lfce[k]),v->f(k+2,j-1),
							penalty(j-1,k+1,ct,data),erg1(k+1,j,k+2,j-1,ct,data)));
						LOG_DEBUG4("[XXXXXX]\tw5[j]: " << rarray);
					}
				}
				//5' dangling end
				if (k>0) {//j==N
					rarray=SUM(rarray,PROD(w5[k],erg4(j,k+1,k,2,ct,data,lfce[k]),
						v->f(k+1,j),penalty(j,k+1,ct,data)));
					LOG_DEBUG4("[XXXXXX]\tw5[j]: " << rarray);

					if ((mod[k+1]||mod[j])) if(inc[ct->numseq[k+2]][ct->numseq[j-1]]&&notgu(k+1,j,ct)&&!(fce->f(k+1,j)&SINGLE)) {
						rarray=SUM(rarray,PROD(w5[k],erg4(j,k+1,k,2,ct,data,lfce[k]),v->f(k+2,j-1),
							penalty(j-1,k+1,ct,data),erg1(k+1,j,k+2,j-1,ct,data)));
						LOG_DEBUG4("[XXXXXX]\tw5[j]: " << rarray);
					}
				}
				//3' dangling end
				if (j!=N) {//k==0
					rarray=SUM(rarray,PROD(w5[k],erg4(j,k+1,j+1,1,ct,data,lfce[j+1]),
						v->f(k+1,j),penalty(j,k+1,ct,data)));
					LOG_DEBUG4("[XXXXXX]\tw5[j]: " << rarray);

					if ((mod[k+1]||mod[j])) if(inc[ct->numseq[k+2]][ct->numseq[j-1]]&&notgu(k+1,j,ct)&&!(fce->f(k+1,j)&SINGLE)) {
						rarray=SUM(rarray,PROD(w5[k],erg4(j,k+1,j+1,1,ct,data,lfce[j+1]),
							v->f(k+2,j-1), penalty(j-1,k+1,ct,data), erg1(k+1,j,k+2,j-1,ct,data)));
						LOG_DEBUG4("[XXXXXX]\tw5[j]: " << rarray);
					}
				}
				//no dangling ends
				else {//k==0 and j==N
					rarray=SUM(rarray,PROD(w5[k], v->f(k+1,j), penalty(j,k+1,ct,data)));
					LOG_DEBUG4("[XXXXXX]\tw5[j]: " << rarray);

					if ((mod[k+1]||mod[j])) if(inc[ct->numseq[k+2]][ct->numseq[j-1]]&&notgu(k+1,j,ct)&&!(fce->f(k+1,j)&SINGLE)) {
						rarray=SUM(rarray,PROD(w5[k], v->f(k+2,j-1),
							penalty(j-1,k+1,ct,data), erg1(k+1,j,k+2,j-1,ct,data)));
						LOG_DEBUG4("[XXXXXX]\tw5[j]: " << rarray);
					}
				}
				#else //!SIMPLEMBLOOP

      			rarray=SUM(rarray,PROD(w5[k], v->f(k+1,j), ptr_to_penalty(j,k+1,ct,data)));
				LOG_DEBUG4("[XO5121]\tw5[j]: " << rarray);

				if ((mod[k+1]||mod[j])) if(ptr_to_can_pair(k+2, j-1,ct,data)&&notgu(k+1,j,ct)&&!(fce->f(k+1,j)&SINGLE)) {
					rarray=SUM(rarray,PROD(w5[k], v->f(k+2,j-1), ptr_to_penalty(j,k+1,ct,data),
						ptr_to_erg1(k+1,j,k+2,j-1,ct,data)));
					LOG_DEBUG4("[XO5121.2]\tw5[j]: " << rarray);
				}

				rarray=SUM(rarray,PROD(w5[k], ptr_to_erg4(j,k+2,k+1,2,ct,data,lfce[k+1]), v->f(k+2,j), ptr_to_penalty(j,k+2,ct,data)));
				LOG_DEBUG4("[0A1GDP]\tw5[j]: " << rarray);

				if((mod[k+2]||mod[j])) if(ptr_to_can_pair(k+3, j-1,ct,data)&&notgu(k+2,j,ct)
					&&!(fce->f(k+2,j)&SINGLE)) {
					rarray=SUM(rarray,PROD(w5[k], ptr_to_erg4(j,k+2,k+1,2,ct,data,lfce[k+1]), v->f(k+3,j-1),
					   ptr_to_penalty(j,k+2,ct,data), ptr_to_erg1(k+2,j,k+3,j-1,ct,data)));
					LOG_DEBUG4("[0A1GDP.2]\tw5[j]: " << rarray);

				}

         		rarray=SUM(rarray,PROD(w5[k], ptr_to_erg4(j-1,k+1,j,1,ct,data,lfce[j]), v->f(k+1,j-1), ptr_to_penalty(j-1,k+1,ct,data)));
				LOG_DEBUG4("[X0HOHA]\tw5[j]: " << rarray);

				if ((mod[k+1]||mod[j-1])) if(ptr_to_can_pair(k+2, j-2,ct,data)&&notgu(k+1,j-1,ct)&&!(fce->f(k+1,j-1)&SINGLE)) {
					rarray=SUM(rarray,PROD(w5[k], ptr_to_erg4(j-1,k+1,j,1,ct,data,lfce[j]), v->f(k+2,j-2),
						ptr_to_penalty(j-1,k+1,ct,data), ptr_to_erg1(k+1,j-1,k+2,j-2,ct,data)));
					LOG_DEBUG4("[X0HOHA.2]\tw5[j]: " << rarray);
				}

				rarray=SUM(rarray,PROD(w5[k], 
// original						data->tstack[ct->numseq[j-1]][ct->numseq[k+2]][ct->numseq[j]][ct->numseq[k+1]],
								ptr_to_tstack(j-1, k+2, j, k+1, ct, data),
								pfchecknp(lfce[j],lfce[k+1]), v->f(k+2,j-1), 
								ptr_to_penalty(j-1,k+2,ct,data)));
				LOG_DEBUG4("[NXM9K1]\tw5[j]: " << rarray);

				if ((mod[k+2]||mod[j-1])) if(ptr_to_can_pair(k+3, j-2,ct,data)&&notgu(k+2,j-1,ct)&&!(fce->f(k+2,j-1)&SINGLE)) {
					rarray=SUM(rarray,PROD(w5[k], 
// original						data->tstack[ct->numseq[j-1]][ct->numseq[k+2]][ct->numseq[j]][ct->numseq[k+1]],
								ptr_to_tstack(j - 1,k + 2, j, k + 1, ct, data),
								pfchecknp(lfce[j],lfce[k+1]), v->f(k+3,j-2), 
								ptr_to_penalty(j-1,k+2,ct,data), ptr_to_erg1(k+2,j-1,k+3,j-2,ct,data)));
					LOG_DEBUG4("[NXM9K1.2]\tw5[j]: " << rarray);
				}

				rarray=SUM(rarray,PROD(w5[k], wca[k+1][j]));
				LOG_DEBUG4("[D753PB]\tw5[j]: " << rarray);

				#endif  //SIMPLEMBLOOP
			}

			w5[j] = rarray;
			LOG_DEBUG3("[0NQO6X]\tw5[j]: " << rarray);
			SET_DEBUG_LEVEL(INFO);
/*
			//check to see if w5 is about to go out of bounds:
			if (w5[j]>PFMAX) {
				rescale(h,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEDOWN);
				twoscaling = twoscaling*SCALEDOWN*SCALEDOWN;
				copyDPArray(*wl,wl_T);
				copyDPArray(*wmbl,wmbl_T);
			}
			else if (w5[j]<PFMIN&&w5[j]>0) {
				rescale(h,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEUP);
				twoscaling=twoscaling*SCALEUP*SCALEUP;
				copyDPArray(*wl,wl_T);
				copyDPArray(*wmbl,wmbl_T);
			} //*/
		}
	}//end if (h<=number-1)

	if (h==number-1) {
		//w3[0] = 0;
		//w3[number+1] = 0;
		for (ii=(number);ii>=(number-minloop);ii--) {    //number+1 ... number-minloop
      		if (lfce[ii]) w3[ii] = (PFPRECISION) ZERO;
			else w3[ii]=PROD(w3[ii+1] SCALING);
			LOG_DEBUG4("[XXXXXX]\tw3[j]: " << w3[ii]);
		}
		//w3[i]=0;
   		for (ii=((number)-minloop-1);ii>=1;ii--) {
      		if (lfce[ii]) rarray = (PFPRECISION) ZERO;

   			else rarray = PROD(w3[ii+1] SCALING);
			LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);

			for (k=((number)+1);k>=(ii+4);k--) {
				#ifdef SIMPLEMBLOOP
				if (ii>1&&k!=N+1) {
					//5' and 3' dangling ends
					rarray=SUM(rarray,PROD(v->f(ii,k-1), erg4(k-1,ii,k,1,ct,data,lfce[k]),
						erg4(k-1,ii,ii-1,2,ct,data,lfce[ii-1]),
						penalty(k-1,ii,ct,data), w3[k]));
					LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);

					if((mod[ii]||mod[k-1]))if(inc[ct->numseq[ii+1]][ct->numseq[k-2]]&&notgu(ii,k-1,ct)&&!(fce->f(ii,k-1)&SINGLE)) {
						rarray=SUM(rarray,PROD(v->f(ii+1,k-2), erg4(k-1,ii,k,1,ct,data,lfce[k]),
						erg4(k-1,ii,ii-1,2,ct,data,lfce[ii-1]),
						penalty(k-1,ii,ct,data), w3[k], erg1(ii,k-1,ii+1,k-2,ct,data)));
						LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);

					}
				}
				else if (k!=N+1) {//i==1
					//3' dangling end
					rarray=SUM(rarray,PROD(v->f(ii,k-1), erg4(k-1,ii,k,1,ct,data,lfce[k]),
						penalty(k-1,ii,ct,data), w3[k]));
					LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);

					if((mod[ii]||mod[k-1]))if(inc[ct->numseq[ii+1]][ct->numseq[k-2]]&&notgu(ii,k-1,ct)&&!(fce->f(ii,k-1)&SINGLE)) {
						rarray=SUM(rarray,PROD(v->f(ii+1,k-2), erg4(k-1,ii,k,1,ct,data,lfce[k]),
						penalty(k-1,ii,ct,data), w3[k], erg1(ii,k-1,ii+1,k-2,ct,data)));
						LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);

					}
				}
				else if (ii>1) {//k==N+1
					//5' dangling end
					rarray=SUM(rarray,PROD(v->f(ii,k-1),
						erg4(k-1,ii,ii-1,2,ct,data,lfce[ii-1]),
						penalty(k-1,ii,ct,data), w3[k]));
					LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);

					if((mod[ii]||mod[k-1]))if(inc[ct->numseq[ii+1]][ct->numseq[k-2]]&&notgu(ii,k-1,ct)&&!(fce->f(ii,k-1)&SINGLE)) {
						rarray=SUM(rarray,PROD(v->f(ii+1,k-2),
						erg4(k-1,ii,ii-1,2,ct,data,lfce[ii-1]),
						penalty(k-1,ii,ct,data), w3[k], erg1(ii,k-1,ii+1,k-2,ct,data)));
						LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);

					}
				}
				else {//ii==1&&k==N+1
					//No dangling ends
					rarray=SUM(rarray,PROD(v->f(ii,k-1),
						penalty(k-1,ii,ct,data), w3[k]));
					LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);

					if((mod[ii]||mod[k-1]))if(inc[ct->numseq[ii+1]][ct->numseq[k-2]]&&notgu(ii,k-1,ct)&&!(fce->f(ii,k-1)&SINGLE)) {
						rarray=SUM(rarray,PROD(v->f(ii+1,k-2), 
						penalty(k-1,ii,ct,data), w3[k], erg1(ii,k-1,ii+1,k-2,ct,data)));
						LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);

					}
				}
				#else //!SIMPLEMBLOOP

      			rarray=SUM(rarray,PROD(v->f(ii,k-1), w3[k], ptr_to_penalty(k-1,ii,ct,data)));
				LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);

				if((mod[ii]||mod[k-1])) if(ptr_to_can_pair(ii+1, k-2,ct,data)&&notgu(ii,k-1,ct)&&!(fce->f(ii,k-1)&SINGLE)) {
					rarray=SUM(rarray,PROD(v->f(ii+1,k-2), w3[k], ptr_to_penalty(k-1,ii,ct,data), ptr_to_erg1(ii,k-1,ii+1,k-2,ct,data)));
					LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);
				}

				rarray=SUM(rarray,PROD(v->f(ii+1,k-1), ptr_to_erg4(k-1,ii+1,ii,2,ct,data,lfce[ii]), ptr_to_penalty(k-1,ii+1,ct,data), w3[k]));
				LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);

				if((mod[ii+1]||mod[k-1])) if(ptr_to_can_pair(ii+2, k-2,ct,data)&&notgu(ii+1,k-1,ct)&&!(fce->f(ii+1,k-1)&SINGLE)) {
					rarray=SUM(rarray,PROD(v->f(ii+2,k-2), ptr_to_erg4(k-1,ii+1,ii,2,ct,data,lfce[ii]),
						ptr_to_penalty(k-1,ii+1,ct,data), w3[k], ptr_to_erg1(ii+1,k-1,ii+2,k-2,ct,data)));
					LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);
				}

				rarray=SUM(rarray,PROD(v->f(ii,k-2), ptr_to_erg4(k-2,ii,k-1,1,ct,data,lfce[k-1]), ptr_to_penalty(k-2,ii,ct,data), w3[k]));
				LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);

				if((mod[ii]||mod[k-2]))if(ptr_to_can_pair(ii+1, k-3,ct,data)&&notgu(ii,k-2,ct)&&!(fce->f(ii,k-2)&SINGLE)) {
					rarray=SUM(rarray,PROD(v->f(ii+1,k-3), ptr_to_erg4(k-2,ii,k-1,1,ct,data,lfce[k-1]),
						ptr_to_penalty(k-2,ii,ct,data), w3[k], ptr_to_erg1(ii,k-2,ii+1,k-3,ct,data)));
					LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);
				}

				if (!lfce[ii]&&!lfce[k-1]) {
					rarray=SUM(rarray,PROD(v->f(ii+1,k-2), 
// original				data->tstack[ct->numseq[k-2]][ct->numseq[ii+1]][ct->numseq[k-1]][ct->numseq[ii]], 
						ptr_to_tstack(k - 2, ii + 1, k - 1, ii, ct, data),
						w3[k], ptr_to_penalty(k-2,ii+1,ct,data)));
					LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);

					if((mod[ii+1]||mod[k-2]))if(ptr_to_can_pair(ii+2,k-3,ct,data)&&notgu(ii+1,k-2,ct)&&!(fce->f(ii+1,k-2)&SINGLE)) {
						rarray=SUM(rarray,PROD(v->f(ii+2,k-3), 
// original					data->tstack[ct->numseq[k-2]][ct->numseq[ii+1]][ct->numseq[k-1]][ct->numseq[ii]], 
							ptr_to_tstack(k - 2, ii + 1, k - 1, ii, ct, data),
							pfchecknp(lfce[k-1],lfce[ii]), w3[k],
							ptr_to_penalty(k-2,ii+1,ct,data), ptr_to_erg1(ii+1,k-2,ii+2,k-3,ct,data)));
						LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);
					}
				}

				//also consider coaxial stacking:
				#ifndef disablecoax //a flag to disable coaxial stacking
				if(!disablecoax){
					for (ip=k+minloop+1;ip<=number+1;ip++) {
						//first consider flush stacking:
						rarray=SUM(rarray,PROD(v->f(ii,k-1), v->f(k,ip-1), w3[ip],
							ptr_to_penalty(ii,k-1,ct,data), ptr_to_penalty(k,ip-1,ct,data),
							ptr_to_ergcoaxflushbases(ii,k-1,k,ip-1,ct,data)));
						LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);

						if(mod[ii]||mod[k-1]||mod[k]||mod[ip-1]) {
							if ((mod[ii]||mod[k-1])&&(mod[k]||mod[ip-1])&&ptr_to_can_pair(ii+1, k-2,ct,data)
								&&ptr_to_can_pair(k+1, ip-2, ct, data)&&notgu(ii,k-1,ct)&&notgu(k,ip-1,ct)
								&&!(fce->f(ii,k-1)&SINGLE)&&!(fce->f(k,ip-1)&SINGLE)) {
								rarray=SUM(rarray,PROD(v->f(ii+1,k-2), v->f(k+1,ip-2), w3[ip],
									ptr_to_penalty(ii,k-1,ct,data), ptr_to_penalty(k,ip-1,ct,data),
									ptr_to_ergcoaxflushbases(ii,k-1,k,ip-1,ct,data),
									ptr_to_erg1(ii,k-1,ii+1,k-2,ct,data), ptr_to_erg1(k,ip-1,k+1,ip-2,ct,data)));
								LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);
							}
							if((mod[ii]||mod[k-1])&&ptr_to_can_pair(ii+1, k-2,ct,data)&&notgu(ii,k-1,ct)&&!(fce->f(ii,k-1)&SINGLE)) {
								rarray=SUM(rarray,PROD(v->f(ii+1,k-2), v->f(k,ip-1), w3[ip],
									ptr_to_penalty(ii,k-1,ct,data), ptr_to_penalty(k,ip-1,ct,data),
									ptr_to_ergcoaxflushbases(ii,k-1,k,ip-1,ct,data),
									ptr_to_erg1(ii,k-1,ii+1,k-2,ct,data)));
								LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);

							}

							if((mod[k]||mod[ip-1])&&ptr_to_can_pair(k+1, ip-2, ct, data)&&notgu(k,ip-1,ct)&&!(fce->f(k,ip-1)&SINGLE)) {
								rarray=SUM(rarray,PROD(v->f(ii,k-1), v->f(k+1,ip-2), w3[ip],
									ptr_to_penalty(ii,k-1,ct,data), ptr_to_penalty(k,ip-1,ct,data), 
									ptr_to_ergcoaxflushbases(ii,k-1,k,ip-1,ct,data),
									ptr_to_erg1(k,ip-1,k+1,ip-2,ct,data)));
								LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);
							}
						}

						//now consider an intervening mismatch:
						if (ii>0&&!lfce[ii]&&!lfce[k-1]) {
							rarray=SUM(rarray,PROD(v->f(ii+1,k-2), v->f(k,ip-1), w3[ip], 
								ptr_to_penalty(ii+1,k-2,ct,data), ptr_to_penalty(k,ip-1,ct,data),
								ptr_to_ergcoaxinterbases1(ii+1,k-2,k,ip-1,ct,data)));
							LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);

							if(mod[ii+1]||mod[k-2]||mod[k]||mod[ip-1]){
								if((mod[ii+1]||mod[k-2])&&(mod[k]||mod[ip-1])&&ptr_to_can_pair(ii+2,k-3,ct,data)
									&&ptr_to_can_pair(k+1, ip-2, ct, data)&&notgu(ii+1,k-2,ct)&&notgu(k,ip-1,ct)
									&&!(fce->f(k,ip-1)&SINGLE)&&!(fce->f(ii+1,k-2)&SINGLE)){
									rarray=SUM(rarray,PROD(v->f(ii+2,k-3), v->f(k+1,ip-2), w3[ip],
										ptr_to_penalty(ii+1,k-2,ct,data), ptr_to_penalty(k,ip-1,ct,data),
										ptr_to_ergcoaxinterbases1(ii+1,k-2,k,ip-1,ct,data),
										ptr_to_erg1(ii+1,k-2,ii+2,k-3,ct,data), ptr_to_erg1(k,ip-1,k+1,ip-2,ct,data)));
									LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);

								}

								if((mod[ii+1]||mod[k-2])&&ptr_to_can_pair(ii+2,k-3,ct,data)&&notgu(ii+1,k-2,ct)&&!(fce->f(ii+1,k-2)&SINGLE)) {
									rarray=SUM(rarray,PROD(v->f(ii+2,k-3), v->f(k,ip-1), w3[ip],
									ptr_to_penalty(ii+1,k-2,ct,data), ptr_to_penalty(k,ip-1,ct,data),
									ptr_to_ergcoaxinterbases1(ii+1,k-2,k,ip-1,ct,data),
									ptr_to_erg1(ii+1,k-2,ii+2,k-3,ct,data)));
									LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);

								}
								if((mod[k]||mod[ip-1])&&ptr_to_can_pair(k+1, ip-2, ct, data)&&notgu(k,ip-1,ct)&&!(fce->f(k,ip-1)&SINGLE)) {
									rarray=SUM(rarray,PROD(v->f(ii+1,k-2), v->f(k+1,ip-2), w3[ip],
										ptr_to_penalty(ii+1,k-2,ct,data), ptr_to_penalty(k,ip-1,ct,data),
										ptr_to_ergcoaxinterbases1(ii+1,k-2,k,ip-1,ct,data),
										ptr_to_erg1(k,ip-1,k+1,ip-2,ct,data)));
									LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);
								}
							}
						}
						if (!lfce[k-1]&&!lfce[ip-1]) {
							rarray=SUM(rarray,PROD(v->f(ii,k-2), v->f(k,ip-2), w3[ip],
								ptr_to_penalty(ii,k-2,ct,data), ptr_to_penalty(k,ip-2,ct,data),
								ptr_to_ergcoaxinterbases2(ii,k-2,k,ip-2,ct,data)));
							LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);

							if (mod[ii]||mod[k-2]||mod[k]||mod[ip-2]) {
								if ((mod[ii]||mod[k-2])&&(mod[k]||mod[ip-2])&&ptr_to_can_pair(ii+1, k-3,ct,data)
									&&ptr_to_can_pair(k+1, ip-3, ct, data)&&notgu(ii,k-2,ct)&&notgu(k,ip-2,ct)
									&&!(fce->f(ii,k-2)&SINGLE)&&!(fce->f(k,ip-2)&SINGLE)) {
									rarray=SUM(rarray,PROD(v->f(ii+1,k-3), v->f(k+1,ip-3), w3[ip],
										ptr_to_penalty(ii,k-2,ct,data), ptr_to_penalty(k,ip-2,ct,data),
										ptr_to_ergcoaxinterbases2(ii,k-2,k,ip-2,ct,data),
										ptr_to_erg1(ii,k-2,ii+1,k-3,ct,data), ptr_to_erg1(k,ip-2,k+1,ip-3,ct,data)));
									LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);
								}

								if ((mod[ii]||mod[k-2])&&ptr_to_can_pair(ii+1, k-3,ct,data)&&notgu(ii,k-2,ct)&&!(fce->f(ii,k-2)&SINGLE)) {
									rarray=SUM(rarray,PROD(v->f(ii+1,k-3), v->f(k,ip-2), w3[ip],
										ptr_to_penalty(ii,k-2,ct,data), ptr_to_penalty(k,ip-2,ct,data),
										ptr_to_ergcoaxinterbases2(ii,k-2,k,ip-2,ct,data),
										ptr_to_erg1(ii,k-2,ii+1,k-3,ct,data)));
									LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);
								}

								if ((mod[k]||mod[ip-2])&&ptr_to_can_pair(k+1, ip-3, ct, data)&&notgu(k,ip-2,ct)&&!(fce->f(k,ip-2)&SINGLE)) {
									rarray=SUM(rarray,PROD(v->f(ii,k-2), v->f(k+1,ip-3), w3[ip],
										ptr_to_penalty(ii,k-2,ct,data), ptr_to_penalty(k,ip-2,ct,data),
										ptr_to_ergcoaxinterbases2(ii,k-2,k,ip-2,ct,data),
										ptr_to_erg1(k,ip-2,k+1,ip-3,ct,data)));
									LOG_DEBUG4("[XXXXXX]\tw3[j]: " << rarray);
								}
							}
						}
					}
				}
				#endif //ifndef disablecoax
				#endif //SIMPLEMBLOOP

			}
			w3[ii] = rarray;
			LOG_DEBUG4("[XXXXXX]\tw3[j] (final): " << rarray);
//*
#ifndef PF_LOG_CALC
			//check to see if w5 is about to go out of bounds:
			if (w3[ii]>PFMAX) {
				rescale(h,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,SCALEDOWN);
				twoscaling=twoscaling*SCALEDOWN*SCALEDOWN;
			}
			else if (w3[ii]<PFMIN&&w3[ii]>0) {
				rescale(h,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,SCALEUP);
				twoscaling = twoscaling*SCALEUP*SCALEUP;
			} 
#endif //*/
		}
	}

	
}

for(ii=0;ii<=number;ii++) {
	delete[] wca[ii];
}
delete[] wca;

//clean up the inc array, which tracked what can pair
for (int letters=0;letters<data->alphabet.size();++letters) delete[] inc[letters];
delete[] inc;

#ifdef timer
timeout << time(NULL)<<"\n";
timeout << time(NULL) - seconds;
timeout.close();
#endif

//////////////////////////
//output V, W, WMB, and W2V:
#if defined (pfdebugmode)
	ofstream foo;
	foo.open("arrays.out");
	foo << "i" << "\t"<<"j"<<"\t"<<"v->f(i,j)"<<"\t"<<"w->f(i,j)"<<"\t"<<"wmb->f(i,j)\twmbl->f(i,j)\twcoax->f(i,j)"<<"\t"<<"wl->f(i,j)"<<"\t"<<"v->f(j,i+number)"<<"\t"<<"w->f(j,i+number)"<<"\t"<<"wmb->f(j,i+number)"<<"\t"<<"wl->f(j,i+number)"<<"\t"<<"wmbl->f(j,i+numer)\twcoax->f(j,i+number)"<<"\n";
	for (j=1;j<=number;j++) {
		for (i=1;i<=j;i++) {
			foo << i << "\t"<<j<<"\t"<<v->f(i,j)<<"\t"<<w->f(i,j)<<"\t"<<wmb->f(i,j)<<"\t"<<wmbl->f(i,j)<<"\t"<<wcoax->f(i,j)<<"\t"<<wl->f(i,j)<<"\t"<<v->f(j,i+number)<<"\t"<<w->f(j,i+number)<<"\t"<<wmb->f(j,i+number)<<"\t"<<wl->f(j,i+number)<<"\t"<<wmbl->f(j,i+number)<<"\t"<<wcoax->f(j,i+number)<<"\n";

		}
	}

	foo <<"\n\n\n";
	foo << "i" << "\t" << "w5[i]" << "\t" << "w3[i]" << "\n";
	for (i=0;i<=number;i++) foo << i << "\t" << w5[i] << "\t" << w3[i] << "\n";

	foo.close();

#endif
//cout << "Q: \t" << w3[1] << endl;
}

//This function cacluates a partition function for the sequence in CT
//If quickQ == true, return the partition function value in Q
//	otherwise, save the partial partition functions in the datafile named save
//If updates on progress are unwanted, set update=NULL
void pfunction(structure* ct,pfdatatable* data, ProgressHandler* update, char* save, bool quickQ, PFPRECISION *Q)
{
int i,j;
bool *lfce,*mod;//[maxbases+1][maxbases+1];
PFPRECISION *w5,*w3;
int number;

#ifdef equiout
	ofstream kout;
	kout.open("k.out");
	kout << "sequence length = "<<ct->GetSequenceLength()<<"\n";
	for (i=1;i<=ct->GetSequenceLength();i++) {
		kout << tobase(ct->numseq[i]);
		if (i%20==0) kout << "\n";
	}
	kout << "\n";
	kout.close();
#endif

//array *v,*w;
//inc is an array that saves time by showing which bases can pair before
//	erg is called

//number is the number of bases being folded
//v[i][j] is the best energy for subsequence i to j when i and j are paired
//	for i<j<n, it is the interior framgment between nucleotides i and j
//	for i<n<j, it is the exterior structure fragmnet that is 5' to j-n and 3' to i
//w[i][j] is the best energy for subsequence i to j

number = (ct->GetSequenceLength());//place the number of bases in an integer

//scaling is a per nucleotide scale factor for which W and V are divided
//This is necessary to keep the doubles from overflowing:

//scaling = 0.2;//this factor assumes about 1 kcal/mol/base

//allocate space for the v and w arrays:
DynProgArray<PFPRECISION> w(number);
DynProgArray<PFPRECISION> v(number);
DynProgArray<PFPRECISION> wmb(number);
DynProgArray<PFPRECISION> wl(number);
DynProgArray<PFPRECISION> wlc(number);
DynProgArray<PFPRECISION> wmbl(number);
DynProgArray<PFPRECISION> wcoax(number);
forceclass fce(number);

if (ct->intermolecular) {
	//take advantage of templating to prevent intramolecular base pairs

	ct->allocatetem();
	for (i=1;i<ct->inter[0];i++) {
		for (j=i+1;j<=ct->inter[2];j++) {
			ct->tem[j][i]=false;

		}
	}
	for (i=ct->inter[2]+1;i<ct->GetSequenceLength();i++) {
		for (j=i+1;j<=ct->GetSequenceLength();j++) {
			ct->tem[j][i]=false;

		}
	}
}

//This code converts the SHAPE array of data to equilibrium constants.  This is
//needed for the partition function.  NOTE, however, there is no going back so a structure
//that has been used for partition functions cannot then be used to predict a structure
//by free energy minimization. This is a compromise for efficiency, but clearly something
//that could cause a problem for programmers.

if (ct->shaped) {
	for (i=1;i<=2*ct->GetSequenceLength();i++) ct->SHAPE[i]=(double) boltzman(ct->SHAPE[i], data->pftemp);

}

//add a second array for intermolecular folding:


lfce = new bool [2*number+1];
mod = new bool [2*number+1];

for (i=0;i<=2*number;i++) {
	lfce[i] = false;
	mod[i] = false;
}

//Register modified nucleotides in the mod array for fast access
for (i=0;i<ct->GetNumberofModified();i++) {
	if (ct->GetModified(i)!=1&&ct->GetModified(i)!=ct->GetSequenceLength()) {
		mod[ct->GetModified(i)]=true;
		mod[ct->GetModified(i)+ct->GetSequenceLength()]=true;
	}
}

w5 = new PFPRECISION [number+1];
w3 = new PFPRECISION [number+2];
//wca = new PFPRECISION *[number+1];

//The next section handles the case where base pairs are not
				//not allowed to form between nucs more distant
				//than ct->GetPairingDistanceLimit()
if (ct->DistanceLimited()) {
	if (!ct->templated) ct->allocatetem();

	for (j=minloop+2;j<=ct->GetSequenceLength();j++) {
		for (i=1;i<j;i++) {
			if (j-i>=ct->GetPairingDistanceLimit()) ct->tem[j][i]=false;
		}
	}
}

calculatepfunction(ct,data,update,save,quickQ,Q,&w,&v,&wmb,&wl,&wlc,&wmbl,&wcoax,&fce,w5,w3,mod,lfce);

if (save!=0) {
	writepfsave(save,ct,w5,w3,&v,&w,&wmb,&wl,&wlc,&wmbl,&wcoax,&fce,mod,lfce,data);
}

if (quickQ) *Q = w5[ct->GetSequenceLength()];

delete[] lfce;
delete[] mod;

delete[] w5;
delete[] w3;

return;
}

template<typename TWriter>
struct DataWriter {
	protected:
	ostream *out = nullptr;
	bool owner = false;
	DataWriter(ostream *os_out, bool owns=true) : out(os_out), owner(owns) { }
	virtual ~DataWriter() { if (owner) delete out; }
	template<typename TData>
	TWriter& write(const TData& data);
};

template<typename TWriter, typename TData>
TWriter& operator<<(TWriter& writer, const TData& data) {
	return writer.write(data);
};

struct BinaryWriter : DataWriter<BinaryWriter> {
	BinaryWriter(ostream& out) : DataWriter<BinaryWriter>(&out, false) { }
	BinaryWriter(const string& file, std::ios_base::openmode mode = std::ios_base::trunc) 
		: DataWriter<BinaryWriter>(new ofstream(file.c_str(), std::ios_base::binary | mode)) { }
	
	template<typename TData>
	BinaryWriter& write(const TData&);
	
	BinaryWriter& write(const char*const data, const size_t size) {
		write(size);
		out->write(data, size);
		return *this;
	}
};

template<typename TData>
BinaryWriter& BinaryWriter::write(const TData& data) {
	out->write(reinterpret_cast<const char*>(data), sizeof(TData));
	return *this;
}

//writepfsave writes a save file with partition function data.
void writepfsave(char *filename, structure *ct,
			 PFPRECISION *w5, PFPRECISION *w3,
			 DynProgArray<PFPRECISION> *v, DynProgArray<PFPRECISION> *w, DynProgArray<PFPRECISION> *wmb, DynProgArray<PFPRECISION> *wl, DynProgArray<PFPRECISION> *wlc, DynProgArray<PFPRECISION> *wmbl, DynProgArray<PFPRECISION> *wcoax,
			 forceclass *fce, bool *mod, bool *lfce, pfdatatable *data) {
	int i,j,k,l,m,n,o,p;
	

	ofstream sav(filename,ios::binary);

	//write the save file information so that the sequence can be re-folded,
		//include thermodynamic data to prevent traceback errors

	//start with a version number of the savefile
	short vers=pfsaveversion;
	write(&sav,&vers); //save a version of the save file

	//start with structure information
	int SequenceLength = ct->GetSequenceLength();
	write(&sav,&(SequenceLength));
	write(&sav,&(ct->intermolecular));
	write(&sav,&data->scaling);

	int constraint;

	//Write out the number of forced pairs
	constraint = ct->GetNumberofPairs();
	write(&sav,&(constraint));
	//Write out the list of forced pairs
	for (i=0;i<ct->GetNumberofPairs();i++) {
		constraint = ct->GetPair5(i);
		write(&sav,&(constraint));
		constraint = ct->GetPair5(i);
		write(&sav,&(constraint));
	}
	for (i=0;i<=ct->GetSequenceLength();i++) {
		write(&sav,&(ct->hnumber[i]));
		sav.write(&(ct->nucs[i]),1);

	}

	for (i=0;i<=2*ct->GetSequenceLength();i++) write(&sav,&(ct->numseq[i]));

	//Write out the number of nucleotides forced double-stranded
	constraint = ct->GetNumberofDoubles();
	write(&sav,&(constraint));

	//Write out the nucleotides forced double stranded
	for (i=0;i<ct->GetNumberofDoubles();i++) {
		constraint = ct->GetDouble(i);
		write(&sav,&(constraint));

	}

	if (ct->intermolecular) {
		for (i=0;i<3;i++) write(&sav,&(ct->inter[i]));

	}

	//Write the number of nucleotides forced single stranded
	constraint = ct->GetNumberofSingles();
	write(&sav,&(constraint));

	//Write the list of nucleotides forced single stranded
	for (i=0;i<ct->GetNumberofSingles();i++) {
		constraint = ct->GetSingle(i);
		write(&sav,&(constraint));

	}

	//Write the number of nucleotides that are modified 
	constraint = ct->GetNumberofModified();
	write(&sav,&(constraint));

	//Write the list of modified nucleotides
	for (i=0;i<ct->GetNumberofModified();i++) {
		constraint = ct->GetModified(i);
		write(&sav,&(constraint));

	}

	//Write the number of Us in GU pairs
	constraint = ct->GetNumberofGU();
	write(&sav,&(constraint));

	//Write the list of Us in GU pairs
	for (i=0;i<ct->GetNumberofGU();i++) {
		constraint = ct->GetGUpair(i);
		write(&sav,&(constraint));

	}
	string label=ct->GetSequenceLabel();
	write(&sav,&(label));

	write(&sav,&(ct->templated));
	if (ct->templated) {
		for (i=0;i<=ct->GetSequenceLength();i++) {
			for (j=0;j<=i;j++) write(&sav,&(ct->tem[i][j]));

		}
	}

	//write the SHAPE data (for pseudo-free energy constraints)
	write(&sav,&(ct->shaped));
	if (ct->shaped) {
		for (i=0;i<=2*ct->GetSequenceLength();i++) write(&sav,&(ct->SHAPE[i]));
		for (i=0;i<=2*ct->GetSequenceLength();i++) write(&sav,&(ct->SHAPEss[i]));

	}

	//now write the array class data for v, w, and wmb:
	for (i=0;i<=ct->GetSequenceLength();i++) {
		write(&sav,&(w3[i]));
		write(&sav,&(w5[i]));
		for (j=0;j<=ct->GetSequenceLength();j++) {
			write(&sav,&(v->dg[i][j+i]));
			write(&sav,&(w->dg[i][j+i]));
			write(&sav,&(wmb->dg[i][j+i]));
			write(&sav,&(wmbl->dg[i][j+i]));
			write(&sav,&(wl->dg[i][j+i]));
			write(&sav,&(wlc->dg[i][j+i]));
			write(&sav,&(wcoax->dg[i][j+i]));
			writesinglechar(&sav,&(fce->dg[i][j]));
		}
	}

	write(&sav,&(w3[ct->GetSequenceLength()+1]));
	for (i=0;i<=2*ct->GetSequenceLength();i++) {
		write(&sav,&(lfce[i]));
		write(&sav,&(mod[i]));

	}

	//write the complete alphabet informtion -- this will be read into a datatable
	write(&sav, &(data->alphabet));
	write(&sav, &(data->pairing));
	write(&sav, &(data->not_pairing));
	write(&sav, &(data->non_interacting));
	write(&sav, &(data->linker));


	//write the alphabet information
	write(&sav, &(data->alphabet));
	write(&sav, &(data->pairing));

	//now write the thermodynamic data:
	write(&sav,&(data->pftemp));
	for (i=0;i<5;i++) write(&sav,&(data->poppen[i]));
	write(&sav,&(data->maxpen));
	for (i=0;i<11;i++) write(&sav,&(data->eparam[i]));
	for (i=0;i<31;i++) {
		write(&sav,&(data->inter[i]));
		write(&sav,&(data->bulge[i]));
		write(&sav,&(data->hairpin[i]));

	}
	for (i=0;i<data->alphabet.size();i++) {
		for (j=0;j<data->alphabet.size();j++) {
			for (k=0;k<data->alphabet.size();k++) {
				for (l=0;l<3;l++) {
					write(&sav,&(data->dangle[i][j][k][l]));
				}
				for (l=0;l<data->alphabet.size();l++) {
					write(&sav,&(data->stack[i][j][k][l]));
					write(&sav,&(data->tstkh[i][j][k][l]));
					write(&sav,&(data->tstki[i][j][k][l]));
					write(&sav,&(data->coax[i][j][k][l]));
					write(&sav,&(data->tstackcoax[i][j][k][l]));
					write(&sav,&(data->coaxstack[i][j][k][l]));
					write(&sav,&(data->tstack[i][j][k][l]));
					write(&sav,&(data->tstkm[i][j][k][l]));
					write(&sav,&(data->tstki23[i][j][k][l]));
					write(&sav,&(data->tstki1n[i][j][k][l]));
					for (m=0;m<data->alphabet.size();m++) {
						for (n=0;n<data->alphabet.size();n++) {
							write(&sav,&(data->iloop11[i][j][k][l][m][n]));
							for (o=0;o<data->alphabet.size();o++) {
								if (data->pairing[i][j]&& data->pairing[n][o]) write(&sav,&(data->iloop21[i][j][k][l][m][n][o]));
								for (p=0;p<data->alphabet.size();p++) {
									if (data->pairing[i][k]&& data->pairing[j][l])
										write(&sav,&(data->iloop22[i][j][k][l][m][n][o][p]));
								}
							}
						}
					}
				}
			}
		}
	}
	write(&sav,&(data->numoftloops));
	for (i=0;i<data->numoftloops;i++) {
		write(&sav,&(data->itloop[i]));
		write(&sav,&(data->tloop[i]));

	}
	write(&sav,&(data->numoftriloops));
	for (i=0;i<data->numoftriloops;i++) {
		write(&sav,&(data->itriloop[i]));
		write(&sav,&(data->triloop[i]));

	}
	write(&sav,&(data->numofhexaloops));
	for (i=0;i<data->numofhexaloops;i++) {
		write(&sav,&(data->ihexaloop[i]));
		write(&sav,&(data->hexaloop[i]));

	}
	write(&sav,&(data->auend));
	write(&sav,&(data->AUappliestoGU));
	write(&sav,&(data->gubonus));
	write(&sav,&(data->cint));
	write(&sav,&(data->cslope));
	write(&sav,&(data->c3));
	write(&sav,&(data->efn2a));
	write(&sav,&(data->efn2b));
	write(&sav,&(data->efn2c));
	write(&sav,&(data->init));
	write(&sav,&(data->mlasym));
	write(&sav,&(data->strain));
	write(&sav,&(data->prelog));
	write(&sav,&(data->singlecbulge));
	write(&sav,&(data->maxintloopsize));
	
	for(int i=0;i<=data->alphabet.size();i++){
		for (int j=0;j<=data->alphabet.size();j++){
			write(&sav,&(data->penalties[i][j]));
		}
	}

	sav.close();


}


//When considering mismatch at the end of a helix, consult this function to check
//	whether the nucs are required to pair
PFPRECISION pfchecknp(bool lfce1,bool lfce2) {
	if (lfce1||lfce2) return ZERO;
	else return ONE;
}

void size4D(pVector4D& matrix, const int sz) {
	matrix.resize(sz);
	for (int i=0;i<sz;++i) {
		matrix[i].resize(sz);
		for (int j=0;j<sz;++j) {
			matrix[i][j].resize(sz);
			for (int k=0;k<sz;++k) matrix[i][j][k].resize(sz);
		}
	}
}

void pfdatatable::allocate_data_tables(const int sz /* =6 */){
	poppen.resize(5);
	eparam.resize(11);
	inter.resize(31);
	bulge.resize(31);
	hairpin.resize(31);

	tloop.resize(maxtloop+1);
	triloop.resize(maxtloop+1);
	hexaloop.resize(maxtloop+1);
	itloop.resize(maxtloop+1);
	itriloop.resize(maxtloop+1);
	ihexaloop.resize(maxtloop+1);

	size4D(stack,sz);
	size4D(tstkh,sz);
	size4D(tstki,sz);
	size4D(coax,sz);
	size4D(tstackcoax,sz);
	size4D(coaxstack,sz);
	size4D(tstack,sz);
	size4D(tstkm,sz);
	size4D(tstki23,sz);
	size4D(tstki1n,sz);

	dangle.resize(sz);
	for(int i=0;i<sz;i++){
		dangle[i].resize(sz);
		for(int j=0;j<sz;j++){
			dangle[i][j].resize(sz);
			for(int k=0;k<sz;k++){
				dangle[i][j][k].resize(3);
			}
		}
	}

	iloop11.resize(sz);
	iloop21.resize(sz);
	iloop22.resize(sz);
	for (int i=0;i<sz;++i) {
		iloop11[i].resize(sz);
		iloop21[i].resize(sz);
		iloop22[i].resize(sz);
		for (int j=0;j<sz;++j) {
			iloop11[i][j].resize(sz);
			iloop21[i][j].resize(sz);
			iloop22[i][j].resize(sz);
			for (int k=0;k<sz;++k) {
				iloop11[i][j][k].resize(sz);
				iloop21[i][j][k].resize(sz);
				iloop22[i][j][k].resize(sz);
				for (int l=0;l<sz;++l) {
					iloop11[i][j][k][l].resize(sz);
					iloop21[i][j][k][l].resize(sz);
					iloop22[i][j][k][l].resize(sz);
					for (int m=0;m<sz;m++) {
						iloop11[i][j][k][l][m].resize(sz);
						iloop21[i][j][k][l][m].resize(sz);
						iloop22[i][j][k][l][m].resize(sz);
						for (int n=0;n<sz;n++) {
							iloop21[i][j][k][l][m][n].resize(sz);
							iloop22[i][j][k][l][m][n].resize(sz);
							for (int o=0;o<sz;o++) {
								iloop22[i][j][k][l][m][n][o].resize(sz);
							}

						}
					}
				}
			}
		}
	}
	
	penalties = new PFPRECISION *[alphabet.size()+1];
	for(int i=0;i<=alphabet.size();i++){
		penalties[i] = new PFPRECISION [alphabet.size()+1];
	}
	
}

pfdatatable::pfdatatable() {
	penalties = NULL;
}

pfdatatable::pfdatatable(datatable *data, const PFPRECISION &initialScaling, const double &Temp) {
	//the partition function datatable needs to be initialized from the datatable
	short i,j,k,l,m,n,o,p;

 #ifdef PF_LOG_CALC
 	// scaling = log(initialScaling);
	scaling = ONE;
 #else
	scaling = initialScaling;
 #endif
 	//store the temperature in the pfdatatable
	pftemp = Temp;	 

	//copy the alphabet information:

	//make the alphabet vector large enough for each letter
	alphabet.resize(data->alphabet.size());
	for (int letter=0;letter<data->alphabet.size();++letter) {
		//alocate space for each exact letter for the current dimension
		alphabet[letter].resize(data->alphabet[letter].size());
		for (int exactletter=0; exactletter<data->alphabet[letter].size(); ++exactletter) {
			//copy each entry
			alphabet[letter][exactletter]=data->alphabet[letter][exactletter];
		}
	}

	baseU=data->baseU;
	baseA = data->baseA;
	

	//Copy the pairing information:
	pairing.resize(alphabet.size());
	for (int letter1=0; letter1<alphabet.size();++letter1) {
		pairing[letter1].resize(alphabet.size());
		for (int letter2=0; letter2<alphabet.size(); ++letter2) {
		//copy the exact contnts stored in disk
			pairing[letter1][letter2]=data->pairing[letter1][letter2];

		}
	}

	//copy also the not pairing, non_interacting, linker, and linker_ints info
	not_pairing.resize(data->not_pairing.size());
	for (int letter=0;letter<not_pairing.size();++letter) not_pairing[letter] = data->not_pairing[letter];

	non_interacting.resize(data->non_interacting.size());
	for (int letter=0;letter<non_interacting.size();++letter) non_interacting[letter] = data->non_interacting[letter];

	linker.resize(data->linker.size());
	for (int letter=0;letter<linker.size();++letter) linker[letter] = data->linker[letter];


	allocate_data_tables(alphabet.size());

	//allocate_data_tables();
	const int sz = alphabet.size();

	for (i=1;i<5;i++) poppen[i]=boltzman(data->poppen[i],pftemp);
	maxpen = boltzman(data->maxpen,pftemp);
	for (i=1;i<11;i++) eparam[i] = boltzman(data->eparam[i],pftemp);
	maxintloopsize = data->eparam[7];
	
	inter.resize(data->inter.size());
	bulge.resize(data->bulge.size());
	hairpin.resize(data->hairpin.size());

	for (i=1;i<data->inter.size();i++) {
		//cout << data->inter[i] << endl;
		// inter[i] = boltzman(data->inter[i],pftemp);
		// cout << inter[i];
		// inter[i] = pow(scaling,(double) i+2);
		// cout << inter[i];
		inter[i] = PROD(boltzman(data->inter[i],pftemp) POWSCALING2(i+2));
	}
	for (i=1;i<data->bulge.size();i++) 
		bulge[i] = PROD(boltzman(data->bulge[i],pftemp) POWSCALING2(i+2));
	
	for (i=1;i<data->hairpin.size();i++) 
		hairpin[i] = PROD(boltzman(data->hairpin[i],pftemp) POWSCALING2(i+2));

	
	for (i=0;i<sz;i++) {
		for (j=0;j<sz;j++) {
			for (k=0;k<sz;k++) {
				for (l=1;l<3;l++) {
					#ifdef SIMPLEMBLOOP
					//In the case of simple multibranch loops, dangles no longer
						//occupy a position in the sequence and do not need a scaling
						//factor
					dangle[i][j][k][l] = boltzman(data->dangle[i][j][k][l],pftemp);
					#else //!SIMPLEMBLOOP
					dangle[i][j][k][l] = PROD(boltzman(data->dangle[i][j][k][l],pftemp) SCALING2);

					#endif
				}
				for (l=0;l<sz;l++) {
					stack[i][j][k][l]=PROD(boltzman(data->stack[i][j][k][l],pftemp) POWSCALING2(2));
					tstkh[i][j][k][l]=boltzman(data->tstkh[i][j][k][l],pftemp);
					tstki[i][j][k][l]=boltzman(data->tstki[i][j][k][l],pftemp);
					//short foo = data->coax[i][j][k][l];
					coax[i][j][k][l]=boltzman(data->coax[i][j][k][l],pftemp);
					tstackcoax[i][j][k][l]=PROD(boltzman(data->tstackcoax[i][j][k][l],pftemp) POWSCALING2(2));
					coaxstack[i][j][k][l] = boltzman(data->coaxstack[i][j][k][l],pftemp);
					tstack[i][j][k][l]=PROD(boltzman(data->tstack[i][j][k][l],pftemp) POWSCALING2(2));
					tstkm[i][j][k][l]=PROD(boltzman(data->tstkm[i][j][k][l],pftemp) POWSCALING2(2));
					tstki23[i][j][k][l]=boltzman(data->tstki23[i][j][k][l],pftemp);
					tstki1n[i][j][k][l]=boltzman(data->tstki1n[i][j][k][l],pftemp);
					for (m=0;m<sz;m++) {
						for (n=0;n<sz;n++) {
							iloop11[i][j][k][l][m][n]=PROD(boltzman(data->iloop11[i][j][k][l][m][n],pftemp) POWSCALING2(4));
							for (o=0;o<sz;o++) {
								iloop21[i][j][k][l][m][n][o]=PROD(boltzman(data->iloop21[i][j][k][l][m][n][o],pftemp) POWSCALING2(5));
								for (p=0;p<sz;p++) {
									iloop22[i][j][k][l][m][n][o][p]=PROD(boltzman(data->iloop22[i][j][k][l][m][n][o][p],pftemp) POWSCALING2(6));
								}
							}
						}
					}
				}
			}
		}
	}
	numoftloops = data->numoftloops;
	for (i=0;i<data->numoftloops;i++) {
		itloop[i]=data->tloop[i][0];
		tloop[i] = PROD(boltzman(data->tloop[i][1],pftemp) POWSCALING2(6));

	}
	numoftriloops=data->numoftriloops;
	for (i=0;i<data->numoftriloops;i++) {
		itriloop[i] = data->triloop[i][0];
		triloop[i] = PROD(boltzman(data->triloop[i][1],pftemp) POWSCALING2(5));
	}
	numofhexaloops=data->numofhexaloops;
	for (i=0;i<data->numofhexaloops;i++) {
		ihexaloop[i]=data->hexaloop[i][0];
		hexaloop[i] = PROD(boltzman(data->hexaloop[i][1],pftemp) POWSCALING2(8));

	}
	auend = boltzman(data->auend,pftemp);
	AUappliestoGU = data->AUappliestoGU;
	gubonus = boltzman(data->gubonus,pftemp);
	cint = boltzman(data->cint,pftemp);
	cslope = boltzman(data->cslope,pftemp);
	c3 = boltzman(data->c3,pftemp);
	efn2a = boltzman(data->efn2a,pftemp);
	efn2b = boltzman(data->efn2b,pftemp);
	efn2c = boltzman(data->efn2c,pftemp);
	init = boltzman(data->init,pftemp);
	mlasym = boltzman(data->mlasym,pftemp);
	strain = boltzman(data->strain,pftemp);
	prelog = data->prelog/conversionfactor;
	singlecbulge = boltzman(data->singlecbulge,pftemp);

	penalties = new PFPRECISION *[alphabet.size()+1];
	

	// Fill end penalties
	for(int i=0;i<=alphabet.size();i++){
		penalties[i] = new PFPRECISION [alphabet.size()+1];
		for (int j=0;j<=alphabet.size();j++){
				if (AUappliestoGU){
					if (i==baseU || j==baseU)
   						penalties[i][j] = auend;
					else
						penalties[i][j] = ONE;
				}
				else {
					if (i == baseA || j == baseA)
						penalties[i][j] = auend;
					else
						penalties[i][j] = ONE;

				}
		}
	}
}


pfdatatable::~pfdatatable() {
	if (penalties != NULL) {
		//delete the penalaties tables

			for (int i = 0; i <= alphabet.size(); i++) {
				delete[] penalties[i];
			}

		delete[] penalties;
	}
}

//calculate the energy of stacked base pairs

//oriented by:
//5' i ip 3'
//   |  |
//   j jp

PFPRECISION erg1(int i,int j,int ip,int jp,structure *ct, pfdatatable *data)
{
		PFPRECISION energy;

		 if ((i==(ct->GetSequenceLength()))||(j==((ct->GetSequenceLength())+1))) {
      		//this is not allowed because n and n+1 are not cavalently attached
			energy = (PFPRECISION) ZERO;
		}
		else {
      		energy = PROD(data->stack[(ct->numseq[i])][(ct->numseq[j])]
				[(ct->numseq[ip])][(ct->numseq[jp])],data->eparam[1]);

				if (ct->shaped) {
					energy=PROD(energy,ct->SHAPE[i]);
					energy=PROD(energy,ct->SHAPE[j]);
					energy=PROD(energy,ct->SHAPE[ip]);
					energy=PROD(energy,ct->SHAPE[jp]);
				}

				if ( ct->experimentalPairBonusExists ) {
				    energy = PROD(energy , ct->EX[i][j] , ct->EX[ip][jp]);

				}
		}

		#ifdef equiout
		ofstream kout;
		kout.open("k.out",ofstream::app);
		kout << "k base pair stack ("<<i<<","<<j<<","<<ip<<","<<jp<<") ="<<energy<< "\n";
		kout.close();
		#endif

		return energy;
}

PFPRECISION multi_erg1(int i, int j, int ip, int jp, structure* ct, pfdatatable* data)
{	// This is eqilibirium constant
	PFPRECISION energy = ONE;
	int N = ct->GetNumberofSequences();
	for (int k = 0; k < N; k++) {
		energy = PROD(energy, erg1(i, j, ip, jp, ct->get_individual_sequence(k), data));
	}
	return energy/N;
}


//calculate the energy of a bulge/internal loop
//where i is paired to j; ip is paired to jp; ip > i; j > jp
PFPRECISION erg2(int i,int j,int ip,int jp,structure *ct, pfdatatable *data,
	char a, char b)
{
	int size,size1,size2, lopsid, count;
	PFPRECISION energy,loginc;
	/* size,size1,size2 = size of a loop
		energy = energy calculated
		loginc = the value of a log used in large hairpin loops
	*/
   	if (((i<=(ct->GetSequenceLength()))&&(ip>(ct->GetSequenceLength())))||((
      	jp<=(ct->GetSequenceLength()))&&(j>(ct->GetSequenceLength())))) {
         //A loop cannot contain the ends of the sequence

         return ZERO;
      }

      size1 = ip-i-1;
		size2 = j - jp - 1;

	if ((a>0)||(b>0)) {
      	if ((a&DUBLE)||(b&DUBLE)) return ZERO;//the loop contains a nuc that
      		//should be double stranded

	}

      //a typical internal or bulge loop:
		//size1 = ip-i-1;
		//size2 = j - jp - 1;
		if (size1==0||size2==0) {//bulge loop

			size = size1+size2;
			if (size==1) {
				count = 1;
				energy = PROD(data->stack[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][ct->numseq[jp]],
						data->bulge[size],PROD(data->eparam[2] POWSCALING(-2)));
				if (size1==1)  {
					//give bonus to C adjacent to single C bulge
					if ((ct->IsNuc(i+1,'C'))&&((ct->IsNuc(i+2,'C'))||(ct->IsNuc(i,'C')))) energy= PROD(energy,data->singlecbulge);
//					if ((ct->IsNuc(i+1,'C')||ct->IsNuc(i+1,'c'))&&((ct->IsNuc(i+2,'C')||ct->IsNuc(i+2,'c'))||(ct->IsNuc(i,'C')||ct->IsNuc(i,'c')))) energy= PROD(energy,data->singlecbulge);
					//if (ct->numseq[i+1]==2&&(ct->numseq[i+2]==2||ct->numseq[i]==2)) energy= energy*data->singlecbulge;

				}
				else {
					//give bonus to C adjacent to single C bulge
					if ((ct->IsNuc(j-1,'C'))&&((ct->IsNuc(j-2,'C'))||(ct->IsNuc(j,'C')))) energy= PROD(energy,data->singlecbulge);
//					if ((ct->IsNuc(j-1,'C')||ct->IsNuc(j-1,'c'))&&((ct->IsNuc(j-2,'C')||ct->IsNuc(j-2,'c'))||(ct->IsNuc(j,'C')||ct->IsNuc(j,'c')))) energy= PROD(energy,data->singlecbulge);
					//if (ct->numseq[j-1]==2&&(ct->numseq[j-2]==2||ct->numseq[j]==2)) energy=energy*data->singlecbulge;

				}
			}
			else if (size>30) {
				loginc = data->prelog*log(size/30.0);
				energy = PROD(DIV(data->bulge[30],PF_EXP_LINEAR(loginc/(RKC*data->pftemp))),data->eparam[2]);
				energy = PROD(energy,penalty(i,j,ct,data),penalty(jp,ip,ct,data) POWSCALING(size-30));
			}
			else {
         		energy = PROD(data->bulge[size],data->eparam[2]);
				energy = PROD(energy,penalty(i,j,ct,data),penalty(jp,ip,ct,data));
			}
		}
		else {//internal loop
			size = size1 + size2;
			lopsid = abs(size1-size2);

			if (size>30) {
				loginc = data->prelog*log(size/30.0);
				if (size1==1||size2==1) {
            		energy = PROD(data->tstki1n[ct->numseq[i]][ct->numseq[j]]
						[ct->numseq[i+1]][ct->numseq[j-1]],
						data->tstki1n[ct->numseq[jp]][ct->numseq[ip]]
						[ct->numseq[jp+1]][ct->numseq[ip-1]],
						DIV(data->inter[30],PF_EXP_LINEAR(loginc/(RKC*data->pftemp))), data->eparam[3],
						max(data->maxpen,
						POWER(data->poppen[min(2,min(size1,size2))],lopsid))

						POWSCALING(size-30));

				}

				else {
					energy = PROD(data->tstki[ct->numseq[i]][ct->numseq[j]]
						[ct->numseq[i+1]][ct->numseq[j-1]],
						data->tstki[ct->numseq[jp]][ct->numseq[ip]]
						[ct->numseq[jp+1]][ct->numseq[ip-1]],
						DIV(data->inter[30], PF_EXP_LINEAR(loginc/(RKC*data->pftemp))),data->eparam[3] ,
						max(data->maxpen,
						POWER(data->poppen[min(2,min(size1,size2))],lopsid)) POWSCALING(size-30));
				}
			}
			else if ((size1==2)&&(size2==2))//2x2 internal loop
			    energy = data->iloop22[ct->numseq[i]][ct->numseq[ip]]
					[ct->numseq[j]][ct->numseq[jp]]
					[ct->numseq[i+1]][ct->numseq[i+2]]
					[ct->numseq[j-1]][ct->numseq[j-2]];

			else if ((size1==1)&&(size2==2)) {//2x1 internal loop
				energy = data->iloop21[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]]
					[ct->numseq[j-1]][ct->numseq[jp+1]][ct->numseq[ip]][ct->numseq[jp]];

			}
			else if ((size1==2)&&(size2==1)) {//1x2 internal loop
				energy = data->iloop21[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]]
					[ct->numseq[ip-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]];

			}

			else if (size==2) //a single mismatch

				energy = data->iloop11[ct->numseq[i]][ct->numseq[i+1]][ct->numseq[ip]]
					[ct->numseq[j]][ct->numseq[j-1]][ct->numseq[jp]];
			else if (size1==1||size2==1) { //this loop is lopsided
         	//this is a 1xn loop:
				energy = PROD(data->tstki1n[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] ,
						data->tstki1n[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] ,
						data->inter[size] , data->eparam[3] ,
					max(data->maxpen,POWER(data->poppen[min(2,min(size1,size2))],lopsid)));
			}

			else if ((size1==2&&size2==3)||(size1==3&&size2==2)) {
			//this is a 2x3 loop
				energy = PROD(data->tstki23[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] ,
					data->tstki23[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] ,
					data->inter[size] , data->eparam[3] ,
					max(data->maxpen,POWER(data->poppen[min(2,min(size1,size2))],lopsid)));

			}
			else {
         		energy = PROD(data->tstki[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] ,
					data->tstki[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] ,
					data->inter[size] , data->eparam[3] ,
					max(data->maxpen,POWER(data->poppen[min(2,min(size1,size2))],lopsid)));
			}
		}
		#ifdef equiout
		ofstream kout;
		kout.open("k.out",ofstream::app);
		kout << "k internal/bulge ("<<i<<","<<j<<","<<ip<<","<<jp<<") ="<<energy<< "\n";
		kout.close();
		#endif

		return energy;
}

//calculate the energy of a bulge/internal loop
//where i is paired to j; ip is paired to jp; ip > i; j > jp
PFPRECISION multi_erg2(int i, int j, int ip, int jp, structure* ct, pfdatatable* data,
	char a, char b) {
	PFPRECISION energy = ONE;
	int N = ct->GetNumberofSequences();
	if (((i <= (ct->GetSequenceLength())) && (ip > (ct->GetSequenceLength()))) || ((
		jp <= (ct->GetSequenceLength())) && (j > (ct->GetSequenceLength())))) {
		//A loop cannot contain the ends of the sequence

		return ZERO;
	}
	for (int k = 0; k < N; k++) {

		structure* ct_i = ct->get_individual_sequence(k);
		int seq_length = ct_i->GetSequenceLength();

		int no_of_nuc_i = ip - i - 1 - ct_i->no_of_gaps_matrix[i][ip];
		int no_of_nuc_j = j - jp - 1 - ct_i->no_of_gaps_matrix[jp][j];

		// if no nucleotides in loop then treat as stack
		if (no_of_nuc_i == 0 && no_of_nuc_j == 0) {
			energy = PROD(energy, erg1(i, j, ip, jp, ct_i, data));
		}

		// if no gaps in internal loop
		else if (no_of_nuc_i == ip - i - 1 && no_of_nuc_j == j - jp - 1) {
			energy = PROD(energy, erg2(i, j, ip, jp, ct_i, data, a, b));
		}

		// if gaps in internal loop
		else if (no_of_nuc_i > 0 || no_of_nuc_j > 0) {
			int i_ = 1;
			int ip_ = 1 + no_of_nuc_i + 1;
			int jp_ = ip_ + 1;
			int j_ = jp_ + no_of_nuc_j + 1;

#ifdef SMP
#pragma omp critical
#endif
			{
				structure* m = remove_gaps(i, j, ip, jp, ct_i, ct->tmp_ct_i);
				//  added 1 to accomodate for the gap introduced by remove_gaps between ip and jp.
				// This will allow sequence to work with alternative bulge loop code in erg2
				energy = PROD(energy, erg2(i_, j_ + 1, ip_, jp_ + 1, m, data, a, b));
				ct->tmp_ct_i->deallocate();
			}
		}

	}

	return energy/N;
}

//calculate the energy of the intorior part of a internal loop
//where i is paired to j; ip is paired to jp; ip > i; j > jp
PFPRECISION erg2in(int i,int j,int ip,int jp,structure *ct, pfdatatable *data,char a,char b)
{
	int size,size1,size2, lopsid;
	PFPRECISION energy;

	if ((a>0)||(b>0)) {
      	if ((a&DUBLE)||(b&DUBLE)) return ZERO;//the loop contains a nuc that
      		//should be double stranded

	}

   	if (((i<=(ct->GetSequenceLength()))&&(ip>(ct->GetSequenceLength())))||((
      	jp<=(ct->GetSequenceLength()))&&(j>(ct->GetSequenceLength())))) {
         //A loop cannot contain the ends of the sequence

         return ZERO;
      }
      size1 = ip-i-1;
	  size2 = j - jp - 1;
      size = size1 + size2;
	  lopsid = abs(size1-size2);
	  if (size1 != 0 && size2 !=0)
	 	{
        energy = 	PROD(data->tstki[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] ,
					 data->eparam[3] ,
					max(data->maxpen,POWER(data->poppen[min(2,min(size1,size2))],lopsid)));
		}

		return energy;
}

//calculate the energy of the exterior part of internal loop
//where i is paired to j; ip is paired to jp; ip > i; j > jp
PFPRECISION erg2ex(int i,int j,int size, structure *ct, pfdatatable *data)
{
	PFPRECISION energy,loginc;
	/* size,size1,size2 = size of a loop
		energy = energy calculated
		loginc = the value of a log used in large hairpin loops
	*/
				if (size>30) {
				loginc = (data->prelog)*log((PFPRECISION ((size))/30.0));
				energy = PROD(data->tstki[ct->numseq[i]][ct->numseq[j]]
						[ct->numseq[i+1]][ct->numseq[j-1]],
						DIV(data->inter[30], PF_EXP_LINEAR(loginc/(RKC*data->pftemp))) POWSCALING(size-30));
			}
						else
         		energy = PROD(data->tstki[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] ,
						 data->inter[size]	);
		return energy;
}

//calculate the energy of a hairpin loop:
#ifndef equiout
PFPRECISION erg3(int i,int j,structure *ct, pfdatatable *data,char dbl) {
#else
PFPRECISION erg3indirect(int i,int j,structure *ct, pfdatatable *data,char dbl) {
#endif
int size,count,key,k;
PFPRECISION energy,loginc;
	/* size,size1,size2 = size of a loop
		energy = energy calculated
		loginc = the value of a log used in large hairpin loops
	*/

	if ((i<=(ct->GetSequenceLength()))&&(j>(ct->GetSequenceLength()))) {
      	//A hairpin cannot contain the ends of the sequence

         return ZERO;
    }

	if (dbl&DUBLE) return ZERO;//the loop contains a base that should be
      										//double stranded

    else if (dbl&INTER) {//intermolecular interaction
      	//intermolecular "hairpin" free energy is that of intermolecular
         //	initiation times the stacked mismatch

         energy =  PROD(data->tstack[ct->numseq[i]][ct->numseq[j]]
         	[ct->numseq[i+1]][ct->numseq[j-1]],penalty(i,j,ct,data) POWSCALING(j-i-1));

		//Add other states that can exist, i.e. 5' and 3' dangles:
		if (ct->numseq[i+1]!=5&&ct->numseq[j-1]!=5) {
			//Add a 3' dangling end

			energy=SUM(energy,PROD(erg4(i,j,i+1,1,ct,data,false) POWSCALING(j-i),penalty(i,j,ct,data)));

			//Add a 5' dangling end

			energy=SUM(energy,PROD(erg4(i,j,j-1,2,ct,data,false) POWSCALING(j-i),penalty(i,j,ct,data)));

			//Add the case where nothing stacks

			energy=SUM(energy,PROD(penalty(i,j,ct,data) POWSCALING(j-i+1)));

		}
		else if (ct->numseq[i+1]!=5||ct->numseq[j-1]!=5) {
			//Add the case where nothing stacks

			energy=SUM(energy,PROD(penalty(i,j,ct,data) POWSCALING(j-i-1)));

		}

         return PROD(data->init,energy);
    }

		size = j-i-1;

		if (size>30) {
			loginc = ((data->prelog)*log((PFPRECISION ((size))/30.0)));

			energy = PROD(data->tstkh[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[i+1]][ct->numseq[j-1]]
				, DIV(data->hairpin[30],PF_EXP_LINEAR(loginc/(RKC*data->pftemp))), data->eparam[4] POWSCALING(size-30));
		}
		else if (size<3) {
      		energy = PROD(data->hairpin[size], data->eparam[4]);
				//if (ct->numseq[i]==4||ct->numseq[j]==4) 
			energy = PROD(energy, penalty(i,j,ct,data));//exp(-.6/(RKC*data->pftemp));
		}
		else if (size==4) {
			int digit_factor = 1;
			key = 0;
			for (count = 0; count < 6; ++count) {
				key += ct->numseq[i+count] * digit_factor;
				digit_factor = digit_factor * data->alphabet.size();
			}
			for (count=0;count<data->numoftloops;count++) {
				if (key==data->itloop[count]) {
					return data->tloop[count];
				}
			}
			energy = PROD(data->tstkh[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[i+1]][ct->numseq[j-1]],
				data->hairpin[size], data->eparam[4]);
		}
		else if (size==3) {
			int digit_factor = 1;
			key = 0;
			for (count = 0; count < 5; ++count) {
				key += ct->numseq[i+count] * digit_factor;
				digit_factor = digit_factor * data->alphabet.size();
			}
			for (count=0;count<data->numoftriloops;count++) {
				if (key==data->itriloop[count]) return data->triloop[count];
			}

			energy =	PROD(data->hairpin[size], data->eparam[4], penalty(i,j,ct,data));
		}
		else if (size==6) {
			int digit_factor = 1;
			key = 0;
			for (count = 0; count < 8; ++count) {
				key += ct->numseq[i+count] * digit_factor;
				digit_factor = digit_factor * data->alphabet.size();
			}
			for (count=0;count<data->numofhexaloops;count++) {
				if (key==data->ihexaloop[count]) {
					return data->hexaloop[count];
				}
			}

			energy = PROD(data->tstkh[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[i+1]][ct->numseq[j-1]],
				data->hairpin[size], data->eparam[4]);
		}

		else {
			energy = PROD(data->tstkh[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[i+1]][ct->numseq[j-1]],
				data->hairpin[size], data->eparam[4]);
		}


		//check for GU closeure preceded by GG
		
	  //if ((ct->IsNuc(i,'G')||ct->IsNuc(i,'g'))&&(ct->IsNuc(j,'U')||ct->IsNuc(j,'u'))) {
      if ((ct->IsNuc(i,'G'))&&(ct->IsNuc(j,'U'))) {
      //if (ct->numseq[i]==3&&ct->numseq[j]==4) {
      	if ((i>2&&i<ct->GetSequenceLength())||(i>ct->GetSequenceLength()+2))
			//if ((ct->IsNuc(i-1,'G')||ct->IsNuc(i-1,'g'))&&(ct->IsNuc(i-2,'G')||ct->IsNuc(i-2,'g'))) {
       		if ((ct->IsNuc(i-1,'G'))&&(ct->IsNuc(i-2,'G'))) {
       		//if (ct->numseq[i-1]==3&&ct->numseq[i-2]==3) {
         		energy = PROD(energy, data->gubonus);

         	}
      }

      //check for an oligo-c loop

      for (k=1;(k<=size);k++) {
//		if (!(ct->IsNuc(i+k,'C')||ct->IsNuc(i+k,'c'))) return energy;//this is not an oligo-C loop
		if (!(ct->IsNuc(i+k,'C'))) return energy;//this is not an oligo-C loop
       	//if (ct->numseq[i+k] != 2) return energy;//this is not an oligo-C loop
      }
      //this is a poly c loop so penalize
      if (size==3) return PROD(energy , data->c3);
      else return PROD(energy, data->cint, POWER(data->cslope,size));

}

PFPRECISION multi_erg3(int i, int j, structure* ct, pfdatatable* data, char dbl) {
	PFPRECISION energy = ONE;
	int N = ct->GetNumberofSequences();
	if ((i <= (ct->GetSequenceLength())) && (j > (ct->GetSequenceLength()))) {
		//A hairpin cannot contain the ends of the sequence
		return ZERO;
	}
	for (int k = 0; k < N; k++) {
		structure* ct_i = ct->get_individual_sequence(k);
		int no_of_nuc = j - i - 1 - ct_i->no_of_gaps_matrix[i][j];
		// if number of nucleotides in hairpin loop < 3
		if (no_of_nuc < 3) { 
			// global variable hairpinlt3k = boltzman(hairpinlt3, pftemp), pftemp=310.15 in defines.h
			energy = PROD(energy, PROD(hairpinlt3k POWSCALING2(no_of_nuc)));
		}
		// if no gaps
		else if (no_of_nuc == j - i - 1) { energy = PROD(energy, erg3(i, j, ct_i, data, dbl)); }
		// if gaps in hairpin loop
		else 
		{
#ifdef SMP
#pragma omp critical
#endif
			{
				int new_i = 1;
				int new_j = 1 + no_of_nuc + 1;
				energy = PROD(energy, erg3(new_i, new_j, remove_gaps(i, j, ct_i, ct->tmp_ct_hp), data, dbl));
				ct->tmp_ct_hp->deallocate();
			}		
		}
	}
	return energy/N;
}




#ifdef equiout
PFPRECISION erg3(int i,int j,structure *ct, pfdatatable *data,char dbl) {
	PFPRECISION energy;

	energy = erg3indirect(i,j,ct, data,dbl,pftemp);
	ofstream kout;
	kout.open("k.out",ofstream::app);
	kout << "k hairpin ("<<i<<","<<j<<") ="<<energy<< "\n";
	kout.close();
	return energy;
}
#endif

//calculate the energy of a dangling end:
PFPRECISION erg4(int i,int j,int ip,int jp,structure *ct, pfdatatable *data, bool lfce)
{
//dangling base
		// jp = 1 => 3' dangle
		// jp = 2 => 5' dangle
      if (lfce) return ZERO;//stacked nuc should be double stranded
	    #ifdef equiout
			ofstream kout;
			kout.open("k.out",ofstream::app);
			kout << "k dangle ("<<i<<","<<j<<","<<ip<<","<<jp<<") ="<<data->dangle[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][jp]<< "\n";
			kout.close();
		#endif

		return data->dangle[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][jp];

}

PFPRECISION multi_erg4(int i, int j, int ip, int jp, structure* ct, pfdatatable* data, bool lfce) {
	PFPRECISION energy = ONE;
	int N = ct->GetNumberofSequences();
	for (int k = 0; k < N; k++) {
		energy = PROD(energy, erg4(i, j, ip, jp, ct->get_individual_sequence(k), data, lfce));
	}
	return energy/N;
}


/*
//this function calculates whether a terminal pair i,j requires the end penalty
PFPRECISION penalty(int i,int j,structure* ct, pfdatatable *data) {
	#ifdef equiout
		ofstream kout;
		kout.open("k.out",ofstream::app);
		if (ct->numseq[i]==4||ct->numseq[j]==4) kout << "k penalty ("<<i<<","<<j<<") ="<<data->auend<< "\n";
		else kout << "k penalty ("<<i<<","<<j<<") ="<<1.0<< "\n";
		kout.close();
	#endif

	if (ct->numseq[i]==data->baseU || ct->numseq[j]==data->baseU)
   	return data->auend;
	else return ONE;//no end penalty
}
/*

inline PFPRECISION penalty(int i,int j,structure* ct, pfdatatable *data) {	
	return data->penalties[ct->numseq[i]][ct->numseq[j]];
}
//*/

//this function calculates whether a terminal pair i,j requires the end penalty
PFPRECISION penalty2(int i,int j, pfdatatable *data) {
   
	if (data->AUappliestoGU&&(std::find(data->alphabet[i].begin(), data->alphabet[i].end(), 'U') != data->alphabet[i].end())) {
		return data->auend;
	}
	else if (data->AUappliestoGU && (std::find(data->alphabet[j].begin(), data->alphabet[j].end(), 'U') != data->alphabet[j].end())) {
		return data->auend;
	}
	else if (!data->AUappliestoGU && (std::find(data->alphabet[i].begin(), data->alphabet[i].end(), 'A') != data->alphabet[i].end())) {
		return data->auend;
	}
	else if (!data->AUappliestoGU && (std::find(data->alphabet[j].begin(), data->alphabet[j].end(), 'A') != data->alphabet[j].end())) {
		return data->auend;
	}
	else return ONE;//no end penalty
	
	//if (i==4||j==4)
   	//return data->auend;
   //else return 1;//no end penalty

}

//this function calculates the free energy for coaxial stacking of pair i-j onto ip-jp
//	require that the backbone continues directly between j and ip without a nucleotide in
//	a canonical pair
//k indicates the intervening nuc in intervening mismatch that is not sandwiched between
//	j and ip
PFPRECISION ergcoaxflushbases(int i, int j, int ip, int jp, structure *ct, pfdatatable *data) {
		//flush stacking
		//remapped 10/11/2013 to match the order of the helical stack table
		return data->coax[ct->numseq[j]][ct->numseq[i]][ct->numseq[ip]][ct->numseq[jp]];

}

PFPRECISION ergcoaxinterbases1(int i, int j, int ip, int jp, structure *ct, pfdatatable *data) {
		//coaxial stacking with an intervening mismatch
			return PROD(data->tstackcoax[ct->numseq[j]][ct->numseq[i]]
				[ct->numseq[j+1]][ct->numseq[i-1]],
				data->coaxstack[ct->numseq[j+1]][ct->numseq[i-1]]
				[ct->numseq[ip]][ct->numseq[jp]]);

}

PFPRECISION ergcoaxinterbases2(int i, int j, int ip, int jp, structure *ct, pfdatatable *data) {
			return PROD(data->tstackcoax[ct->numseq[jp]][ct->numseq[ip]]
				[ct->numseq[jp+1]][ct->numseq[ip-1]],
				data->coaxstack[ct->numseq[j]][ct->numseq[i]]
				[ct->numseq[j+1]][ct->numseq[jp+1]]);
}

PFPRECISION multi_ergcoaxflushbases(int i, int j, int ip, int jp, structure* ct, pfdatatable* data) { 
	PFPRECISION energy = ONE;
	int N = ct->GetNumberofSequences();
	for (int k = 0; k < N; k++) {
		energy = PROD(energy, ergcoaxflushbases(i, j, ip, jp, ct->get_individual_sequence(k), data));
	}
	return energy / N;
}
PFPRECISION multi_ergcoaxinterbases1(int i, int j, int ip, int jp, structure* ct, pfdatatable* data) { 
	PFPRECISION energy = ONE;
	int N = ct->GetNumberofSequences();
	for (int k = 0; k < N; k++) {
		energy = PROD(energy, ergcoaxinterbases1(i, j, ip, jp, ct->get_individual_sequence(k), data));
	}
	return energy / N;
}
PFPRECISION multi_ergcoaxinterbases2(int i, int j, int ip, int jp, structure* ct, pfdatatable* data) { 
	PFPRECISION energy = ONE;
	int N = ct->GetNumberofSequences();
	for (int k = 0; k < N; k++) {
		energy = PROD(energy, ergcoaxinterbases2(i, j, ip, jp, ct->get_individual_sequence(k), data));
	}
	return energy / N;
}

PFPRECISION tstkm(int a, int b, int c, int d, structure* ct, pfdatatable* data) {
	return data->tstkm[ct->numseq[a]][ct->numseq[b]][ct->numseq[c]][ct->numseq[d]];
}

PFPRECISION multi_tstkm(int a, int b, int c, int d, structure* ct, pfdatatable* data) {
	PFPRECISION energy = ONE;
	int N = ct->GetNumberofSequences();
	for (int k = 0; k < N; k++) {
		energy = PROD(energy, tstkm(a, b, c, d, ct->get_individual_sequence(k), data));
	}
	return energy / N;
}

PFPRECISION tstack(int a, int b, int c, int d, structure* ct, pfdatatable* data) {
	return data->tstack[ct->numseq[a]][ct->numseq[b]][ct->numseq[c]][ct->numseq[d]];
}

PFPRECISION multi_tstack(int a, int b, int c, int d, structure* ct, pfdatatable* data) {
	PFPRECISION energy = ONE;
	int N = ct->number_of_sequences;
	for (int k = 0; k < N; k++) {
		energy = PROD(energy, tstack(a, b, c, d, ct->get_individual_sequence(k), data));
	}
	return energy / N;
}

PFPRECISION eparam(int a, structure* ct, pfdatatable* data) {
	return data->eparam[a];
}

PFPRECISION multi_eparam(int a, structure* ct, pfdatatable* data) {
	return POWER(data->eparam[a], ct->number_of_sequences)/ ct->number_of_sequences;
}

void readpfsave(const char *filename, structure *ct,
			 PFPRECISION *w5, PFPRECISION *w3,
			 DynProgArray<PFPRECISION> *v, DynProgArray<PFPRECISION> *w, DynProgArray<PFPRECISION> *wmb, DynProgArray<PFPRECISION> *wl, DynProgArray<PFPRECISION> *wlc, DynProgArray<PFPRECISION> *wmbl, DynProgArray<PFPRECISION> *wcoax,
			 forceclass *fce, PFPRECISION *scaling, bool *mod, bool *lfce, pfdatatable *data, datatable *data2) {
	 int i,j,k,l,m,n,o,p;
	ifstream sav(filename,ios::binary);
	
	short vers;

	

	//Associate the datatable with the structure
	ct->SetThermodynamicDataTable(data2);

	//read the save file

	//read the file version first
	read(&sav,&vers);

	//start with structure information
	int SequenceLength;
	read(&sav,&(SequenceLength));
	//ct->numofbases = SequenceLength;
	read(&sav,&(ct->intermolecular));
	read(&sav,scaling);

	data->scaling=*scaling;

	int constraint,constraint2,numberofconstraints;

	//Read information about the forced pairs
	read(&sav,&(numberofconstraints));
	for (i=0;i<numberofconstraints;i++) {
		read(&sav,&(constraint));
		read(&sav,&(constraint2));

		ct->AddPair(constraint,constraint2);
	}
	for (i=0;i<=ct->GetSequenceLength();i++) {
		read(&sav,&(ct->hnumber[i]));
		sav.read(&(ct->nucs[i]),1);

	}

	for (i=0;i<=2*ct->GetSequenceLength();i++) read(&sav,&(ct->numseq[i]));

	//Read information about nucleotides forced to be double stranded.
	read(&sav,&(numberofconstraints));
	for (i=0;i<numberofconstraints;i++) {
		read(&sav,&(constraint));

		ct->AddDouble(constraint);

	}
	if (ct->intermolecular) {
		for (i=0;i<3;i++) read(&sav,&(ct->inter[i]));

	}

	//Read information about nucleotides not allowed to pair
	read(&sav,&(numberofconstraints));
	for (i=0;i<numberofconstraints;i++) {
		read(&sav,&(constraint));

		ct->AddSingle(constraint);

	}

	//Read information about nucleotides that are accessible to chemical modification:
	read(&sav,&(numberofconstraints));
	for (i=0;i<numberofconstraints;i++) {
		read(&sav,&(constraint));

		ct->AddModified(constraint);

	}

	//Read information about Us in GU pairs:
	read(&sav,&(numberofconstraints));
	for (i=0;i<numberofconstraints;i++) {
		read(&sav,&(constraint));

		ct->AddGUPair(constraint);

	}

	string label;
	read(&sav,&(label));
	ct->SetSequenceLabel(label);

	read(&sav,&(ct->templated));
	if (ct->templated) {
		ct->allocatetem();
		for (i=0;i<=ct->GetSequenceLength();i++) {
			for (j=0;j<=i;j++) read(&sav,&(ct->tem[i][j]));

		}
	}

	read(&sav,&(ct->shaped));
	if (ct->shaped) {
		ct->SHAPE = new double [2*ct->GetSequenceLength()+1];
		for (i=0;i<=2*ct->GetSequenceLength();i++) read(&sav,&(ct->SHAPE[i]));
		ct->SHAPEss = new double [2 * ct->GetSequenceLength() + 1];
		for (i=0;i<=2*ct->GetSequenceLength();i++) read(&sav,&(ct->SHAPEss[i]));

	}

	//now read the array class data for v, w, and wmb:
	for (i=0;i<=ct->GetSequenceLength();i++) {
		read(&sav,&(w3[i]));
		read(&sav,&(w5[i]));
		for (j=0;j<=ct->GetSequenceLength();j++) {
			read(&sav,&(v->dg[i][j+i]));
			read(&sav,&(w->dg[i][j+i]));
			read(&sav,&(wmb->dg[i][j+i]));
			read(&sav,&(wmbl->dg[i][j+i]));
			read(&sav,&(wl->dg[i][j+i]));
			read(&sav,&(wlc->dg[i][j+i]));
			read(&sav,&(wcoax->dg[i][j+i]));
			readsinglechar(&sav,&(fce->dg[i][j]));
		}
	}

	read(&sav,&(w3[ct->GetSequenceLength()+1]));
	for (i=0;i<=2*ct->GetSequenceLength();i++) {
		read(&sav,&(lfce[i]));
		read(&sav,&(mod[i]));

	}

	//read the alphabet in data2
	read(&sav, &(data2->alphabet));
	read(&sav, &(data2->pairing));
	read(&sav, &(data2->not_pairing));
	read(&sav, &(data2->non_interacting));
	read(&sav, &(data2->linker));
	/*
	for (i=0; i<data2->alphabet.size(); i++){
		for (j=0; j<data2->alphabet[i].size(); j++){
			data2->alphabet_map[j<data2->alphabet[i][j]]=i;
			ct->data->alphabet_map[j<data2->alphabet[i][j]]=i;
		}
	}
//*/
	data2->LinkerInts.resize(data2->alphabet.size());
	//build the LinkerInts array:
	for (int li=0;li<data2->LinkerInts.size();++li) {
		data2->LinkerInts[li]=false;

	}
	
	for (int li=0;li<data2->linker.size();++li) {

		int n = data2->basetonum(data2->linker[li]);
		data2->LinkerInts[data2->basetonum(data2->linker[li])]=true;

	}

	//read the alphabet
	read(&sav, &(data->alphabet));
	read(&sav, &(data->pairing));

	if(data->alphabet.size()==0){
		cerr<<"WARNING: no alphabet was read. Is the specification file formatted correctly?\n";
	}

	

	data->allocate_data_tables(data->alphabet.size());

	//now read the thermodynamic data:
	read(&sav,&(data->pftemp));
	for (i=0;i<5;i++) read(&sav,&(data->poppen[i]));
	read(&sav,&(data->maxpen));
	for (i=0;i<11;i++) read(&sav,&(data->eparam[i]));
	for (i=0;i<31;i++) {
		read(&sav,&(data->inter[i]));
		read(&sav,&(data->bulge[i]));
		read(&sav,&(data->hairpin[i]));

	}
	for (i=0;i<data->alphabet.size();i++) {
		for (j=0;j<data->alphabet.size();j++) {
			for (k=0;k<data->alphabet.size();k++) {
				for (l=0;l<3;l++) {
					read(&sav,&(data->dangle[i][j][k][l]));
				}
				for (l=0;l<data->alphabet.size();l++) {
					read(&sav,&(data->stack[i][j][k][l]));
					read(&sav,&(data->tstkh[i][j][k][l]));
					read(&sav,&(data->tstki[i][j][k][l]));
					read(&sav,&(data->coax[i][j][k][l]));
					read(&sav,&(data->tstackcoax[i][j][k][l]));
					read(&sav,&(data->coaxstack[i][j][k][l]));
					read(&sav,&(data->tstack[i][j][k][l]));
					read(&sav,&(data->tstkm[i][j][k][l]));
					read(&sav,&(data->tstki23[i][j][k][l]));
					read(&sav,&(data->tstki1n[i][j][k][l]));
					for (m=0;m<data->alphabet.size();m++) {
						for (n=0;n<data->alphabet.size();n++) {
							read(&sav,&(data->iloop11[i][j][k][l][m][n]));
							for (o=0;o<data->alphabet.size();o++) {
								if (data->pairing[i][j]&& data->pairing[n][o]) read(&sav,&(data->iloop21[i][j][k][l][m][n][o]));
								for (p=0;p<data->alphabet.size();p++) {
									if (data->pairing[i][k]&& data->pairing[j][l])
										read(&sav,&(data->iloop22[i][j][k][l][m][n][o][p]));
								}
							}
						}
					}
				}
			}
		}
	}
	read(&sav,&(data->numoftloops));
	for (i=0;i<data->numoftloops;i++) {
		read(&sav,&(data->itloop[i]));
		read(&sav,&(data->tloop[i]));

	}
	read(&sav,&(data->numoftriloops));
	for (i=0;i<data->numoftriloops;i++) {
		read(&sav,&(data->itriloop[i]));
		read(&sav,&(data->triloop[i]));

	}
	read(&sav,&(data->numofhexaloops));
	for (i=0;i<data->numofhexaloops;i++) {
		read(&sav,&(data->ihexaloop[i]));
		read(&sav,&(data->hexaloop[i]));

	}
	read(&sav,&(data->auend));
	read(&sav,&(data->AUappliestoGU));
	read(&sav,&(data->gubonus));
	read(&sav,&(data->cint));
	read(&sav,&(data->cslope));
	read(&sav,&(data->c3));
	read(&sav,&(data->efn2a));
	read(&sav,&(data->efn2b));
	read(&sav,&(data->efn2c));
	read(&sav,&(data->init));
	read(&sav,&(data->mlasym));
	read(&sav,&(data->strain));
	read(&sav,&(data->prelog));
	read(&sav,&(data->singlecbulge));
	read(&sav,&(data->maxintloopsize));
	
	for(int i=0;i<=data->alphabet.size();i++){
		for (int j=0;j<=data->alphabet.size();j++){
			read(&sav,&(data->penalties[i][j]));
		}
	}

	sav.close();


}

//return the pairing probability of the i=j pair, where i<j.
double calculateprobability(int i, int j, DynProgArray<PFPRECISION> *v, PFPRECISION *w5, structure *ct, pfdatatable *data, bool *lfce, bool *mod, PFPRECISION scaling, forceclass *fce) {
	PFPRECISION interior, exterior;
	bool before,after;
	//int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
	//{0,1,0,1,0,0},{0,0,0,0,0,0}};
	bool adjacentgu;
	
	if (!mod[i]&&!mod[j]) {
		if (ct->constant==NULL) return TO_LINEAR(DIV(PROD(v->f(i,j), v->f(j,i+ct->GetSequenceLength())), PROD(w5[ct->GetSequenceLength()] POWSCALING2(2))));
		else {
			//constant is being used.
			if (ct->constant[j][i]<TO_XLOG(EPSILON)) return 0.0;
			return TO_LINEAR(DIV(PROD(v->f(i,j),v->f(j,i+ct->GetSequenceLength())),PROD(w5[ct->GetSequenceLength()] POWSCALING2(2),ct->constant[j][i])));
		}
	}
	else {
		if (!(fce->f(i,j)&SINGLE)) {
			before =0;
			if ((i>1&&j<(2*ct->GetSequenceLength())&&j!=ct->GetSequenceLength())) {
				if ((j>ct->GetSequenceLength()&&((i-j+ct->GetSequenceLength())>minloop+2))||j<ct->GetSequenceLength()) {
					before = data->pairing[ct->numseq[i-1]][ct->numseq[j+1]];
				}
			}

			//after = 0 if a stacked pair cannot form 3' to i
			if ((((j-i)>minloop+2)&&(j<=ct->GetSequenceLength())||(j>ct->GetSequenceLength()+1))&&(i!=ct->GetSequenceLength())) {
				after = data->pairing[ct->numseq[i+1]][ct->numseq[j-1]];

			}
			else after = 0;

			adjacentgu = false;
			//check for preceding or following GU or whether the pair itself is gu
			if (ct->numseq[i+1]==3&&ct->numseq[j-1]==4) adjacentgu = true;
			else if (ct->numseq[i+1]==4&&ct->numseq[j-1]==3) adjacentgu = true;
			else if (ct->numseq[i]==3&&ct->numseq[j]==4) adjacentgu = true;
			else if (ct->numseq[i]==4&&ct->numseq[j]==3) adjacentgu = true;
			else if (i-1>0&&j+1<=ct->GetSequenceLength()) {
				if (ct->numseq[i-1]==3&&ct->numseq[j+1]==4) adjacentgu = true;
				else if (ct->numseq[i-1]==4&&ct->numseq[j+1]==3) adjacentgu = true;

			}

			//if there are no stackable pairs to i.j then don't allow a pair i,j
			if ((before!=0)||(after!=0)) {
				if (i+1<j-1&&!adjacentgu) interior = PROD(erg1(i,j,i+1,j-1,ct,data),v->f(i+1,j-1));
				else interior = (PFPRECISION) ZERO;
				if (j+1<=ct->GetSequenceLength()&&!adjacentgu) exterior = PROD(erg1(j,i+ct->GetSequenceLength(),j+1,i+ct->GetSequenceLength()-1,ct,data), v->f(j+1,i+ct->GetSequenceLength()-1));
				else exterior = (PFPRECISION) ZERO;
				return TO_LINEAR(DIV(DIFF(PROD(SUM(v->f(i,j),interior),SUM(v->f(j,i+ct->GetSequenceLength()),exterior)),PROD(interior,exterior)),PROD(w5[ct->GetSequenceLength()] POWSCALING2(2))));
			}
			else return 0.0;

		}
		else return 0.0;

	}
}

//function to rescale all arrays when partition function calculation is headed out of bounds

void rescale(int currenth, structure *ct, pfdatatable *data, DynProgArray<PFPRECISION> *v, DynProgArray<PFPRECISION> *w, DynProgArray<PFPRECISION> *wl, DynProgArray<PFPRECISION> *wcoax,
			 DynProgArray<PFPRECISION> *wmb,DynProgArray<PFPRECISION> *wmbl, PFPRECISION *w5, PFPRECISION *w3, PFPRECISION **wca, 
			 PFPRECISION rescalefactor) {
//	int d,dp,ii,jj,hh,nucs,index,lowi,highi,number;
	int d, nucs, index, number;
	PFPRECISION multiplier;
	cout << "RESCALE factor: " << rescalefactor << endl;
	number=ct->GetSequenceLength();
	for (int h=0;h<=currenth;h++){
		d=(h<=(number-1))?h:(h-number+1);
		int start;
		int end;
		if (h<=(number-1)) {
			start = 1;
			end = number-h;
		}
		else {
			start = 2*number-h;
			end = number;
		}
		for (int locali=start;locali<=end;locali++){
			int localj=locali+d;

	#ifdef pfdebugmode
		ofstream ufout;
		ufout.open(pfdebugpath,ios::app);
		ufout<<"rescale factor = "<<rescalefactor<<"\n";
		ufout.close();
	#endif
		//rescale v,w,wl,wcoax,wmb,wmbl,wca
			nucs = localj-locali+1; //this work even jj> number
			multiplier = pow(rescalefactor,(PFPRECISION) nucs);

		#ifdef oldpfdebugmodenotuptodate
			//look for positions that will underflow

			if (multiplier<0) {
				//values will get smaller...
				if (v->f(ii,jj)>0) {
					if (log10(v->f(ii,jj))+log10(multiplier)<log10(PFMIN)) {
					ufout.open(pfdebugpath,ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" v= "<<v->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (w->f(ii,jj)>0) {
					if (log10(w->f(ii,jj))+log10(multiplier)<log10(PFMIN)) {
					ufout.open(pfdebugpath,ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" w= "<<w->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wl->f(ii,jj)>0) {
					if (log10(wl->f(ii,jj))+log10(multiplier)<log10(PFMIN)) {
					ufout.open(pfdebugpath,ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wl= "<<wl->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wcoax->f(ii,jj)>0) {
					if (log10(wcoax->f(ii,jj))+log10(multiplier)<log10(PFMIN)) {
					ufout.open(pfdebugpath,ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wcoax= "<<wcoax->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wmb->f(ii,jj)>0) {
					if (log10(wmb->f(ii,jj))+log10(multiplier)<log10(PFMIN)) {
					ufout.open(pfdebugpath,ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wmb= "<<wmb->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wmbl->f(ii,jj)>0) {
					if (log10(wmbl->f(ii,jj))+log10(multiplier)<log10(PFMIN)) {
					ufout.open(pfdebugpath,ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wmbl= "<<wmbl->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
			}

			if (multiplier>0) {
				//values will get smaller...
				if (v->f(ii,jj)>0) {
					if (log10(v->f(ii,jj))+log10(multiplier)>log10(PFMAX)) {
					ufout.open(pfdebugpath,ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" v= "<<v->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (w->f(ii,jj)>0) {
					if (log10(w->f(ii,jj))+log10(multiplier)>log10(PFMAX)) {
					ufout.open(pfdebugpath,ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" w= "<<w->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wl->f(ii,jj)>0) {
					if (log10(wl->f(ii,jj))+log10(multiplier)>log10(PFMAX)) {
					ufout.open(pfdebugpath,ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wl= "<<wl->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wcoax->f(ii,jj)>0) {
					if (log10(wcoax->f(ii,jj))+log10(multiplier)>log10(PFMAX)) {
					ufout.open(pfdebugpath,ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wcoax= "<<wcoax->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wmb->f(ii,jj)>0) {
					if (log10(wmb->f(ii,jj))+log10(multiplier)>log10(PFMAX)) {
					ufout.open(pfdebugpath,ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wmb= "<<wmb->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wmbl->f(ii,jj)>0) {
					if (log10(wmbl->f(ii,jj))+log10(multiplier)>log10(PFMAX)) {
					ufout.open(pfdebugpath,ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wmbl= "<<wmbl->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
			}

		#endif //debug
			//rescale v,w,wl,wcoax,wmb,wmbl,wca
			v->f(locali,localj)=v->f(locali,localj)*multiplier;
			w->f(locali,localj)=w->f(locali,localj)*multiplier;
			wl->f(locali,localj)=wl->f(locali,localj)*multiplier;
			wcoax->f(locali,localj)=wcoax->f(locali,localj)*multiplier;
			wmb->f(locali,localj)=wmb->f(locali,localj)*multiplier;
			wmbl->f(locali,localj)=wmbl->f(locali,localj)*multiplier;
			if (localj<=number) wca[locali][localj]=wca[locali][localj]*multiplier;
			if (locali==1 && localj<=number) {
				//rescale w5
				w5[localj]=w5[localj]*pow(rescalefactor,(PFPRECISION) localj);

				if (localj==number) {
					//rescale w3
					for (index=1;index<=number;index++) w3[index]=w3[index]*pow(rescalefactor,(PFPRECISION) (number-index+1));

				}
			}
		}//end for locali 
	}//end for h

	//rescale datatable
	data->rescaledatatable(rescalefactor);
}

void rescale(int currenth, structure *ct, pfdatatable *data, DynProgArray<PFPRECISION> *v, DynProgArray<PFPRECISION> *w, DynProgArray<PFPRECISION> *wl, DynProgArray<PFPRECISION> *wcoax,
			 DynProgArray<PFPRECISION> *wmb,DynProgArray<PFPRECISION> *wmbl, PFPRECISION *w5, PFPRECISION *w3, PFPRECISION **wca, 
			 PFPRECISION **curE, PFPRECISION **prevE, PFPRECISION rescalefactor) {
	//int d,dp,ii,jj,hh,nucs,index,lowi,highi,number;
	int d, dp, ii, nucs, index, number;
	double multiplier;
	// cout << "RESCALE #" << ++TotalRescaleCount << " factor: " << rescalefactor << endl;
	number=ct->GetSequenceLength();
	for (int h=0;h<=currenth;h++){
		d=(h<=(number-1))?h:(h-number+1);
		int start;
		int end;
		if (h<=(number-1)) {
			start = 1;
			end = number-h;
		}
		else {
			start = 2*number-h;
			end = number;
		}
		for (int locali=start;locali<=end;locali++){
			int localj=locali+d;

	#ifdef pfdebugmode
		ofstream ufout;
		ufout.open(pfdebugpath,ios::app);
		ufout<<"rescale factor = "<<rescalefactor<<"\n";
		ufout.close();
	#endif
		//rescale v,w,wl,wcoax,wmb,wmbl,wca
			nucs = localj-locali+1; //this work even jj> number
			multiplier = pow(rescalefactor,(PFPRECISION) nucs);

		#ifdef oldpfdebugmodenotuptodate
			//look for positions that will underflow

			if (multiplier<0) {
				//values will get smaller...
				if (v->f(ii,jj)>0) {
					if (log10(v->f(ii,jj))+log10(multiplier)<log10(PFMIN)) {
					ufout.open(pfdebugpath,ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" v= "<<v->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (w->f(ii,jj)>0) {
					if (log10(w->f(ii,jj))+log10(multiplier)<log10(PFMIN)) {
					ufout.open(pfdebugpath,ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" w= "<<w->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wl->f(ii,jj)>0) {
					if (log10(wl->f(ii,jj))+log10(multiplier)<log10(PFMIN)) {
					ufout.open(pfdebugpath,ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wl= "<<wl->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wcoax->f(ii,jj)>0) {
					if (log10(wcoax->f(ii,jj))+log10(multiplier)<log10(PFMIN)) {
					ufout.open(pfdebugpath,ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wcoax= "<<wcoax->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wmb->f(ii,jj)>0) {
					if (log10(wmb->f(ii,jj))+log10(multiplier)<log10(PFMIN)) {
					ufout.open(pfdebugpath,ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wmb= "<<wmb->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wmbl->f(ii,jj)>0) {
					if (log10(wmbl->f(ii,jj))+log10(multiplier)<log10(PFMIN)) {
					ufout.open(pfdebugpath,ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wmbl= "<<wmbl->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
			}

			if (multiplier>0) {
				//values will get smaller...
				if (v->f(ii,jj)>0) {
					if (log10(v->f(ii,jj))+log10(multiplier)>log10(PFMAX)) {
					ufout.open(pfdebugpath,ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" v= "<<v->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (w->f(ii,jj)>0) {
					if (log10(w->f(ii,jj))+log10(multiplier)>log10(PFMAX)) {
					ufout.open(pfdebugpath,ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" w= "<<w->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wl->f(ii,jj)>0) {
					if (log10(wl->f(ii,jj))+log10(multiplier)>log10(PFMAX)) {
					ufout.open(pfdebugpath,ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wl= "<<wl->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wcoax->f(ii,jj)>0) {
					if (log10(wcoax->f(ii,jj))+log10(multiplier)>log10(PFMAX)) {
					ufout.open(pfdebugpath,ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wcoax= "<<wcoax->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wmb->f(ii,jj)>0) {
					if (log10(wmb->f(ii,jj))+log10(multiplier)>log10(PFMAX)) {
					ufout.open(pfdebugpath,ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wmb= "<<wmb->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wmbl->f(ii,jj)>0) {
					if (log10(wmbl->f(ii,jj))+log10(multiplier)>log10(PFMAX)) {
					ufout.open(pfdebugpath,ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wmbl= "<<wmbl->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
			}

		#endif //debug
			//rescale v,w,wl,wcoax,wmb,wmbl,wca
			v->f(locali,localj)=v->f(locali,localj)*multiplier;
			w->f(locali,localj)=w->f(locali,localj)*multiplier;
			wl->f(locali,localj)=wl->f(locali,localj)*multiplier;
			wcoax->f(locali,localj)=wcoax->f(locali,localj)*multiplier;
			wmb->f(locali,localj)=wmb->f(locali,localj)*multiplier;
			wmbl->f(locali,localj)=wmbl->f(locali,localj)*multiplier;
			if (localj<=number) wca[locali][localj]=wca[locali][localj]*multiplier;
			if (locali==1 && localj<=number) {
				//rescale w5
				w5[localj]=w5[localj]*pow(rescalefactor,(PFPRECISION) localj);

				if (localj==number) {
					//rescale w3
					for (index=1;index<=number;index++) w3[index]=w3[index]*pow(rescalefactor,(PFPRECISION) (number-index+1));

				}
			}
		}//end for locali 
	}//end for h

	if (curE!=NULL) {
			//curE and preE are not used in SMP mode:
		//rescale curE and prevE
		for (ii=((currenth<=(number-2))?1:(2*number-currenth-1));ii<=((currenth<=(number-2))?(number-currenth):number);ii++)
		for (dp=1;dp<=d-1;dp++){
			if (ii<number) {
				curE[dp][ii] = curE[dp][ii]*pow(rescalefactor,(PFPRECISION)(dp+1));
				prevE[dp][ii+1] = prevE[dp][ii+1]*pow(rescalefactor,(PFPRECISION)(dp+1));
			}
		}
	}

	//rescale datatable
	data->rescaledatatable(rescalefactor);
}
//*/

//function to rescale all arrays when partition function calculation is headed out
	//of bounds when claculating w3
//void rescaleatw3(int ii, structure *ct, pfdatatable *data, DynProgArray<PFPRECISION> *v, DynProgArray<PFPRECISION> *w, DynProgArray<PFPRECISION> *wl, DynProgArray<PFPRECISION> *wcoax,
//				 DynProgArray<PFPRECISION> *wmb,DynProgArray<PFPRECISION> *wmbl, PFPRECISION *w5, PFPRECISION *w3, double rescalefactor) {
//}

//function to rescale all arrays when partition function calculation is headed out
	//of bounds when claculating w5
//void rescaleatw5(int jj,structure *ct, pfdatatable *data, DynProgArray<PFPRECISION> *v, DynProgArray<PFPRECISION> *w, DynProgArray<PFPRECISION> *wl, DynProgArray<PFPRECISION> *wcoax,
//				 DynProgArray<PFPRECISION> *wmb,DynProgArray<PFPRECISION> *wmbl, PFPRECISION *w5, PFPRECISION *w3, double rescalefactor) {
	//rescale the previously filled arrays
//	rescale(1, jj, ct, data, v, w, wl, wcoax, wmb, wmbl, w5, w3, rescalefactor);

//	w5[jj]=w5[jj]*pow(rescalefactor,(double) jj);

//}



//rescale the entries in datatable
void pfdatatable::rescaledatatable(const PFPRECISION &rescalefactor) {
	scaling=scaling*rescalefactor;
	int i,j,k,l,m,n,o,p;

	for (i=0;i<31;i++) {
		inter[i] = inter[i]*pow(rescalefactor,i+2);
		bulge[i] = bulge[i]*pow(rescalefactor,i+2);
		hairpin[i] = hairpin[i]*pow(rescalefactor,i+2);

	}
	const int sz = alphabet.size();
	for (i=0;i<sz;i++) {
		for (j=0;j<sz;j++) {
			for (k=0;k<sz;k++) {
				for (l=0;l<3;l++) {
					dangle[i][j][k][l] = dangle[i][j][k][l]*rescalefactor;
				}
				for (l=0;l<sz;l++) {
					stack[i][j][k][l]=stack[i][j][k][l]*pow(rescalefactor,2);

					tstackcoax[i][j][k][l]=tstackcoax[i][j][k][l]*pow(rescalefactor,2);

					tstack[i][j][k][l]=tstack[i][j][k][l]*pow(rescalefactor,2);
					tstkm[i][j][k][l]=tstkm[i][j][k][l]*pow(rescalefactor,2);

					for (m=0;m<sz;m++) {
						for (n=0;n<sz;n++) {
							iloop11[i][j][k][l][m][n]=iloop11[i][j][k][l][m][n]*pow(rescalefactor,4);
							for (o=0;o<sz;o++) {
								iloop21[i][j][k][l][m][n][o]=iloop21[i][j][k][l][m][n][o]*pow(rescalefactor,5);
								for (p=0;p<sz;p++) {
									iloop22[i][j][k][l][m][n][o][p]=iloop22[i][j][k][l][m][n][o][p]*pow(rescalefactor,6);
								}
							}
						}
					}
				}
			}
		}
	}

	for (i=0;i<numoftloops;i++) {
		tloop[i] = tloop[i]*pow(rescalefactor,6);

	}

	for (i=0;i<numoftriloops;i++) {
		triloop[i] = triloop[i]*pow(rescalefactor,5);
	}

	for (i=0;i<numofhexaloops;i++) {
		hexaloop[i] = hexaloop[i]*pow(rescalefactor,8);
	}
}



//thresh-structure builds a structure containing all base pairs above the probability thresh (expressed as a fraction from 0 to 1).
//Note that thresh must be 0.5 or larger for the resulting structure to be a valid secondary structure.
void thresh_structure(structure *ct, char *pfsfile, double thresh) {
	int i,j;
	short vers;
	PFPRECISION *w5, *w3, scaling;
	DynProgArray<PFPRECISION> *v, *w,*wmb,*wmbl,*wl,*wlc, *wcoax;
	forceclass *fce;
	bool *mod,*lfce;
	pfdatatable *data;
	datatable *data2;

	//allocate the ct file by reading the save file:
	ifstream sav(pfsfile,ios::binary);

	read(&sav,&(vers));//read the version of the save file
		//right now there is no infrastructure to indicate the wrong version is being read.
		//This should be changed in the future...

	int SequenceLength;
	read(&sav,&(SequenceLength));
	//ct->numofbases = SequenceLength;

	sav.close();

	//allocate everything

	//array = new PFPRECISION *[ct.GetSequenceLength()+1];
	//for (i=1;i<=ct.GetSequenceLength();i++) {
	//	array[i] = new PFPRECISION [i+1];
	//}
	//for (i=1;i<=ct.GetSequenceLength();i++) {
	//	for (j=1;j<=i;j++) {
	//		array[i][j]=0;
	//	}
	//}

	ct->allocate(SequenceLength);

	w = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
	v = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
	wmb = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
	wmbl = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
	wcoax = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
	wl = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
	wlc = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
	fce = new forceclass(ct->GetSequenceLength());

	w5 = new PFPRECISION [ct->GetSequenceLength()+1];
	w3 = new PFPRECISION [ct->GetSequenceLength()+2];

	lfce = new bool [2*ct->GetSequenceLength()+1];
    mod = new bool [2*ct->GetSequenceLength()+1];

	data = new pfdatatable();
	data2 = new datatable();

	//load all the data from the pfsavefile:
	readpfsave(pfsfile, ct, w5, w3,v, w, wmb,wl,wlc, wmbl, wcoax, fce,&scaling,mod,lfce,data,data2);

	//reset the base pairing info:
	//for (i=1;i<=ct->GetSequenceLength();i++) ct->basepr[1][i]=0;

	ct->AddStructure();
	//fill array with the values for the plot:
	for (i=1;i<ct->GetSequenceLength();i++) {
		for (j=i+1;j<=ct->GetSequenceLength();j++) {
			if(calculateprobability(i,j,v,w5,ct,data,lfce,mod,scaling,fce)>thresh) {
				ct->SetPair(i,j);

			}
		}
	}

	//now build the structure

	delete w;
	delete v;
	delete wmb;
	delete fce;
	delete[] w5;
	delete[] w3;
	//for (i=1;i<=ct.GetSequenceLength();i++) {
	//	delete[] array[i];
	//}
	//delete[] array;
	delete[] lfce;
	delete[] mod;
	delete data;
	delete data2;
	delete wmbl;
	delete wl;
	delete wcoax;

}

bool can_pair(int i, int j, structure* ct, pfdatatable* data) {
	return data->pairing[ct->numseq[i]][ct->numseq[j]];
}

bool multi_can_pair(int i, int j, structure* ct, pfdatatable* data) {
	return ct->pairing_matrix[i][j];
}

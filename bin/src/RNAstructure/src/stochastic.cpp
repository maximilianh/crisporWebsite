
#define ENABLE_DEBUG_LOGS
#include "debug_logging.h"
#include "stochastic.h"
#include "random.h"
#include "pfunction.h"
#include "stackclass.h"
#include <iostream>
#include <cstdlib>
#include <cstring>
using namespace std;

// A flag to force the calculation of all possible configurations
// Comment out to stop each loop when the roll value is exceeded
#define force_continue

//Define 1 for the purposes of taking an inverse
//static PFPRECISION ONE=1; //include code for extended double if neseded

//register a base pair between two nucleotides
inline void regbp(structure *ct, int structurenumber, short i, short j) {
	ct->SetPair(i,j,structurenumber);

}


int stochastictraceback(DynProgArray<PFPRECISION> *w,DynProgArray<PFPRECISION> *wmb,DynProgArray<PFPRECISION> *wmbl,DynProgArray<PFPRECISION> *wcoax,DynProgArray<PFPRECISION> *wl,DynProgArray<PFPRECISION> *wlc,DynProgArray<PFPRECISION> *v,	forceclass *fce, PFPRECISION *w3,PFPRECISION *w5,PFPRECISION scaling, bool *lfce, bool *mod, pfdatatable *data, int numberofstructures,
	structure *ct, int randomseed, ProgressHandler *progress) {

	PFPRECISION twoscaling; //scalinginv,
//	int number, d, ll;
	int number;
	// int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
	// {0,1,0,1,0,0},{0,0,0,0,0,0}};
	bool** inc;
	//scalinginv = ONE/data->scaling;
	twoscaling=data->scaling*data->scaling;
	int tracebackerror=0;

	inc = new bool *[data->alphabet.size()];
	for (int letters=0; letters<data->alphabet.size();++letters) {
		inc[letters]=new bool[data->alphabet.size()];
		for (int letterpairs=0;letterpairs<data->alphabet.size();++letterpairs) {
			inc[letters][letterpairs]=data->pairing[letters][letterpairs];
		}
	}
	
	SET_DEBUG_LEVEL(INFO);
	// SET_DEBUG_LEVEL(DEBUG)

	//allocation of space for the sampled structures
	//must be done OUTSIDE the main loop (which might
	//be run in parallel)
	for(int i = 1; i <= numberofstructures; i++){
		ct->AddStructure();
		ct->SetCtLabel(ct->GetSequenceLabel(), i);
	}

	PFPRECISION cumulative_max = TO_XLOG(1 + 1e-5);
	PFPRECISION cumulative_min = TO_XLOG(1 - 1e-5);
	PFPRECISION max_cumulative = 0;

	if (ct->constant!=NULL) cout << "compensating for constant" << endl;

	#ifdef SMP
	#pragma omp parallel for
	#endif

	for (number = 1; number <= numberofstructures; number++) {
		randomnumber rand;
		rand.seed(randomseed+number);
		PFPRECISION roll;
		PFPRECISION cumulative, denominator;
		PFPRECISION locale, accumulator, tmp_cumulative;


		stackclass stack;
		integersize dummy1;
		short dummy2;
		bool found;
//		short switchcase,i,j,k,ip,jp,d;
		short switchcase, i, j, k, ip, jp;

		if (progress!=NULL) {
			progress->update((int) (((float) 100*number)/((float) numberofstructures)));
		}
		//start by putting the whole fragment on the stack:
		// [Debug_mod]
		stack.push(1,ct->GetSequenceLength(),0,0,0);
		// stack.push(619, 622, 0,0,0);
		// SET_DEBUG_LEVEL(DEBUG4)

		while (stack.pull(&i,&j,&switchcase,&dummy1,&dummy2)) {

			roll=TO_XLOG(rand.roll());
			cumulative = ZERO;
			found = false;

			switch(switchcase) {
				case 0: //switchcase=0, dealing with w5 fragment
					// if (i==619 && j==622 && only_log_once){
					// 	SET_DEBUG_LEVEL(DEBUG4)
					// 	only_log_once = false;
					// }
					LOG_DEBUG("switchcase: " << switchcase << "\ti,j: " << i << ", " << j);

					//Try adding a nucleotide to existing w5 fragment:
					if (j==0) {
						found = true;

						#ifndef force_continue
							break;
						#endif
					}

					denominator = w5[j];

					LOG_DEBUG4("denominator: " << denominator);

					// Scale the roll probability
					roll = PROD(roll,denominator);

					if (!lfce[j]&&!found) {
						cumulative = PROD(w5[j-1] SCALING);
						if (cumulative > roll) {
							LOG_DEBUG2("Pushing to stack");
							stack.push(1,j-1,0,0,0);
							found=true;

							#ifndef force_continue
								break;
							#endif

						}
					}

					LOG_DEBUG3("[332V5Q]\tcumulative: " << cumulative);
					
					
					
					for (k=0;k<=j-4;k++) {
						// accumulator+=cumulative;
						// cumulative = 0;

						#ifdef SIMPLEMBLOOP
						// cumulative+=(w5[k]*erg4(j-1,k+1,j,1,ct,data,lfce[j])*v->f(k+1,j-1)*penalty(j-1,k+1,ct,data));
						cumulative=SUM(cumulative,PROD(w5[k], erg4(j-1,k+1,j,1,ct,data,lfce[j]), v->f(k+1,j-1), penalty(j-1,k+1,ct,data)));

						LOG_DEBUG4("[XXXXXX]\tcumulative: " << cumulative);

						if (cumulative>roll&&!found) {
							LOG_DEBUG2("Pushing to stack");
							stack.push(1,k,0,0,0);
							stack.push(k+1,j-1,1,0,0);
							found = true;

							#ifndef force_continue
								break;
							#endif
						}

						if ((mod[k+1]||mod[j-1])) if(inc[ct->numseq[k+2]][ct->numseq[j-2]]&&notgu(k+1,j-1,ct)&&!(fce->f(k+1,j-1)&SINGLE)) {

							// cumulative+=(w5[k]*erg4(j-1,k+1,j,1,ct,data,lfce[j])*v->f(k+2,j-2)
							// 	*penalty(j-1,k+1,ct,data)*erg1(k+1,j-1,k+2,j-2,ct,data));
							cumulative=SUM(sumulative,PROD(w5[k],erg4(j-1,k+1,j,1,ct,data,lfce[j]),v->f(k+2,j-2),
								penalty(j-1,k+1,ct,data),erg1(k+1,j-1,k+2,j-2,ct,data)));

							LOG_DEBUG4("[XXXXXX]\tcumulative: " << cumulative);

							if (cumulative>roll&&!found) {
								LOG_DEBUG2("Pushing to stack");
								stack.push(1,k,0,0,0);
								stack.push(k+2,j-2,1,0,0);
								regbp(ct,number,k+1,j-1);
								found=true;

								#ifndef force_continue
									break;
								#endif
							}
						}

						#else  //!SIMPLEMBLOOP
						// cumulative +=(w5[k]*v->f(k+1,j)*penalty(j,k+1,ct,data));
						cumulative =SUM(cumulative,PROD(w5[k], v->f(k+1,j), penalty(j,k+1,ct,data)));

						LOG_DEBUG4("[XO5121]\tcumulative: " << cumulative);

						if (cumulative > roll&&!found) {
							LOG_DEBUG2("Pushing to stack");
							stack.push(1,k,0,0,0);
							stack.push(k+1,j,1,0,0);
							found = true;

							#ifndef force_continue
								break;
							#endif
						}

						if ((mod[k+1]||mod[j])) if(inc[ct->numseq[k+2]][ct->numseq[j-1]]&&notgu(k+1,j,ct)&&!(fce->f(k+1,j)&SINGLE)) {
							// cumulative +=(w5[k]*v->f(k+2,j-1)*penalty(j,k+1,ct,data)
							// 	*erg1(k+1,j,k+2,j-1,ct,data));
							cumulative =SUM(cumulative, PROD(w5[k], v->f(k+2,j-1), penalty(j,k+1,ct,data), 
								erg1(k+1,j,k+2,j-1,ct,data)));

							LOG_DEBUG4("[XO5121.2]\tcumulative: " << cumulative);

							if (cumulative > roll&&!found) {
								LOG_DEBUG2("Pushing to stack");
								stack.push(1,k,0,0,0);
								stack.push(k+2,j-1,1,0,0);
								regbp(ct,number,k+1,j);
								found=true;

								#ifndef force_continue
									break;
								#endif
							}
						}

						// cumulative+=(w5[k]*erg4(j,k+2,k+1,2,ct,data,lfce[k+1])*v->f(k+2,j)*penalty(j,k+2,ct,data));
						cumulative=SUM(cumulative, PROD(w5[k], erg4(j,k+2,k+1,2,ct,data,lfce[k+1]), v->f(k+2,j), penalty(j,k+2,ct,data)));

						LOG_DEBUG4("[0A1GDP]\tcumulative: " << cumulative);

						if (cumulative > roll&&!found) {
							LOG_DEBUG2("Pushing to stack");
							stack.push(1,k,0,0,0);
							stack.push(k+2,j,1,0,0);
							found=true;

							#ifndef force_continue
								break;
							#endif
						}

						if ((mod[k+2]||mod[j])) if(inc[ct->numseq[k+3]][ct->numseq[j-1]]&&notgu(k+1,j-1,ct)&&!(fce->f(k+1,j-1)&SINGLE)) {
							// cumulative += (w5[k]*erg4(j,k+2,k+1,2,ct,data,lfce[j])*v->f(k+3,j-1)
							// 	*penalty(j,k+2,ct,data)*erg1(k+2,j,k+3,j-1,ct,data));
							cumulative = SUM(cumulative, PROD(w5[k], erg4(j,k+2,k+1,2,ct,data,lfce[j]), v->f(k+3,j-1),
								penalty(j,k+2,ct,data), erg1(k+2,j,k+3,j-1,ct,data)));

							LOG_DEBUG4("[0A1GDP.2]\tcumulative: " << cumulative);

							if (cumulative>roll&&!found) {
								LOG_DEBUG2("Pushing to stack");
								stack.push(1,k,0,0,0);
								stack.push(k+2,j,1,0,0);
								regbp(ct,number,k+3,j-1);
								found=true;

								#ifndef force_continue
									break;
								#endif
							}
						}

						// cumulative+=(w5[k]*erg4(j-1,k+1,j,1,ct,data,lfce[j])*v->f(k+1,j-1)*penalty(j-1,k+1,ct,data));
						cumulative = SUM(cumulative, PROD(w5[k], erg4(j-1,k+1,j,1,ct,data,lfce[j]), v->f(k+1,j-1), penalty(j-1,k+1,ct,data)));

						LOG_DEBUG4("[X0HOHA]\tcumulative: " << cumulative);

						if (cumulative>roll&&!found) {
							LOG_DEBUG2("Pushing to stack");
							stack.push(1,k,0,0,0);
							stack.push(k+1,j-1,1,0,0);
							found = true;

							#ifndef force_continue
								break;
							#endif
						}

						if ((mod[k+1]||mod[j-1])) if(inc[ct->numseq[k+2]][ct->numseq[j-2]]&&notgu(k+1,j-1,ct)&&!(fce->f(k+1,j-1)&SINGLE)) {
							// cumulative+=(w5[k]*erg4(j-1,k+1,j,1,ct,data,lfce[j])*v->f(k+2,j-2)
							// 	*penalty(j-1,k+1,ct,data)*erg1(k+1,j-1,k+2,j-2,ct,data));
							cumulative = SUM(cumulative, PROD(w5[k], erg4(j-1,k+1,j,1,ct,data,lfce[j]), v->f(k+2,j-2),
								penalty(j-1,k+1,ct,data), erg1(k+1,j-1,k+2,j-2,ct,data)));

							LOG_DEBUG4("[X0HOHA.2]\tcumulative: " << cumulative);

							if (cumulative>roll&&!found) {
								LOG_DEBUG2("Pushing to stack");
								stack.push(1,k,0,0,0);
								stack.push(k+2,j-2,1,0,0);
								regbp(ct,number,k+1,j-1);
								found=true;

								#ifndef force_continue
									break;
								#endif
							}
						}

						// cumulative+=(w5[k]*data->tstack[ct->numseq[j-1]][ct->numseq[k+2]][ct->numseq[j]][ct->numseq[k+1]]
						// 			*pfchecknp(lfce[j],lfce[k+1]) * v->f(k+2,j-1)*penalty(j-1,k+2,ct,data));
						cumulative = SUM(cumulative, PROD(w5[k], data->tstack[ct->numseq[j-1]][ct->numseq[k+2]][ct->numseq[j]][ct->numseq[k+1]],
									pfchecknp(lfce[j],lfce[k+1]), v->f(k+2,j-1), penalty(j-1,k+2,ct,data)));

						LOG_DEBUG4("[NXM9K1]\tcumulative: " << accumulator+cumulative);

						if (cumulative>roll&&!found) {
							LOG_DEBUG2("Pushing to stack");
							stack.push(1,k,0,0,0);
							stack.push(k+2,j-1,1,0,0);
							found=true;

							#ifndef force_continue
								break;
							#endif
						}

						if ((mod[k+2]||mod[j-1])) if(inc[ct->numseq[k+3]][ct->numseq[j-2]]&&notgu(k+2,j-1,ct)&&!(fce->f(k+2,j-1)&SINGLE)) {
							// cumulative+=(w5[k]*data->tstack[ct->numseq[j-1]][ct->numseq[k+2]][ct->numseq[j]][ct->numseq[k+1]]
							// 		*pfchecknp(lfce[j],lfce[k+1]) * v->f(k+3,j-2)*penalty(j-1,k+2,ct,data)*erg1(k+2,j-1,k+3,j-2,ct,data));
							cumulative = SUM(cumulative, PROD(w5[k], data->tstack[ct->numseq[j-1]][ct->numseq[k+2]][ct->numseq[j]][ct->numseq[k+1]],
									pfchecknp(lfce[j],lfce[k+1]), v->f(k+3,j-2), penalty(j-1,k+2,ct,data), erg1(k+2,j-1,k+3,j-2,ct,data)));

							LOG_DEBUG4("[NXM9K1.2]\tcumulative: " << cumulative);

							if (cumulative>roll&&!found) {
								LOG_DEBUG2("Pushing to stack");
								stack.push(1,k,0,0,0);
								stack.push(k+3,j-2,1,0,0);
								regbp(ct,number,k+2,j-1);
								found=true;

								#ifndef force_continue
									break;
								#endif
							}
						}

                        // The next section is equivalent to localrarray+=w5[k]*wca[k+1][j] in partition
						
                         //recheck all the coaxial stacking possibilities:
						i = k+1;
                        if (((j-i)>(2*minloop+2))) {

							locale=ZERO;
							accumulator = ZERO;

							

							for (ip=i+minloop+1;ip<j-minloop-1;ip++) {

								//first consider flush stacking
								accumulator = SUM(accumulator, PROD(w5[k], v->f(i,ip), v->f(ip+1,j), penalty(i,ip,ct,data),
									penalty(ip+1,j,ct,data), ergcoaxflushbases(i,ip,ip+1,j,ct,data)));

								LOG_TRACE("[TYG2WQ.1]\tcumulative: " << SUM(accumulator, cumulative));

								if (SUM(accumulator, cumulative, locale) > roll&&!found) {
									LOG_DEBUG2("Pushing to stack");
									stack.push(i,ip,1,0,0);
									stack.push(ip+1,j,1,0,0);
									found=true;
									stack.push(1,k,0,0,0);

									#ifndef force_continue
										break;
									#endif
								}

								if ((mod[i]||mod[ip]||mod[ip+1]||mod[j])) {

									if ((mod[i]||mod[ip])&&(mod[ip+1]||mod[j])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]
										&&inc[ct->numseq[ip+2]][ct->numseq[j-1]]&&notgu(i,ip,ct)&&notgu(ip+1,j,ct)
										&&!(fce->f(ip+1,j)&SINGLE)&&!(fce->f(i,ip)&SINGLE)) {

										// accumulator+=w5[k]*v->f(i+1,ip-1)*v->f(ip+2,j-1)*penalty(i,ip,ct,data)
										// 	*penalty(ip+1,j,ct,data)*ergcoaxflushbases(i,ip,ip+1,j,ct,data)
										// 	*erg1(i,ip,i+1,ip-1,ct,data)*erg1(ip+1,j,ip+2,j-1,ct,data);
										accumulator= SUM(accumulator, PROD(w5[k], v->f(i+1,ip-1), v->f(ip+2,j-1), penalty(i,ip,ct,data),
											penalty(ip+1,j,ct,data), ergcoaxflushbases(i,ip,ip+1,j,ct,data),
											erg1(i,ip,i+1,ip-1,ct,data), erg1(ip+1,j,ip+2,j-1,ct,data)));

										LOG_TRACE("[XXXXXX]\tcumulative: " << SUM(accumulator, cumulative));

										if (!found&&SUM(accumulator, cumulative, locale)>roll) {
											LOG_DEBUG2("Pushing to stack");
											stack.push(i+1,ip-1,1,0,0);
											stack.push(ip+2,j-1,1,0,0);
											found=true;
											regbp(ct,number,i,ip);
											regbp(ct,number,ip+1,j);
											stack.push(1,k,0,0,0);

											#ifndef force_continue
												break;
											#endif
										}
									}

									if ((mod[i]||mod[ip])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]&&notgu(i,ip,ct)&&!(fce->f(i,ip)&SINGLE)) {
										// accumulator+=w5[k]*v->f(i+1,ip-1)*v->f(ip+1,j)*penalty(i,ip,ct,data)
										// 	*penalty(ip+1,j,ct,data)*ergcoaxflushbases(i,ip,ip+1,j,ct,data)
										// 	*erg1(i,ip,i+1,ip-1,ct,data);
										accumulator = SUM(accumulator, PROD(w5[k], v->f(i+1,ip-1), v->f(ip+1,j), penalty(i,ip,ct,data),
											penalty(ip+1,j,ct,data), ergcoaxflushbases(i,ip,ip+1,j,ct,data), erg1(i,ip,i+1,ip-1,ct,data)));

										LOG_TRACE("[XXXXXX]\tcumulative: " << SUM(accumulator, cumulative));

										if (!found&&SUM(accumulator, cumulative, locale)>roll) {
											LOG_DEBUG2("Pushing to stack");
											stack.push(i+1,ip-1,1,0,0);
											stack.push(ip+1,j,1,0,0);
											found=true;
											regbp(ct,number,i,ip);
											stack.push(1,k,0,0,0);

											#ifndef force_continue
												break;
											#endif
										}
									}

									if ((mod[ip+1]||mod[j])&&inc[ct->numseq[ip+2]][ct->numseq[j-1]]&&notgu(ip+1,j,ct)&&!(fce->f(ip+1,j)&SINGLE)) {
										// accumulator+=w5[k]*v->f(i,ip)*v->f(ip+2,j-1)*penalty(i,ip,ct,data)
										// 	*penalty(ip+1,j,ct,data)*ergcoaxflushbases(i,ip,ip+1,j,ct,data)
										// 	*erg1(ip+1,j,ip+2,j-1,ct,data);
										accumulator = SUM(accumulator, PROD(w5[k], v->f(i,ip), v->f(ip+2,j-1), penalty(i,ip,ct,data),
											penalty(ip+1,j,ct,data), ergcoaxflushbases(i,ip,ip+1,j,ct,data),
											erg1(ip+1,j,ip+2,j-1,ct,data)));

										LOG_TRACE("[XXXXXX]\tcumulative: " << SUM(accumulator, cumulative));

										if (!found&&SUM(accumulator, cumulative, locale)>roll) {
											LOG_DEBUG2("Pushing to stack");
											stack.push(i,ip,1,0,0);
											stack.push(ip+2,j-1,1,0,0);
											found=true;
											regbp(ct,number,ip+1,j);
											stack.push(1,k,0,0,0);

											#ifndef force_continue
												break;
											#endif
										}
									}
								} // if ((mod[i]||mod[ip]||mod[ip+1]||mod[j]))

								if (!lfce[ip+1]&&!lfce[j]) {
									//now consider an intervening mismatch
									// locale+=w5[k]*v->f(i,ip)*v->f(ip+2,j-1)*penalty(i,ip,ct,data)
									// 	*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data);
									locale=SUM(locale, PROD(w5[k], v->f(i,ip), v->f(ip+2,j-1), penalty(i,ip,ct,data),
										penalty(ip+2,j-1,ct,data), ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)));

									LOG_TRACE("[XXXXXX]\tcumulative: " << SUM(accumulator, cumulative));

									if (SUM(accumulator, cumulative, locale)>roll&&!found) {
										LOG_DEBUG2("Pushing to stack");
										stack.push(i,ip,1,0,0);
										stack.push(ip+2,j-1,1,0,0);
										found=true;
										stack.push(1,k,0,0,0);

										#ifndef force_continue
											break;
										#endif
									}

									if (mod[i]||mod[ip]||mod[ip+2]||mod[j-1]) {
										if ((mod[i]||mod[ip])&&(mod[ip+2]||mod[j-1])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]
											&&inc[ct->numseq[ip+3]][ct->numseq[j-2]]&&notgu(i,ip,ct)&&notgu(ip+2,j-1,ct)
												&&!(fce->f(i,ip)&SINGLE)&&!(fce->f(ip+2,j-1)&SINGLE)) {

											// locale+=w5[k]*v->f(i+1,ip-1)*v->f(ip+3,j-2)*penalty(i,ip,ct,data)
											// 	*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)
											// 	*erg1(i,ip,i+1,ip-1,ct,data)*erg1(ip+2,j-1,ip+3,j-2,ct,data);
											locale = SUM(locale, PROD(w5[k], v->f(i+1,ip-1), v->f(ip+3,j-2), penalty(i,ip,ct,data),
												penalty(ip+2,j-1,ct,data), ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data),
												erg1(i,ip,i+1,ip-1,ct,data), erg1(ip+2,j-1,ip+3,j-2,ct,data)));

											LOG_TRACE("[XXXXXX]\tcumulative: " << SUM(accumulator, cumulative));

											if (!found&&SUM(accumulator, cumulative, locale)>roll) {
												LOG_DEBUG2("Pushing to stack");
												stack.push(i+1,ip-1,1,0,0);
												stack.push(ip+3,j-2,1,0,0);
												regbp(ct,number,i,ip);
												regbp(ct,number,ip+2,j-1);
												found=true;
												stack.push(1,k,0,0,0);

												#ifndef force_continue
													break;
												#endif
											}
										}

										if ((mod[i]||mod[ip])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]&&notgu(i,ip,ct)&&!(fce->f(i,ip)&SINGLE)) {

											// locale+=w5[k]*v->f(i+1,ip-1)*v->f(ip+2,j-1)*penalty(i,ip,ct,data)
											// 	*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)
											// 	*erg1(i,ip,i+1,ip-1,ct,data);
											locale = SUM(locale, PROD(w5[k], v->f(i+1,ip-1), v->f(ip+2,j-1), penalty(i,ip,ct,data),
												penalty(ip+2,j-1,ct,data), ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data),
												erg1(i,ip,i+1,ip-1,ct,data)));

											LOG_TRACE("[XXXXXX]\tcumulative: " << SUM(accumulator, cumulative));

											if(!found&&SUM(accumulator, cumulative, locale)>roll) {
												LOG_DEBUG2("Pushing to stack");
												stack.push(i+1,ip-1,1,0,0);
												stack.push(ip+2,j-1,1,0,0);
												regbp(ct,number,i,ip);
												found=true;
												stack.push(1,k,0,0,0);

												#ifndef force_continue
													break;
												#endif
											}
										}

										if ((mod[ip+2]||mod[j-1])&&inc[ct->numseq[ip+3]][ct->numseq[j-2]]&&notgu(ip+2,j-1,ct)&&!(fce->f(ip+2,j-1)&SINGLE)) {

											// locale+=w5[k]*v->f(i,ip)*v->f(ip+3,j-2)*penalty(i,ip,ct,data)
											// 	*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)
											// 	*erg1(ip+2,j-1,ip+3,j-2,ct,data);
											locale = SUM(locale, PROD(w5[k], v->f(i,ip), v->f(ip+3,j-2), penalty(i,ip,ct,data),
												penalty(ip+2,j-1,ct,data), ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data),
												erg1(ip+2,j-1,ip+3,j-2,ct,data)));

											LOG_TRACE("[XXXXXX]\tcumulative: " << SUM(accumulator, cumulative));

											if (!found&&SUM(accumulator, cumulative, locale)>roll) {
												LOG_DEBUG2("Pushing to stack");
												stack.push(i,ip,1,0,0);
												stack.push(ip+3,j-2,1,0,0);
												regbp(ct,number,ip+2,j-1);
												found=true;
												stack.push(1,k,0,0,0);

												#ifndef force_continue
													break;
												#endif
											}
										}
									} // if (mod[i]||mod[ip]||mod[ip+2]||mod[j-1])
								} // if (!lfce[ip+1]&&!lfce[j])

								if(!lfce[i]&&!lfce[ip+1]) {
									// locale+=w5[k]*v->f(i+1,ip)*v->f(ip+2,j)*penalty(i+1,ip,ct,data)
									// 	*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data);
									locale = SUM(locale, PROD(w5[k], v->f(i+1,ip), v->f(ip+2,j), penalty(i+1,ip,ct,data),
										penalty(ip+2,j,ct,data), ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)));

									LOG_TRACE("[XXXXXX]\tcumulative: " << SUM(accumulator, cumulative));

									if (!found&&SUM(accumulator, cumulative, locale)>roll) {
										LOG_DEBUG2("Pushing to stack");
										stack.push(i+1,ip,1,0,0);
										stack.push(ip+2,j,1,0,0);
										stack.push(1,k,0,0,0);
										found=true;

										#ifndef force_continue
											break;
										#endif
									}

									if (mod[i+1]||mod[ip]||mod[ip+2]||mod[j]) {
										if ((mod[i+1]||mod[ip])&&(mod[ip+2]||mod[j])&&inc[ct->numseq[i+2]][ct->numseq[ip-1]]
											&&inc[ct->numseq[ip+3]][ct->numseq[j-1]]&&notgu(i+1,ip,ct)&&notgu(ip+2,j,ct)
											&&!(fce->f(i+1,ip)&SINGLE)&&!(fce->f(ip+2,j)&SINGLE)	) {

											// locale+=w5[k]*v->f(i+2,ip-1)*v->f(ip+3,j-1)*penalty(i+1,ip,ct,data)
											// 	*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)
											// 	*erg1(i+1,ip,i+2,ip-1,ct,data)*erg1(ip+2,j,ip+3,j-1,ct,data);
											locale = SUM(locale, PROD(w5[k], v->f(i+2,ip-1), v->f(ip+3,j-1), penalty(i+1,ip,ct,data),
												penalty(ip+2,j,ct,data), ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data),
												erg1(i+1,ip,i+2,ip-1,ct,data), erg1(ip+2,j,ip+3,j-1,ct,data)));

											LOG_TRACE("[XXXXXX]\tcumulative: " << SUM(accumulator, cumulative));

											if (!found&&SUM(accumulator, cumulative, locale)>roll) {
												LOG_DEBUG2("Pushing to stack");
												stack.push(i+2,ip-1,1,0,0);
												stack.push(ip+3,j-1,1,0,0);
												regbp(ct,number,i+1,ip);
												regbp(ct,number,ip+2,j);
												found=true;
												stack.push(1,k,0,0,0);

												#ifndef force_continue
													break;
												#endif
											}
										}

										if ((mod[i+1]||mod[ip])&&inc[ct->numseq[i+2]][ct->numseq[ip-1]]&&notgu(i+1,ip,ct)&&!(fce->f(i+1,ip)&SINGLE)) {

											// locale+=w5[k]*v->f(i+2,ip-1)*v->f(ip+2,j)*penalty(i+1,ip,ct,data)
											// 	*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)
											// 	*erg1(i+1,ip,i+2,ip-1,ct,data);
											locale = SUM(locale, PROD(w5[k], v->f(i+2,ip-1), v->f(ip+2,j), penalty(i+1,ip,ct,data),
												penalty(ip+2,j,ct,data), ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data),
												erg1(i+1,ip,i+2,ip-1,ct,data)));

											LOG_TRACE("[XXXXXX]\tcumulative: " << SUM(accumulator, cumulative));

											if (!found&&SUM(accumulator, cumulative, locale)>roll) {
												LOG_DEBUG2("Pushing to stack");
												stack.push(i+2,ip-1,1,0,0);
												stack.push(ip+2,j,1,0,0);
												regbp(ct,number,i+1,ip);
												stack.push(1,k,0,0,0);

												found=true;

												#ifndef force_continue
													break;
												#endif
											}
										}

										if ((mod[ip+2]||mod[j])&&inc[ct->numseq[ip+3]][ct->numseq[j-1]]&&notgu(ip+2,j,ct)&&!(fce->f(ip+2,j)&SINGLE)) {
											// locale+=w5[k]*v->f(i+1,ip)*v->f(ip+3,j-1)*penalty(i+1,ip,ct,data)
											// 	*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)
											// 	*erg1(ip+2,j,ip+3,j-1,ct,data);
											locale = SUM(locale, PROD(w5[k], v->f(i+1,ip), v->f(ip+3,j-1), penalty(i+1,ip,ct,data),
												penalty(ip+2,j,ct,data), ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data),
												erg1(ip+2,j,ip+3,j-1,ct,data)));

											LOG_TRACE("[XXXXXX]\tcumulative: " << SUM(accumulator, cumulative));

											if (!found&&SUM(accumulator, cumulative, locale)>roll) {
												LOG_DEBUG2("Pushing to stack");
												stack.push(i+1,ip,1,0,0);
												stack.push(ip+3,j-1,1,0,0);
												stack.push(1,k,0,0,0);

												regbp(ct,number,ip+2,j);
												found=true;

												#ifndef force_continue
													break;
												#endif
											}
										}
									} // if (mod[i+1]||mod[ip]||mod[ip+2]||mod[j])
								} // if(!lfce[i]&&!lfce[ip+1])
							} // for (ip=i+minloop+1;ip<j-minloop-1;ip++)
							cumulative=SUM(cumulative, accumulator, locale);

							LOG_DEBUG4("[D753PB]\tcumulative: " << cumulative);

							#endif //SIMPLEMBLOOP
						} // if (((j-i)>(2*minloop+2)))
                    } // for (k=0;k<=j-4;k++)

					LOG_DEBUG4("[0NQO6X]\tcumulative: " << cumulative);

					if (!found) {
						cout << "Traceback error at w5\n";
						LOG_ERROR("Traceback errors at w5\t" << "i:  " << i << "\tj:  " << j << "\tCase:  " << switchcase << "\tCumulative:  " << DIV(cumulative,denominator));
						tracebackerror=14;
					}

					break;


				case 1: //switchcase=1, dealing with a v fragment

					LOG_DEBUG("switchcase: " << switchcase << "\ti,j: " << i << ", " << j);
					regbp(ct,number,i,j);

					//check to see if constant is used.
					//If it is, then v->f(i,j) was multiplied by ct->constant[j][i] and
					//this influence needs to be removed.
					if (ct->constant!=NULL) denominator = DIV(v->f(i,j), ct->constant[j][i]);
					else denominator = v->f(i,j);
					LOG_DEBUG2("denominator: " << denominator);

					// Scale roll probability
					roll = PROD(roll, denominator);
					LOG_DEBUG2("[Case 1]\troll: " << roll);

					//try closing a hairpin
					//[X5Y34F]
					cumulative = SUM(cumulative, erg3(i,j,ct,data,fce->f(i,j)));

					LOG_DEBUG3("[X5Y34F]\tcumulative: " << cumulative);

					if (cumulative>roll) {
						LOG_DEBUG2("Terminating branch");
						found = true; //nothing to put on the stack

						#ifndef force_continue
							break;
						#endif
					}

					//try stacking on a previous pair
					//[NDYIAA]
					if (!mod[i]&&!mod[j]) {
						// cumulative+=erg1(i,j,i+1,j-1,ct,data)*v->f(i+1,j-1);
						cumulative = SUM(cumulative, PROD(erg1(i,j,i+1,j-1,ct,data), v->f(i+1,j-1)));


						LOG_DEBUG3("[NDYIAA]\tcumulative: " << cumulative);

						if (!found&&cumulative>roll) {
							LOG_DEBUG2("Pushing to stack");
							stack.push(i+1,j-1,1,0,0);
							found=true;

							#ifndef force_continue
								break;
							#endif
						}
					}
					else {
						if ((ct->numseq[i]==3&&ct->numseq[j]==4)||(ct->numseq[i]==4&&ct->numseq[j]==3)) {
							// cumulative+=erg1(i,j,i+1,j-1,ct,data)*v->f(i+1,j-1);
							cumulative = SUM(cumulative, PROD(erg1(i,j,i+1,j-1,ct,data), v->f(i+1,j-1)));

							LOG_DEBUG3("[NDYIAA.2]\tcumulative: " << cumulative);

							if (cumulative>roll&&!found) {
								LOG_DEBUG2("Pushing to stack");
								stack.push(i+1,j-1,1,0,0);
								found=true;

								#ifndef force_continue
									break;
								#endif
							}
						}
						else if ((ct->numseq[i+1]==3&&ct->numseq[j-1]==4)||(ct->numseq[i+1]==4&&ct->numseq[j-1]==3)) {

							// cumulative+=erg1(i,j,i+1,j-1,ct,data)*v->f(i+1,j-1);
							cumulative = SUM(cumulative, PROD(erg1(i,j,i+1,j-1,ct,data), v->f(i+1,j-1)));

							LOG_DEBUG3("[NDYIAA.3]\tcumulative: " << cumulative);

							if (cumulative>roll&&!found) {
								LOG_DEBUG2("Pushing to stack");
								stack.push(i+1,j-1,1,0,0);
								found=true;

								#ifndef force_continue
									break;
								#endif
							}
						}
						else if (i-1>0) {
							if ((ct->numseq[i-1]==3&&ct->numseq[j+1]==4)||(ct->numseq[i-1]==4&&ct->numseq[j+1]==3)) {

								// cumulative+=erg1(i,j,i+1,j-1,ct,data)*v->f(i+1,j-1);
								cumulative = SUM(cumulative, PROD(erg1(i,j,i+1,j-1,ct,data), v->f(i+1,j-1)));

								LOG_DEBUG3("[NDYIAA.4]\tcumulative: " << cumulative);

								if (cumulative>roll&&!found) {
									LOG_DEBUG2("Pushing to stack");
									stack.push(i+1,j-1,1,0,0);
									found=true;

									#ifndef force_continue
										break;
									#endif
								}
							}
						}
					}
					// Check possible internal loops
					// [FQGH8G]
					/*
					// This follows O(N^3) code
					d = j-i;
					if ((d-1)>=(minloop+3)) { //[NIK7BX] If statement not in stochastic
						for (int dp=d-3;dp>=(minloop+1);dp--) {
							ll=d-dp-2;

							//calculate every ip,jp when ll <=5: 0x1,0x2,0x3,0x4,0x5,1x1,1x2,1x3,1x4,2x2,2x3
							if(ll>=1&&ll<=5) {
								for (ip=i+1;ip<=j-1-dp;ip++){
	  								jp=ip+dp;
									if (inc[ct->numseq[ip]][ct->numseq[jp]])
	  								if ( (ip<=number&&jp>number) || (ip<=number&&j<=number) ) {
   										cumulative=SUM(cumulative, PROD(erg2(i,j,ip,jp,ct,data,fce->f(i,ip), fce->f(jp,j)), v->f(ip,jp)));
										LOG_TRACE("[FQGH8G]\tcumulative_(v(i,j)): " << cumulative);
										
										//locali or localj is modified
										if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]&&!(fce->f(ip,jp)&SINGLE)){
		  									cumulative+=SUM(cumultive, PROD(erg2(i,j,ip,jp,ct,data,fce->f(i,ip),fce->f(jp,j)),
                  								v->f(ip+1,jp-1), erg1(ip,jp,ip+1,jp-1,ct,data)));
											LOG_TRACE("[FQGH8G.2]\tcumulative_(v(i,j)): " << cumulative);
										}

      								}
								}
							}
							//when size >=6 and <=30;

							else if (ll>=6&&ll<=maxinter)
							{
								//interior energy prefilled (after sub2:) + exterior engergy (the function for erg2in and erg2ex is stored in rna_library_inter.cpp
								cumulative+=curE[dp][i] * erg2ex(i,j,ll,ct,data);
								LOG_TRACE("[XXXXXX]\tcumulative_(v(i,j)): " << cumulative);

							//considering loop 1x(n-1) (bl=1)and 0 x n(bulge)(bl=0) as stacking bonus on 1x(n-1) and bulge is not allowed
								for (int bl=0;bl<=1;bl++) {
									ip=i+1+bl;
									jp=ip+dp;
									
									if (inc[ct->numseq[ip]][ct->numseq[jp]])
									if (abs(ip-i+jp-j)<=maxasym) {
		  								cumulative+=erg2(i,j,ip,jp,ct,data,fce->f(i,ip),fce->f(jp,j)) * v->f(ip,jp);
										LOG_TRACE("[XXXXXX]\tcumulative_(v(i,j)): " << cumulative);

										//locali or localj is modified
		   								if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]&&!(fce->f(ip,jp)&SINGLE)){
											cumulative+=erg2(i,j,ip,jp,ct,data,fce->f(i,ip),fce->f(jp,j)) *
                  								v->f(ip+1,jp-1) * erg1(ip,jp,ip+1,jp-1,ct,data);
											LOG_TRACE("[XXXXXX]\tcumulative_(v(i,j)): " << cumulative);
										}
									}
								//interior energy prefilled (after sub2:) + exterior engergy (the function for erg2in and erg2ex is stored in rna_library_inter.cpp
//								cumulative+=(curE[dp][locali]) *erg2ex(locali,localj,ll,ct,data);

							//considering loop 1x(n-1) (bl=1)and 0 x n(bulge)(bl=0) as stacking bonus on 1x(n-1) and bulge is not allowed
									jp=j-1-bl;
									ip=jp-dp;
									
									
									if (inc[ct->numseq[ip]][ct->numseq[jp]])
									if (abs(ip-i+jp-j)<=maxasym) {
										cumulative+=erg2(i,j,ip,jp,ct,data,fce->f(i,ip),fce->f(jp,j)) * v->f(ip,jp);
										LOG_TRACE("[XXXXXX]\tcumulative_(v(i,j)): " << cumulative);
			   							if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]&&!(fce->f(ip,jp)&SINGLE)){
											cumulative+=erg2(i,j,ip,jp,ct,data,fce->f(i,ip),fce->f(jp,j)) *
                  								v->f(ip+1,jp-1) * erg1(ip,jp,ip+1,jp-1,ct,data);
											LOG_TRACE("[XXXXXX]\tcumulative_(v(i,j)): " << cumulative);
										}
									}
								}
							}
						}
                    } /*/
					// This follows O(N^4) code (used in partition-smp)
					if ((j-i-1)>=(minloop+3)) {
						int maxsize;
						maxsize = min(data->maxintloopsize,j-i-minloop-3);//interior fragment

						for (int size=1;size<=maxsize;++size) {
							for (ip=i+1,jp=j-size-1;ip<=i+size+1;++ip,++jp) {
									{
										// cumulative+=erg2(i,j,ip,jp,ct,data,fce->f(i,ip), fce->f(jp,j)) * v->f(ip,jp);
										cumulative = SUM(cumulative, PROD(erg2(i,j,ip,jp,ct,data,fce->f(i,ip), fce->f(jp,j)), v->f(ip,jp)));

										LOG_TRACE("[FQGH8G]\tcumulative: " << cumulative);

										if (!found&&cumulative>roll) {
											LOG_DEBUG2("Pushing to stack");
											stack.push(ip,jp,1,0,0);
											found=true;

											#ifndef force_continue
												break;
											#endif
										}

										if ((mod[ip]||mod[jp])) if (inc[ct->numseq[ip]][ct->numseq[jp]]&&notgu(ip,jp,ct)&&!(fce->f(ip,jp)&SINGLE)) {
											//i or j is modified
											// cumulative+=erg2(i,j,ip,jp,ct,data,fce->f(i,ip),fce->f(jp,j))*
											// 	v->f(ip+1,jp-1)*erg1(ip,jp,ip+1,jp-1,ct,data);
											cumulative = SUM(cumulative, PROD(erg2(i,j,ip,jp,ct,data,fce->f(i,ip),fce->f(jp,j)),
												v->f(ip+1,jp-1), erg1(ip,jp,ip+1,jp-1,ct,data)));
											LOG_TRACE("[FQGH8G.2]\tcumulative: " << cumulative);

											if (cumulative>roll&&!found) {
												LOG_DEBUG2("Pushing to stack");									
												regbp(ct,number,ip,jp);
												stack.push(ip+1,jp-1,1,0,0);
												found=true;

												#ifndef force_continue
													break;
												#endif
											}
										}
									}
							}//for size
						}//for localip
					}//if ((localj-locali-1)>=(minloop+3))||(localj>(number))
					//*/
					LOG_DEBUG3("[RSYVZ0]\tcumulative: " << cumulative);
					//consider the multiloop closed by i,j
					//[7TI8Q8]
					if ((j-i)>(2*minloop+4)) {

						#ifdef SIMPLEMBLOOP
						//i+1 dangles
						// cumulative+=erg4(i,j,i+1,1,ct,data,lfce[i+1])*penalty(i,j,ct,data)*
            			// 	wmb->f(i+2,j-1)* data->eparam[5] * data->eparam[6] * data->eparam[10]*twoscaling;
						cumulative = SUM(cumulative, PROD(erg4(i,j,i+1,1,ct,data,lfce[i+1]), penalty(i,j,ct,data), 
            				wmb->f(i+2,j-1),  data->eparam[5] ,  data->eparam[6] ,  data->eparam[10] TWOSCALING));

						LOG_DEBUG4("[7TI8Q8]\tcumulative: " << cumulative);

						if (!found&&cumulative>roll) {
							LOG_DEBUG2("Pushing to stack");
							stack.push(i+2,j-1,3,0,0);
							found=true;

							#ifndef force_continue
								break;
							#endif
						}
						#else  //!SIMPLEMBLOOP

						//no dangling ends on i-j pair: [M1RT7B]
						// cumulative+=wmb->f(i+1,j-1)*data->eparam[5]*data->eparam[10]
            			// 	*penalty(i,j,ct,data)*twoscaling;
						cumulative = SUM(cumulative, PROD(wmb->f(i+1,j-1), data->eparam[5], data->eparam[10],
            				penalty(i,j,ct,data) TWOSCALING));

						LOG_DEBUG4("[M1RT7B]\tcumulative: " << cumulative);

						if (cumulative>roll&&!found) {
							LOG_DEBUG2("Pushing to stack");
							stack.push(i+1,j-1,3,0,0);
							found=true;

							#ifndef force_continue
								break;
							#endif
						}

						//i+1 dangles [L23JLU]
						// cumulative+=erg4(i,j,i+1,1,ct,data,lfce[i+1])*penalty(i,j,ct,data)*
            			// 	wmb->f(i+2,j-1)* data->eparam[5] * data->eparam[6] * data->eparam[10]*twoscaling;
						cumulative = SUM(cumulative, PROD(erg4(i,j,i+1,1,ct,data,lfce[i+1]), penalty(i,j,ct,data), 
            				wmb->f(i+2,j-1), data->eparam[5], data->eparam[6], data->eparam[10] TWOSCALING));

						LOG_DEBUG4("[L23JLU]\tcumulative: " << cumulative);

						if (!found&&cumulative>roll) {
							LOG_DEBUG2("Pushing to stack");
							stack.push(i+2,j-1,3,0,0);
							found=true;

							#ifndef force_continue
								break;
							#endif
						}

						//j-1 dangles [C4BBSC]
						// cumulative+=erg4(i,j,j-1,2,ct,data,lfce[j-1]) * penalty(i,j,ct,data) *
            			// 	wmb->f(i+1,j-2) * data->eparam[5] * data->eparam[6] * data->eparam[10]*twoscaling;
						cumulative = SUM(cumulative, PROD(erg4(i,j,j-1,2,ct,data,lfce[j-1]), penalty(i,j,ct,data),
            				wmb->f(i+1,j-2), data->eparam[5], data->eparam[6], data->eparam[10] TWOSCALING));

						LOG_DEBUG4("[C4BBSC]\tcumulative: " << cumulative);

						if (!found&&cumulative>roll) {
							LOG_DEBUG2("Pushing to stack");
							stack.push(i+1,j-2,3,0,0);
							found=true;

							#ifndef force_continue
								break;
							#endif
						}

						//both i+1 and j-1 dangle
            			// cumulative+=data->tstkm[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]]*
						// 		pfchecknp(lfce[i+1],lfce[j-1])*
						// 		wmb->f(i+2,j-2) * data->eparam[5] * data->eparam[6] * data->eparam[6]* data->eparam[10]
						// 		*penalty(i,j,ct,data)*twoscaling;
            			cumulative = SUM(cumulative, PROD(data->tstkm[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]],
								pfchecknp(lfce[i+1],lfce[j-1]),
								wmb->f(i+2,j-2), data->eparam[5], data->eparam[6], data->eparam[6], data->eparam[10],
								penalty(i,j,ct,data) TWOSCALING));

						LOG_DEBUG4("[LXCS43]\tcumulative: " << cumulative);

						if (!found&&cumulative>roll) {
							LOG_DEBUG2("Pushing to stack");
							stack.push(i+2,j-2,3,0,0);
							found=true;

							#ifndef force_continue
								break;
							#endif
						}

						#ifndef disablecoax //a flag to turn off coaxial stacking

						//consider the coaxial stacking of a helix from i to j onto helix i+1 or i+2 to ip:
						//
						for (ip=i+1;(ip<j);ip++) {
							// LOG_DEBUG4("[Coax1.start]\tip: " << ip);

							//if (data->pairing[ct->numseq[i+1]][ct->numseq[ip]]){
							if (inc[ct->numseq[i+1]][ct->numseq[ip]]){
								//first consider flush stacking [SMWF7W] (w and wmb are considered in a single equation in partition, two equations in stochastic)
								// cumulative+=penalty(i,j,ct,data)*v->f(i+1,ip)*
								// 	penalty(i+1,ip,ct,data)*data->eparam[5]
								// 	*data->eparam[10]*data->eparam[10]*(w->f(ip+1,j-1))*ergcoaxflushbases(j,i,i+1,ip,ct,data)
								// 	*twoscaling;
								cumulative = SUM(cumulative, PROD(penalty(i,j,ct,data), v->f(i+1,ip), 
									penalty(i+1,ip,ct,data), data->eparam[5]
									, data->eparam[10], data->eparam[10], (w->f(ip+1,j-1)), ergcoaxflushbases(j,i,i+1,ip,ct,data)
									TWOSCALING));

								LOG_TRACE("[SMWF7W.1]\tcumulative: " << cumulative);

								if (!found&&cumulative>roll) {
									LOG_DEBUG2("Pushing to stack");
									stack.push(ip+1,j-1,4,0,0);
									stack.push(i+1,ip,1,0,0);
									found=true;

									#ifndef force_continue
										break;
									#endif
								}

								// cumulative+=penalty(i,j,ct,data)*v->f(i+1,ip)*
								// 	penalty(i+1,ip,ct,data)*data->eparam[5]
								// 	*data->eparam[10]*data->eparam[10]*(wmb->f(ip+1,j-1))*ergcoaxflushbases(j,i,i+1,ip,ct,data)
								// 	*twoscaling;
								cumulative = SUM(cumulative, PROD(penalty(i,j,ct,data), v->f(i+1,ip), 
									penalty(i+1,ip,ct,data), data->eparam[5],
									data->eparam[10], data->eparam[10], (wmb->f(ip+1,j-1)), ergcoaxflushbases(j,i,i+1,ip,ct,data)
									TWOSCALING));

								LOG_DEBUG4("[SMWF7W.2]\tcumulative: " << cumulative);

								if (!found&&cumulative>roll) {
									LOG_DEBUG2("Pushing to stack");
									stack.push(ip+1,j-1,3,0,0);
									stack.push(i+1,ip,1,0,0);
									found=true;

									#ifndef force_continue
										break;
									#endif
								}

								if((mod[i+1]||mod[ip])) if (inc[ct->numseq[i+2]][ct->numseq[ip-1]]&&notgu(i+1,ip,ct)&&!(fce->f(i+1,ip)&SINGLE)) {
									// cumulative+=penalty(i,j,ct,data)*v->f(i+2,ip-1)*
									// 	penalty(i+1,ip,ct,data)*data->eparam[5]
									// 	*data->eparam[10]*data->eparam[10]*(w->f(ip+1,j-1))*ergcoaxflushbases(j,i,i+1,ip,ct,data)
									// 	*erg1(i+1,ip,i+2,ip-1,ct,data)*twoscaling;
									cumulative = SUM(cumulative, PROD(penalty(i,j,ct,data), v->f(i+2,ip-1), 
										penalty(i+1,ip,ct,data), data->eparam[5]
										, data->eparam[10], data->eparam[10], (w->f(ip+1,j-1)), ergcoaxflushbases(j,i,i+1,ip,ct,data)
										, erg1(i+1,ip,i+2,ip-1,ct,data) TWOSCALING));

									LOG_TRACE("[SMWF7W.3]\tcumulative: " << cumulative);

									if (!found&&cumulative>roll) {
										LOG_DEBUG2("Pushing to stack");
										stack.push(ip+1,j-1,4,0,0);
										stack.push(i+2,ip-1,1,0,0);
										regbp(ct,number,i+1,ip);
										found=true;

										#ifndef force_continue
											break;
										#endif
									}

									// cumulative+=penalty(i,j,ct,data)*v->f(i+2,ip-1)*
									// 	penalty(i+1,ip,ct,data)*data->eparam[5]
									// 	*data->eparam[10]*data->eparam[10]*(wmb->f(ip+1,j-1))*ergcoaxflushbases(j,i,i+1,ip,ct,data)
									// 	*erg1(i+1,ip,i+2,ip-1,ct,data)*twoscaling;
									cumulative = SUM(cumulative, PROD(penalty(i,j,ct,data), v->f(i+2,ip-1), 
										penalty(i+1,ip,ct,data), data->eparam[5],
										data->eparam[10], data->eparam[10], (wmb->f(ip+1,j-1)), ergcoaxflushbases(j,i,i+1,ip,ct,data),
										erg1(i+1,ip,i+2,ip-1,ct,data) TWOSCALING));

									LOG_DEBUG4("[SMWF7W.4]\tcumulative: " << cumulative);

									if (!found&&cumulative>roll) {
										LOG_DEBUG2("Pushing to stack");
										stack.push(ip+1,j-1,3,0,0);
										stack.push(i+2,ip-1,1,0,0);
										regbp(ct,number,i+1,ip);
										found=true;

										#ifndef force_continue
											break;
										#endif
									}
								}
							} // if (inc[ct->numseq[i+1]][ct->numseq[ip]])

							//if (data->pairing[ct->numseq[i+2]][ct->numseq[ip]]){
							if (inc[ct->numseq[i+2]][ct->numseq[ip]]){
								//Now calculate ca stacking with intervening mismatch [60M14C]
								if ((ip+2<j-1)) {
									// cumulative+=penalty(i,j,ct,data)*v->f(i+2,ip)*
									// 	penalty(i+2,ip,ct,data)*data->eparam[5]
									// 	*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
									// 	(w->f(ip+2,j-1))
									// 	*ergcoaxinterbases2(j,i,i+2,ip,ct,data)*twoscaling*pfchecknp(lfce[i+1],lfce[ip+1]);
									cumulative = SUM(cumulative, PROD(penalty(i,j,ct,data), v->f(i+2,ip), 
										penalty(i+2,ip,ct,data), data->eparam[5],
										data->eparam[6], data->eparam[6], data->eparam[10], data->eparam[10], 
										(w->f(ip+2,j-1)),
										ergcoaxinterbases2(j,i,i+2,ip,ct,data), pfchecknp(lfce[i+1],lfce[ip+1]) TWOSCALING));

									LOG_TRACE("[60M14C.1]\tcumulative: " << cumulative);

									if (!found&&cumulative>roll) {
										LOG_DEBUG2("Pushing to stack");
										stack.push(ip+2,j-1,4,0,0);
										stack.push(i+2,ip,1,0,0);
										found=true;

										#ifndef force_continue
											break;
										#endif
									}

									// cumulative+=penalty(i,j,ct,data)*v->f(i+2,ip)*
									// 	penalty(i+2,ip,ct,data)*data->eparam[5]
									// 	*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
									// 	(wmb->f(ip+2,j-1))
									// 	*ergcoaxinterbases2(j,i,i+2,ip,ct,data)*twoscaling*pfchecknp(lfce[i+1],lfce[ip+1]);
									cumulative = SUM(cumulative, PROD(penalty(i,j,ct,data), v->f(i+2,ip),
										penalty(i+2,ip,ct,data), data->eparam[5],
										data->eparam[6], data->eparam[6], data->eparam[10], data->eparam[10], 
										(wmb->f(ip+2,j-1)),
										ergcoaxinterbases2(j,i,i+2,ip,ct,data), pfchecknp(lfce[i+1],lfce[ip+1]) TWOSCALING));

									LOG_DEBUG4("[60M14C.2]\tcumulative: " << cumulative);

									if (!found&&cumulative>roll) {
										LOG_DEBUG2("Pushing to stack");
										stack.push(ip+2,j-1,3,0,0);
										stack.push(i+2,ip,1,0,0);
										found=true;

										#ifndef force_continue
											break;
										#endif

									}

									if((mod[i+2]||mod[ip])) if (inc[ct->numseq[i+3]][ct->numseq[ip-1]]&&notgu(i+2,ip,ct)
										&&!(fce->f(i+2,ip)&SINGLE)) {

										// cumulative+=penalty(i,j,ct,data)*v->f(i+3,ip-1)*
										// 	penalty(i+2,ip,ct,data)*data->eparam[5]
										// 	*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
										// 	(w->f(ip+2,j-1))
										// 	*ergcoaxinterbases2(j,i,i+2,ip,ct,data)
										// 	*erg1(i+2,ip,i+3,ip-1,ct,data)*twoscaling*pfchecknp(lfce[i+1],lfce[ip+1]);
										cumulative = SUM(cumulative, PROD(penalty(i,j,ct,data), v->f(i+3,ip-1),
											penalty(i+2,ip,ct,data), data->eparam[5],
											data->eparam[6], data->eparam[6], data->eparam[10], data->eparam[10],
											(w->f(ip+2,j-1)),
											ergcoaxinterbases2(j,i,i+2,ip,ct,data),
											erg1(i+2,ip,i+3,ip-1,ct,data), pfchecknp(lfce[i+1],lfce[ip+1]) TWOSCALING));

										LOG_TRACE("[60M14C.3]\tcumulative: " << cumulative);

										if (!found&&cumulative>roll) {
											LOG_DEBUG2("Pushing to stack");
											stack.push(ip+2,j-1,4,0,0);
											stack.push(i+3,ip-1,1,0,0);
											regbp(ct,number,i+2,ip);
											found=true;

											#ifndef force_continue
												break;
											#endif
										}

										// cumulative+=penalty(i,j,ct,data)*v->f(i+3,ip-1)*
										// 	penalty(i+2,ip,ct,data)*data->eparam[5]
										// 	*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
										// 	(wmb->f(ip+2,j-1))
										// 	*ergcoaxinterbases2(j,i,i+2,ip,ct,data)
										// 	*erg1(i+2,ip,i+3,ip-1,ct,data)*twoscaling*pfchecknp(lfce[i+1],lfce[ip+1]);
										cumulative = SUM(cumulative, PROD(penalty(i,j,ct,data), v->f(i+3,ip-1),
											penalty(i+2,ip,ct,data), data->eparam[5],
											data->eparam[6], data->eparam[6], data->eparam[10], data->eparam[10],
											(wmb->f(ip+2,j-1)),
											ergcoaxinterbases2(j,i,i+2,ip,ct,data),
											erg1(i+2,ip,i+3,ip-1,ct,data), pfchecknp(lfce[i+1],lfce[ip+1]) TWOSCALING));

										LOG_DEBUG4("[60M14C.4]\tcumulative: " << cumulative);

										if (!found&&cumulative>roll) {
											LOG_DEBUG2("Pushing to stack");
											stack.push(ip+2,j-1,3,0,0);
											stack.push(i+3,ip-1,1,0,0);
											regbp(ct,number,i+2,ip);
											found=true;

											#ifndef force_continue
												break;
											#endif
										}
									}
								} // if ((ip+2<j-1))

								// [AWFSPY]
								if (ip+1<j-2) {
									// cumulative+=penalty(i,j,ct,data)*v->f(i+2,ip)*
									// 	penalty(i+2,ip,ct,data)*data->eparam[5]
									// 	*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
									// 	(w->f(ip+1,j-2))
									// 	*ergcoaxinterbases1(j,i,i+2,ip,ct,data)*twoscaling
									// 	*pfchecknp(lfce[i+1],lfce[j-1]);
									cumulative = SUM(cumulative, PROD(penalty(i,j,ct,data), v->f(i+2,ip),
										penalty(i+2,ip,ct,data), data->eparam[5],
										data->eparam[6], data->eparam[6], data->eparam[10], data->eparam[10],
										(w->f(ip+1,j-2)),
										ergcoaxinterbases1(j,i,i+2,ip,ct,data),
										pfchecknp(lfce[i+1],lfce[j-1]) TWOSCALING));

									LOG_TRACE("[AWFSPY.1]\tcumulative: " << cumulative);

									if (!found&&cumulative>roll) {
										LOG_DEBUG2("Pushing to stack");
										stack.push(ip+1,j-2,4,0,0);
										stack.push(i+2,ip,1,0,0);
										found=true;

										#ifndef force_continue
											break;
										#endif
									}

									// cumulative+=penalty(i,j,ct,data)*v->f(i+2,ip)*
									// 	penalty(i+2,ip,ct,data)*data->eparam[5]
									// 	*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
									// 	(wmb->f(ip+1,j-2))
									// 	*ergcoaxinterbases1(j,i,i+2,ip,ct,data)*twoscaling
									// 	*pfchecknp(lfce[i+1],lfce[j-1]);
									cumulative = SUM(cumulative, PROD(penalty(i,j,ct,data), v->f(i+2,ip),
										penalty(i+2,ip,ct,data), data->eparam[5],
										data->eparam[6], data->eparam[6], data->eparam[10], data->eparam[10],
										(wmb->f(ip+1,j-2)),
										ergcoaxinterbases1(j,i,i+2,ip,ct,data),
										pfchecknp(lfce[i+1],lfce[j-1]) TWOSCALING));

									LOG_DEBUG4("[AWFSPY.2]\tcumulative: " << cumulative);

									if (!found&&cumulative>roll) {
										LOG_DEBUG2("Pushing to stack");
										stack.push(ip+1,j-2,3,0,0);
										stack.push(i+2,ip,1,0,0);

										found=true;

										#ifndef force_continue
											break;
										#endif
									}

									if((mod[i+2]||mod[ip])) if (inc[ct->numseq[i+3]][ct->numseq[ip-1]]&&notgu(i+2,ip,ct)
										&&!(fce->f(i+2,ip)&SINGLE)) {

										// cumulative+=penalty(i,j,ct,data)*v->f(i+3,ip-1)*
										// 	penalty(i+2,ip,ct,data)*data->eparam[5]
										// 	*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
										// 	(w->f(ip+1,j-2))
										// 	*ergcoaxinterbases1(j,i,i+2,ip,ct,data)
										// 	*erg1(i+2,ip,i+3,ip-1,ct,data)*twoscaling*pfchecknp(lfce[i+1],lfce[j-1]);
										cumulative = SUM(cumulative, PROD(penalty(i,j,ct,data), v->f(i+3,ip-1), 
											penalty(i+2,ip,ct,data), data->eparam[5],
											data->eparam[6], data->eparam[6], data->eparam[10], data->eparam[10], 
											(w->f(ip+1,j-2)),
											ergcoaxinterbases1(j,i,i+2,ip,ct,data),
											erg1(i+2,ip,i+3,ip-1,ct,data), pfchecknp(lfce[i+1],lfce[j-1]) TWOSCALING));

										LOG_TRACE("[AWFSPY.3]\tcumulative: " << cumulative);

										if (!found&&cumulative>roll) {
											LOG_DEBUG2("Pushing to stack");
											stack.push(ip+1,j-2,4,0,0);
											stack.push(i+3,ip-1,1,0,0);
											regbp(ct,number,i+2,ip);
											found=true;

											#ifndef force_continue
												break;
											#endif
										}

										// cumulative+=penalty(i,j,ct,data)*v->f(i+3,ip-1)*
										// 	penalty(i+2,ip,ct,data)*data->eparam[5]
										// 	*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
										// 	(wmb->f(ip+1,j-2))
										// 	*ergcoaxinterbases1(j,i,i+2,ip,ct,data)
										// 	*erg1(i+2,ip,i+3,ip-1,ct,data)*twoscaling*pfchecknp(lfce[i+1],lfce[j-1]);
										cumulative = SUM(cumulative, PROD(penalty(i,j,ct,data), v->f(i+3,ip-1), 
											penalty(i+2,ip,ct,data), data->eparam[5],
											data->eparam[6], data->eparam[6], data->eparam[10], data->eparam[10], 
											(wmb->f(ip+1,j-2)),
											ergcoaxinterbases1(j,i,i+2,ip,ct,data),
											erg1(i+2,ip,i+3,ip-1,ct,data), pfchecknp(lfce[i+1],lfce[j-1]) TWOSCALING));

										LOG_DEBUG4("[AWFSPY.4]\tcumulative: " << cumulative);

										if (!found&&cumulative>roll) {
											LOG_DEBUG2("Pushing to stack");
											stack.push(ip+1,j-2,3,0,0);
											stack.push(i+3,ip-1,1,0,0);
											regbp(ct,number,i+2,ip);
											found=true;

											#ifndef force_continue
												break;
											#endif
										}
									}
								}
							} //if (inc[ct->numseq[i+2]][ct->numseq[ip]])
						} //for (ip=i+1;(ip<j);ip++)

						LOG_DEBUG4("[Coax1]\tcumulative: " << cumulative);

						//consider the coaxial stacking of a helix from i to j onto helix ip to j-2 or j-1:
						for (ip=j-1;ip>i;ip--) {
							LOG_DEBUG4("[Coax2.start]\tip: " << ip);

							//first consider flush stacking [T0FCK6]
							// cumulative+=penalty(i,j,ct,data)*v->f(ip,j-1)*
							// 	penalty(j-1,ip,ct,data)*data->eparam[5]
							// 	*data->eparam[10]*data->eparam[10]*(w->f(i+1,ip-1))*ergcoaxflushbases(ip,j-1,j,i,ct,data)
							// 	*twoscaling;
							cumulative = SUM(cumulative, PROD(penalty(i,j,ct,data), v->f(ip,j-1),
								penalty(j-1,ip,ct,data), data->eparam[5],
								data->eparam[10], data->eparam[10], (w->f(i+1,ip-1)), ergcoaxflushbases(ip,j-1,j,i,ct,data)
								TWOSCALING));

							LOG_TRACE("[T0FCK6.1]\tcumulative: " << cumulative);

							if (!found&&cumulative>roll) {
								LOG_DEBUG2("Pushing to stack");
								stack.push(ip,j-1,1,0,0);
								stack.push(i+1,ip-1,4,0,0);
								found=true;

								#ifndef force_continue
									break;
								#endif
							}

							//[GPL3QX]
							// cumulative+=penalty(i,j,ct,data)*v->f(ip,j-1)*
							// 	penalty(j-1,ip,ct,data)*data->eparam[5]
							// 	*data->eparam[10]*data->eparam[10]*(wmb->f(i+1,ip-1))*ergcoaxflushbases(ip,j-1,j,i,ct,data)
							// 	*twoscaling;
							cumulative = SUM(cumulative, PROD(penalty(i,j,ct,data), v->f(ip,j-1),
								penalty(j-1,ip,ct,data), data->eparam[5],
								data->eparam[10], data->eparam[10], (wmb->f(i+1,ip-1)), ergcoaxflushbases(ip,j-1,j,i,ct,data)
								TWOSCALING));

							LOG_DEBUG4("[T0FCK6.2]\tcumulative: " << cumulative);

							if (!found&&cumulative>roll) {
								LOG_DEBUG2("Pushing to stack");
								stack.push(ip,j-1,1,0,0);
								stack.push(i+1,ip-1,3,0,0);
								found=true;

								#ifndef force_continue
									break;
								#endif
							}

							if((mod[ip]||mod[j-1])) if(inc[ct->numseq[ip+1]][ct->numseq[j-2]]&&notgu(ip,j-1,ct)&&!(fce->f(ip,j-1)&SINGLE)) {

								// cumulative+=penalty(i,j,ct,data)*v->f(ip+1,j-2)*
								// 	penalty(j-1,ip,ct,data)*data->eparam[5]
								// 	*data->eparam[10]*data->eparam[10]*(w->f(i+1,ip-1))
								// 	*ergcoaxflushbases(ip,j-1,j,i,ct,data)
								// 	*erg1(ip,j-1,ip+1,j-2,ct,data)*twoscaling;
								cumulative = SUM(cumulative, PROD(penalty(i,j,ct,data), v->f(ip+1,j-2), 
									penalty(j-1,ip,ct,data), data->eparam[5],
									data->eparam[10], data->eparam[10], (w->f(i+1,ip-1)),
									ergcoaxflushbases(ip,j-1,j,i,ct,data),
									erg1(ip,j-1,ip+1,j-2,ct,data) TWOSCALING));

								LOG_TRACE("[T0FCK6.3]\tcumulative: " << cumulative);

								if (!found&&cumulative>roll) {
									LOG_DEBUG2("Pushing to stack");
									stack.push(ip+1,j-2,1,0,0);
									stack.push(i+1,ip-1,4,0,0);
									regbp(ct,number,ip,j-1);
									found=true;

									#ifndef force_continue
										break;
									#endif
								}

								// cumulative+=penalty(i,j,ct,data)*v->f(ip+1,j-2)*
								// 	penalty(j-1,ip,ct,data)*data->eparam[5]
								// 	*data->eparam[10]*data->eparam[10]*(wmb->f(i+1,ip-1))
								// 	*ergcoaxflushbases(ip,j-1,j,i,ct,data)
								// 	*erg1(ip,j-1,ip+1,j-2,ct,data)*twoscaling;
								cumulative = SUM(cumulative, PROD(penalty(i,j,ct,data), v->f(ip+1,j-2), 
									penalty(j-1,ip,ct,data), data->eparam[5],
									data->eparam[10], data->eparam[10], (wmb->f(i+1,ip-1)),
									ergcoaxflushbases(ip,j-1,j,i,ct,data),
									erg1(ip,j-1,ip+1,j-2,ct,data) TWOSCALING));

								LOG_DEBUG4("[T0FCK6.4]\tcumulative: " << cumulative);

								if (!found&&cumulative>roll) {
									LOG_DEBUG2("Pushing to stack");
									stack.push(ip+1,j-2,1,0,0);
									stack.push(i+1,ip-1,3,0,0);
									regbp(ct,number,ip,j-1);
									found=true;

									#ifndef force_continue
										break;
									#endif
								}
							}

							//now consider an intervening nuc [332V5Q]
							if (ip-2>i+1) {
								// cumulative+=penalty(i,j,ct,data)*v->f(ip,j-2)*
								// 	penalty(j-2,ip,ct,data)*data->eparam[5]
								// 	*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
								// 	(w->f(i+1,ip-2))
								// 	*ergcoaxinterbases1(ip,j-2,j,i,ct,data)*twoscaling*pfchecknp(lfce[j-1],lfce[ip-1]);
								cumulative = SUM(cumulative, PROD(penalty(i,j,ct,data), v->f(ip,j-2), 
									penalty(j-2,ip,ct,data), data->eparam[5],
									data->eparam[6], data->eparam[6], data->eparam[10], data->eparam[10], 
									(w->f(i+1,ip-2)),
									ergcoaxinterbases1(ip,j-2,j,i,ct,data), pfchecknp(lfce[j-1],lfce[ip-1]) TWOSCALING));

								LOG_TRACE("[332V5Q.1]\tcumulative: " << cumulative);

								if (!found&&cumulative>roll) {
									LOG_DEBUG2("Pushing to stack");
									stack.push(ip,j-2,1,0,0);
									stack.push(i+1,ip-2,4,0,0);
									found=true;

									#ifndef force_continue
										break;
									#endif
								}

								// cumulative+=penalty(i,j,ct,data)*v->f(ip,j-2)*
								// 	penalty(j-2,ip,ct,data)*data->eparam[5]
								// 	*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
								// 	(wmb->f(i+1,ip-2))
								// 	*ergcoaxinterbases1(ip,j-2,j,i,ct,data)*twoscaling*pfchecknp(lfce[j-1],lfce[ip-1]);
								cumulative= SUM(cumulative, PROD(penalty(i,j,ct,data), v->f(ip,j-2), 
									penalty(j-2,ip,ct,data), data->eparam[5]
									, data->eparam[6], data->eparam[6], data->eparam[10], data->eparam[10], 
									(wmb->f(i+1,ip-2))
									, ergcoaxinterbases1(ip,j-2,j,i,ct,data), pfchecknp(lfce[j-1],lfce[ip-1]) TWOSCALING));

								LOG_DEBUG4("[332V5Q.2]\tcumulative: " << cumulative);
								//LOG_DEBUG4("[332V5Q.2]\tw(i_1,ip-2): " << w->f(i+1,ip-2))
								//LOG_DEBUG4("[332V5Q.2]\twmb(i_1,ip-2): " << wmb->f(i+1,ip-2))

								if (!found&&cumulative>roll) {
									LOG_DEBUG2("Pushing to stack");
									stack.push(ip,j-2,1,0,0);
									stack.push(i+1,ip-2,3,0,0);
									found=true;

									#ifndef force_continue
										break;
									#endif
								}

								if((mod[ip]||mod[j-2])) if(inc[ct->numseq[ip+1]][ct->numseq[j-3]]&&notgu(ip,j-2,ct)&&!(fce->f(ip,j-2)&SINGLE)) {

									// cumulative+=penalty(i,j,ct,data)*v->f(ip+1,j-3)*
									// 	penalty(j-2,ip,ct,data)*data->eparam[5]
									// 	*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
									// 	(w->f(i+1,ip-2))
									// 	*ergcoaxinterbases1(ip,j-2,j,i,ct,data)
									// 	*erg1(ip,j-2,ip+1,j-3,ct,data)*twoscaling*pfchecknp(lfce[j-1],lfce[ip-1]);
									cumulative = SUM(cumulative, PROD(penalty(i,j,ct,data), v->f(ip+1,j-3), 
										penalty(j-2,ip,ct,data), data->eparam[5],
										data->eparam[6], data->eparam[6], data->eparam[10], data->eparam[10], 
										(w->f(i+1,ip-2)),
										ergcoaxinterbases1(ip,j-2,j,i,ct,data),
										erg1(ip,j-2,ip+1,j-3,ct,data), pfchecknp(lfce[j-1],lfce[ip-1]) TWOSCALING));

									LOG_TRACE("[332V5Q.3]\tcumulative: " << cumulative);

									if (!found&&cumulative>roll) {
										LOG_DEBUG2("Pushing to stack");
										stack.push(ip+1,j-3,1,0,0);
										stack.push(i+1,ip-2,4,0,0);
										regbp(ct,number,ip,j-2);
										found=true;

										#ifndef force_continue
											break;
										#endif
									}

									// cumulative+=penalty(i,j,ct,data)*v->f(ip+1,j-3)*
									// 	penalty(j-2,ip,ct,data)*data->eparam[5]
									// 	*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
									// 	(wmb->f(i+1,ip-2))
									// 	*ergcoaxinterbases1(ip,j-2,j,i,ct,data)
									// 	*erg1(ip,j-2,ip+1,j-3,ct,data)*twoscaling*pfchecknp(lfce[j-1],lfce[ip-1]);
									cumulative = SUM(cumulative, PROD(penalty(i,j,ct,data), v->f(ip+1,j-3), 
										penalty(j-2,ip,ct,data), data->eparam[5],
										data->eparam[6], data->eparam[6], data->eparam[10], data->eparam[10], 
										(wmb->f(i+1,ip-2)),
										ergcoaxinterbases1(ip,j-2,j,i,ct,data),
										erg1(ip,j-2,ip+1,j-3,ct,data), pfchecknp(lfce[j-1],lfce[ip-1]) TWOSCALING));

									LOG_DEBUG4("[332V5Q.4]\tcumulative: " << cumulative);

									if (!found&&cumulative>roll) {
										LOG_DEBUG2("Pushing to stack");
										stack.push(ip+1,j-3,1,0,0);
										stack.push(i+1,ip-2,3,0,0);
										regbp(ct,number,ip,j-2);
										found=true;

										#ifndef force_continue
											break;
										#endif
									}
								}
							}

							if ((ip-1>i+2)) {
								// [Q4LM6F]
								// cumulative+=penalty(i,j,ct,data)*v->f(ip,j-2)*
								// 	penalty(j-2,ip,ct,data)*data->eparam[5]
								// 	*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
								// 	(w->f(i+2,ip-1))
								// 	*ergcoaxinterbases2(ip,j-2,j,i,ct,data)*twoscaling*pfchecknp(lfce[j-1],lfce[i+1]);
								cumulative = SUM(cumulative, PROD(penalty(i,j,ct,data), v->f(ip,j-2), 
									penalty(j-2,ip,ct,data), data->eparam[5],
									data->eparam[6], data->eparam[6], data->eparam[10], data->eparam[10], 
									(w->f(i+2,ip-1)),
									ergcoaxinterbases2(ip,j-2,j,i,ct,data), pfchecknp(lfce[j-1],lfce[i+1]) TWOSCALING));

								LOG_TRACE("[Q4LM6F.1]\tcumulative: " << cumulative);

								if (!found&&cumulative>roll) {
									LOG_DEBUG2("Pushing to stack");
									stack.push(ip,j-2,1,0,0);
									stack.push(i+2,ip-1,4,0,0);
									found=true;

									#ifndef force_continue
										break;
									#endif
								}

								//[DM5IZM]
								// cumulative+=penalty(i,j,ct,data)*v->f(ip,j-2)*
								// 	penalty(j-2,ip,ct,data)*data->eparam[5]
								// 	*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
								// 	(wmb->f(i+2,ip-1))
								// 	*ergcoaxinterbases2(ip,j-2,j,i,ct,data)*twoscaling*pfchecknp(lfce[j-1],lfce[i+1]);
								cumulative = SUM(cumulative, PROD(penalty(i,j,ct,data), v->f(ip,j-2), 
									penalty(j-2,ip,ct,data), data->eparam[5],
									data->eparam[6], data->eparam[6], data->eparam[10], data->eparam[10], 
									(wmb->f(i+2,ip-1)),
									ergcoaxinterbases2(ip,j-2,j,i,ct,data), pfchecknp(lfce[j-1],lfce[i+1]) TWOSCALING));

								LOG_DEBUG4("[Q4LM6F.2]\tcumulative: " << cumulative);

								if (!found&&cumulative>roll) {
									LOG_DEBUG2("Pushing to stack");
									stack.push(ip,j-2,1,0,0);
									stack.push(i+2,ip-1,3,0,0);
									found=true;

									#ifndef force_continue
										break;
									#endif
								}

								if((mod[ip]||mod[j-2])) if(inc[ct->numseq[ip+1]][ct->numseq[j-3]]&&notgu(ip,j-2,ct)
									&&!(fce->f(ip,j-2)&SINGLE)) {

									// cumulative+=penalty(i,j,ct,data)*v->f(ip+1,j-3)*
									// 	penalty(j-2,ip,ct,data)*data->eparam[5]
									// 	*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
									// 	(w->f(i+2,ip-1))
									// 	*ergcoaxinterbases2(ip,j-2,j,i,ct,data)
									// 	*erg1(ip,j-2,ip+1,j-3,ct,data)*twoscaling*pfchecknp(lfce[j-1],lfce[i+1]);
									cumulative = SUM(cumulative, PROD(penalty(i,j,ct,data), v->f(ip+1,j-3), 
										penalty(j-2,ip,ct,data), data->eparam[5],
										data->eparam[6], data->eparam[6], data->eparam[10], data->eparam[10], 
										(w->f(i+2,ip-1)),
										ergcoaxinterbases2(ip,j-2,j,i,ct,data),
										erg1(ip,j-2,ip+1,j-3,ct,data), pfchecknp(lfce[j-1],lfce[i+1]) TWOSCALING));

									LOG_TRACE("[Q4LM6F.3]\tcumulative: " << cumulative);

									if (!found&&cumulative>roll) {
										LOG_DEBUG2("Pushing to stack");
										stack.push(ip+1,j-3,1,0,0);
										stack.push(i+2,ip-1,4,0,0);
										regbp(ct,number,ip,j-2);
										found=true;

										#ifndef force_continue
											break;
										#endif
									}

									// cumulative+=penalty(i,j,ct,data)*v->f(ip+1,j-3)*
									// 	penalty(j-2,ip,ct,data)*data->eparam[5]
									// 	*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
									// 	(wmb->f(i+2,ip-1))
									// 	*ergcoaxinterbases2(ip,j-2,j,i,ct,data)
									// 	*erg1(ip,j-2,ip+1,j-3,ct,data)*twoscaling*pfchecknp(lfce[j-1],lfce[i+1]);
									cumulative = SUM(cumulative, PROD(penalty(i,j,ct,data), v->f(ip+1,j-3), 
										penalty(j-2,ip,ct,data), data->eparam[5],
										data->eparam[6], data->eparam[6], data->eparam[10], data->eparam[10], 
										(wmb->f(i+2,ip-1)),
										ergcoaxinterbases2(ip,j-2,j,i,ct,data),
										erg1(ip,j-2,ip+1,j-3,ct,data), pfchecknp(lfce[j-1],lfce[i+1]) TWOSCALING));

									LOG_DEBUG4("[Q4LM6F.4]\tcumulative: " << cumulative);

									if (!found&&cumulative>roll) {
										LOG_DEBUG2("Pushing to stack");
										stack.push(ip+1,j-3,1,0,0);
										stack.push(i+2,ip-1,3,0,0);
										regbp(ct,number,ip,j-2);
										found=true;

										#ifndef force_continue
											break;
										#endif
									}
								}
							}
						} // for (ip=j-1;ip>i;ip--)

						LOG_DEBUG4("[Coax2]\tcumulative: " << cumulative);

						#endif  //disable coax stacking
						#endif  //SIMPLEMBLOOP
					} // if ((j-i)>(2*minloop+4))

					if (!found) {
						cout << "Traceback error in v!\n";
						LOG_ERROR("Traceback errors at v\t" << "i:  " << i << "\tj:  " << j << "\tCase:  " << switchcase << "\tCumulative:  " << cumulative/denominator);
						tracebackerror=14;
					}

					break;

				case 2: //switchcase = 2, dealing with a wcoax fragment
					// if (i==619 && j==622 && only_log_once){
					// 	SET_DEBUG_LEVEL(DEBUG4)
					// 	only_log_once = false;
					// }

					LOG_DEBUG("switchcase: " << switchcase << "\ti,j: " << i << ", " << j);
					#ifndef SIMPLEMBLOOP
					#ifndef disablecoax
					
					denominator = wcoax->f(i,j);
					LOG_DEBUG2("denominator: " << denominator);

					roll = PROD(roll, denominator);
					LOG_DEBUG2("roll: " << roll);

					locale = ZERO;

					for (ip=i+minloop+1;ip<j-minloop-1;ip++) {
						//first consider flush stacking
						// cumulative+=v->f(i,ip)*v->f(ip+1,j)*penalty(i,ip,ct,data)
						// 	*penalty(ip+1,j,ct,data)*ergcoaxflushbases(i,ip,ip+1,j,ct,data)*data->eparam[10]*data->eparam[10];
						cumulative = SUM(cumulative, PROD(v->f(i,ip), v->f(ip+1,j), penalty(i,ip,ct,data),
							penalty(ip+1,j,ct,data), ergcoaxflushbases(i,ip,ip+1,j,ct,data), data->eparam[10], data->eparam[10]));

						LOG_DEBUG4("[TYG2WQ]\tcumulative: " << cumulative);

						if (SUM(cumulative,locale) > roll&&!found) {
							LOG_DEBUG2("Coax(wcoax)\tip:  " << ip);
							LOG_DEBUG2("Pushing to stack");
							stack.push(i,ip,1,0,0);
							stack.push(ip+1,j,1,0,0);
							found=true;

							#ifndef force_continue
								break;
							#endif
						}

						if ((mod[i]||mod[ip]||mod[ip+1]||mod[j])) {
							if ((mod[i]||mod[ip])&&(mod[ip+1]||mod[j])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]
								&&inc[ct->numseq[ip+2]][ct->numseq[j-1]]&&notgu(i,ip,ct)&&notgu(ip+1,j,ct)
								&&!(fce->f(ip+1,j)&SINGLE)&&!(fce->f(i,ip)&SINGLE)) {

								// cumulative+=v->f(i+1,ip-1)*v->f(ip+2,j-1)*penalty(i,ip,ct,data)
								// 	*penalty(ip+1,j,ct,data)*ergcoaxflushbases(i,ip,ip+1,j,ct,data)
								// 	*erg1(i,ip,i+1,ip-1,ct,data)*erg1(ip+1,j,ip+2,j-1,ct,data)*data->eparam[10]*data->eparam[10];
								cumulative = SUM(cumulative, PROD(v->f(i+1,ip-1), v->f(ip+2,j-1), penalty(i,ip,ct,data),
									penalty(ip+1,j,ct,data), ergcoaxflushbases(i,ip,ip+1,j,ct,data),
									erg1(i,ip,i+1,ip-1,ct,data), erg1(ip+1,j,ip+2,j-1,ct,data), data->eparam[10], data->eparam[10]));

								LOG_DEBUG4("[TYG2WQ.2]\tcumulative: " << cumulative);

								if (!found&&SUM(cumulative,locale)>roll) {
									LOG_DEBUG2("Coax(wcoax)\tip:  " << ip);
									LOG_DEBUG2("Pushing to stack");
									stack.push(i+1,ip-1,1,0,0);
									stack.push(ip+2,j-1,1,0,0);
									found=true;
									regbp(ct,number,i,ip);
									regbp(ct,number,ip+1,j);

									#ifndef force_continue
										break;
									#endif
								}
							}

							if ((mod[i]||mod[ip])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]&&notgu(i,ip,ct)&&!(fce->f(i,ip)&SINGLE)) {

								// cumulative+=v->f(i+1,ip-1)*v->f(ip+1,j)*penalty(i,ip,ct,data)
								// 	*penalty(ip+1,j,ct,data)*ergcoaxflushbases(i,ip,ip+1,j,ct,data)
								// 	*erg1(i,ip,i+1,ip-1,ct,data)*data->eparam[10]*data->eparam[10];
								cumulative = SUM(cumulative, PROD(v->f(i+1,ip-1), v->f(ip+1,j), penalty(i,ip,ct,data),
									penalty(ip+1,j,ct,data), ergcoaxflushbases(i,ip,ip+1,j,ct,data),
									erg1(i,ip,i+1,ip-1,ct,data), data->eparam[10], data->eparam[10]));

								LOG_DEBUG4("[TYG2WQ.3]\tcumulative: " << cumulative);

								if (!found&&SUM(cumulative,locale)>roll) {
									LOG_DEBUG2("Coax(wcoax)\tip:  " << ip);
									LOG_DEBUG2("Pushing to stack");
									stack.push(i+1,ip-1,1,0,0);
									stack.push(ip+1,j,1,0,0);
									found=true;
									regbp(ct,number,i,ip);

									#ifndef force_continue
										break;
									#endif
								}
							}

							if ((mod[ip+1]||mod[j])&&inc[ct->numseq[ip+2]][ct->numseq[j-1]]&&notgu(ip+1,j,ct)&&!(fce->f(ip+1,j)&SINGLE)) {

								// cumulative+=v->f(i,ip)*v->f(ip+2,j-1)*penalty(i,ip,ct,data)
								// 	*penalty(ip+1,j,ct,data)*ergcoaxflushbases(i,ip,ip+1,j,ct,data)
								// 	*erg1(ip+1,j,ip+2,j-1,ct,data)*data->eparam[10]*data->eparam[10];
								cumulative = SUM(cumulative, PROD(v->f(i,ip), v->f(ip+2,j-1), penalty(i,ip,ct,data),
									penalty(ip+1,j,ct,data), ergcoaxflushbases(i,ip,ip+1,j,ct,data),
									erg1(ip+1,j,ip+2,j-1,ct,data), data->eparam[10], data->eparam[10]));

								LOG_DEBUG4("[TYG2WQ.4]\tcumulative: " << cumulative);

								if (!found&&SUM(cumulative,locale)>roll) {
									LOG_DEBUG2("Coax(wcoax)\tip:  " << ip);
									LOG_DEBUG2("Pushing to stack");
									stack.push(i,ip,1,0,0);
									stack.push(ip+2,j-1,1,0,0);
									found=true;
									regbp(ct,number,ip+1,j);

									#ifndef force_continue
										break;
									#endif
								}
							}
						}

						if (!lfce[ip+1]&&!lfce[j]) {
							//now consider an intervening mismatch //[NGXFM9]
							// locale+=v->f(i,ip)*v->f(ip+2,j-1)*penalty(i,ip,ct,data)
							// 	*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)*data->eparam[10]*data->eparam[10]*data->eparam[6]*data->eparam[6];
							locale = SUM(locale, PROD(v->f(i,ip), v->f(ip+2,j-1), penalty(i,ip,ct,data),
								penalty(ip+2,j-1,ct,data), ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data), data->eparam[10], data->eparam[10], data->eparam[6], data->eparam[6]));

							LOG_DEBUG4("[XXXXXX]\tlocale: " << locale);

							if ((SUM(cumulative,locale))>roll&&!found) {
								LOG_DEBUG2("Coax(wcoax)\tip:  " << ip);
								LOG_DEBUG2("Increment:  " << v->f(i,ip)*v->f(ip+2,j-1)*penalty(i,ip,ct,data)*penalty(ip+2,j-1,ct,data)
										   *ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)*data->eparam[10]*data->eparam[10]*data->eparam[6]*data->eparam[6];);
								LOG_DEBUG2("Pushing to stack");
								stack.push(i,ip,1,0,0);
								stack.push(ip+2,j-1,1,0,0);
								found=true;

								#ifndef force_continue
									//cumulative += locale;
									break;
								#endif
							}

							if (mod[i]||mod[ip]||mod[ip+2]||mod[j-1]) {
								if ((mod[i]||mod[ip])&&(mod[ip+2]||mod[j-1])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]
									&&inc[ct->numseq[ip+3]][ct->numseq[j-2]]&&notgu(i,ip,ct)&&notgu(ip+2,j-1,ct)
										&&!(fce->f(i,ip)&SINGLE)&&!(fce->f(ip+2,j-1)&SINGLE)) {

									// locale+=v->f(i+1,ip-1)*v->f(ip+3,j-2)*penalty(i,ip,ct,data)
									// 	*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)
									// 	*erg1(i,ip,i+1,ip-1,ct,data)*erg1(ip+2,j-1,ip+3,j-2,ct,data)*data->eparam[10]*data->eparam[10]*data->eparam[6]*data->eparam[6];
									locale = SUM(locale, PROD(v->f(i+1,ip-1), v->f(ip+3,j-2), penalty(i,ip,ct,data),
										penalty(ip+2,j-1,ct,data), ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data),
										erg1(i,ip,i+1,ip-1,ct,data), erg1(ip+2,j-1,ip+3,j-2,ct,data), data->eparam[10], data->eparam[10], data->eparam[6], data->eparam[6]));

									LOG_DEBUG4("[XXXXXX]\tlocale: " << locale);

									if (!found&&(SUM(cumulative,locale))>roll) {
										LOG_DEBUG2("Coax(wcoax)\tip:  " << ip);
										LOG_DEBUG2("Pushing to stack");
										stack.push(i+1,ip-1,1,0,0);
										stack.push(ip+3,j-2,1,0,0);
										regbp(ct,number,i,ip);
										regbp(ct,number,ip+2,j-1);
										found=true;

										#ifndef force_continue
											//cumulative += locale;
											break;
										#endif
									}
								}

								if ((mod[i]||mod[ip])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]&&notgu(i,ip,ct)&&!(fce->f(i,ip)&SINGLE)) {

									// locale+=v->f(i+1,ip-1)*v->f(ip+2,j-1)*penalty(i,ip,ct,data)
									// 	*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)
									// 	*erg1(i,ip,i+1,ip-1,ct,data)*data->eparam[10]*data->eparam[10]*data->eparam[6]*data->eparam[6];
									locale = SUM(locale, PROD(v->f(i+1,ip-1), v->f(ip+2,j-1), penalty(i,ip,ct,data),
										penalty(ip+2,j-1,ct,data), ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data),
										erg1(i,ip,i+1,ip-1,ct,data), data->eparam[10], data->eparam[10], data->eparam[6], data->eparam[6]));

									LOG_DEBUG4("[XXXXXX]\tlocale: " << locale);

									if(!found&&(SUM(cumulative,locale))>roll) {
										LOG_DEBUG2("Coax(wcoax)\tip:  " << ip);
										LOG_DEBUG2("Pushing to stack");
										stack.push(i+1,ip-1,1,0,0);
										stack.push(ip+2,j-1,1,0,0);
										regbp(ct,number,i,ip);
										found=true;

										#ifndef force_continue
											//cumulative += locale;
											break;
										#endif
									}
								}

								if ((mod[ip+2]||mod[j-1])&&inc[ct->numseq[ip+3]][ct->numseq[j-2]]&&notgu(ip+2,j-1,ct)&&!(fce->f(ip+2,j-1)&SINGLE)) {

									// locale+=v->f(i,ip)*v->f(ip+3,j-2)*penalty(i,ip,ct,data)
									// 	*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)
									// 	*erg1(ip+2,j-1,ip+3,j-2,ct,data)*data->eparam[10]*data->eparam[10]*data->eparam[6]*data->eparam[6];
									locale = SUM(locale, PROD(v->f(i,ip), v->f(ip+3,j-2), penalty(i,ip,ct,data),
										penalty(ip+2,j-1,ct,data), ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data),
										erg1(ip+2,j-1,ip+3,j-2,ct,data), data->eparam[10], data->eparam[10], data->eparam[6], data->eparam[6]));

									LOG_DEBUG4("[XXXXXX]\tlocale: " << locale);

									if (!found&&(SUM(cumulative,locale))>roll) {
										LOG_DEBUG2("Coax(wcoax)\tip:  " << ip);
										LOG_DEBUG2("Pushing to stack");
										stack.push(i,ip,1,0,0);
										stack.push(ip+3,j-2,1,0,0);
										regbp(ct,number,ip+2,j-1);
										found=true;

										#ifndef force_continue
											//cumulative += locale;
											break;
										#endif
									}
								}
							}
						}

						if(!lfce[i]&&!lfce[ip+1]) {
							// locale+=v->f(i+1,ip)*v->f(ip+2,j)*penalty(i+1,ip,ct,data)
							// 	*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)*data->eparam[10]*data->eparam[10]*data->eparam[6]*data->eparam[6];
							locale = SUM(locale, PROD(v->f(i+1,ip), v->f(ip+2,j), penalty(i+1,ip,ct,data),
								penalty(ip+2,j,ct,data), ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data), data->eparam[10], data->eparam[10], data->eparam[6], data->eparam[6]));

							LOG_DEBUG4("[XXXXXX]\tlocale: " << locale);

							if (!found&&(SUM(cumulative,locale))>roll) {
								LOG_DEBUG2("Coax(wcoax)\tip:  " << ip);
								LOG_DEBUG2("Pushing to stack");
								stack.push(i+1,ip,1,0,0);
								stack.push(ip+2,j,1,0,0);
								found=true;

								#ifndef force_continue
									//cumulative += locale;
									break;
								#endif
							}

							if (mod[i+1]||mod[ip]||mod[ip+2]||mod[j]) {
								if ((mod[i+1]||mod[ip])&&(mod[ip+2]||mod[j])&&inc[ct->numseq[i+2]][ct->numseq[ip-1]]
									&&inc[ct->numseq[ip+3]][ct->numseq[j-1]]&&notgu(i+1,ip,ct)&&notgu(ip+2,j,ct)
									&&!(fce->f(i+1,ip)&SINGLE)&&!(fce->f(ip+2,j)&SINGLE)	) {

									// locale+=v->f(i+2,ip-1)*v->f(ip+3,j-1)*penalty(i+1,ip,ct,data)
									// 	*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)
									// 	*erg1(i+1,ip,i+2,ip-1,ct,data)*erg1(ip+2,j,ip+3,j-1,ct,data)*data->eparam[10]*data->eparam[10]*data->eparam[6]*data->eparam[6];
									locale = SUM(locale, PROD(v->f(i+2,ip-1), v->f(ip+3,j-1), penalty(i+1,ip,ct,data)
										, penalty(ip+2,j,ct,data), ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)
										, erg1(i+1,ip,i+2,ip-1,ct,data), erg1(ip+2,j,ip+3,j-1,ct,data), data->eparam[10], data->eparam[10], data->eparam[6], data->eparam[6]));

									LOG_DEBUG4("[XXXXXX]\tlocale: " << locale);

									if (!found&&SUM(cumulative,locale)>roll) {
										LOG_DEBUG2("Coax(wcoax)\tip:  " << ip);
										LOG_DEBUG2("Pushing to stack");
										stack.push(i+2,ip-1,1,0,0);
										stack.push(ip+3,j-1,1,0,0);
										regbp(ct,number,i+1,ip);
										regbp(ct,number,ip+2,j);
										found=true;

										#ifndef force_continue
											//cumulative += locale;
											break;
										#endif
									}
								}

								if ((mod[i+1]||mod[ip])&&inc[ct->numseq[i+2]][ct->numseq[ip-1]]&&notgu(i+1,ip,ct)&&!(fce->f(i+1,ip)&SINGLE)) {

									// locale+=v->f(i+2,ip-1)*v->f(ip+2,j)*penalty(i+1,ip,ct,data)
									// 	*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)
									// 	*erg1(i+1,ip,i+2,ip-1,ct,data)*data->eparam[10]*data->eparam[10]*data->eparam[6]*data->eparam[6];
									locale = SUM(locale, PROD(v->f(i+2,ip-1), v->f(ip+2,j), penalty(i+1,ip,ct,data),
										penalty(ip+2,j,ct,data), ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data),
										erg1(i+1,ip,i+2,ip-1,ct,data), data->eparam[10], data->eparam[10], data->eparam[6], data->eparam[6]));

									LOG_DEBUG4("[XXXXXX]\tlocale: " << locale);

									if (!found&&SUM(cumulative,locale)>roll) {
										LOG_DEBUG2("Coax(wcoax)\tip:  " << ip);
										LOG_DEBUG2("Pushing to stack");
										stack.push(i+2,ip-1,1,0,0);
										stack.push(ip+2,j,1,0,0);
										regbp(ct,number,i+1,ip);

										found=true;

										#ifndef force_continue
											//cumulative += locale;
											break;
										#endif
									}
								}

								if ((mod[ip+2]||mod[j])&&inc[ct->numseq[ip+3]][ct->numseq[j-1]]&&notgu(ip+2,j,ct)&&!(fce->f(ip+2,j)&SINGLE)) {

									// locale+=v->f(i+1,ip)*v->f(ip+3,j-1)*penalty(i+1,ip,ct,data)
									// 	*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)
									// 	*erg1(ip+2,j,ip+3,j-1,ct,data)*data->eparam[10]*data->eparam[10]*data->eparam[6]*data->eparam[6];
									locale = SUM(locale, PROD(v->f(i+1,ip), v->f(ip+3,j-1), penalty(i+1,ip,ct,data),
										penalty(ip+2,j,ct,data), ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data),
										erg1(ip+2,j,ip+3,j-1,ct,data), data->eparam[10], data->eparam[10], data->eparam[6], data->eparam[6]));

									LOG_DEBUG4("[XXXXXX]\tlocale: " << locale);

									if (!found&&SUM(cumulative,locale)>roll) {
										LOG_DEBUG2("Coax(wcoax)\tip:  " << ip);
										LOG_DEBUG2("Pushing to stack");
										stack.push(i+1,ip,1,0,0);
										stack.push(ip+3,j-1,1,0,0);

										regbp(ct,number,ip+2,j);
										found=true;

										#ifndef force_continue
											//cumulative += locale;
											break;
										#endif
									}

								}
							}
						} // if(!lfce[i]&&!lfce[ip+1])
					} //for (ip=i+minloop+1;ip<j-minloop-1;ip++)

					cumulative = SUM(cumulative, locale);

					if (!found) {

						cout << "Traceback error at wcoax/n";
						LOG_ERROR("Traceback errors at wcoax\t" << "i:  " << i << "\tj:  " << j << "\tCase:  " << switchcase << "\tCumulative:  " << cumulative/denominator << "\troll:  " << roll/denominator);
						tracebackerror=14;

					}

					#endif
					#endif //SIMPLEMBLOOP

					LOG_DEBUG4("denominator: " << denominator);
					break;

				case 3: //switchcase = 3, dealing with wmb fragment
					// if (i==619 && j==622 && only_log_once){
					// 	SET_DEBUG_LEVEL(DEBUG4)
					// 	only_log_once = false;
					// }

					LOG_DEBUG("switchcase: " << switchcase << "\ti,j: " << i << ", " << j);
					denominator = wmb->f(i,j);
					roll = PROD(roll, denominator);

					cumulative = wmbl->f(i,j);
					LOG_DEBUG4("[XXXXXX]\tcumulative: " << cumulative);

					if (cumulative>roll) {
						LOG_DEBUG2("Pushing to stack");
						found=true;
						stack.push(i,j,5,0,0);

						#ifndef force_continue
							break;
						#endif
					}

					if (!lfce[j]) {
						// cumulative+=wmb->f(i,j-1)*data->eparam[6]*data->scaling;
						cumulative=SUM(cumulative, PROD(wmb->f(i,j-1), data->eparam[6] SCALING));
						LOG_DEBUG4("[XXXXXX]\tcumulative: " << cumulative);

						if (!found&&cumulative>roll) {
							LOG_DEBUG2("Pushing to stack");
							found=true;
							stack.push(i,j-1,3,0,0);

							#ifndef force_continue
								break;
							#endif
						}
					}

					if (!found) {
						cout << "Traceback error at wmb\n";
						LOG_ERROR("Traceback errors at wmb\t" << "i:  " << i << "\tj:  " << j << "\tCase:  " << switchcase << "\tCumulative:  " << cumulative/denominator);
						tracebackerror=14;
					}

					break;

				case 4: //switchcase = 4, dealing with a w fragment
					// if (i==619 && j==622 && only_log_once){
					// 	SET_DEBUG_LEVEL(DEBUG4)
					// 	only_log_once = false;
					// }

					LOG_DEBUG("switchcase: " << switchcase << "\ti,j: " << i << ", " << j);
					denominator = w->f(i,j);
					roll = PROD(roll, denominator);

					cumulative = wl->f(i,j);

					LOG_DEBUG4("[XXXXXX]\tcumulative: " << cumulative);

					if (cumulative>roll) {
						LOG_DEBUG2("Pushing to stack");
						stack.push(i,j,6,0,0);
						found=true;

						#ifndef force_continue
							break;
						#endif
					}

					if (!lfce[j]) {
               			// cumulative+= w->f(i,j-1) * data->eparam[6]*data->scaling;
               			cumulative = SUM(cumulative, PROD(w->f(i,j-1), data->eparam[6] SCALING));

						LOG_DEBUG4("[XXXXXX]\tcumulative: " << cumulative);

						if (cumulative>roll&&!found) {
							LOG_DEBUG2("Pushing to stack");
							stack.push(i,j-1,4,0,0);
							found=true;

							#ifndef force_continue
								break;
							#endif
						}
					}

					if (!found) {
						cout << "Traceback error at w\n";
						LOG_ERROR("Traceback errors at w\t" << "i:  " << i << "\tj:  " << j << "\tCase:  " << switchcase << "\tCumulative:  " << cumulative/denominator);
						tracebackerror=14;
					}

					break;

				case 5:  //switchcase = 5, wmbl fragment
					LOG_DEBUG("switchcase: " << switchcase << "\ti,j: " << i << ", " << j);
					denominator = wmbl->f(i,j);
					roll = PROD(roll, denominator);

					cumulative=wcoax->f(i,j);

					LOG_DEBUG3("[304YK0]\tcumulative: " << cumulative);

					if (cumulative>roll&&!found) {
						LOG_DEBUG2("Pushing to stack");
						stack.push(i,j,2,0,0);

						found=true;

						#ifndef force_continue
							break;
						#endif
					}

					for (k=i+1;k<j;k++) {
						tmp_cumulative = cumulative;
						if (!lfce[i]) {
							cumulative=SUM(cumulative, PROD(SUM(wlc->f(i,k), wcoax->f(i,k)), SUM(wl->f(k+1,j), wmbl->f(k+1,j))));
							LOG_DEBUG4("[0LNQ9E]\tcumulative: " << cumulative);

							if (cumulative>roll&&!found) {							
								// tmp_cumulative+=wlc->f(i,k)*(wl->f(k+1,j));
								tmp_cumulative=SUM(tmp_cumulative, PROD(wlc->f(i,k), wl->f(k+1,j)));

								if (tmp_cumulative>roll&&!found) {
									LOG_DEBUG2("Pushing to stack");
									stack.push(i,k,7,0,0);
									stack.push(k+1,j,6,0,0);
									found=true;

									#ifndef force_continue
										break;
									#endif
								}

								tmp_cumulative=SUM(tmp_cumulative, PROD(wcoax->f(i,k), wl->f(k+1,j)));

								if (tmp_cumulative>roll&&!found) {
									LOG_DEBUG2("Pushing to stack");
									stack.push(i,k,2,0,0);
									stack.push(k+1,j,6,0,0);
									found=true;

									#ifndef force_continue
										break;
									#endif
								}

								tmp_cumulative=SUM(tmp_cumulative, PROD(wlc->f(i,k), wmbl->f(k+1,j)));

								if (tmp_cumulative>roll&&!found) {
									LOG_DEBUG2("Pushing to stack");
									stack.push(i,k,7,0,0);
									stack.push(k+1,j,5,0,0);
									found=true;

									#ifndef force_continue
										break;
									#endif
								}

								// cumulative+=(wcoax->f(i,k))*(wmbl->f(k+1,j));

								// if (cumulative>roll&&!found) {
								else{
									if (!found){
									LOG_DEBUG2("Pushing to stack");
									stack.push(i,k,2,0,0);
									stack.push(k+1,j,5,0,0);
									found=true;

									#ifndef force_continue
										break;
									#endif
									}
								}
							}
						}
						else {
							// cumulative+=(wl->f(i,k)+ wcoax->f(i,k))* (wl->f(j,k+1)+ wmbl->f(j,k+1));
							cumulative=SUM(cumulative, PROD(SUM(wl->f(i,k)+ wcoax->f(i,k))* SUM(wl->f(j,k+1)+ wmbl->f(j,k+1))));
							LOG_DEBUG4("[NQ84MW]\tcumulative: " << cumulative);

							if (cumulative>roll&&!found) {
								tmp_cumulative=SUM(tmp_cumulative, PROD(wl->f(i,k), wl->f(k+1,j)));

								if (tmp_cumulative>roll&&!found) {
									LOG_DEBUG2("Pushing to stack");
									stack.push(i,k,6,0,0);
									stack.push(k+1,j,6,0,0);
									found=true;

									#ifndef force_continue
										break;
									#endif
								}

								tmp_cumulative+=SUM(tmp_cumulative, PROD(wcoax->f(i,k), wl->f(k+1,j)));

								if (tmp_cumulative>roll&&!found) {
									LOG_DEBUG2("Pushing to stack");
									stack.push(i,k,2,0,0);
									stack.push(k+1,j,6,0,0);
									found=true;

									#ifndef force_continue
										break;
									#endif
								}

								tmp_cumulative=SUM(tmp_cumulative, PROD(wl->f(i,k), wmbl->f(k+1,j)));

								if (tmp_cumulative>roll&&!found) {
									LOG_DEBUG2("Pushing to stack");
									stack.push(i,k,6,0,0);
									stack.push(k+1,j,5,0,0);
									found=true;

									#ifndef force_continue
										break;
									#endif
								}

								// cumulative+=(wcoax->f(i,k))*(wmbl->f(k+1,j));

								// if (cumulative>roll&&!found) {
								else{
									if (!found){
									LOG_DEBUG2("Pushing to stack");
									stack.push(i,k,2,0,0);
									stack.push(k+1,j,5,0,0);
									found=true;

									#ifndef force_continue
										break;
									#endif
									}
								}
							}
						}
					} //for (k=i+1;k<j;k++)

					if (!lfce[i]) {
						// cumulative+=wmbl->f(i+1,j)*data->eparam[6]*data->scaling; //[572PKS]
						cumulative=SUM(cumulative, PROD(wmbl->f(i+1,j), data->eparam[6] SCALING)); //[572PKS]

						LOG_DEBUG4("[NYT5QE]\tcumulative: " << cumulative);

						if (cumulative>roll&&!found) {
							LOG_DEBUG2("Pushing to stack");
							stack.push(i+1,j,5,0,0);
							found=true;

							#ifndef force_continue
								break;
							#endif
						}
					}

					if (!found) {
						cout << "Traceback error at wmbl\n";
						LOG_ERROR("Traceback error at wmbl\t" << "i:  " << i << "\tj:  " << j << "\tCase:  " << switchcase << "\tCumulative:  " << cumulative/denominator);
						tracebackerror=14;
					}

					break;

				case 6: //switchcase 6, wl fragment
					// if (i==619 && j==622 && only_log_once){
					// 	SET_DEBUG_LEVEL(DEBUG4)
					// 	only_log_once = false;
					// }

					LOG_DEBUG("switchcase: " << switchcase << "\ti,j: " << i << ", " << j);
					denominator = wl->f(i,j);
					roll = PROD(roll, denominator);

					#ifdef SIMPLEMBLOOP
					cumulative = PROD(v->f(i,j-1), data->eparam[10], data->eparam[6],
         				erg4(j-1,i,j,1,ct,data,lfce[j]), penalty(i,j-1,ct,data));

					LOG_DEBUG4("[XXXXXX]\tcumulative: " << cumulative);

					if (cumulative>roll&&!found) {
						LOG_DEBUG2("Pushing to stack");
						stack.push(i,j-1,1,0,0);
						found=true;

						#ifndef force_continue
							break;
						#endif
					}

					if ((mod[i]||mod[j-1])) if(inc[ct->numseq[i+1]][ct->numseq[j-2]]&&notgu(i,j-1,ct)&&!(fce->f(i,j-1)&SINGLE)) {

						// cumulative+= v->f(i+1,j-2) * data->eparam[10] * data->eparam[6] *
         				// 	erg4(j-1,i,j,1,ct,data,lfce[j])*penalty(i,j-1,ct,data)
						// 	*erg1(i,j-1,i+1,j-2,ct,data);
						cumulative = SUM(cumulative, PROD(v->f(i+1,j-2), data->eparam[10], data->eparam[6],
         					erg4(j-1,i,j,1,ct,data,lfce[j]), penalty(i,j-1,ct,data),
							erg1(i,j-1,i+1,j-2,ct,data)));

						LOG_DEBUG4("[XXXXXX]\tcumulative: " << cumulative);

						if (cumulative>roll&&!found) {
							LOG_DEBUG2("Pushing to stack");
							stack.push(i+1,j-2,1,0,0);
							regbp(ct,number,i,j-1);
							found=true;

							#ifndef force_continue
								break;
							#endif
						}
					}
					#else //!SIMPLEMBLOOP

					// cumulative= data->eparam[10]*v->f(i,j)*penalty(j,i,ct,data);
					cumulative= PROD(data->eparam[10], v->f(i,j), penalty(j,i,ct,data));

					LOG_DEBUG4("[89AGXL]\tcumulative: " << cumulative);

					if (cumulative>roll&&!found) {
						LOG_DEBUG2("Pushing to stack");
						found=true;
						stack.push(i,j,1,0,0);

						#ifndef force_continue
							break;
						#endif
					}

					if ((mod[i]||mod[j])) if (inc[ct->numseq[i+1]][ct->numseq[j-1]]&&notgu(i,j,ct)) {

						// cumulative+= data->eparam[10]*v->f(i+1,j-1)*penalty(j,i,ct,data)*erg1(i,j,i+1,j-1,ct,data);
						cumulative = SUM(cumulative, PROD(data->eparam[10], v->f(i+1,j-1), penalty(j,i,ct,data), erg1(i,j,i+1,j-1,ct,data)));

						LOG_DEBUG4("[89AGXL.2]\tcumulative: " << cumulative);

						if (!found&&cumulative>roll) {
							LOG_DEBUG2("Pushing to stack");
							stack.push(i+1,j-1,1,0,0);
							regbp(ct,number,i,j);
							found = true;

							#ifndef force_continue
								break;
							#endif
						}
					}
					//[6ZD0LC]
					// cumulative+= v->f(i+1,j)*data->eparam[10]*data->eparam[6]*
         			// 	erg4(j,i+1,i,2,ct,data,lfce[i])*penalty(i+1,j,ct,data);
					cumulative = SUM(cumulative, PROD(v->f(i+1,j), data->eparam[10], data->eparam[6], 
         				erg4(j,i+1,i,2,ct,data,lfce[i]), penalty(i+1,j,ct,data)));

					LOG_DEBUG4("[G98UOW]\tcumulative: " << cumulative);

					if (!found&&cumulative>roll) {
						LOG_DEBUG2("Pushing to stack");
						stack.push(i+1,j,1,0,0);
						found=true;

						#ifndef force_continue
							break;
						#endif
					}

					if ((mod[i+1]||mod[j])) if(inc[ct->numseq[i+2]][ct->numseq[j-1]]&&notgu(i+1,j,ct)&&!(fce->f(i+1,j)&SINGLE)) {

						// cumulative+= v->f(i+2,j-1) * data->eparam[10] *data->eparam[6] *
         				// 	erg4(j,i+1,i,2,ct,data,lfce[i])*penalty(i+1,j,ct,data)
						// 	*erg1(i+1,j,i+2,j-1,ct,data);
						cumulative = SUM(cumulative, PROD(v->f(i+2,j-1), data->eparam[10], data->eparam[6], 
         					erg4(j,i+1,i,2,ct,data,lfce[i]), penalty(i+1,j,ct,data), erg1(i+1,j,i+2,j-1,ct,data)));

						LOG_DEBUG4("[G98UOW.2]\tcumulative: " << cumulative);

						if (!found&&cumulative>roll) {
							LOG_DEBUG2("Pushing to stack");
							stack.push(i+2,j-1,1,0,0);
							found=true;
							regbp(ct,number,i+1,j);

							#ifndef force_continue
								break;
							#endif
						}
					}

					// cumulative+= v->f(i,j-1)* data->eparam[10] * data->eparam[6] *
         			// 	erg4(j-1,i,j,1,ct,data,lfce[j])*penalty(i,j-1,ct,data);
					cumulative = SUM(cumulative, PROD(v->f(i,j-1), data->eparam[10], data->eparam[6],
         				erg4(j-1,i,j,1,ct,data,lfce[j]), penalty(i,j-1,ct,data)));

					LOG_DEBUG4("[YBN850\tcumulative: " << cumulative);

					if (cumulative>roll&&!found) {
						LOG_DEBUG2("Pushing to stack");
						stack.push(i,j-1,1,0,0);
						found=true;

						#ifndef force_continue
							break;
						#endif
					}

					if ((mod[i]||mod[j-1])) if(inc[ct->numseq[i+1]][ct->numseq[j-2]]&&notgu(i,j-1,ct)&&!(fce->f(i,j-1)&SINGLE)) {

						// cumulative+= v->f(i+1,j-2) * data->eparam[10] * data->eparam[6] *
         				// 	erg4(j-1,i,j,1,ct,data,lfce[j])*penalty(i,j-1,ct,data)
						// 	*erg1(i,j-1,i+1,j-2,ct,data);
						cumulative = SUM(cumulative, PROD(v->f(i+1,j-2), data->eparam[10], data->eparam[6],
         					erg4(j-1,i,j,1,ct,data,lfce[j]), penalty(i,j-1,ct,data), erg1(i,j-1,i+1,j-2,ct,data)));

						LOG_DEBUG4("[YBN850.2]\tcumulative: " << cumulative);

						if (cumulative>roll&&!found) {
							LOG_DEBUG2("Pushing to stack");
							stack.push(i+1,j-2,1,0,0);
							regbp(ct,number,i,j-1);
							found=true;

							#ifndef force_continue
								break;
							#endif
						}
					}

					if (j!=1&&!lfce[i]&&!lfce[j]) {
         				// cumulative+= v->f(i+1,j-1) *data->eparam[10] * (data->eparam[6]*data->eparam[6]) * // [LX1N9X]
						// 	data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
						// 	*penalty(j-1,i+1,ct,data);
         				cumulative = SUM(cumulative, PROD(v->f(i+1,j-1), data->eparam[10], data->eparam[6], data->eparam[6], // [LX1N9X]
							data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]], penalty(j-1,i+1,ct,data)));

						LOG_DEBUG4("[LX1N9X]\tcumulative: " << cumulative);

						if (cumulative>roll&&!found) {
							LOG_DEBUG2("Pushing to stack");
							stack.push(i+1,j-1,1,0,0);
							found=true;

							#ifndef force_continue
								break;
							#endif
						}

						if ((mod[i+1]||mod[j-1])) if((j-2>0)&&!(fce->f(i+1,j-1)&SINGLE)) {
							if(inc[ct->numseq[i+2]][ct->numseq[j-2]]&&notgu(i+1,j-1,ct)) {

								// cumulative+= v->f(i+2,j-2) * data->eparam[10] * (data->eparam[6]*data->eparam[6]) *
         						// 	data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
								// 	*penalty(j-1,i+1,ct,data)*erg1(i+1,j-1,i+2,j-2,ct,data);
								cumulative = SUM(cumulative, PROD(v->f(i+2,j-2), data->eparam[10], data->eparam[6], data->eparam[6], 
         							data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]],
									penalty(j-1,i+1,ct,data), erg1(i+1,j-1,i+2,j-2,ct,data)));

								LOG_DEBUG4("[LX1N9X.2]\tcumulative: " << cumulative);

								if (cumulative>roll&&!found) {
									LOG_DEBUG2("Pushing to stack");
									stack.push(i+2,j-2,1,0,0);
									regbp(ct,number,i+1,j-1);
									found=true;

									#ifndef force_continue
										break;
									#endif
								}
							}
						}
					}
					#endif  //SIMPLEMBLOOP

					if (!lfce[i]) {
         				// cumulative+=  wl->f(i+1,j)*data->eparam[6]*data->scaling;
         				cumulative = SUM(cumulative, PROD(wl->f(i+1,j), data->eparam[6] SCALING));

						LOG_DEBUG4("[YJDEPA]\tcumulative: " << cumulative);

						if (cumulative>roll&&!found) {
							LOG_DEBUG2("Pushing to stack");
							stack.push(i+1,j,6,0,0);
							found=true;

							#ifndef force_continue
								break;
							#endif
						}
					}

					if (!found) {
						cout << "Traceback errors at wl\n";
						LOG_ERROR("Traceback errors at wl\t" << "i:  " << i << "\tj:  " << j << "\tCase:  " << switchcase << "\tCumulative:  " << cumulative/denominator);
						tracebackerror=14;
					}

					LOG_DEBUG4("denominator: " << denominator);

					break;

				case 7: //switchcase 7, a wl seed (helix only)
					// if (i==619 && j==622 && only_log_once){
					// 	SET_DEBUG_LEVEL(DEBUG4)
					// 	only_log_once = false;
					// }

					LOG_DEBUG("switchcase: " << switchcase << "\ti,j: " << i << ", " << j);

					// denominator=(wl->f(i,j)-wl->f(i+1,j)*data->eparam[6]*data->scaling);
					denominator=wlc->f(i,j);
					roll = PROD(roll, denominator);
					#ifdef SIMPLEMBLOOP

					// cumulative+= v->f(i,j-1)* data->eparam[10] * data->eparam[6] *
         			// 	erg4(j-1,i,j,1,ct,data,lfce[j])*penalty(i,j-1,ct,data);
					cumulative = SUM(cumulative, PROD(v->f(i,j-1), data->eparam[10], data->eparam[6], 
         				erg4(j-1,i,j,1,ct,data,lfce[j]), penalty(i,j-1,ct,data)));

					LOG_DEBUG4("[XXXXXX]\tcumulative: " << cumulative);

					if (cumulative>roll&&!found) {
						LOG_DEBUG2("Pushing to stack");
						stack.push(i,j-1,1,0,0);
						found=true;

						#ifndef force_continue
							break;
						#endif
					}

					if ((mod[i]||mod[j-1])) if(inc[ct->numseq[i+1]][ct->numseq[j-2]]&&notgu(i,j-1,ct)&&!(fce->f(i,j-1)&SINGLE)) {

						// cumulative+= v->f(i+1,j-2) * data->eparam[10] * data->eparam[6] *
         				// 	erg4(j-1,i,j,1,ct,data,lfce[j])*penalty(i,j-1,ct,data)
						// 	*erg1(i,j-1,i+1,j-2,ct,data);
						cumulative = SUM(cumulative, PROD(v->f(i+1,j-2), data->eparam[10], data->eparam[6],
         					erg4(j-1,i,j,1,ct,data,lfce[j]), penalty(i,j-1,ct,data), erg1(i,j-1,i+1,j-2,ct,data)));

						LOG_DEBUG4("[XXXXXX]\tcumulative: " << cumulative);

						if (cumulative>roll&&!found) {
							LOG_DEBUG2("Pushing to stack");
							stack.push(i+1,j-2,1,0,0);
							regbp(ct,number,i,j-1);
							found=true;

							#ifndef force_continue
								break;
							#endif
						}
					}
					#else   //!SIMPLEMBLOOP

					// cumulative+= data->eparam[10]*v->f(i,j)*penalty(j,i,ct,data);
					cumulative = SUM(cumulative, PROD(data->eparam[10], v->f(i,j), penalty(j,i,ct,data)));

					LOG_DEBUG4("[XXXXXX]\tcumulative: " << cumulative);

					if (cumulative>roll&&!found) {
						LOG_DEBUG2("Pushing to stack");
						found=true;
						stack.push(i,j,1,0,0);

						#ifndef force_continue
							break;
						#endif
					}

					if ((mod[i]||mod[j])) if (inc[ct->numseq[i+1]][ct->numseq[j-1]]&&notgu(i,j,ct)) {

						// cumulative+= data->eparam[10]*v->f(i+1,j-1)*penalty(j,i,ct,data)*erg1(i,j,i+1,j-1,ct,data);
						cumulative = SUM(cumulative, PROD(data->eparam[10], v->f(i+1,j-1), penalty(j,i,ct,data), erg1(i,j,i+1,j-1,ct,data)));

						LOG_DEBUG4("[XXXXXX]\tcumulative: " << cumulative);

						if (!found&&cumulative>roll) {
							LOG_DEBUG2("Pushing to stack");
							stack.push(i+1,j-1,1,0,0);
							regbp(ct,number,i,j);
							found = true;

							#ifndef force_continue
								break;
							#endif
						}
					}

					//[J1C6M9]
					// cumulative+= v->f(i+1,j)*data->eparam[10]*data->eparam[6]*
         			// 	erg4(j,i+1,i,2,ct,data,lfce[i])*penalty(i+1,j,ct,data);
					cumulative = SUM(cumulative, PROD(v->f(i+1,j), data->eparam[10], data->eparam[6], 
         				erg4(j,i+1,i,2,ct,data,lfce[i]), penalty(i+1,j,ct,data)));

					LOG_DEBUG4("[XXXXXX]\tcumulative: " << cumulative);

					if (!found&&cumulative>roll) {
						LOG_DEBUG2("Pushing to stack");
						stack.push(i+1,j,1,0,0);
						found=true;

						#ifndef force_continue
							break;
						#endif
					}

					if ((mod[i+1]||mod[j])) if(inc[ct->numseq[i+2]][ct->numseq[j-1]]&&notgu(i+1,j,ct)&&!(fce->f(i+1,j)&SINGLE)) {

						// cumulative+= v->f(i+2,j-1) * data->eparam[10] *data->eparam[6] *
         				// 	erg4(j,i+1,i,2,ct,data,lfce[i])*penalty(i+1,j,ct,data)
						// 	*erg1(i+1,j,i+2,j-1,ct,data);
						cumulative = SUM(cumulative, PROD(v->f(i+2,j-1), data->eparam[10], data->eparam[6],
         					erg4(j,i+1,i,2,ct,data,lfce[i]), penalty(i+1,j,ct,data),
							erg1(i+1,j,i+2,j-1,ct,data)));

						LOG_DEBUG4("[XXXXXX]\tcumulative: " << cumulative);

						if (!found&&cumulative>roll) {
							LOG_DEBUG2("Pushing to stack");
							stack.push(i+2,j-1,1,0,0);
							found=true;
							regbp(ct,number,i+1,j);

							#ifndef force_continue
								break;
							#endif
						}
					}

					// [LX1N9X]
					// cumulative+= v->f(i,j-1)* data->eparam[10] * data->eparam[6] *
         			// 	erg4(j-1,i,j,1,ct,data,lfce[j])*penalty(i,j-1,ct,data);
					cumulative = SUM(cumulative, PROD(v->f(i,j-1), data->eparam[10], data->eparam[6],
         				erg4(j-1,i,j,1,ct,data,lfce[j]), penalty(i,j-1,ct,data)));

					LOG_DEBUG4("[XXXXXX]\tcumulative: " << cumulative);

					if (cumulative>roll&&!found) {
						LOG_DEBUG2("Pushing to stack");
						stack.push(i,j-1,1,0,0);
						found=true;

						#ifndef force_continue
							break;
						#endif
					}

					if ((mod[i]||mod[j-1])) if(inc[ct->numseq[i+1]][ct->numseq[j-2]]&&notgu(i,j-1,ct)&&!(fce->f(i,j-1)&SINGLE)) {

						// cumulative+= v->f(i+1,j-2) * data->eparam[10] * data->eparam[6] *
         				// 	erg4(j-1,i,j,1,ct,data,lfce[j])*penalty(i,j-1,ct,data)
						// 	*erg1(i,j-1,i+1,j-2,ct,data);
						cumulative = SUM(cumulative, PROD(v->f(i+1,j-2), data->eparam[10], data->eparam[6],
         					erg4(j-1,i,j,1,ct,data,lfce[j]), penalty(i,j-1,ct,data), erg1(i,j-1,i+1,j-2,ct,data)));

						LOG_DEBUG4("[XXXXXX]\tcumulative: " << cumulative);

						if (cumulative>roll&&!found) {
							LOG_DEBUG2("Pushing to stack");
							stack.push(i+1,j-2,1,0,0);
							regbp(ct,number,i,j-1);
							found=true;

							#ifndef force_continue
								break;
							#endif
						}
					}

					if (j!=1&&!lfce[i]&&!lfce[j]) {
         				// cumulative+= v->f(i+1,j-1) *data->eparam[10] * (data->eparam[6]*data->eparam[6]) *  //[WERF93]
         				// data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
						// *penalty(j-1,i+1,ct,data);
         				cumulative = SUM(cumulative, PROD(v->f(i+1,j-1), data->eparam[10], data->eparam[6], data->eparam[6],
							data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]],
							penalty(j-1,i+1,ct,data)));

						LOG_DEBUG4("[XXXXXX]\tcumulative: " << cumulative);

						if (cumulative>roll&&!found) {
							LOG_DEBUG2("Pushing to stack");
							stack.push(i+1,j-1,1,0,0);
							found=true;

							#ifndef force_continue
								break;
							#endif
						}

						if ((mod[i+1]||mod[j-1])) if((j-2>0)&&!(fce->f(i+1,j-1)&SINGLE)) {
							if(inc[ct->numseq[i+2]][ct->numseq[j-2]]&&notgu(i+1,j-1,ct)) {

								// cumulative+= v->f(i+2,j-2) * data->eparam[10] * (data->eparam[6]*data->eparam[6]) *
         						// 	data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
								// 	*penalty(j-1,i+1,ct,data)*erg1(i+1,j-1,i+2,j-2,ct,data);
								cumulative = SUM(cumulative, PROD(v->f(i+2,j-2), data->eparam[10], data->eparam[6], data->eparam[6],
         							data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]],
									penalty(j-1,i+1,ct,data), erg1(i+1,j-1,i+2,j-2,ct,data)));

								LOG_DEBUG4("[XXXXXX]\tcumulative: " << cumulative);

								if (cumulative>roll&&!found) {
									LOG_DEBUG2("Pushing to stack");
									stack.push(i+2,j-2,1,0,0);
									regbp(ct,number,i+1,j-1);
									found=true;

									#ifndef force_continue
										break;
									#endif
								}

							}
						}
					}
					#endif //SIMPLEMBLOOP
					if (!found) {
						cout << "Traceback error at wl seed\n";
						LOG_ERROR("Traceback errors at wl seed\t" << "i:  " << i << "\tj:  " << j << "\tCase:  " << switchcase << "\tCumulative:  " << cumulative/denominator);
						tracebackerror=14;
					}

					break;
			} //switch(switchcase)

			LOG_DEBUG("denominator: " << denominator);

			// if (DIV(cumulative,denominator)>cumulative_max || (DIV(cumulative,denominator) < cumulative_min && switchcase!=0 && j!=0)) {
			if (cumulative>PROD(cumulative_max, denominator)) {
				// cout << "Over 1 probability error\n";
				//max_cumulative = max(max_cumulative,cumulative/denominator);
				LOG_WARN("i:  " << i << "\tj:  " << j << "\tCase:  " << switchcase << "\tCumulative:  " << DIV(cumulative,denominator));
				tracebackerror=21;
			}

			if (!found) {
				cout << "Overall traceback error\n";
				// LOG_ERROR("Overall traceback error\t" << "i:  " << i << "\tj:  " << j << "\tCase:  " << switchcase << "\tCumulative:  " << cumulative/denominator);
				tracebackerror=14;
			}

			SET_DEBUG_LEVEL(INFO);
			// SET_DEBUG_LEVEL(DEBUG2)
		} //while (stack.pull(&i,&j,&switchcase,&dummy1,&dummy2))


	} //for (number = 1; number <= numberofstructures; number++)

	//DEBUG_LOG(INFO, "\tmax_cumulative: " << max_cumulative);
	return tracebackerror;
}


//Stochastic() starts stochastic traceback from the point of reading a save file from disk.
//Stochastic() calls stochastictraceback() to do the actual work of sampling.
void stochastic(structure *ct, char *savefilename, int numberofstructures, int randomseed, ProgressHandler *progress) {


	short vers;

	DynProgArray<PFPRECISION> *w,*wmb,*wmbl,*wcoax,*wl,*wlc,*v;
	forceclass *fce;
	PFPRECISION *w3,*w5,scaling;
	bool *lfce,*mod;
	pfdatatable *data;
	datatable *data2;


	ifstream sav(savefilename,ios::binary);

	read(&sav,&vers);//read the save file version
		//right now there is no infrastructure to indicate an error in file version
		//that should be added at some point

	int sequencelength;
	read(&sav,&(sequencelength));

	sav.close();
	//allocate everything

	data = new pfdatatable;
	data2 = new datatable;

	ct->allocate(sequencelength);

	w = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
	v = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
	wmb = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
	fce = new forceclass(ct->GetSequenceLength());
	wl = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
	wlc = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
	wcoax = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
	wmbl = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());

	w5 = new PFPRECISION [ct->GetSequenceLength()+1];
	w3 = new PFPRECISION [ct->GetSequenceLength()+2];

	lfce = new bool [2*ct->GetSequenceLength()+1];
	mod = new bool [2*ct->GetSequenceLength()+1];




	//load all the data from the pfsavefile:
	readpfsave(savefilename, ct, w5, w3,v, w, wmb,wl, wlc, wmbl, wcoax, fce,&scaling,mod,lfce,data,data2);

	data->scaling = scaling;

	//now that partition function data was read from disk, do the sampling:
	stochastictraceback(w,wmb,wmbl,wcoax,wl,wlc,v,
		fce, w3,w5,scaling, lfce, mod, data, numberofstructures,
		ct, randomseed, progress);


	//delete everything
	delete data;
	delete data2;
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


}

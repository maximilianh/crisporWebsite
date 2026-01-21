/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005,2006
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 * Contributors: Arif Harmanci and Gaurav Sharma, 2006, 2007
 */

#include "dynalign.h"
#include <math.h>
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>

#include "algorithm.h"
#include "DynProgArray.h"
#include "defines.h"
#include "dynalignarray.h"
#include "dynalignheap.h"
#include "dynalignstackclass.h"
#include "forceclass.h"
#include "stackclass.h"
#include "../src/phmm/phmm_aln.h"
#include "../src/phmm/structure/structure_object.h"


#ifdef _WINDOWS_GUI
	#include "../Windows_interface_deprecated/platform.h"
#endif //_WINDOWS_GUI
#include "rna_library.h"
#include "structure.h"

#include "TProgressDialog.h"
#include "varray.h"
#include "wendarray.h"
//#include "align_services.h"

#ifdef COMPILE_SMP
#include "observingtextprogressbar.h"
#include "rankconsumer.h"
#include "rankmanager.h"
#include "rankproducer.h"
#endif //COMPILE_SMP
using namespace std;

//#define timer //flag to turn on a timer that writes a file
#undef timer

//maximum size of unpaired nucs on each side of an internal loop
#define maxloop 20

//specify the maximum size of an internal loop
#define maxinternal 30

#ifdef CHECK_ARRAY
extern string seq1, seq2;
#endif

#define gap_branch(x) iintercept+islope*x

//Determine an alignment envelope of nucleotide pairs that are allowed to be aligned, uses HMM
void calculate_coinc_probs_env(structure *ct1, structure *ct2, bool **allowed_alignments, short **forcealign) {

	vector<char>* seq1_nucs = new vector<char>();
    for(int i = 1; i <= ct1->GetSequenceLength(); i++)
    {
            seq1_nucs->push_back(ct1->nucs[i]);
    } // i loop.

    // Sequence2 nucleotides.
    vector<char>* seq2_nucs = new vector<char>();
    for(int j = 1; j <= ct2->GetSequenceLength(); j++)
    {
            seq2_nucs->push_back(ct2->nucs[j]);
    } // i loop.

    // Allocate the structures.
    t_structure* str1 = new t_structure("seq1", seq1_nucs);
    t_structure* str2 = new t_structure("seq2", seq2_nucs);

    // Allocate the structures.
    t_phmm_aln* phmm_aln = new t_phmm_aln(str1, str2);

    // Set the constraints.
    if(forcealign != NULL)
    {
            int* aln_constraints = new int[ct1->GetSequenceLength() + 2];
            for(int i = 1; i <= ct1->GetSequenceLength(); i++)
            {
                    aln_constraints[i] = forcealign[0][i];
            } // i loop.
            phmm_aln->set_aln_constraints(aln_constraints);
            delete [] aln_constraints;
    }

    // Get the alignment envelope.
    t_aln_env_result* aln_env_result = phmm_aln->compute_alignment_envelope(PROB_ALN_ENV, 0);

    // Set the allowed alignments.
    for(int i = 1; i <= ct1->GetSequenceLength(); i++)
    {
            for(int j = 1; j <= ct2->GetSequenceLength(); j++)
            {
                    if(aln_env_result->low_limits[i] <= j &&
                            aln_env_result->high_limits[i] >= j)
                    {
                            allowed_alignments[i][j] = true;
                    }
                    else
                    {
                            allowed_alignments[i][j] = false;
                    }
            } // j loop
    } // i loop.

// Free memory
    phmm_aln->free_aln_env_result(aln_env_result);
    delete(phmm_aln);
    delete(seq1_nucs);
    delete(seq2_nucs);
    delete(str1);
    delete(str2);

}

short int edangle5noforce(int i,int j,int ip,structure* ct,datatable* data) {
  return data->dangle[ct->numseq[j]][ct->numseq[i]][ct->numseq[ip]][2];
}

short int edangle3noforce(int i,int j,int ip,structure* ct,datatable* data) {
  return data->dangle[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][1];
}

short int edangle5force(int i,int j,int ip,structure* ct,datatable* data) {
  if (ct->fcedbl[ip]) {
    return DYNALIGN_INFINITY;
  } else {
    return data->dangle[ct->numseq[j]][ct->numseq[i]][ct->numseq[ip]][2];
  }
}

short int edangle3force(int i,int j,int ip,structure* ct,datatable* data) {
  if (ct->fcedbl[ip]) {
    return DYNALIGN_INFINITY;
  } else {
    return data->dangle[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][1];
  }
}

#ifndef DYNALIGN_II
#else

short int branch_1(DynProgArray<integersize> *single1_w, dynalignarray *w,int i, int j, int k, int l,int N,int c,short *lowend, short *highend, short iintercept, short islope){

  short int en1=DYNALIGN_INFINITY;

  if(c<N){
    if(k>=lowend[c+1]&&k<=highend[c+1]&&l>=lowend[j]&&l<=highend[j]){en1=min(en1,w->f(c+1,j,k,l)+single1_w->f(i,c)+gap_branch(c-i+1));/*cerr<<"c-i+1="<<c-i+1<<"\n";cerr<<"lala"<<gap_branch(c-i+1)<<"\n";*/}
    //  if(i-2==91&&j+1==203&&k-2==109&&l+2==195&&c==178)cerr<<"c<N 1 en1 "<<en1<<"\n";
    if(j<=N&&l>=lowend[c]&&l<=highend[c]&&k>=lowend[i]&&k<=highend[i])en1=min(en1,w->f(i,c,k,l)+single1_w->f(c+1,j)+gap_branch(j-c));
    //    if(i-2==91&&j+1==203&&k-2==109&&l+2==195&&c==178)cerr<<"c<N 2 en1 "<<en1<<"\n";
  }
  else if(l>=lowend[c]&&l<=highend[c]&&k>=lowend[i]&&k<=highend[i])en1=min(en1,w->f(i,c,k,l)+single1_w->f(c+1,j)+gap_branch(j-c));
  //  if(i-2==91&&j+1==203&&k-2==109&&l+2==195&&c==178)cerr<<"c>=N en1 "<<en1<<"\n";
  return en1;

}

short int branch_2(DynProgArray<integersize> *single2_w, dynalignarray *w,int i, int j, int k, int l,int N2,int d,short *lowend, short *highend , short iintercept, short islope){
  
  short int en1=DYNALIGN_INFINITY;
  if(d<N2){
    if(d+1>=lowend[i]&&d+1<=highend[i]&&l>=lowend[j]&&l<=highend[j])en1=min(en1,w->f(i,j,d+1,l)+single2_w->f(k,d)+gap_branch(d-k+1));
    //  if(i==28&&j==45&&k==28&&l==60&&en1==-132)cerr<<"1 d "<<d<<" en1 "<<en1<<"\n";
    // if(i==28&&j==45&&k==28&&l==60&&en1==-132)cerr<<"lowend[i] "<<lowend[i]<<" highend[i] "<<highend[i]<<" lowend[j] "<<lowend[j]<<" highend[j] "<<highend[j]<<"\n";
    if(l<=N2&&d>=lowend[j]&&d<=highend[j]&&k>=lowend[i]&&k<=highend[i])en1=min(en1,w->f(i,j,k,d)+single2_w->f(d+1,l)+gap_branch(l-d));
    // if(i==28&&j==45&&k==28&&l==60&&en1==-132)cerr<<"2 d "<<d<<" en1 "<<en1<<"\n";
    //  if(i==28&&j==45&&k==28&&l==60&&en1==-132)cerr<<"lowend[i] "<<lowend[i]<<" highend[i] "<<highend[i]<<" lowend[j] "<<lowend[j]<<" highend[j] "<<highend[j]<<"\n";

  }
  else if(d>=lowend[j]&&d<=highend[j]&&k>=lowend[i]&&k<=highend[i])en1=min(en1,w->f(i,j,k,d)+single2_w->f(d+1,l)+gap_branch(l-d));
  //  if(i==28&&j==45&&k==28&&l==60&&en1==-132)cerr<<"3 d "<<d<<" en1 "<<en1<<"\n";
  //  if(i==28&&j==45&&k==28&&l==60&&en1==-132)cerr<<"lowend[i] "<<lowend[i]<<" highend[i] "<<highend[i]<<" lowend[j] "<<lowend[j]<<" highend[j] "<<highend[j]<<"\n";

  return en1;


}

void trace_branch_1(int i,int j,int k,int l,int c,bool& found,structure *ct1,datatable *data,DynProgArray<integersize> *single1_v, DynProgArray<integersize> *single1_w, DynProgArray<integersize> *single1_wmb, 
		    bool *single1_lfce, forceclass *single1_fce, bool *single1_mod,dynalignarray *w,dynalignstackclass& stack,short *lowend,short *highend,short int energy,short int addition , short iintercept, short islope){
//  if(i-2==91&&j+1==203&&k-2==109&&l+2==195&&c==178)cerr<<"traceback energy "<<energy<<" addition "<<addition<<"\n";
  	if(c<ct1->GetSequenceLength()){
	  if(!found&&k>=lowend[c+1]&&k<=highend[c+1]&&l>=lowend[j]&&l<=highend[j]&&energy==w->f(c+1,j,k,l)+single1_w->f(i,c)+gap_branch(c-i+1)+addition){
	    //	    cout<<1<<' '<<i<<' '<<c<<"\n";
	    trace(ct1,data,i,c,single1_v,single1_w,single1_wmb,NULL,NULL,single1_lfce,single1_fce,NULL,NULL,single1_mod,NULL,single1_w->f(i,c),0,0);
	    found=true;
	    stack.push(c+1,j,k,l,w->f(c+1,j,k,l));
	  }
	  if(!found&&j<=ct1->GetSequenceLength()&&l>=lowend[c]&&l<=highend[c]&&k>=lowend[i]&&k<=highend[i]&&energy==w->f(i,c,k,l)+single1_w->f(c+1,j)+gap_branch(j-c)+addition){
	    //  cout<<1<<' '<<c+1<<' '<<j<<"\n";
	    trace(ct1,data,c+1,j,single1_v,single1_w,single1_wmb,NULL,NULL,single1_lfce,single1_fce,NULL,NULL,single1_mod,NULL,single1_w->f(c+1,j),0,0);
	    found=true;
	    stack.push(i,c,k,l,w->f(i,c,k,l));
	  }
	}
	else if(!found&&l>=lowend[c]&&l<=highend[c]&&k>=lowend[i]&&k<=highend[i]&&energy==w->f(i,c,k,l)+single1_w->f(c+1,j)+gap_branch(j-c)+addition){
	  //	  cout<<1<<' '<<c+1<<' '<<j<<"\n";
	  trace(ct1,data,c+1-ct1->GetSequenceLength(),j-ct1->GetSequenceLength(),single1_v,single1_w,single1_wmb,NULL,NULL,single1_lfce,single1_fce,NULL,NULL,single1_mod,NULL,single1_w->f(c+1,j),0,0);
	  found=true;
	  stack.push(i,c,k,l,w->f(i,c,k,l));
	}
}

void trace_branch_2(int i,int j,int k,int l,int d,bool& found,structure *ct2,datatable *data,DynProgArray<integersize> *single2_v, DynProgArray<integersize> *single2_w, DynProgArray<integersize> *single2_wmb, 
		    bool *single2_lfce, forceclass *single2_fce, bool *single2_mod,dynalignarray *w,dynalignstackclass& stack,short *lowend,short *highend,short int energy,short int addition , short iintercept, short islope){
  	if(d<ct2->GetSequenceLength()){
	  //	  if(i==28&&j==45&&k==28&&l==60&&energy==-132&&d==45&&d+1>=lowend[i]&&d+1<=highend[i]&&l>=lowend[j]&&l<=highend[j])cerr<<"d "<<d<<" 1 trace "<<w->f(i,j,d+1,l)+single2_w->f(k,d)+gap_branch(d-k+1)+addition<<"\n";
	  if(!found&&d+1>=lowend[i]&&d+1<=highend[i]&&l>=lowend[j]&&l<=highend[j]&&energy==w->f(i,j,d+1,l)+single2_w->f(k,d)+gap_branch(d-k+1)+addition){
	    //  cout<<2<<' '<<k<<' '<<d<<"\n";
	    trace(ct2,data,k,d,single2_v,single2_w,single2_wmb,NULL,NULL,single2_lfce,single2_fce,NULL,NULL,single2_mod,NULL,single2_w->f(k,d),0,0);
	    found=true;
	    stack.push(i,j,d+1,l,w->f(i,j,d+1,l));
	  }
	  //	  if(i==28&&j==45&&k==28&&l==60&&energy==-132&&d==45&&l<=ct2->GetSequenceLength()&&d>=lowend[j]&&d<=highend[j]&&k>=lowend[i]&&k<=lowend[i])cerr<<"d "<<d<<" 2 trace "<<w->f(i,j,k,d)+single2_w->f(d+1,l)+gap_branch(l-d)+addition<<"\n";
	  //	  if(i==28&&j==45&&k==28&&l==60&&energy==-132&&d==45)cerr<<"trace lowend[i] "<<lowend[i]<<" highend[i] "<<highend[i]<<" lowend[j] "<<lowend[j]<<" highend[j] "<<highend[j]<<" w "<<w->f(i,j,k,d)<<"\n";

	  if(!found&&l<=ct2->GetSequenceLength()&&d>=lowend[j]&&d<=highend[j]&&k>=lowend[i]&&k<=highend[i]&&energy==w->f(i,j,k,d)+single2_w->f(d+1,l)+gap_branch(l-d)+addition){
	    //  if(i==28&&j==45&&k==28&&l==60&&energy==-132&&d==45/*&&l<=ct2->GetSequenceLength()&&d>=lowend[j]&&d<=highend[j]&&k>=lowend[i]&&k<=lowend[i]*/)cerr<<"found!\n";
  //  cout<<2<<' '<<d+1<<' '<<l<<"\n";
	    trace(ct2,data,d+1,l,single2_v,single2_w,single2_wmb,NULL,NULL,single2_lfce,single2_fce,NULL,NULL,single2_mod,NULL,single2_w->f(d+1,l),0,0);
	    found=true;
	    stack.push(i,j,k,d,w->f(i,j,k,d));
	  }
	}
	else 
	  {
	    //   if(i==28&&j==45&&k==28&&l==60&&energy==-132&&d==45&&d>=lowend[j]&&d<=highend[j]&&k>=lowend[i]&&k<=highend[i])cerr<<"d "<<d<<" 3 trace "<<w->f(i,j,k,d)+single2_w->f(d+1,l)+gap_branch(l-d)+addition<<"\n";
if(!found&&d>=lowend[j]&&d<=highend[j]&&k>=lowend[i]&&k<=highend[i]&&energy==w->f(i,j,k,d)+single2_w->f(d+1,l)+gap_branch(l-d)+addition){
	  //	  cout<<2<<' '<<d+1<<' '<<l<<"\n";
	  trace(ct2,data,d+1-ct2->GetSequenceLength(),l-ct2->GetSequenceLength(),single2_v,single2_w,single2_wmb,NULL,NULL,single2_lfce,single2_fce,NULL,NULL,single2_mod,NULL,single2_w->f(d+1,l),0,0);
	  found=true;
	  stack.push(i,j,k,d,w->f(i,j,k,d));
	}
	  }
}
//report to Dave9 can i directly use i,j,k,l? or deref?
#endif

// Return the most 5' nucleotide in sequence 2 that can be aligned to
// nucleotide i in sequence 1, using the M constraint:
short lowlimit(short i, short M, short N1, short N2) {
  if (i<=N1) {
    return i*N2/N1-M;
  } else {
    return (i-N1)*N2/N1+N2-M;
  }
}
//overload to use a bool array of allowed alignments:
short lowlimit(short i, bool **allowed_alignments, short N1, short N2) {
	int j;

	if (i==0) return 0;

	if (i<=N1) {

		for (j=1;j<=N2;j++) {
			if (allowed_alignments[i][j]) {
				return j;
			}

		}
		return min(i,N2);//return something for program execution
	}
	else {
		for (j=1;j<=N2;j++) {
			if (allowed_alignments[i-N1][j]) return j+N2;
		}
		return min(i+N1,2*N2);
	}



}

// Return the most 3' nucleotide in sequence 2 that can be aligned to
// nucleotide i in sequence 1, using the M constraint:
short highlimit(short i, short M, short N1, short N2) {
  if (i<=N1) {
    return i*N2/N1+M;
  } else {
    return (i-N1)*N2/N1+N2+M;
  }
}
//overload to use a bool array of allowed alignments:
short highlimit(short i, bool **allowed_alignments, short N1, short N2) {
	int j;

    if (i==0) return N2;

	if (i<=N1) {
		for (j=N2;j>0;j--) {
			if (allowed_alignments[i][j]) {
				return j;
			}
		}
		return min(i,N2);//return something for program execution
	}
	else {
		for (j=N2;j>0;j--) {
			if (allowed_alignments[i-N1][j]) return j+N2;
		}
		return min(i+N1,2*N2);//return something for program execution

	}
}


//short reference(short i, short k, short N, short N2, short M) {
//  return k-i*N2/N+M;
//}


void dynalignfceunpaired(structure *ct,char **fce,int nopair) {

  int i;

  for (i=nopair+1;i<nopair+(ct->GetSequenceLength());++i) {
    fce[jref(nopair,i,ct->GetSequenceLength())][iref(nopair,i,ct->GetSequenceLength())]=fce[jref(nopair,i,ct->GetSequenceLength())][iref(nopair,i,ct->GetSequenceLength())]|SINGLE;
  }
  for (i=1;i<nopair;++i) {
    fce[nopair][i]=fce[nopair][i]|SINGLE;
  }
  for (i=nopair+1;i<=ct->GetSequenceLength();++i) {
    fce[jref(i,nopair+ct->GetSequenceLength(),ct->GetSequenceLength())][iref(i,nopair+ct->GetSequenceLength(),ct->GetSequenceLength())]=
      fce[jref(i,nopair+ct->GetSequenceLength(),ct->GetSequenceLength())][iref(i,nopair+ct->GetSequenceLength(),ct->GetSequenceLength())]|SINGLE;
  }

}

void dynforcedbl(int dbl,structure* ct,char **fce,bool *lineardbl) {
  int i,j;

  lineardbl[dbl]=true;
  lineardbl[dbl+ct->GetSequenceLength()]=true;

  for(i=dbl+1;i<=ct->GetSequenceLength();++i) {
    for (j=1;j<dbl;++j) {
      fce[i][j] = fce[i][j]|DUBLE;
    }
  }
  for(j=(dbl+(ct->GetSequenceLength())-1);j>ct->GetSequenceLength();j--) {
    for (i=dbl+1;i<=ct->GetSequenceLength();++i) {
      fce[jref(i,j,ct->GetSequenceLength())][iref(i,j,ct->GetSequenceLength())] = fce[jref(i,j,ct->GetSequenceLength())][iref(i,j,ct->GetSequenceLength())]|DUBLE;
    }
  }
}

void dynforcepair(int x,int y,structure *ct,char **v) {
  int i,j;
  //v->f(x,y) = v->f(x,y)|PAIR;
  //v->f(y,x+ct->GetSequenceLength())=v->f(y,x+ct->GetSequenceLength())|PAIR;
  for (i=y+1;i<=x-1+ct->GetSequenceLength();++i) {
    v[jref(x,i,ct->GetSequenceLength())][iref(x,i,ct->GetSequenceLength())] = v[jref(x,i,ct->GetSequenceLength())][iref(x,i,ct->GetSequenceLength())]|NOPAIR;
  }
  for (i=x;i<=y-1;++i) {
    v[jref(x,i,ct->GetSequenceLength())][iref(x,i,ct->GetSequenceLength())] = v[jref(x,i,ct->GetSequenceLength())][iref(x,i,ct->GetSequenceLength())]|NOPAIR;
  }
  for (i=1;i<=x-1;++i) {
    v[jref(i,y,ct->GetSequenceLength())][iref(i,y,ct->GetSequenceLength())] = v[jref(i,y,ct->GetSequenceLength())][iref(i,y,ct->GetSequenceLength())]|NOPAIR;
  }
  for (i=x+1;i<=y;++i) {
    v[jref(i,y,ct->GetSequenceLength())][iref(i,y,ct->GetSequenceLength())] = v[jref(i,y,ct->GetSequenceLength())][iref(i,y,ct->GetSequenceLength())]|NOPAIR;
  }
  for (i=1;i<=x-1;++i) {
    v[jref(i,x,ct->GetSequenceLength())][iref(i,x,ct->GetSequenceLength())] = v[jref(i,x,ct->GetSequenceLength())][iref(i,x,ct->GetSequenceLength())]|NOPAIR;
  }
  for (i=y+1;i<=ct->GetSequenceLength();++i) {
    v[jref(i,y+ct->GetSequenceLength(),ct->GetSequenceLength())][iref(i,y+ct->GetSequenceLength(),ct->GetSequenceLength())]=v[jref(i,y+ct->GetSequenceLength(),ct->GetSequenceLength())][iref(i,y+ct->GetSequenceLength(),ct->GetSequenceLength())]|NOPAIR;
  }
  for (i=y;i<=x-1+(ct->GetSequenceLength());++i) {
    v[jref(y,i,ct->GetSequenceLength())][iref(y,i,ct->GetSequenceLength())] = v[jref(y,i,ct->GetSequenceLength())][iref(y,i,ct->GetSequenceLength())]|NOPAIR;
  }
  for (i=(ct->GetSequenceLength())+x+1;i<=(ct->GetSequenceLength())+y-1;++i) {
    v[jref(y,i,ct->GetSequenceLength())][iref(y,i,ct->GetSequenceLength())] = v[jref(y,i,ct->GetSequenceLength())][iref(y,i,ct->GetSequenceLength())]|NOPAIR;
  }
  for (i=x+1;i<=y-1;++i) {
    v[jref(i,x+ct->GetSequenceLength(),ct->GetSequenceLength())][iref(i,x+ct->GetSequenceLength(),ct->GetSequenceLength())] = v[jref(i,x+ct->GetSequenceLength(),ct->GetSequenceLength())][iref(i,x+ct->GetSequenceLength(),ct->GetSequenceLength())]|NOPAIR;
  }
  for (i=y+1;i<=ct->GetSequenceLength();++i) {
    v[jref(i,x+ct->GetSequenceLength(),ct->GetSequenceLength())][iref(i,x+ct->GetSequenceLength(),ct->GetSequenceLength())] = v[jref(i,x+ct->GetSequenceLength(),ct->GetSequenceLength())][iref(i,x+ct->GetSequenceLength(),ct->GetSequenceLength())]|NOPAIR;
  }
  for (i=1;i<=x-1;++i) {
    for (j = x+1;j<=y-1;++j){
      v[jref(i,j,ct->GetSequenceLength())][iref(i,j,ct->GetSequenceLength())] = v[jref(i,j,ct->GetSequenceLength())][iref(i,j,ct->GetSequenceLength())] |NOPAIR;
    }
  }
  for (i=x+1;i<=y-1;++i) {
    for (j=y+1;j<=(ct->GetSequenceLength())+x-1;++j) {
      v[jref(i,j,ct->GetSequenceLength())][iref(i,j,ct->GetSequenceLength())] = v[jref(i,j,ct->GetSequenceLength())][iref(i,j,ct->GetSequenceLength())]|NOPAIR;
    }
  }
  for (i=y+1;i<=ct->GetSequenceLength();++i) {
    for (j=(ct->GetSequenceLength())+x+1;j<=(ct->GetSequenceLength())+y-1;++j) {
      v[jref(i,j,ct->GetSequenceLength())][iref(i,j,ct->GetSequenceLength())] = v[jref(i,j,ct->GetSequenceLength())][iref(i,j,ct->GetSequenceLength())]|NOPAIR;
    }
  }
}

void dynforcepairg(int x,structure *ct,char **fce) {
  int i,j;


  for (j=x+1;j<x+ct->GetSequenceLength();++j) {
    if (ct->numseq[j]!=3) fce[jref(x,j,ct->GetSequenceLength())][iref(x,j,ct->GetSequenceLength())]=fce[jref(x,j,ct->GetSequenceLength())][iref(x,j,ct->GetSequenceLength())]|NOPAIR;
  }

  for (j=x+ct->GetSequenceLength()+1;j<2*ct->GetSequenceLength();++j) {
    if (ct->numseq[j]!=3) fce[jref(x+ct->GetSequenceLength(),j,ct->GetSequenceLength())][iref(x+ct->GetSequenceLength(),j,ct->GetSequenceLength())]=fce[jref(x+ct->GetSequenceLength(),j,ct->GetSequenceLength())][iref(x+ct->GetSequenceLength(),j,ct->GetSequenceLength())]|NOPAIR;
  }

  for (i=x-1;i>0;--i) {
    if (ct->numseq[i]!=3) fce[jref(i,x,ct->GetSequenceLength())][iref(i,x,ct->GetSequenceLength())]=fce[jref(i,x,ct->GetSequenceLength())][iref(i,x,ct->GetSequenceLength())]|NOPAIR;
  }

  for (i=x+ct->GetSequenceLength()-1;i>x;--i) {
    if (ct->numseq[i]!=3) fce[jref(i,x+ct->GetSequenceLength(),ct->GetSequenceLength())][iref(i,x+ct->GetSequenceLength(),ct->GetSequenceLength())]=fce[jref(i,x+ct->GetSequenceLength(),ct->GetSequenceLength())][iref(i,x+ct->GetSequenceLength(),ct->GetSequenceLength())]|NOPAIR;
  }
}

void dynalignforce(structure *ct1, structure *ct2/*, short ****fce*/, char **fce1, char **fce2, bool *dbl1,
                   bool *dbl2, bool *mod1, bool *mod2) {
  int count;
  //fill the fce array for a dynalign calculation using the following key:

  //SINGLE applies to any fce(i,j) s.t. i or j should be single stranded
  //PAIR applies to any fce(i,j) where i is paired to j
  //NOPAIR applies to any fce(i,j) where either i or j is paired to
  //    another nucleotide or i and j are forbidden to pair
  //DUBLE applies to any fce(i.j) where an nuc, k, i<k<j is double stranded
  //INTER applies to any fce(i,j) such that some k, i<k<j, is the virtual linker
  //    used for intermolecular folding

  //refer to fce[j][i][a][b]
  //where a = k-i+maxsep and i+maxsep>=k>=i-maxsep and i<j
  //and for l<N and j<N2 b = l-j+maxsep and j+maxsep>=l>=j-maxsep
  //for l>N and j>N2, b = l-j+maxsep+N-N2 so that maxsep is "reset" beyond the ends of the sequences when they are of different lengths


  //start with unpaired nucs in seq1
  for (count=0;count<ct1->GetNumberofSingles();++count) {
    //dynalignforceunpaired1(ct1, ct2, fce, maxseparation, ct1->GetSingle(count]);
      dynalignfceunpaired(ct1,fce1,ct1->GetSingle(count));

  }

  //now encode unpaired nucs in sequence 2
  for (count=0;count<ct2->GetNumberofSingles();++count) {
    //dynalignforceunpaired2(ct1,ct2,fce,maxseparation,ct2->GetSingle(count]);
      dynalignfceunpaired(ct2,fce2,ct2->GetSingle(count));

  }

  //force nucleotides double-stranded in seq1
  for(count=0;count<ct1->GetNumberofDoubles();++count) {
      dynforcedbl(ct1->GetDouble(count),ct1,fce1,dbl1);
  }

  //force nucleotides double-stranded in seq2
  for(count=0;count<ct2->GetNumberofDoubles();++count) {
      dynforcedbl(ct2->GetDouble(count),ct2,fce2,dbl2);
  }


  //now handle pairs that are required:
  for (count = 0; count <ct1->GetNumberofPairs(); ++count) {

    dynforcepair(ct1->GetPair5(count),ct1->GetPair3(count),ct1,fce1);
    dynforcedbl(ct1->GetPair5(count),ct1,fce1,dbl1);
    dynforcedbl(ct1->GetPair3(count),ct1,fce1,dbl1);


  }

  for (count = 0; count <ct2->GetNumberofPairs();++count) {
    dynforcepair(ct2->GetPair5(count),ct2->GetPair3(count),ct2,fce2);
    dynforcedbl(ct2->GetPair5(count),ct2,fce2,dbl2);
    dynforcedbl(ct2->GetPair3(count),ct2,fce2,dbl2);

  }


  //now handle FMN cleavage (U in GU pair)
  for (count=0;count<ct1->GetNumberofGU();++count) {
      dynforcedbl(ct1->GetGUpair(count),ct1,fce1,dbl1);
    dynforcepairg(ct1->GetGUpair(count),ct1,fce1);

  }

  for (count=0;count<ct2->GetNumberofGU();++count) {
      dynforcedbl(ct2->GetGUpair(count),ct2,fce2,dbl2);
    dynforcepairg(ct2->GetGUpair(count),ct2,fce2);

  }

  //now handle prohibited base pairs
  for (count=0;count<ct1->GetNumberofForbiddenPairs();++count) {
    fce1[jref(ct1->GetForbiddenPair5(count),ct1->GetForbiddenPair3(count),ct1->GetSequenceLength())][iref(ct1->GetForbiddenPair5(count),ct1->GetForbiddenPair3(count),ct1->GetSequenceLength())]=
      fce1[jref(ct1->GetForbiddenPair5(count),ct1->GetForbiddenPair3(count),ct1->GetSequenceLength())][iref(ct1->GetForbiddenPair5(count),ct1->GetForbiddenPair3(count),ct1->GetSequenceLength())]|NOPAIR;

    fce1[jref(ct1->GetForbiddenPair3(count),ct1->GetForbiddenPair5(count)+ct1->GetSequenceLength(),ct1->GetSequenceLength())][iref(ct1->GetForbiddenPair3(count),ct1->GetForbiddenPair5(count)+ct1->GetSequenceLength(),ct1->GetSequenceLength())]=
      fce1[jref(ct1->GetForbiddenPair3(count),ct1->GetForbiddenPair5(count)+ct1->GetSequenceLength(),ct1->GetSequenceLength())][iref(ct1->GetForbiddenPair3(count),ct1->GetForbiddenPair5(count)+ct1->GetSequenceLength(),ct1->GetSequenceLength())]|NOPAIR;
  }

  for (count=0;count<ct2->GetNumberofForbiddenPairs();++count) {
    fce2[jref(ct2->GetForbiddenPair5(count),ct2->GetForbiddenPair3(count),ct2->GetSequenceLength())][iref(ct2->GetForbiddenPair5(count),ct2->GetForbiddenPair3(count),ct2->GetSequenceLength())]=
      fce1[jref(ct2->GetForbiddenPair5(count),ct2->GetForbiddenPair3(count),ct2->GetSequenceLength())][iref(ct2->GetForbiddenPair5(count),ct2->GetForbiddenPair3(count),ct2->GetSequenceLength())]|NOPAIR;

    fce2[jref(ct2->GetForbiddenPair3(count),ct2->GetForbiddenPair5(count)+ct2->GetSequenceLength(),ct2->GetSequenceLength())][iref(ct2->GetForbiddenPair3(count),ct2->GetForbiddenPair5(count)+ct2->GetSequenceLength(),ct2->GetSequenceLength())]=
      fce1[jref(ct2->GetForbiddenPair3(count),ct2->GetForbiddenPair5(count)+ct2->GetSequenceLength(),ct2->GetSequenceLength())][iref(ct2->GetForbiddenPair3(count),ct2->GetForbiddenPair5(count)+ct2->GetSequenceLength(),ct2->GetSequenceLength())]|NOPAIR;
  }

  //now handle chemical modification
  for (count=0;count<ct1->GetNumberofModified();++count) {

    if (ct1->GetModified(count)!=1&&ct1->GetModified(count)!=ct1->GetSequenceLength()) {
      mod1[ct1->GetModified(count)]=true;
      mod1[ct1->GetModified(count)+ct1->GetSequenceLength()]=true;
    }
  }
  for (count=0;count<ct2->GetNumberofModified();++count) {

    if (ct2->GetModified(count)!=1&&ct2->GetModified(count)!=ct2->GetSequenceLength()) {
      mod2[ct2->GetModified(count)]=true;
      mod2[ct2->GetModified(count)+ct2->GetSequenceLength()]=true;
    }
  }
}


//maxseparation is the M parameter that limits the alignment space
// maxseparation < 0 implies that the alignment constraints are instead coming from allowed_alignments
//return an int that indicates whether an error occurred.  0 = no error
#ifdef DYNALIGN_II
int dynalign(structure *ct1, structure *ct2, short **alignment,
	     short int maxseparation, short int islope, short int iintercept, short int gapincrease, datatable *data,
	     short maxtracebacks, short window, short awindow, short percentsort, short **forcealign,int max_elongation,
	     bool **allowed_alignments, ProgressHandler *progress, const char *Savefile, bool energyonly,
	     bool local, bool forced, short int numProcessors)
#else
int dynalign(structure *ct1, structure *ct2, short **alignment,
	     short int maxseparation, short int gapincrease, datatable *data,
	     bool singleinsert, short maxtracebacks, short window, short awindow, short percentsort, short **forcealign,
	     bool **allowed_alignments, ProgressHandler *progress, const char *Savefile, bool energyonly,
	     bool local, bool forced, short int numProcessors)
#endif

 {
  
	//cout << "maxtrace:\t" << maxtracebacks << endl;
	//cout << "bpwin:\t" << window << endl;
	//cout << "awin:\t" << awindow << endl;
	//cout << "percent:\t" << percentsort << endl;
	//cout << "imaxseparation:\t" << maxseparation << endl;
	//cout << "gap:\t" << gapincrease << endl;
	//cout << "singleinsert:\t" << singleinsert << endl;
	//cout << "energyonly:\t" << energyonly << endl;
	//cout << "forced:\t" << forced << endl;
	//cout << "local:\t" << local << endl;
	//cout << "numProcessors:\t" << numProcessors << endl;


  short int i,j,k,l,a,b,c,d,/*maxsep,*/gap,N,N2,Ndiff;
  short ip,kp,jp,lp,dp;
  integersize en1;
  short lowest;
  short int I;
  int flag;
  bool ikincrement,jldecrement,jldecrement2;
  integersize imin, jmax, kmin, lmax, wval;
  int error;
  //  int stacking_elongation_1,stacking_elongation_2;

//#ifdef COMPILE_SMP
//  bool removeprogress;
//#endif


  dynalignarray *w;
  varray *v;

  wendarray *w5,*w3;

  vector< vector<bool> > inc = data->pairing;

  dynalignarray *vmod;


  short *lowend,*highend;//store the limits calculated from lowlimit and highlimit

  // Select the version of edangle to use

  short int (*edangle5)(int i, int j, int ip, structure *ct, datatable *data) = edangle5noforce;
  short int (*edangle3)(int i, int j, int ip, structure *ct, datatable *data) = edangle3noforce;

  if (forced) {
    edangle5 = &edangle5force;
    edangle3 = &edangle3force;
  }

  // Local symbols edangle3 and edangle5 are now function pointers to
  // actual functions edangle3noforce and edangle5noforce or
  // edangle3force and edangle5force, respectively, if force is true

#ifdef timer
#include <time.h>
  ofstream timeout;
  int seconds;
  char timerstring[100];
  char timelength[10];
  strcpy(timerstring,"time_pf_");
  sprintf(timelength,"%i",ct1->GetSequenceLength());
  strcat(timerstring,timelength);
  strcat(timerstring,".out");

  timeout.open(timerstring);
  timeout<<time(NULL)<<"\n";
  seconds = time(NULL);
#endif //timer


  //here is a 4-d array to store data about folding constraints that involve two sequences
  //dynalignarray *fce;
  //short ****fce;

  //fce1 and fce2 contain configuration constraints for single sequences
  char **fce1;
  char **fce2;

  //dbl1 and dbl2 flag nucleotides that must be double-stranded in seq1 and seq2, respectively
  //false means a nucleotide is not forced double
  //true indicates a nucleotide is forced double
  bool *dbl1,*dbl2;

  //mod indicates whether a nucleotide is chemically modified
  //mod1 for seq1 and mod2 or seq2
  bool *mod1,*mod2;


  //the following are used for chemical modfication cases
  bool modification;
  bool alignmentforced;


  short crit;
  crit = DYNALIGN_INFINITY;

  //store the number of bases in ct1 and ct2 in register shorts
  N = ct1->GetSequenceLength();//the length of sequence 1
  N2 = ct2->GetSequenceLength();//the length of sequence 2
  Ndiff = N-N2;//the difference in sequence lengths
  //maxsep = maxseparation;



  lowest = DYNALIGN_INFINITY;





  //double up the sequences as part of suboptimal traceback support:
  for (i=1;i<=N;i++) {
    ct1->numseq[(N)+i] = ct1->numseq[i];
  }
  for (i=1;i<ct2->GetSequenceLength();i++) ct2->numseq[ct2->GetSequenceLength()+i]=ct2->numseq[i];


  //fill the low and highend arrays:
  //allocate space:
  lowend = new short [2*N];
  highend = new short [2*N];

  if (maxseparation>=0) {
	//Using the traditional M parameter to constrain the alignment space
	for (i=0;i<2*N;++i) {
		lowend[i]=lowlimit(i,maxseparation,N,N2);
		highend[i]=highlimit(i,maxseparation,N,N2);

		//printf("%d -> (%d, %d)\n", i, lowend[i], highend[i]);
	}
  }
  else {
	//allowed_alignments must be defined to constrain the alignment space
    //need to determine lowend and highend
	  for (i=0;i<2*N;++i) {
		lowend[i]=lowlimit(i,allowed_alignments,N,N2);
		highend[i]=highlimit(i,allowed_alignments,N,N2);

		//printf("%d -> (%d, %d)\n", i, lowend[i], highend[i]);
	  }


  }

  //Allocate the 2-d and 4-d arrays for storing mfe's for subfragments:
  w5 = new wendarray(N,N2,lowend,highend);
  w3 = new wendarray(N,N2,lowend,highend);

  v = new varray(N,N2,lowend,highend,ct1->tem,energyonly);
  w = new dynalignarray(N,N2,lowend,highend,energyonly);
  //report to Dave1
  
  // structure *ct=new structure();
  // openseq(ct,argv[1]);

#ifdef DYNALIGN_II
  int single1_vmin=0;
  int single2_vmin=0;
  integersize *single1_w5,*single1_w3,*single2_w5,*single2_w3;
  bool *single1_lfce,*single2_lfce,*single1_mod,*single2_mod;
  // DynProgArray<integersize> *single1_w2,*single1_wmb2,*single2_w2,*single2_wmb2;
  
  //allocate everything
  DynProgArray<integersize> single1_w(ct1->GetSequenceLength()),single2_w(ct2->GetSequenceLength());
  DynProgArray<integersize> single1_v(ct1->GetSequenceLength()),single2_v(ct2->GetSequenceLength());
  DynProgArray<integersize> single1_wmb(ct1->GetSequenceLength()),single2_wmb(ct2->GetSequenceLength());
  DynProgArray<integersize> single1_we(ct1->GetSequenceLength(),0),single2_we(ct2->GetSequenceLength(),0);
  forceclass single1_fce(ct1->GetSequenceLength()),single2_fce(ct2->GetSequenceLength());


  single1_lfce = new bool [2*ct1->GetSequenceLength()+1];
  single2_lfce = new bool [2*ct2->GetSequenceLength()+1];
  single1_mod = new bool [2*ct1->GetSequenceLength()+1];
  single2_mod = new bool [2*ct2->GetSequenceLength()+1];

  for (size_t i=0;i<=2*ct1->GetSequenceLength();i++) {
    single1_lfce[i] = false;
    single1_mod[i] = false;
  }
  for (size_t i=0;i<=2*ct2->GetSequenceLength();i++) {
    single2_lfce[i] = false;
    single2_mod[i] = false;
  }
 
  single1_w5 = new integersize [ct1->GetSequenceLength()+1];
  single1_w3 = new integersize [ct1->GetSequenceLength()+2];
  

  for (size_t i=0;i<=ct1->GetSequenceLength();i++) {
    single1_w5[i] = 0;
    single1_w3[i] = 0;

  }
  single1_w3[ct1->GetSequenceLength()+1] = 0;


  single2_w5 = new integersize [ct2->GetSequenceLength()+1];
  single2_w3 = new integersize [ct2->GetSequenceLength()+2];
  

  for (size_t i=0;i<=ct2->GetSequenceLength();i++) {
    single2_w5[i] = 0;
    single2_w3[i] = 0;

  }
  single2_w3[ct2->GetSequenceLength()+1] = 0;

  /*
  if (ct->intermolecular) {
    w2 = new DynProgArray<integersize>(ct->GetSequenceLength());
    wmb2 = new DynProgArray<integersize>(ct->GetSequenceLength());



    }*/
  
  force(ct1,&single1_fce,single1_lfce);
  force(ct2,&single2_fce,single2_lfce);

  for (i=1;i<=ct1->GetNumberofModified();i++) {
      if (ct1->GetModified(i)>1&&ct1->GetModified(i)<ct1->GetSequenceLength()) {
          single1_mod[ct1->GetModified(i)]=true;
          single1_mod[ct1->GetModified(i)+ct1->GetSequenceLength()]=true;
    }
  }
  
  for (i=1;i<=ct2->GetNumberofModified();i++) {
      if (ct2->GetModified(i)>1&&ct2->GetModified(i)<ct2->GetSequenceLength()) {
          single2_mod[ct2->GetModified(i)]=true;
          single2_mod[ct2->GetModified(i)+ct2->GetSequenceLength()]=true;
    }
  }

  

  fill(ct1,single1_v,single1_w,single1_wmb,single1_fce,single1_vmin,single1_lfce,single1_mod,single1_w5,single1_w3,0,data,NULL,NULL,&single1_we);
  fill(ct2,single2_v,single2_w,single2_wmb,single2_fce,single2_vmin,single2_lfce,single2_mod,single2_w5,single2_w3,0,data,NULL,NULL,&single2_we);
 
#else
#endif
  //report to Dave1
  if (forced) {

    //allocate the fce array for constraints
    //The following definitions are bitwise applied to the fce array:
    //SINGLE applies to any fce(i,j) s.t. i or j should be single stranded
    //PAIR applies to any fce(i,j) where i is paired to j
    //NOPAIR applies to any fce(i,j) where either i or j is paired to
    //  another nucleotide or i and j are forbidden to pair
    //DUBLE applies to any fce(i.j) where an nuc, k, i<k<j is double stranded
    //INTER applies to any fce(i,j) such that some k, i<k<j, is the virtual linker
    //  used for intermolecular folding
    //INTER applies to any fce(i,i) s.t. i is the center of the virtual linker
    //The above terms are defined in define.h

    //note that INTER is not currently supported in Dynalign
    //fce = new dynalignarray(N,N2,maxseparation);


    //allocate space for single sequence configuration in fce1 and fce2

    fce1 = new char *[2*N+1];

    for (i=0;i<=2*N;++i) {
      if (i<=N) {
        fce1[i] = new char [i+1];
        for (j=0;j<=i;++j) fce1[i][j]=0;
      }
      else {
        fce1[i] = new char [2*N-i+1];
        for (j=0;j<=2*N-i;++j) fce1[i][j]=0;
      }
    }

    fce2 = new char *[2*N2+1];

    for (i=0;i<=2*N2;++i) {
      if (i<=N2) {
        fce2[i] = new char [i+1];
        for (j=0;j<=i;++j) fce2[i][j]=0;
      }
      else {
        fce2[i] = new char [2*N2-i+1];
        for (j=0;j<=2*N2-i;++j) fce2[i][j]=0;
      }
    }

    dbl1 = new bool [2*N];
    mod1=new bool [2*N];
    for (i=0;i<2*N;++i) {
      dbl1[i]=false;
      mod1[i]=false;
    }
    dbl2 = new bool [2*N2];
    mod2 = new bool [2*N2];
    for (i=0;i<2*N2;++i) {
      dbl2[i]=false;
      mod2[i]=false;
    }

    //store the locations of dbl1 and dbl2 in ct1 and ct2, respectively
    //this facilitates checks done in functions edangle5 and edangle3
    ct1->fcedbl=dbl1;
    ct2->fcedbl=dbl2;

    //now assign values to fce according to the above key using the function dynalignforce
    dynalignforce(ct1,ct2/*,fce*/,fce1,fce2,dbl1,dbl2,mod1,mod2);

    if (ct1->GetNumberofModified()>0||ct2->GetNumberofModified()>0) {
      modification = true;
      //For chemical modification, a second v array, vmod is required
      vmod = new dynalignarray(N,N2,lowend,highend);

    }
    else {
      modification =false;
      vmod = NULL;
    }

    alignmentforced = (forcealign != NULL);

  } else {
    vmod = NULL;
	dbl1=NULL;
	dbl2=NULL;
	alignmentforced=false;
	fce1=NULL;fce2=NULL;
	modification=false;
	mod1=NULL;
	mod2=NULL;
  }

  //maxsep = maxseparation; //place maxseparation into a register int
  gap = gapincrease;//place the gap penalty in a register short


  // look for alignment of every i paired to j in sequence 1 to k
  // paired to l in sequence 2

  //for (size=minloop;size<=ct1->GetSequenceLength();size++) {

#ifndef COMPILE_SMP
  for (j = minloop; j <= N; ++j) {
    if (progress != NULL) {
      if (!energyonly) progress->update((100 * j) / (2*N));
	  else progress->update((100 * j) / (N));
    if (progress->canceled()) break; //exit if the user canceled the operation.
    }

    for (i = min(N,j-1); i >= 1; --i) {
      kmin = max(lowend[i],1);
      for (k = min(highend[i], N2); k >= kmin; k--) {
        lmax = min(highend[j], ct2->GetSequenceLength());
        for (l = max(max(lowend[j],1),k+minloop); l <= lmax; ++l) {
          if (((j - i) > minloop) && ((l - k) > minloop)) {
	    // cerr<<"infinite!\n"<<i<<" "<<j<<" "<<k<<" "<<j<<"\n";
//              if(fce1 == NULL) cerr << "fce1null" << endl;
#ifdef DYNALIGN_II        
            dynalignstep(ct1, ct2, data,
                         v, w, w5, w3,
                         vmod, &single1_w,&single2_w,mod1, mod2, modification,
                         fce1, fce2,
                         alignmentforced, forcealign,
                         dbl1, dbl2,max_elongation,
                         i, j, k, l, N, N2,
                         lowend, highend, islope, iintercept, gap,
                         forced,local,
                         edangle5, edangle3);
#else
                dynalignstep(ct1, ct2, data,
                         v, w, w5, w3,
                         vmod, mod1, mod2, modification,
                         fce1, fce2,
                         alignmentforced, forcealign,
                         dbl1, dbl2,
                         i, j, k, l, N, N2,
                         lowend, highend, gap, singleinsert,
                         forced,local,
                         edangle5, edangle3);
#endif
          }
        }
      }
    }
  }
 
#else //COMPILE_SMP defined
  pthread_t consumerids[numProcessors];
  pthread_t producerid;

  // If progress is null, progressObserver will be silent.
  ObservingTextProgressBar *progressObserver = new ObservingTextProgressBar(progress);
  progressObserver->setMaximum(getNumInternalRanks(N, N2)+getNumExternalRanks(N, N2));

  RankManager internalManager(getNumInternalRanks(N, N2),
                              numProcessors);

  internalManager.subscribe(progressObserver);
#ifdef DYNALIGN_II
  rankconsumerargs rca = {&internalManager, &single1_w, &single2_w, max_elongation, islope, iintercept, data, vmod, w, ct1, ct2, v, w5, w3,
                          mod1, mod2, dbl1, dbl2, lowend, highend,
                          fce1, fce2, forcealign,
                          edangle5, edangle3,
                          gap, modification,
                          alignmentforced,  forced,local};
#else
   rankconsumerargs rca = {&internalManager, data, vmod, w, ct1, ct2, v, w5, w3,
                          mod1, mod2, dbl1, dbl2, lowend, highend,
                          fce1, fce2, forcealign,
                          edangle5, edangle3,
                          gap, modification,
                          alignmentforced, singleinsert, forced,local};
#endif
  rankproducerargs rpa = {&internalManager, lowend, highend,
                          N, N2, minloop};


  for (int i = 0; i < numProcessors; i++) {
    pthread_create(&(consumerids[i]), NULL, &rankconsumer, static_cast<void *>(&rca));
  }

  pthread_create(&producerid, NULL, &internalrankproducer, static_cast<void *>(&rpa));
  pthread_join(producerid, NULL);

  for (int i = 0; i < numProcessors; i++) {
    pthread_join(consumerids[i], NULL);
  }
#endif //COMPILE_SMP

  //calculate w3 and w5
  //w3[c][d] = lowest free energy of fragment c->N and d->N2
  //w5[c][d] = lowest free energy of fragment 1->c and 1->d
  //initialize w5[i][k]:
  if (local) {
    for (ip=0;ip<=N;++ip){

      //kp = highend[ip];
      //kp = lowend[j];
      //kp = lowend[ip];
      for (kp=0;kp<=N2;++kp){
        w5->f(ip,kp)=0;//w5[ip][kp] = 0;
      }
    }
  } else {

    for (ip=0;ip<=N;ip++){

		for (kp = 0; kp <= N2; ++kp) {

			w5->f(ip,kp) = gap* (abs(ip - kp));
		}
	}
  }

  for (ip=1;ip<=N;++ip) {
    //cp = min(N2,highend[ip]);//min(N2,ip+maxsep);
    for (kp=1/*max(1,ip-maxsep)*/;kp<=N2;++kp) {
      //ap = kp-ip+maxsep;

      //although this is really a decrement...
      if(ip>1)  ikincrement = kp-1 <= highend[ip-1] && kp-1 >= lowend[ip-1];

      if (
          (!forced ||
           (!dbl1[ip] && !dbl2[kp]))) {
        en1=w5->f(ip-1,kp-1);//w5[ip-1][ap]; //adding one nuc on both that is unpaired and unstacked
      } else {
        en1 = DYNALIGN_INFINITY;
      }
//      if(ip==131&&kp==116)cerr<<"w5->f(i,k) "<<en1<<"\n";

      if (
          (!forced||
           (!dbl1[ip]))) {
        en1 = min(en1,w5->f(ip-1,kp)/*w5[ip-1][ap+1]*/+gap);//add a gapped nuc, ip
      }
      //    if(ip==131&&kp==116)cerr<<"w5->f(i,k) "<<en1<<"\n";

      if (
          (!forced ||
           (!dbl2[kp]))) {
        en1 = min(en1,w5->f(ip,kp-1)/*w5[ip][ap-1]*/+gap);//add a gapped nuc, kp (aka ap)
      }
//      if(ip==131&&kp==116)cerr<<"w5->f(i,k) "<<en1<<"\n";

      for (jp=0;jp+minloop<ip;++jp) {
        //dp = N2;
        dp = kp-minloop-1;//min(dp,kp-minloop-1);
        for (lp=0;lp<=dp;++lp) {


          //check whether mine[i][a] is split so that the lowest free energy is
          //an exterior fragment from 1 to j and a helix from j+1 to i;

          //check all possible alignments (in index a(k) and b(l))

          //must consider whether an unpaired nucleotide is stacked
          //is stacked onto each of the nucleotides in a helix
          //There are 16 cases:
          //consider them in this order (0 unstacked, 1 stacked)
          //		j+1	i	l+1	k
          //1		0	0	0	0
          //2		0	0	0	1
          //3		0	0	1	0
          //4		0	0	1	1
          //5		0	1	0	0
          //6		0	1	0	1
          //7		0	1	1	0
          //8		0	1	1	1
          //9		1	0	0	0
          //10		1	0	0	1
          //11		1	0	1	0
          //12		1	0	1	1
          //13		1	1	0	0
          //14		1	1	0	1
          //15		1	1	1	0
          //16		1	1	1	1

          //note that for these exterior loops:
          //j<i
          //l<k

          //although this is an increment...
          jldecrement  = lp+1 >= lowend[jp+1] && lp+1 <= highend[jp+1];
          //although this is an increment...
          jldecrement2 = lp+2 >= lowend[jp+2] && lp+2 <= highend[jp+2];
		  bool ipkp = kp <= highend[ip]&&kp >= lowend[ip];
		  //bool jplp = lp <= highend[jp]&&lp >= lowend[jp];


          //no stacking
          //case 1: 0	0	0	0
          if(jldecrement&&ipkp) en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+1,ip,lp+1,kp)/*v[ip][jp+1][bp][ap]*/+
                                    penalty(ip,jp+1,ct1,data)+penalty(kp,lp+1,ct2,data));
//	  if(ip==131&&kp==116)cerr<<"w5->f(i,k) 1 "<<en1<<" jp "<<jp<<" lp "<<lp<<"\n";

          //case 6: 0	1	0	1


          if(ip>1&&ikincrement&&jldecrement) en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+1,ip-1,lp+1,kp-1)/*v[ip-1][jp+1][bp][ap]*/+edangle3(ip-1,jp+1,ip,ct1,data)+
                                                 edangle3(kp-1,lp+1,kp,ct2,data)+
                                                 penalty(ip-1,jp+1,ct1,data)+penalty(kp-1,lp+1,ct2,data));
//	  if(ip==131&&kp==116)cerr<<"w5->f(i,k) 6 "<<en1<<" jp "<<jp<<" lp "<<lp<<"\n";

          //case 11	1	0	1	0
          if (jldecrement2&&ipkp) en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+2,ip,lp+2,kp)/*v[ip][jp+2][bp][ap]*/+edangle5(jp+2,ip,jp+1,ct1,data)+
                                      edangle5(lp+2,kp,lp+1,ct2,data)+penalty(ip,jp+2,ct1,data)+penalty(lp+2,kp,ct2,data));
//	  if(ip==131&&kp==116)cerr<<"w5->f(i,k) 11 "<<en1<<" jp "<<jp<<" lp "<<lp<<"\n";


          //case 16	1	1	1	1
          if (ip>1&&ikincrement&&jldecrement2) en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+2,ip-1,lp+2,kp-1)/*v[ip-1][jp+2][bp][ap]*/+edangle5(jp+2,ip-1,jp+1,ct1,data)+
                                                   edangle3(ip-1,jp+2,ip,ct1,data)+edangle5(lp+2,kp-1,lp+1,ct2,data)+
                                                   edangle3(kp-1,lp+2,kp,ct2,data)+penalty(ip-1,jp+2,ct1,data)+penalty(lp+2,kp-1,ct2,data));
//	  if(ip==131&&kp==116)cerr<<"w5->f(i,k) 16 "<<en1<<" jp "<<jp<<" lp "<<lp<<"\n";



          if (kp-1>=lowend[ip]&&kp-1<=highend[ip]) {
            //case 2: 0	0	0	1
            if (jldecrement) en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+1,ip,lp+1,kp-1)/*v[ip][jp+1][bp][ap-1]*/+edangle3(kp-1,lp+1,kp,ct2,data)+
                                       penalty(ip,jp+1,ct1,data)+penalty(kp-1,lp+1,ct2,data)+gap);
//	  if(ip==131&&kp==116)cerr<<"w5->f(i,k) 2 "<<en1<<" jp "<<jp<<" lp "<<lp<<"\n";

            //case 12	1	0	1	1
            if (jldecrement2) en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+2,ip,lp+2,kp-1)/*v[ip][jp+2][bp][ap-1]*/+edangle5(jp+2,ip,jp+1,ct1,data)+
                                        edangle5(lp+2,kp-1,lp+1,ct2,data) + edangle3(kp-1,lp+2,kp,ct2,data)+
                                        penalty(ip,jp+2,ct1,data)+penalty(kp-1,lp+2,ct2,data)+gap);
//	  if(ip==131&&kp==116)cerr<<"w5->f(i,k) 12 "<<en1<<" jp "<<jp<<" lp "<<lp<<"\n";

            if (lp+2<=highend[jp+1]/*bp+1<2*maxsep+2*/&&lp+2>=lowend[jp+1]) {
              //case 4
              en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+1,ip,lp+2,kp-1)/*v[ip][jp+1][bp+1][ap-1]*/+edangle3(kp-1,lp+2,kp,ct2,data)+
                        edangle5(lp+2,kp-1,lp+1,ct2,data)+penalty(ip,jp+1,ct1,data)+penalty(kp-1,lp+2,ct2,data)+2*gap);
            }
          }
//	  if(ip==131&&kp==116)cerr<<"w5->f(i,k) 4 "<<en1<<" jp "<<jp<<" lp "<<lp<<"\n";

          if (lp+2<=highend[jp+1]&&lp+2>=lowend[jp+1]) {
            //case 3: 0	0	1	0
            if(ipkp) en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+1,ip,lp+2,kp)/*v[ip][jp+1][bp+1][ap]*/+edangle5(lp+2,kp,lp+1,ct2,data)+
                      penalty(ip,jp+1,ct1,data)+penalty(lp+2,kp,ct2,data)+gap);
//	  if(ip==131&&kp==116)cerr<<"w5->f(i,k) 3 "<<en1<<" jp "<<jp<<" lp "<<lp<<"\n";

            //case 8:		0	1	1	1
            if(ip>1&&ikincrement) en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+1,ip-1,lp+2,kp-1)/*v[ip-1][jp+1][bp+1][ap]*/+edangle3(ip-1,jp+1,ip,ct1,data)+
                                      edangle5(lp+2,kp-1,lp+1,ct2,data)+edangle3(kp-1,lp+2,kp,ct2,data)+
                                      penalty(ip-1,jp+1,ct1,data)+penalty(lp+2,kp-1,ct2,data)+gap);
//	  if(ip==131&&kp==116)cerr<<"w5->f(i,k) 8 "<<en1<<" jp "<<jp<<" lp "<<lp<<"\n";

            if (ip>1&&kp<=highend[ip-1]&&kp>=lowend[ip-1]) {
              //case 7		0	1	1	0
              en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+1,ip-1,lp+2,kp)/*v[ip-1][jp+1][bp+1][ap+1]*/+edangle5(lp+2,kp,lp+1,ct2,data)+
                        edangle3(ip-1,jp+1,ip,ct1,data)+penalty(ip-1,jp+1,ct1,data)+penalty(lp+2,kp,ct2,data)+2*gap);
//	  if(ip==131&&kp==116)cerr<<"w5->f(i,k) 7 "<<en1<<" jp "<<jp<<" lp "<<lp<<"\n";

            }
          }

          if (ip>1&&kp<=highend[ip-1]&&kp>=lowend[ip-1]) {
            //case5: 0	1	0	0
            if (jldecrement) en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+1,ip-1,lp+1,kp)/*v[ip-1][jp+1][bp][ap+1]*/+edangle3(ip-1,jp+1,ip,ct1,data)+
                                       penalty(ip-1,jp+1,ct1,data)+penalty(kp,lp+1,ct2,data)+gap);
//	  if(ip==131&&kp==116)cerr<<"w5->f(i,k) 5 "<<en1<<" jp "<<jp<<" lp "<<lp<<"\n";

            //case 15	1	1	1	0
            if (jldecrement2) en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+2,ip-1,lp+2,kp)/*v[ip-1][jp+2][bp][ap+1]*/ + edangle5(jp+2,ip-1,jp+1,ct1,data)+
                                        edangle3(ip-1,jp+2,ip,ct1,data)+edangle5(lp+2,kp,lp+1,ct2,data)+
                                        penalty(ip-1,jp+2,ct1,data)+penalty(lp+2,kp,ct2,data)+gap);
//	  if(ip==131&&kp==116)cerr<<"w5->f(i,k) 15 "<<en1<<" jp "<<jp<<" lp "<<lp<<"\n";

            if (lp+1>=lowend[jp+2]&&lp+1<=highend[jp+2]/*bp>0*/) {
              //case 13	1	1	0	0
              en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+2,ip-1,lp+1,kp)/*v[ip-1][jp+2][bp-1][ap+1]*/+edangle5(jp+2,ip-1,jp+1,ct1,data)+
                        edangle3(ip-1,jp+2,ip,ct1,data)+
                        penalty(ip-1,jp+2,ct1,data)+penalty(kp,lp+1,ct2,data)+2*gap);

            }
//	  if(ip==131&&kp==116)cerr<<"w5->f(i,k) 13 "<<en1<<" jp "<<jp<<" lp "<<lp<<"\n";


          }
          if (lp+1>=lowend[jp+2]&&lp+1<=highend[jp+2]/*bp>0*/) {
            //case 9		1	0	0	0
            if(ipkp) en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+2,ip,lp+1,kp)/*v[ip][jp+2][bp-1][ap]*/+edangle5(jp+2,ip,jp+1,ct1,data)+
                      penalty(ip,jp+2,ct1,data)+penalty(kp,lp+1,ct2,data)+gap);

//	  if(ip==131&&kp==116)cerr<<"w5->f(i,k) 9 "<<en1<<" jp "<<jp<<" lp "<<lp<<"\n";

            //case 14	1	1	0	1
            if (ip>1&&ikincrement) en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+2,ip-1,lp+1,kp-1)/*v[ip-1][jp+2][bp-1][ap]*/+edangle5(jp+2,ip-1,jp+1,ct1,data)+
                                       edangle3(ip-1,jp+2,ip,ct1,data)+edangle3(kp-1,lp+1,kp,ct2,data)+
                                       penalty(ip-1,jp+2,ct1,data)+penalty(kp-1,lp+1,ct2,data)+gap);
//	  if(ip==131&&kp==116)cerr<<"w5->f(i,k) 14 "<<en1<<" jp "<<jp<<" lp "<<lp<<"\n";

            if (kp-1>=lowend[ip]&&kp-1<=highend[ip]) {
              //case 10	1	0	0	1
              en1 = min(en1,w5->f(jp,lp)/*w5[jp][bp]*/+v->f(jp+2,ip,lp+1,kp-1)/*v[ip][jp+2][bp-1][ap-1]*/+edangle5(jp+2,ip,jp+1,ct1,data)+
                        edangle3(kp-1,lp+1,kp,ct2,data)+penalty(ip,jp+2,ct1,data)+penalty(kp-1,lp+1,ct2,data)+2*gap);
            }
//	    	  if(ip==131&&kp==116)cerr<<"w5->f(i,k) 10 "<<en1<<" jp "<<jp<<" lp "<<lp<<"\n";

          }
        }
      }
      //report to Dave2
//      if(ip==131&&kp==116)cerr<<"w5->f(i,k) "<<en1<<"\n";

#ifdef DYNALIGN_II   

      for (jp=0;jp+minloop<ip;++jp) {
	en1=min(en1,w5->f(jp,kp)+single1_we.f(jp+1,ip)+gap_branch(ip-jp));
//	if(ip==131&&kp==116)cerr<<"w5->f(i,k) single1 "<<" jp "<<jp<<" "<<en1<<"\n";
      }
	
      for (lp=0;lp+minloop<kp;++lp) {
	en1=min(en1,w5->f(ip,lp)+single2_we.f(lp+1,kp)+gap_branch(kp-lp));
//	if(ip==131&&kp==116)cerr<<"w5->f(i,k) single2 "<<" lp "<<lp<<" "<<en1<<"\n";
      }
#else
#endif
      w5->f(ip,kp) /*w5[ip][ap]*/ = min(en1,w5->f(ip,kp));
//      if(ip==131&&kp==116)cerr<<"w5->f(i,k) "<<w5->f(ip,kp)<<"\n";
      



      if (local) {
        if (en1<lowest) {
          lowest = en1;
        }
      } else if (/*ip==N&&*/abs((N-ip)-(N2-kp))*gap+en1/*en1+gap*abs(ap-maxsep)*/<lowest) {
        lowest = abs((N-ip)-(N2-kp))*gap+en1;
      }
    }
  }

  if (!energyonly) { // calculate w3 if proceeding with extra
                     // calculations beyond energy
    //initialize w3[i][k]:
    if (local) {
      //local alignment calculation
	  //for (kp=N2;kp>=0;kp--){
      //    w3->f(N2+1,kp) = 0;
      //}
		for (kp=N2+1;kp>=0;kp--){
			w3->f(N+1,kp) = 0;
		}
		//for (ip=N+1;ip>0;ip--){//Unallocated space!
		//	w3->f(ip,N2+1) = 0;
		//}
      for (ip=N;ip>0;ip--){

        for (kp=N2;kp>=0;kp--){
          w3->f(ip,kp)=0;
        }
      }
    }
    else {
      //global alignment calculation
		for (kp=N2+1;kp>=0;kp--){
			w3->f(N+1,kp) = gap*(abs(N2-kp+1));
		}

		//uncommented section, DHM 11/13/2010
		for (ip=N+1;ip>0;ip--){ //Now allocated space!
			w3->f(ip,N2+1) = gap*(abs(N-ip+1));
		}
      for (ip=N;ip>0;ip--){
        for (kp=N2;kp>=0;kp--){
          w3->f(ip,kp) = gap*(abs((N-ip)-(N2-kp)));
        }
      }
    }

    for (ip=N;ip>=1;ip--) {
      //cp = max(1,lowend[ip]);
      for (kp=N2;kp>=1;--kp) {

        ikincrement =  kp+1 <= highend[ip+1] && kp+1 >= lowend[ip+1];

        if(kp+1 <= N2 &&
           (!forced ||
            (!dbl1[ip] && !dbl2[kp]))) {
          en1 = w3->f(ip+1,kp+1);//adding one nuc on both that is unpaired and unstacked
        } else {
          en1 = DYNALIGN_INFINITY;
        }

        if (
            (!forced ||
             (!dbl1[ip]))) {
          en1 = min(en1,w3->f(ip+1,kp)/*w3[ip+1][ap-1]*/+gap);//add a gapped nuc, ip
        }

        if (kp+1<=N2 /*ap+1<2*maxsep+2*/ &&
            (!forced ||
             (!dbl2[kp]))) {
          en1 = min(en1,w3->f(ip,kp+1)/*w3[ip][ap+1]*/+gap);//add a gapped nuc, kp (aka ap)
        }

		//if (kp>=lowend[ip]&&kp<=highend[ip]) {
        for (jp=N;jp>=ip+minloop;jp--) {
          dp = max(1/*jp-maxsep-1*/,kp+minloop);//report to Dave3 is this necessary?
          for (lp=N2;lp>=dp;lp--) {
            //bp = lp-jp+maxsep;

			  bool ipkp,jplp;
			  ipkp = kp>=lowend[ip]&&kp<=highend[ip];
			  jplp = lp>=lowend[jp]&&lp<=highend[jp];

            //check whether mine[i][a] is split so that the lowest free energy is
            //an exterior fragment from 1 to j and a helix from j+1 to i;

            //check all possible alignments (in index a(k) and b(l))

            //must consider whether an unpaired nucleotide is stacked
            //is stacked onto each of the nucleotides in a helix
            //There are 16 cases:
            //consider them in this order (0 unstacked, 1 stacked)
            //		jp+1	ip	lp+1	k
            //1		0		0	0		0
            //2		0		0	0		1
            //3		0		0	1		0
            //4		0		0	1		1
            //5		0		1	0		0
            //6		0		1	0		1
            //7		0		1	1		0
            //8		0		1	1		1
            //9		1		0	0		0
            //10	1		0	0		1
            //11	1		0	1		0
            //12	1		0	1		1
            //13	1		1	0		0
            //14	1		1	0		1
            //15	1		1	1		0
            //16	1		1	1		1

			  //if (lp+1>=lowend[jp+1]) {

				  //if ((lp+1<=highend[jp+1])) {
						//This code is no longer needed, DHM 11/13/2010
						//if (lp+1==N2+1) {
							//This is a case where there is no extension on sequence 2, so handle that correctly:
						//	if (local) wval = 0;
						//	else {
						//		wval = gap*(N-jp);

						//	}

						//}
						//else
			  wval=w3->f(jp+1,lp+1);
				  //}

				  //else wval = DYNALIGN_INFINITY;//In the future, this could probably be a condition around these cases...


				  jldecrement = lp-1 <= highend[jp-1] && lp-1 >= lowend[jp-1];

				  //no stacking
				  //case 1: 0		0	0		0
				  if (ipkp&&jplp) en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip,jp,kp,lp)/*v[jp][ip][ap][bp]*/+
							penalty(ip,jp,ct1,data)+penalty(kp,lp,ct2,data));

				  //case 6: 0		1	0		1


				  if (ikincrement&&jplp) en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip+1,jp,kp+1,lp)/*v[jp][ip+1][ap][bp]*/+edangle5(ip+1,jp,ip,ct1,data)+
											 edangle5(kp+1,lp,kp,ct2,data)+
											 penalty(ip+1,jp,ct1,data)+penalty(kp+1,lp,ct2,data));

				  //case 11	1	0	1	0
				  if (jldecrement&&ipkp) en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip,jp-1,kp,lp-1)/*v[jp-1][ip][ap][bp]*/+edangle3(jp-1,ip,jp,ct1,data)+
											 edangle3(lp-1,kp,lp,ct2,data)+penalty(ip,jp-1,ct1,data)+penalty(lp-1,kp,ct2,data));


				  //case 16	1	1	1	1
				  if(ikincrement&&jldecrement) en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip+1,jp-1,kp+1,lp-1)/*v[jp-1][ip+1][ap][bp]*/+edangle3(jp-1,ip+1,jp,ct1,data)+
														 edangle5(ip+1,jp-1,ip,ct1,data)+edangle3(lp-1,kp+1,lp,ct2,data)+
														 edangle5(kp+1,lp-1,kp,ct2,data)+penalty(ip+1,jp-1,ct1,data)+penalty(lp-1,kp+1,ct2,data));




				  if (kp+1<=highend[ip]&&kp+1>=lowend[ip]) {
					//case 2: 0		0	0		1
					if (jplp) en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip,jp,kp+1,lp)/*v[jp][ip][ap+1][bp]*/+edangle5(kp+1,lp,kp,ct2,data)+
							  penalty(ip,jp,ct1,data)+penalty(kp+1,lp,ct2,data)+gap);

					//case 12	1	0	1	1
					if (jldecrement) en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip,jp-1,kp+1,lp-1)/*v[jp-1][ip][ap+1][bp]*/+edangle3(jp-1,ip,jp,ct1,data)+
								   edangle3(lp-1,kp+1,lp,ct2,data) + edangle5(kp+1,lp-1,kp,ct2,data)+//report to Dave Nov.6 change kp->kp+1
								   			   penalty(ip,jp-1,ct1,data)+penalty(kp+1,lp-1,ct2,data)+gap);

					if (lp-1>=lowend[jp]&&lp-1<=highend[jp]) {
					  //case 4 0		0	1		1
					  en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip,jp,kp+1,lp-1)/*v[jp][ip][ap+1][bp-1]*/+edangle5(kp+1,lp-1,kp,ct2,data)+
								edangle3(lp-1,kp+1,lp,ct2,data)+penalty(ip,jp,ct1,data)+penalty(kp+1,lp-1,ct2,data)+2*gap);

					}

				  }


				  if (lp-1>=lowend[jp]&&lp-1<=highend[jp]) {
					//case 3: 0		0	1		0
					if (ipkp) en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip,jp,kp,lp-1)/*v[jp][ip][ap][bp-1]*/+edangle3(lp-1,kp,lp,ct2,data)+
							  penalty(ip,jp,ct1,data)+penalty(lp-1,kp,ct2,data)+gap);

					//case 8:		0	1	1	1
					if (ikincrement) en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip+1,jp,kp+1,lp-1)/*v[jp][ip+1][ap][bp-1]*/+edangle5(ip+1,jp,ip,ct1,data)+
											   edangle3(lp-1,kp+1,lp,ct2,data)+edangle5(kp+1,lp-1,kp,ct2,data)+
											   penalty(ip+1,jp,ct1,data)+penalty(lp-1,kp+1,ct2,data)+gap);

					if (kp>=lowend[ip+1]&&kp<=highend[ip+1]) {
					  //case 7		0	1	1	0
					  en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip+1,jp,kp,lp-1)/*v[jp][ip+1][ap-1][bp-1]*/+edangle3(lp-1,kp,lp,ct2,data)+
								edangle5(ip+1,jp,ip,ct1,data)+penalty(ip+1,jp,ct1,data)+penalty(lp-1,kp,ct2,data)+2*gap);

					}


				  }


				  if (kp>=lowend[ip+1]&&kp<=highend[ip+1]) {
					//case5: 0 1 0 0
					if (jplp) en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip+1,jp,kp,lp)/*v[jp][ip+1][ap-1][bp]*/+edangle5(ip+1,jp,ip,ct1,data)+
							  penalty(ip+1,jp,ct1,data)+penalty(kp,lp,ct2,data)+gap);

					//case 15	1	1	1	0
					if (jldecrement) en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip+1,jp-1,kp,lp-1)/*v[jp-1][ip+1][ap-1][bp]*/ + edangle3(jp-1,ip+1,jp,ct1,data)+
											   edangle5(ip+1,jp-1,ip,ct1,data)+edangle3(lp-1,kp,lp,ct2,data)+
											   penalty(ip+1,jp-1,ct1,data)+penalty(lp-1,kp,ct2,data)+gap);

					if (lp<=highend[jp-1]&&lp>=lowend[jp-1]) {
					  //case 13	1	1	0	0
					  en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip+1,jp-1,kp,lp)/*v[jp-1][ip+1][ap-1][bp+1]*/+edangle3(jp-1,ip+1,jp,ct1,data)+
								edangle5(ip+1,jp-1,ip,ct1,data)+
								penalty(ip+1,jp-1,ct1,data)+penalty(kp,lp,ct2,data)+2*gap);

					}


				  }
				  if (lp<=highend[jp-1]&&lp>=lowend[jp-1]) {
					//case 9		1	0	0	0
					if (ipkp) en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip,jp-1,kp,lp)/*v[jp-1][ip][ap][bp+1]*/+edangle3(jp-1,ip,jp,ct1,data)+
							  penalty(ip,jp-1,ct1,data)+penalty(kp,lp,ct2,data)+gap);


					//case 14	1	1	0	1
					if (ikincrement) en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip+1,jp-1,kp+1,lp)/*v[jp-1][ip+1][ap][bp+1]*/+edangle3(jp-1,ip+1,jp,ct1,data)+
											   edangle5(ip+1,jp-1,ip,ct1,data)+edangle5(kp+1,lp,kp,ct2,data)+
											   penalty(ip+1,jp-1,ct1,data)+penalty(kp+1,lp,ct2,data)+gap);


					if (kp+1<=highend[ip]&&kp+1>=lowend[ip]/*ap<2*maxsep+2*/) {
					  //case 10	1	0	0	1
					  en1 = min(en1,wval/*w3[jp+1][bp]*/+v->f(ip,jp-1,kp+1,lp)/*v[jp-1][ip][ap+1][bp+1]*/+edangle3(jp-1,ip,jp,ct1,data)+
								edangle5(kp+1,lp,kp,ct2,data)+penalty(ip,jp-1,ct1,data)+penalty(kp+1,lp,ct2,data)+2*gap);

					}
				  }
				}
          //}
        }
	//report to Dave4
#ifdef DYNALIGN_II

	for (jp=N;jp>=ip+minloop;jp--) {
	  en1=min(en1,single1_we.f(ip,jp)+w3->f(jp+1,kp)+gap_branch(jp-ip+1));
	}

	for (lp=N2;lp>=kp+minloop;lp--) {
	  en1=min(en1,single2_we.f(kp,lp)+w3->f(ip,lp+1)+gap_branch(lp-kp+1));
	}

#else
#endif

	//is it plausible to just insert a stem in one side?
        w3->f(ip,kp)/*w3[ip][ap]*/ = min(en1,w3->f(ip,kp));
      }
	  //}
    }

    // w3 is now filled, proceed with extra filling operations on v
    // and w

#ifndef COMPILE_SMP
    jmax = 2*N - 1;
    for (; j <= jmax; ++j) {
      if (progress != NULL) {
        progress->update((100 * j) / (2*N));
        if (progress->canceled()) break;
      }

      imin = j - N + 1;
      for (i = min(N, j - 1); i >= imin; --i) {

        kmin = max(lowend[i], 1);
        for (k = min(highend[i], N2); k >= kmin; k--) {

          lmax = min(highend[j],  N2+k);//2 * N2 - 1);
          for (l = max(lowend[j], N2+1); l <= lmax; ++l) {
            //fprintf(stderr, "('v', %d, %d, %d, %d),\n", i, j, k, l);
	    //	    cerr<<i<<" "<<j<<" "<<k<<" "<<l<<" over!\n";
//              if(fce1 == NULL) {cerr << "fce1null" << endl;}
#ifdef DYNALIGN_II
            dynalignstep(ct1, ct2, data,
                         v, w, w5, w3,
                         vmod, &single1_w,&single2_w,mod1, mod2, modification,
                         fce1, fce2,
                         alignmentforced, forcealign,
                         dbl1, dbl2,max_elongation,
                         i, j, k, l,
                         N, N2, lowend, highend, islope, iintercept, gap,
                         forced,local,
                         edangle5, edangle3);
#else
               dynalignstep(ct1, ct2, data,
                         v, w, w5, w3,
                         vmod, mod1, mod2, modification,
                         fce1, fce2,
                         alignmentforced, forcealign,
                         dbl1, dbl2,
                         i, j, k, l, N, N2,
                         lowend, highend, gap, singleinsert,
                         forced,local,
                         edangle5, edangle3);
#endif
          }
        }
      }
    }

    if (progress != NULL && !progress->canceled()) {
      progress->update(100);
    }
#else
    RankManager externalManager(getNumExternalRanks(N, N2),
                                numProcessors);

    externalManager.subscribe(progressObserver);

    rca.manager = &externalManager;
    rpa.manager = &externalManager;

    for (int i = 0; i < numProcessors; i++) {
      pthread_create(&(consumerids[i]), NULL, &rankconsumer, static_cast<void *>(&rca));
    }

    pthread_create(&producerid, NULL, &externalrankproducer, static_cast<void *>(&rpa));
    pthread_join(producerid, NULL);

    for (int i = 0; i < numProcessors; i++) {
      pthread_join(consumerids[i], NULL);
    }
#endif

  }

  //Removed DHM 11/11/2010 -- This is too complicated and not needed.
  //See the end of W5 above for the old, more elegant solution.
  //Find lowest, the lowest folding free energy change structure:
  //if (!local) {
	//  for (i=1;i<=ct1->GetSequenceLength();i++) {
	//	for (k=max(lowend[i],1);k<=highend[i];k++) {
	//	  if (k<N2&&i<N) {
	//		  if (k+1>=lowend[i+1]&&k+1<=highend[i+1]) {
	//			if ((w5->f(i,k)  + w3->f(i+1,k+1) ) < lowest) {
	//			  lowest = w5->f(i,k) + w3->f(i+1,k+1);
	//			}
	//		  }
	//	  }
	//	  else if (k==N2&&i==N) {

	//		if ((w5->f(i,k) ) < lowest) {
	//			  lowest = w5->f(i,k) ;
	//		}

	//	  }
	//	  else if (k==N2) {
	//		if ((w5->f(i,k)+ gap*(N-i) ) < lowest) {
	//			  lowest = w5->f(i,k) + gap*(N-i);
	//		}
	//	  }
	//	  else {
			//i==N
	//		if ((w5->f(i,k)+ gap*(N2-k) ) < lowest) {
	//			  lowest = w5->f(i,k) + gap*(N2-k);
	//		}
	//	  }

	//	}
	//  }
  //}

  if (!(progress&&progress->canceled())) {

  if (Savefile!=NULL) {
    ofstream sav(Savefile,ios::binary);



    //write the save file information so that the sequence can be re-folded,
    //include thermodynamic data to prevent traceback errors


	//Write the save file version
	int version = dynalignsavefileversion;
	write(&sav, &version);
	


    //write a flag as to whether vmod will be needed
    //also include whether the data needed for suboptimal structure prediction has been included
    //or not.  For backwards compatibility, use one digit.
    //0 - suboptimal, no modifications
    //1 - suboptimal, with modifications
    //2 - optimal only, no modifications
    //3 - optimal only, with modifications

    if (ct1->GetNumberofModified()>0||ct2->GetNumberofModified()>0) {
      if (energyonly) flag=3;
      else flag = 1;

    }
    else {
      if (energyonly) flag = 2;
      else flag = 0;

    }
    write(&sav,&flag);


    //write the 3 parameters needed for array allocation
  
	int localint = ct1->GetSequenceLength();
    write(&sav,&(localint));
	localint = ct2->GetSequenceLength();
    write(&sav,&(localint));

    write(&sav,&maxseparation);
  
    //start with structure information
	writestructuresave(&sav,ct1);

	writestructuresave(&sav,ct2);

	//write the datatable
	write(&sav, data);

    //now write the folding data:

    write(&sav,&gap);
    write(&sav,&lowest);
#ifdef DYNALIGN_II
    write(&sav,&max_elongation);
#else
    write(&sav,&singleinsert);
#endif

#ifdef DYNALIGN_II
    write(&sav,&islope);
    write(&sav,&iintercept);
#else
#endif
 

	//write allowed_alignments if maxseparation < 0
	if (maxseparation < 0) {
		for (i=0;i<=ct1->GetSequenceLength();++i) {
			for (j=0;j<=ct2->GetSequenceLength();++j) {
				write(&sav,&(allowed_alignments[i][j]));
			}
		}
	}

    short int kmax, lmax;
    for (i=0;i<=N;++i) {
      if (energyonly) {
        I = N;
      } else {
        I = i+N-1;
      }

      kmax = highend[i];//limit(i, maxsep, N, N2);
      for (j=i;j<=I;++j) {
        lmax = highend[j];//limit(j, maxsep, N, N2);
        for (k = lowend[i]/*limit(i, maxsep, N, N2)*/; k <= kmax; k++) {
          for (l = lowend[j]/*limit(j, maxsep, N, N2)*/; l <= lmax; l++) {
            if (j > N) {
              b = i;
              a = j-N;
            } else {
              b = j;
              a = i;
            }

            if (ct1->tem[b][a]) {
              write(&sav,&(v->f(i,j,k,l)));
            }

            write(&sav,&(w->f(i,j,k,l)));

            if (modification) {
              write(&sav,&(vmod->f(i,j,k,l)));
            }
          }
        }
      }
    }

  
    for (i=0;i<=N+1;++i) {
      for (j=0;j<=N2+1/*2*maxsep+2*/;++j) {
        if (!energyonly) {
          write(&sav,&(w3->f(i,j))/*w3[i][j]*/);
        }
        write(&sav,&(w5->f(i,j))/*w5[i][j]*/);
      }
    }
    if (!energyonly) {
      for (j=0/*lowend[N+1]/*0*/;j<=N2+1/*highend[N+1]/*2*maxsep+2*/;++j) {
        write(&sav,&(w3->f(N+1,j))/*w3[N+1][j]*/);
      }
    }
 
    if (local) {
      flag=1;
    }
    else flag=0;

    write(&sav,&flag);
  
    //for (i=1;i<=2*ct1->GetSequenceLength();i++) {
    //	write(&sav,&pair[0][i]);

    //}
    //for (i=1;i<=2*ct2->GetSequenceLength();i++) {
    //	write(&sav,&pair[1][i]);

    //}
#ifdef DYNALIGN_II
	//write the save file information so that the sequence can be re-folded,
		//include thermodynamic data to prevent traceback errors
		short vers=safiversion;
		write(&sav,&vers); //save a version of the save file 
		//start with structure information
		int sequencelength=ct1->GetSequenceLength();
		write(&sav,&(sequencelength));
		write(&sav,&(ct1->intermolecular));

		int pairs = ct1->GetNumberofPairs();
		write(&sav,&(pairs));
		for (i=0;i<ct1->GetNumberofPairs();i++) {
			pairs=ct1->GetPair5(i);
			write(&sav,&(pairs));
			pairs=ct1->GetPair3(i);
			write(&sav,&(pairs));
		}

		pairs = ct1->GetNumberofForbiddenPairs();
		write(&sav,&(pairs));
		for (i=0;i<ct1->GetNumberofForbiddenPairs();i++) {
			pairs = ct1->GetForbiddenPair5(i);
			write(&sav,&(pairs));
			pairs = ct1->GetForbiddenPair3(i);
			write(&sav,&(pairs));
		}
		for (i=0;i<=ct1->GetSequenceLength();i++) {
		
			write(&sav,&(ct1->hnumber[i]));
			sav.write(&(ct1->nucs[i]),1);

		}

		for (i=0;i<=2*ct1->GetSequenceLength();i++) write(&sav,&(ct1->numseq[i]));
	
		int doubles = ct1->GetNumberofDoubles();
		write(&sav,&(doubles));
		for (i=0;i<ct1->GetNumberofDoubles();i++) {
			doubles = ct1->GetDouble(i);
			write(&sav,&(doubles));
		}

	
		int singles = ct1->GetNumberofSingles();
		write(&sav,&(singles));
		for (i=0;i<ct1->GetNumberofSingles();i++) {
			singles = ct1->GetSingle(i);	
			write(&sav,&(singles));
		}

		int modified;
		modified = ct1->GetNumberofModified();
		write(&sav,&(modified));
		for (i=0;i<ct1->GetNumberofModified();i++) {
			modified = ct1->GetModified(i);
			write(&sav,&(modified));
		}
	
		modified = ct1->GetNumberofGU();
		write(&sav,&(modified));
		for (i=0;i<ct1->GetNumberofGU();i++) {
			modified = ct1->GetGUpair(i);
			write(&sav,&(modified));
		}
	
		string label=ct1->GetSequenceLabel();
		write(&sav,&(label));

		write(&sav,&(ct1->templated));
		if (ct1->templated) {
			for (i=0;i<=ct1->GetSequenceLength();i++) {
				for (j=0;j<=i;j++) write(&sav,&(ct1->tem[i][j]));	

			}

		}

		//write the SHAPE data (for pseudo-free energy constraints)
		write(&sav,&(ct1->shaped));
		if (ct1->shaped) {
			for (i=0;i<=2*ct1->GetSequenceLength();i++) write(&sav,&(ct1->SHAPE[i]));

		}

	
		//now write the array class data for v, w, and wmb:
		for (i=0;i<=ct1->GetSequenceLength();i++) {
                    for (j=0;j<=ct1->GetSequenceLength();j++) {
				write(&sav,&(single1_v.dg[i][j+i]));
				write(&sav,&(single1_w.dg[i][j+i]));
				write(&sav,&(single1_wmb.dg[i][j+i]));
                                write(&sav,&(single1_we.dg[i][j+i]));
				writesinglechar(&sav,&(single1_fce.dg[i][j]));
			}	
		}

                for (i=0;i<=2*ct1->GetSequenceLength();i++) {
			write(&sav,&(single1_lfce[i]));
			write(&sav,&(single1_mod[i]));
		}

		write(&sav,&single1_vmin);



	//write the save file information so that the sequence can be re-folded,
		//include thermodynamic data to prevent traceback errors
//		vers=safiversion;
//		write(&sav,&vers); //save a version of the save file 
		//start with structure information
		sequencelength=ct2->GetSequenceLength();
		write(&sav,&(sequencelength));
		write(&sav,&(ct2->intermolecular));

		pairs = ct2->GetNumberofPairs();
		write(&sav,&(pairs));
		for (i=0;i<ct2->GetNumberofPairs();i++) {
			pairs=ct2->GetPair5(i);
			write(&sav,&(pairs));
			pairs=ct2->GetPair3(i);
			write(&sav,&(pairs));
		}

		pairs = ct2->GetNumberofForbiddenPairs();
		write(&sav,&(pairs));
		for (i=0;i<ct2->GetNumberofForbiddenPairs();i++) {
			pairs = ct2->GetForbiddenPair5(i);
			write(&sav,&(pairs));
			pairs = ct2->GetForbiddenPair3(i);
			write(&sav,&(pairs));
		}
		for (i=0;i<=ct2->GetSequenceLength();i++) {
		
			write(&sav,&(ct2->hnumber[i]));
			sav.write(&(ct2->nucs[i]),1);

		}

		for (i=0;i<=2*ct2->GetSequenceLength();i++) write(&sav,&(ct2->numseq[i]));
	
		doubles = ct2->GetNumberofDoubles();
		write(&sav,&(doubles));
		for (i=0;i<ct2->GetNumberofDoubles();i++) {
			doubles = ct2->GetDouble(i);
			write(&sav,&(doubles));
		}

		singles = ct2->GetNumberofSingles();
		write(&sav,&(singles));
		for (i=0;i<ct2->GetNumberofSingles();i++) {
			singles = ct2->GetSingle(i);	
			write(&sav,&(singles));
		}

		modified = ct2->GetNumberofModified();
		write(&sav,&(modified));
		for (i=0;i<ct2->GetNumberofModified();i++) {
			modified = ct2->GetModified(i);
			write(&sav,&(modified));
		}
	
		modified = ct2->GetNumberofGU();
		write(&sav,&(modified));
		for (i=0;i<ct2->GetNumberofGU();i++) {
			modified = ct2->GetGUpair(i);
			write(&sav,&(modified));
		}
	
		label=ct2->GetSequenceLabel();
		write(&sav,&(label));

		write(&sav,&(ct2->templated));
		if (ct2->templated) {
			for (i=0;i<=ct2->GetSequenceLength();i++) {
				for (j=0;j<=i;j++) write(&sav,&(ct2->tem[i][j]));	

			}

		}

		//write the SHAPE data (for pseudo-free energy constraints)
		write(&sav,&(ct2->shaped));
		if (ct2->shaped) {
			for (i=0;i<=2*ct2->GetSequenceLength();i++) write(&sav,&(ct2->SHAPE[i]));

		}

	
		//now write the array class data for v, w, and wmb:
		for (i=0;i<=ct2->GetSequenceLength();i++) {
			for (j=0;j<=ct2->GetSequenceLength();j++) {
				write(&sav,&(single2_v.dg[i][j+i]));
				write(&sav,&(single2_w.dg[i][j+i]));
				write(&sav,&(single2_wmb.dg[i][j+i]));
				writesinglechar(&sav,&(single2_fce.dg[i][j]));

			}	
		}

		for (i=0;i<=2*ct2->GetSequenceLength();i++) {
			write(&sav,&(single2_lfce[i]));
			write(&sav,&(single2_mod[i]));
		}

		write(&sav,&single2_vmin);

#else
#endif

    sav.close();
  }

  if (energyonly) {
	  //traceback the lowest free energy structure



	//set the alignment to zero
	for (int index=0;index<=ct1->GetSequenceLength();index++) alignment[0][index]=0;

	//calculate the total energy
	ct1->SetEnergy(1,lowest);
	ct2->SetEnergy(1,lowest); 

	ct1->AddStructure();
	ct2->AddStructure();

#ifdef DYNALIGN_II
	error = dyntrace(1, N, 1, N2, ct1, ct2, 1, alignment[0],
                         w, v, w3, w5, lowend, highend,data, islope, iintercept, gap, vmod, &single1_w,&single2_w,&single1_v,&single2_v,&single1_wmb,&single2_wmb,&single1_we,&single2_we,single1_lfce,single2_lfce,single1_mod,single2_mod,&single1_fce,&single2_fce,
			 local,max_elongation,mod1, mod2, modification,
                         fce1, fce2,
                         alignmentforced, forcealign,forced,true);
#else
       	error = dyntrace(1, N, 1, N2, ct1, ct2, 1, alignment[0],
		  w, v, w3, w5, lowend, highend,data, gap, vmod, local, true);
#endif
  } else {
#ifdef DYNALIGN_II
    error = dyntraceback(maxtracebacks,window,awindow,percentsort,
			 v,w,w3,w5,&single1_w,&single2_w,&single1_v,&single2_v,&single1_wmb,&single2_wmb,&single1_we,&single2_we,
			 single1_lfce,single2_lfce,single1_mod,single2_mod,&single1_fce,&single2_fce,
			 ct1,ct2,
			 alignment,lowend,highend,islope,iintercept,gapincrease,data,
			 lowest,vmod,
			 local,max_elongation,mod1, mod2, modification,
                         fce1, fce2,
			    alignmentforced, forcealign,forced);
#else
    error = dyntraceback(maxtracebacks,window,awindow,percentsort,
                         v,w,w3,w5,
                         ct1,ct2,
                         alignment,lowend,highend,gapincrease,data,
                         singleinsert,lowest,vmod,
                         local);
#endif

  }
  // cerr<<"dyntraceback over!\n";
  // zane: sanity check for the arrays
  #ifdef CHECK_ARRAY
  int k0=0, l0=0, lowestij, mfe = w5->f(N, N2);
  bool flag_out=true;
  // check W5 array is equal to the lowest free energy
  if (lowest != mfe)
      cerr << "lowest and W5(N1, N2) are different:"
           << "    lowest: " << lowest
           << "    W5(N1,N2)" << mfe << endl;
  // check the equivalence of w5 and v arrays
  for (i=1;i <= ct1->GetSequenceLength();i++) {
    j = ct1->basepr[1][i];

    if (j > i) {
      lowestij = DYNALIGN_INFINITY;
      for (k=max(lowend[i],1);k<=(min(ct2->GetSequenceLength(),highend[i]));k++)
	for (l=max(lowend[j],k);l<=(min(ct2->GetSequenceLength(),highend[j]));l++) {
	  int sumofv = v->f(i,j,k,l)+v->f(j,i+ct1->GetSequenceLength(),l,k+ct2->GetSequenceLength());
	  if (sumofv <= lowestij) { //find lowest energy for each predicted base pair
	    lowestij = sumofv;
	    k0=k;
	    l0=l;
	  }
	}
      if ( lowestij != mfe ){
	if (flag_out){
          cerr << "Conflicts between w5 and v arrays: "
               << seq1 << " " << seq2
               << endl;
	  cerr << "    i\tj\tk\tl\tMFE_with_i-j(v)\tMFE(w5)\n";
	  flag_out = false;
	}
	cerr << "    " << i << "\t" << j << '\t' << k0 << '\t' << l0 << '\t' << lowestij << '\t' << mfe << endl;
      }
    }
  }

  // check w3(1,1) == w5(N, N2) ?
  if ( w5->f(N, N2) != w3->f(1, 1) )
    cerr << "Conflicts between w5 and w3 arrays: "
         << seq1 << " " << seq2 << "\n"
         << "    w5(N, N2): " << w5->f(N, N2) << "\n"
         << "    w3(1, 1):  " << w3->f(1, 1)
         << endl;


  #endif
  } // if (!progress->canceled())
  delete w;
  delete v;


  delete w3;
  delete w5;



  if (forced) {




    for (i=0;i<2*N;++i) {
      delete[] fce1[i];
    }

    delete[] fce1;



    for (i=0;i<2*N2;++i) {
      delete[] fce2[i];
    }

    delete[] fce2;
    delete[] dbl1;
    delete[] dbl2;

    if (ct1->GetNumberofModified()>0||ct2->GetNumberofModified()>0) {
      //delete vmod
      delete vmod;


    }
  }

  delete[] lowend;
  delete[] highend;

#ifdef timer
  timeout << time(NULL)<<"\n";
  timeout << time(NULL) - seconds;
  timeout.close();
#endif


#ifdef COMPILE_SMP
  //If progress was allocated in this function, remove it.
  //if (removeprogress) delete progress;
  delete progressObserver;
#endif
  //return the errorcode, which was returned during traceback
  // cerr<<"dynalign over!\n";
  return error;


}



//calculate a single point in the v and w arrays - allowing constraints if force is true
#ifdef DYNALIGN_II
void dynalignstep(structure *ct1, structure *ct2, datatable *data,
                  varray *v, dynalignarray *w, wendarray *w5, wendarray *w3,
                  dynalignarray *vmod,   DynProgArray<integersize> *single1_w, DynProgArray<integersize> *single2_w,bool *mod1, bool *mod2, bool modification,
                  char **fce1, char **fce2, bool alignmentforced, short **forcealign,
                  bool *dbl1, bool *dbl2, int max_elongation,
                  int i, int j, int k, int l, int N, int N2,
                  short *lowend, short *highend, int islope, int iintercept, int gap,
                  bool force,bool local,
                  short int (*edangle5)(int i, int j, int ip, structure *ct, datatable *data),
                  short int (*edangle3)(int i, int j, int ip, structure *ct, datatable *data)
		  ) 
#else
void dynalignstep(structure *ct1, structure *ct2, datatable *data,
                  varray *v, dynalignarray *w, wendarray *w5, wendarray *w3,
                  dynalignarray *vmod, bool *mod1, bool *mod2, bool modification,
                  char **fce1, char **fce2, bool alignmentforced, short **forcealign,
                  bool *dbl1, bool *dbl2,
                  int i, int j, int k, int l, int N, int N2,
                  short *lowend, short *highend, int gap, bool singleinsert,
                  bool force,bool local,
                  short int (*edangle5)(int i, int j, int ip, structure *ct, datatable *data),
                  short int (*edangle3)(int i, int j, int ip, structure *ct, datatable *data))
#endif
{

  short int c,d,e,f;
  integersize idangle,jdangle,ldangle,kdangle,ijdangle,kldangle;
  bool ikincrement,jldecrement,ikincrement2,jldecrement2;
  short int startd,endd,starte,ende,startf,endf;

  vector< vector<bool> > inc = data->pairing;

  integersize einternal[maxloop+1][maxloop+1]; //precalculate internal loop free energies
  integersize einternal2[maxloop+1][maxloop+1]; //precalculate internal loop free energies

  integersize en1,en2,en3;
#ifdef DYNALIGN_II
  integersize en1_temp;
#else
#endif

  //use a and b when testing whether a recursion is valid
  //a = k-i+maxsep;
  //if (l<=N2) b = l-j+maxsep;
  //else b = l-j+maxsep+Ndiff;

  //calculate some values for bounds checking
  /*kmin = lowlimit(i,maxsep,N,N2);
    kmax = highlimit(i,maxsep,N,N2);

    lmin = lowlimit(j,maxsep,N,N2);
    lmax = highlimit(j,maxsep,N,N2);*/


  //filter out isolated base pairs
  if (i+1<N&&j-1!=N&&i+1<j-1) e = inc[ct1->numseq[i+1]][ct1->numseq[j-1]];
  else e = 0;
  if (i-1>0&&j+1<2*N&&j!=N) e+=inc[ct1->numseq[i-1]][ct1->numseq[j+1]];


  if (k+1<ct2->GetSequenceLength()&&l-1!=ct2->GetSequenceLength()&&k+1<l-1) f = inc[ct2->numseq[k+1]][ct2->numseq[l-1]];
  else f = 0;
  if (k-1>0&&l+1<2*N2&&l!=N2) f += inc[ct2->numseq[k-1]][ct2->numseq[l+1]];
  //  cerr<<"i "<<i<<" j "<<j<<" k "<<k<<" l "<<l<<"\n";

  //declare some variables needed for processing calculations with constraints
  short g;
  integersize en1m,en2m,en3m;
  bool helicalstack1,helicalstack2,stuck;
  //report to Dave 17
//  cerr << "force" << force << endl;

  if (force) {
    //if i,j,k, or l should be single stranded, don't allow the pair
    if (fce1[jref(i,j,N)][iref(i,j,N)]&SINGLE||fce2[jref(k,l,N2)][iref(k,l,N2)]&SINGLE||
        fce1[jref(i,j,N)][iref(i,j,N)]&NOPAIR||fce2[jref(k,l,N2)][iref(k,l,N2)]&NOPAIR) {
      e=0;
      f=0;
    }
    if (alignmentforced) {
      //if alignment is being forced, be sure that the pairs satisfy the requirement

      if (forcealign[0][i]>0&&forcealign[0][i]!=k) {
        e=0;
        f=0;

      }
      else if (forcealign[0][j]>0&&forcealign[0][j]!=l) {
        e=0;
        f=0;

      }
      else if (forcealign[1][k]>0&&forcealign[1][k]!=i) {
        e=0;
        f=0;

      }
      else if (forcealign[1][l]>0&&forcealign[1][l]!=j) {
        e=0;
        f=0;

      }
      else {
        //scan through to make sure that the base pair isn't precluding an alignment
        for (g=i+1;g<j;g++) {
          if (forcealign[0][g]>0) {

            if (forcealign[0][g]<k||forcealign[0][g]>l) {
              e=0;
              f=0;
            }
          }
        }

        for (g=k+1;g<l;g++) {
          if (forcealign[1][g]>0) {

            if (forcealign[1][g]<i||forcealign[1][g]>j) {
              e=0;
              f=0;
            }
          }
        }
      }
    }
  }

  // Use template information to decide whether to proceed with
  // calculating v(i,j,k,l)
  if (ct1->templated) {
    if (j<=N) {
      c=j;
      d=i;
    }
    else {
      d=j-N;
      c=i;
    }

    if (!ct1->tem[c][d]) {
      //forbind this pair
      e = 0;
    }
  }
  if (ct2->templated) {
    if (l<=N2) {
      c=l;
      d=k;
    }
    else {
      d=l-N2;
      c=k;
    }

    if (!ct2->tem[c][d]) {
      //forbind this pair
      f = 0;
    }
  }

  //test whether it is safe to increment or decrement the aligned nucleotides


  ikincrement = k+1 >= lowend[i+1] && k+1 <= highend[i+1];
  jldecrement = l-1 <= highend[j-1] && l-1 >= lowend[j-1];

  ikincrement2 = k+2 >= lowend[i+2] && k+2 <= highend[i+2];
  jldecrement2 = l-2 <= highend[j-2] && l-2 >= lowend[j-2];


  //Fill v, the best energy with i-j paired, k-l paired,
  //and the pair of basepairs aligned:
  if (inc[ct1->numseq[i]][ct1->numseq[j]]&&
      inc[ct2->numseq[k]][ct2->numseq[l]]&&e&&f)
    {
      //Now consider hairpins for both
      if (j<=N) v->f(i,j,k,l) =
                    erg3(i, j, ct1, data, force ? fce1[jref(i,j,N)][iref(i,j,N)] : 0) +
                    erg3(k, l, ct2, data, force ? fce2[jref(k,l,N2)][iref(k,l,N2)] : 0) +
                    gap*abs(j-i-l+k);
      //if(i==54&&j==61&&k==66&&l==74)cerr<<"v->f(i,j,k,l) after"<<v->f(i,j,k,l)<<"\n";
      //   if(i==197&&j==234&&k==204&&l==226)cerr<<"v->f(i,j,k,l) "<<v->f(i,j,k,l)<<"\n";
      
        en1 = DYNALIGN_INFINITY;

		if (j > N ) {//&& jldecrement && (ikincrement):This needs to be more tailored for each case

			//Consider the exterior loop closed by i-j (actually j-N paired to i):
          //must consider whether an unpaired nucleotide is stacked
          //is stacked onto each of the nucleotides in a helix
          //There are 16 cases:
          //consider them in this order (0 unstacked, 1 stacked)
          //		i	j	k	l
          //1		0	0	0	0
          //2		0	0	0	1
          //3		0	0	1	0
          //4		0	0	1	1
          //5		0	1	0	0
          //6		0	1	0	1
          //7		0	1	1	0
          //8		0	1	1	1
          //9		1	0	0	0
          //10		1	0	0	1
          //11		1	0	1	0
          //12		1	0	1	1
          //13		1	1	0	0
          //14		1	1	0	1
          //15		1	1	1	0
          //16		1	1	1	1

          //to save time, only do the function calls once:
          idangle = edangle3(i,j,i+1,ct1,data);
          jdangle = edangle5(j,i,j-1,ct1,data);
          if (k<N2) kdangle = edangle3(k,l,k+1,ct2,data);
          else kdangle = DYNALIGN_INFINITY;
          if (l>N2+1) ldangle = edangle5(l,k,l-1,ct2,data);
          else ldangle=DYNALIGN_INFINITY;

		  //DHM: Depracated
		  //allow an exception for jldecrement2 if this reaches the end of the sequence
		  //if (j-2==N) jldecrement2=true;


          //case 1, no stacks
		  //if (jldecrement) {
		  /*if (i<N&&k<N2)*/ 
		  en1 = min(en1,w5->f(j-1-N,l-1-N2)+w3->f(i+1,k+1));
		  
		  //else if (k==N2&&local) en1 = min(en1,w5->f(j-1-N,l-1-N2));
		  //else if (k==N2) en1 = min(en1,w5->f(j-1-N,l-1-N2)+gap*(N-i));

		  //}

          //case 16, all stack
          if ((j-2-N>=0)&&(l-2-N2>=0)) if ((i<N)&&(j>N+1)&&(k<N2)) en1 = min(en1,w5->f(j-2-N,l-2-N2)+w3->f(i+2,k+2)/*w5[j-2-N][b]+w3[i+2][a]*/+jdangle+ldangle
                                                                    +idangle+kdangle);

          //case 6, j and l stacked
		  if (j-2-N>=0&&l-2-N2>=0) if (j>N+1) {
			//if (ikincrement) 
			  en1 = min(en1,w5->f(j-2-N,l-2-N2)+w3->f(i+1,k+1)/*w5[j-2-N][b]+w3[i+1][a]*/+jdangle+ldangle);
			//else if (k==N2&&local)  en1 = min(en1,w5->f(j-2-N,l-2-N2)+jdangle+ldangle);
			//else if (k==N2) en1 = min(en1,w5->f(j-2-N,l-2-N2)+jdangle+ldangle+gap*(N-i));

		  }

          //case 11, i and k stacked
          if ((i<N)&&(k<N2)) en1 = min(en1,w5->f(j-1-N,l-1-N2)+w3->f(i+2,k+2)/*w5[j-1-N][b]+w3[i+2][a]*/+idangle+kdangle);

          if (l-2-N2>=0) {
            //case 2, l stacked
            /*if (i+1<=N&&k+1<=N2)*/ en1 = min(en1,w5->f(j-1-N,l-2-N2)+w3->f(i+1,k+1)/*w5[j-1-N][b-1]+w3[i+1][a]*/+ldangle+gap);
			//else if (k==N2&&local) en1 = min(en1,w5->f(j-1-N,l-2-N2)+ldangle+gap);
			//else if (k==N2) en1 = min(en1,w5->f(j-1-N,l-2-N2)+ldangle+gap+gap*(N-i));

            //case 12, i, k, and l stacked
            if ((i<N)&&(k<N2)) en1 = min(en1,w5->f(j-1-N,l-2-N2)+w3->f(i+2,k+2)/*w5[j-1-N][b-1]+w3[i+2][a]*/+ldangle+idangle
                                               +kdangle+gap);

            if ((k<N2)) {
              //case 4, l and k stacked
              en1 = min(en1,w5->f(j-1-N,l-2-N2)+w3->f(i+1,k+2)/*w5[j-1-N][b-1]+w3[i+1][a+1]*/+ldangle+kdangle+2*gap);
            }

            if ((i<N)) {
              //case 10, i and l stacked
              /*if (k+1<=N2)*/ en1 = min(en1,w5->f(j-1-N,l-2-N2)+w3->f(i+2,k+1)/*w5[j-1-N][b-1]+w3[i+2][a-1]*/+ldangle+idangle+2*gap);
			  //else if (k==N2 && local) en1 = min(en1,w5->f(j-1-N,l-2-N2)+ldangle+idangle+2*gap);
			  //else if (k==N2) en1 = min(en1,w5->f(j-1-N,l-2-N2)+ldangle+idangle+2*gap+gap*(N-i-1));
            }



          }

          if ((k<N2)) {
            //case 3, k alone stacked
            en1 = min(en1,w5->f(j-1-N,l-1-N2)+w3->f(i+1,k+2)/*w5[j-1-N][b]+w3[i+1][a+1]*/+kdangle+gap);

			if (j-2-N>=0) {
            //case 8, j, k, and l stacked
				if (l-1-N2>0) if ((j>N+1)) en1 = min(en1,w5->f(j-2-N,l-2-N2)+w3->f(i+1,k+2)/*w5[j-2-N][b]+w3[i+1][a+1]*/+kdangle
													 +jdangle+ldangle+gap);

				//if (l-1-N2<=highend[j-2-N]&&l-1-N2>=lowend[j-2-N]/*b<2*maxsep+2*/) {
				  //case 7, j and k stacked
				  if ((j>N+1)) en1 = min(en1,w5->f(j-2-N,l-1-N2)+w3->f(i+1,k+2)/*w5[j-2-N][b+1]+w3[i+1][a+1]*/+kdangle+jdangle+2*gap);
				//}
			}
          }

		  if (j-2-N>=0) {
			  //if (l-1-N2<=highend[j-2-N]&&l-1-N2>=lowend[j-2-N]) {
				//case 5, j stacked
				  if ((j>N+1)) {

					  /*if (i+1<N&&k+1<N2)*/ en1 = min(en1,w5->f(j-2-N,l-1-N2)+w3->f(i+1,k+1)/*w5[j-2-N][b+1]+w3[i+1][a]*/+jdangle+gap);
						//else if (k==N2&&local) en1 = min(en1,w5->f(j-2-N,l-1-N2)+jdangle+gap);
						//else if (k==N2) en1 = min(en1,w5->f(j-2-N,l-1-N2)+jdangle+gap+gap*(N-i));
				  }
				//case 15, i, j, and k stacked
				if ((i<N)&&(j>N+1)&&k+1<=N2) en1 = min(en1,w5->f(j-2-N,l-1-N2)+w3->f(i+2,k+2)/*w5[j-2-N][b+1]+w3[i+2][a]*/+jdangle
															+idangle+kdangle+gap);

				if ((i<N)&&(j>N+1)) {
				  //case 13, j and i stacked
					//if (k+1<N2) {
						en1 = min(en1,w5->f(j-2-N,l-1-N2)+w3->f(i+2,k+1)/*w5[j-2-N][b+1]+w3[i+2][a-1]*/+jdangle+idangle+2*gap);
					//}
					//else if (k==N2&&local) en1 = min(en1,w5->f(j-2-N,l-1-N2)+jdangle+idangle+2*gap);
					//else if (k==N2)en1 = min(en1,w5->f(j-2-N,l-1-N2)+jdangle+idangle+2*gap+gap*(N-i-1));
				}
			  //}
		  }

          //if (k+1>=lowend[i+2]) {
            //case 9, i stacked
			  if ((i<N)) {
				  /*if (k+1<=N2)*/ en1 = min(en1,w5->f(j-1-N,l-1-N2)+w3->f(i+2,k+1)/*w5[j-1-N][b]+w3[i+2][a-1]*/+idangle+gap);
				  //else if (k==N2&&local) en1 = min(en1,w5->f(j-1-N,l-1-N2)+idangle+gap);
				  //else if (k==N2) en1 = min(en1,w5->f(j-1-N,l-1-N2)+idangle+gap+gap*(N-i-1));
			  }
            //case 14, i, j, and l stacked
			  if ((i<N)&&(j>N+1)&&(j-2-N>=0)&&(l-2-N2>0) ) {
				  /*if (k+1<=N2)*/ en1 = min(en1,w5->f(j-2-N,l-2-N2)+w3->f(i+2,k+1)/*w5[j-2-N][b]+w3[i+2][a-1]*/+idangle
                                                        +jdangle+ldangle+gap);
				  //else if (k==N2&&local) en1 = min(en1,w5->f(j-2-N,l-2-N2)+idangle
                  //                                      +jdangle+ldangle+gap);
				 // else if (k==N2) en1 = min(en1,w5->f(j-2-N,l-2-N2)+idangle
                 //                                       +jdangle+ldangle+gap+gap*(N-i-1));

			  }
          //}


			//DHM:Depracated
		  //reset jl decrement2 for remainder of code... i.e. undo the exception above
		  //if (j-2==N) jldecrement2 = l-2 <= highend[j-2] && l-2 >= lowend[j-2];

		}

//    if(i==91&&j==203&&k==109&&l==195)cerr<<"v->f(i,j,k,l) "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<"\n";

       if (i!=N&&j!=N+1&&k!=N2&&l!=N2+1&&ikincrement&&jldecrement) {
          //Now consider multiloop:
          //junction closed by i-j pair aligned with a and b
          //calculate the free energy of 2 fragments merged:
	  //report to Dave7 Dec.26 why do we need ikincrement and jldecrement?
          idangle = edangle3(i,j,i+1,ct1,data);
          jdangle = edangle5(j,i,j-1,ct1,data);
          kdangle = edangle3(k,l,k+1,ct2,data);
          ldangle = edangle5(l,k,l-1,ct2,data);


          for (c=i+minloop+1;c<j-minloop;++c) {
            //if (c>N) kp = Ndiff;
            //else kp=0;
            startd=max(k+minloop+1,lowend[c]);
            endd=min(highend[c],l-minloop);
            for (d=startd/*c-maxsep-kp*/;
                 d<endd/*d<l-minloop&&d<=c+maxsep-kp*/;++d) {

              //e = d-c+maxsep+kp;

              if (((c<N)&&(d<N2)||((c>N)&&(d>N2)))&&d+1>=lowend[c+1]&&d+1<=highend[c+1]) {


                //must consider whether an unpaired nucleotide
                //is stacked onto each of the nucleotides in a helix
                //There are 16 cases:
                //consider rthem in this order (0 unstacked, 1 stacked)
                //		i	j	k	l
                //1		0	0	0	0
                //2		0	0	0	1
                //3		0	0	1	0
                //4		0	0	1	1
                //5		0	1	0	0
                //6		0	1	0	1
                //7		0	1	1	0
                //8		0	1	1	1
                //9		1	0	0	0
                //10	1	0	0	1
                //11	1	0	1	0
                //12	1	0	1	1
                //13	1	1	0	0
                //14	1	1	0	1
                //15	1	1	1	0
                //16	1	1	1	1



                //case 1 - no stacks
                en1 = min(en1,w->f(i+1,c,k+1,d)+w->f(c+1,j-1,d+1,l-1)/*w[c][iref(i+1,c,N)][a][e]*/
                          /*w[jref(c+1,j-1,N)][iref(c+1,j-1,N)][e][b]*/+2*data->eparam[5]+2*data->eparam[10]);


                //case 16 - all four stacked

                if (((i+1)!=N)&&((j-1)!=N)&&((k+1)!=N2)&&((l-1)!=N2)&&ikincrement2&&jldecrement2) {
                  en1 = min(en1,w->f(i+2,c,k+2,d)+w->f(c+1,j-2,d+1,l-2)/*w[c][iref(i+2,c,N)][a][e]
                                                                         +w[jref(c+1,j-2,N)][iref(c+1,j-2,N)][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                            +4*data->eparam[6]+ldangle
                            +kdangle
                            +jdangle
                            +idangle);
                }


                //case 6 - j and l stacked:
                if (((j-1)!=N)&&((l-1)!=N2)&&jldecrement2) {
                  en1 = min(en1,w->f(i+1,c,k+1,d)+w->f(c+1,j-2,d+1,l-2)/*w[c][iref(i+1,c,N)][a][e]
                                                                         +w[jref(c+1,j-2,N)][iref(c+1,j-2,N)][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                            +2*data->eparam[6]+jdangle
                            +ldangle);
                }


                //case 11 - i and k stacked
                if (((i+1)!=N)&&((k+1)!=N2)&&ikincrement2) {
                  en1 = min(en1,w->f(i+2,c,k+2,d)+w->f(c+1,j-1,d+1,l-1)/*w[c][iref(i+2,c,N)][a][e]
                                                                         +w[jref(c+1,j-1,N)][iref(c+1,j-1,N)][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                            +2*data->eparam[6]+kdangle+
                            idangle);
                }



                if (l-2>=lowend[j-1]&&l-2<=highend[j-1]&&((l-1)!=N2)) {
                  //case 2 - stack on l
                  en1 = min(en1,w->f(i+1,c,k+1,d)+w->f(c+1,j-1,d+1,l-2)/*w[c][iref(i+1,c,N)][a][e]
                                                                         +w[jref(c+1,j-1,N)][iref(c+1,j-1,N)][e][b-1]*/+2*data->eparam[5]+2*data->eparam[10]
                            +data->eparam[6]+ldangle+gap);

                  if ((i+1)!=N) {
                    //case 12 - i, k, and l stacked
                    if ((k+1)!=N2&&ikincrement2) {
                      en1 = min(en1,w->f(i+2,c,k+2,d)+w->f(c+1,j-1,d+1,l-2)/*w[c][iref(i+2,c,N)][a][e]
                                                                             +w[jref(c+1,j-1,N)][iref(c+1,j-1,N)][e][b-1]*/+2*data->eparam[5]+2*data->eparam[10]
                                +3*data->eparam[6]+ldangle
                                +idangle+kdangle+gap);
                    }


                    if (k+1>=lowend[i+2]&&k+1<=highend[i+2]) {
                      //case 10 - l and i stacked
                      en1 = min(en1,w->f(i+2,c,k+1,d)+w->f(c+1,j-1,d+1,l-2)/*w[c][iref(i+2,c,N)][a-1][e]
                                                                             +w[jref(c+1,j-1,N)][iref(c+1,j-1,N)][e][b-1]*/+2*data->eparam[5]+2*data->eparam[10]
                                +2*data->eparam[6]+ldangle
                                +idangle+2*gap);


                    }
                  }
                  if (k+2<=highend[i+1]&&k+2>=lowend[i+1]&&((k+1)!=N2)) {
                    //case 4 - k and l stacked
                    en1 = min(en1,w->f(i+1,c,k+2,d)+w->f(c+1,j-1,d+1,l-2)/*w[c][iref(i+1,c,N)][a+1][e]
                                                                           +w[jref(c+1,j-1,N)][iref(c+1,j-1,N)][e][b-1]*/+2*data->eparam[5]+2*data->eparam[10]
                              +2*data->eparam[6]+ldangle
                              +kdangle+2*gap);
                  }
                }
                if (k+1>=lowend[i+2]&&k+1<=highend[i+2]&&(i+1!=N)) {

                  //case 9 - i stacked
                  en1 = min(en1,w->f(i+2,c,k+1,d)+w->f(c+1,j-1,d+1,l-1)/*w[c][iref(i+2,c,N)][a-1][e]
                                                                         +w[jref(c+1,j-1,N)][iref(c+1,j-1,N)][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                            +data->eparam[6]+idangle+gap);

                  if (j-1!=N) {
                    //case 14 - i, j, and l stacked
                    if (l-1!=N2&&jldecrement2) {
                      en1 = min(en1,w->f(i+2,c,k+1,d)+w->f(c+1,j-2,d+1,l-2)/*w[c][iref(i+2,c,N)][a-1][e]
                                                                             +w[jref(c+1,j-2,N)][iref(c+1,j-2,N)][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                                +3*data->eparam[6]+ldangle
                                +jdangle
                                +idangle+gap);
                    }

                    if (l-1>=lowend[j-2]&&l-1<=highend[j-2]/*b+1<2*maxsep+2*/) {
                      //case 13 - i and j stacked
                      en1 = min(en1,w->f(i+2,c,k+1,d)+w->f(c+1,j-2,d+1,l-1)/*w[c][iref(i+2,c,N)][a-1][e]
                                                                             +w[jref(c+1,j-2,N)][iref(c+1,j-2,N)][e][b+1]*/+2*data->eparam[5]+2*data->eparam[10]
                                +2*data->eparam[6]+jdangle
                                +idangle+2*gap);
                    }
                  }
                }

                if (k+2<=highend[i+1]&&k+2>=lowend[i+1]/*(a+1<2*maxsep+2)*/&&(k+1!=N2)) {
                  //case 3 - stack on k
                  en1 = min(en1,w->f(i+1,c,k+2,d)+w->f(c+1,j-1,d+1,l-1)/*w[c][iref(i+1,c,N)][a+1][e]*/+2*data->eparam[10]
                            +/*w[jref(c+1,j-1,N)][iref(c+1,j-1,N)][e][b]*/+2*data->eparam[5]
                            +data->eparam[6]+kdangle+gap);

                  if (j-1!=N) {
                    //case 8 - j, k, and l stacked
                    if (l-1!=N2&&jldecrement2) {
                      en1 = min(en1,w->f(i+1,c,k+2,d)+w->f(c+1,j-2,d+1,l-2)/*w[c][iref(i+1,c,N)][a+1][e]
                                                                             +w[jref(c+1,j-2,N)][iref(c+1,j-2,N)][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                                +3*data->eparam[6]+jdangle
                                +ldangle+kdangle+gap);
                    }
                    if (l-1<=highend[j-2]&&l-1>=lowend[j-2]) {
                      //case 7 - j and k stacked
                      en1 = min(en1,w->f(i+1,c,k+2,d)+w->f(c+1,j-2,d+1,l-1)/*w[c][iref(i+1,c,N)][a+1][e]
                                                                             +w[jref(c+1,j-2,N)][iref(c+1,j-2,N)][e][b+1]*/+2*data->eparam[5]+2*data->eparam[10]
                                +2*data->eparam[6]+jdangle
                                +kdangle+2*gap);
                    }
                  }


                }
                if (l-1<=highend[j-2]&&l-1>=lowend[j-2]&&(j-1!=N)) {
                  //case 15 - i, j, and k stacked
                  if ((i+1!=N)&&(k+1!=N2)&&ikincrement2) {
                    en1 = min(en1,w->f(i+2,c,k+2,d)+w->f(c+1,j-2,d+1,l-1)/*w[c][iref(i+2,c,N)][a][e]
                                                                           +w[jref(c+1,j-2,N)][iref(c+1,j-2,N)][e][b+1]*/+2*data->eparam[5]+2*data->eparam[10]
                              +3*data->eparam[6]+kdangle
                              +jdangle
                              +idangle+gap);
                  }

                  //case 5 - j stacked
                  en1 = min(en1,w->f(i+1,c,k+1,d)+w->f(c+1,j-2,d+1,l-1)/*w[c][iref(i+1,c,N)][a][e]
                                                                         +w[jref(c+1,j-2,N)][iref(c+1,j-2,N)][e][b+1]*/+2*data->eparam[5]+2*data->eparam[10]
                            +data->eparam[6]+jdangle+gap);
                }
              }
            }
          }
        }
   
//       if(i==91&&j==203&&k==109&&l==195)cerr<<"v->f(i,j,k,l) "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<"\n";
#ifdef DYNALIGN_II

       if (i!=N&&j!=N+1&&k!=N2&&l!=N2+1) {
	 
	 idangle = edangle3(i,j,i+1,ct1,data);
	 jdangle = edangle5(j,i,j-1,ct1,data);
	 kdangle = edangle3(k,l,k+1,ct2,data);
	 ldangle = edangle5(l,k,l-1,ct2,data);
	
	 for (c=i+1+minloop;c+1<j-1-minloop&&c!=N;++c) {
	   //case1 0 0 0 0
        en1=min(en1,branch_1(single1_w,w,i+1,j-1,k+1,l-1,N,c,lowend,highend, iintercept, islope)+2*data->eparam[5]+2*data->eparam[10]);
//	      if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 1 1 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" c "<<c<<"\n";
	   //case2 0 0 0 1
	    en1=min(en1,branch_1(single1_w,w,i+1,j-1,k+1,l-2,N,c,lowend,highend, iintercept, islope)+ldangle+gap+2*data->eparam[5]+2*data->eparam[10]
		    +data->eparam[6]);
//	      if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 1 2 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" c "<<c<<"\n";
	    //case3 0 0 1 0
	     en1=min(en1,branch_1(single1_w,w,i+1,j-1,k+2,l-1,N,c,lowend,highend, iintercept, islope)+kdangle+gap+2*data->eparam[5]+2*data->eparam[10]
		     +data->eparam[6]);
//	         if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 1 3 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" c "<<c<<"\n";
	     //case4 0 0 1 1
	     en1=min(en1,branch_1(single1_w,w,i+1,j-1,k+2,l-2,N,c,lowend,highend, iintercept, islope)+kdangle+ldangle+2*gap+2*data->eparam[5]+2*data->eparam[10]
		     +2*data->eparam[6]);
//	         if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 1 4 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" c "<<c<<"\n";
	     //case5 0 1 0 0
	     en1=min(en1,branch_1(single1_w,w,i+1,j-2,k+1,l-1,N,c,lowend,highend, iintercept, islope)+jdangle+gap+2*data->eparam[5]+2*data->eparam[10]
		     +data->eparam[6]);
//	          if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 1 5 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" c "<<c<<"\n";
	     //case6 0 1 0 1
	     en1=min(en1,branch_1(single1_w,w,i+1,j-2,k+1,l-2,N,c,lowend,highend, iintercept, islope)+jdangle+ldangle+2*data->eparam[5]+2*data->eparam[10]
		     +2*data->eparam[6]);
//	          if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 1 6 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" c "<<c<<"\n";
	     //case7 0 1 1 0
	      en1=min(en1,branch_1(single1_w,w,i+1,j-2,k+2,l-1,N,c,lowend,highend, iintercept, islope)+jdangle+kdangle+2*gap+2*data->eparam[5]+2*data->eparam[10]
		      +2*data->eparam[6]);
//	            if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 1 7 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" c "<<c<<"\n";
	      //case8 0 1 1 1
	      en1=min(en1,branch_1(single1_w,w,i+1,j-2,k+2,l-2,N,c,lowend,highend, iintercept, islope)+jdangle+kdangle+ldangle+gap+2*data->eparam[5]+2*data->eparam[10]
		       +3*data->eparam[6]);
//	      if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 1 8 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" c "<<c<<"\n";
	      //case9 1 0 0 0
	      en1=min(en1,branch_1(single1_w,w,i+2,j-1,k+1,l-1,N,c,lowend,highend, iintercept, islope)+idangle+gap+2*data->eparam[5]+2*data->eparam[10]
		      +data->eparam[6]);
//	           if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 1 9 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" c "<<c<<"\n";
	      //case10 1 0 0 1
	      en1=min(en1,branch_1(single1_w,w,i+2,j-1,k+1,l-2,N,c,lowend,highend, iintercept, islope)+idangle+ldangle+2*gap+2*data->eparam[5]+2*data->eparam[10]
		      +2*data->eparam[6]);
//	           if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 1 10 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" c "<<c<<"\n";
	      //case11 1 0 1 0
	      en1=min(en1,branch_1(single1_w,w,i+2,j-1,k+2,l-1,N,c,lowend,highend, iintercept, islope)+idangle+kdangle+2*data->eparam[5]+2*data->eparam[10]
		      +2*data->eparam[6]);
//	      if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 1 11 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" c "<<c<<"\n";
	      //case12 1 0 1 1
	      en1=min(en1,branch_1(single1_w,w,i+2,j-1,k+2,l-2,N,c,lowend,highend, iintercept, islope)+idangle+kdangle+ldangle+gap+2*data->eparam[5]+2*data->eparam[10]
		      +3*data->eparam[6]);
//	            if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 1 12 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" c "<<c<<"\n";
	      //case13 1 1 0 0
	      en1=min(en1,branch_1(single1_w,w,i+2,j-2,k+1,l-1,N,c,lowend,highend, iintercept, islope)+idangle+jdangle+2*gap+2*data->eparam[5]+2*data->eparam[10]
	            +2*data->eparam[6]);
//	           if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 1 13 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" c "<<c<<"\n";
	      //case14 1 1 0 1
	      en1=min(en1,branch_1(single1_w,w,i+2,j-2,k+1,l-2,N,c,lowend,highend, iintercept, islope)+idangle+jdangle+ldangle+gap+2*data->eparam[5]+2*data->eparam[10]
		      +3*data->eparam[6]);
//	          if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 1 14 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" c "<<c<<"\n";
	      //case15 1 1 1 0
	      en1=min(en1,branch_1(single1_w,w,i+2,j-2,k+2,l-1,N,c,lowend,highend, iintercept, islope)+idangle+jdangle+kdangle+gap+2*data->eparam[5]+2*data->eparam[10]
		      +3*data->eparam[6]);
//	      if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 1 15 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" c "<<c<<"\n";
	      //case16 1 1 1 1
	      en1=min(en1,branch_1(single1_w,w,i+2,j-2,k+2,l-2,N,c,lowend,highend, iintercept, islope)+idangle+jdangle+kdangle+ldangle+2*data->eparam[5]+2*data->eparam[10]
		      +4*data->eparam[6]);
//	            if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 1 16 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" c "<<c<<"\n";
	 }
   

	 for (d=k+1+minloop;d+1<l-1-minloop&&d!=N2;++d) {
	   	   //case1 0 0 0 0
	    en1=min(en1,branch_2(single2_w,w,i+1,j-1,k+1,l-1,N2,d,lowend,highend, iintercept, islope)+2*data->eparam[5]+2*data->eparam[10]);
//	        if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 2 1 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" d "<<d<<"\n";
	   //case2 0 0 0 1
	    en1=min(en1,branch_2(single2_w,w,i+1,j-1,k+1,l-2,N2,d,lowend,highend, iintercept, islope)+ldangle+gap+2*data->eparam[5]+2*data->eparam[10]
		    +data->eparam[6]);
//	           if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 2 2 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" d "<<d<<"\n";
	    //case3 0 0 1 0
	    en1=min(en1,branch_2(single2_w,w,i+1,j-1,k+2,l-1,N2,d,lowend,highend, iintercept, islope)+kdangle+gap+2*data->eparam[5]+2*data->eparam[10]
		     +data->eparam[6]);
//	      if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 2 3 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" d "<<d<<"\n";
	    //case4 0 0 1 1
	    en1=min(en1,branch_2(single2_w,w,i+1,j-1,k+2,l-2,N2,d,lowend,highend, iintercept, islope)+kdangle+ldangle+2*gap+2*data->eparam[5]+2*data->eparam[10]
		     +2*data->eparam[6]);
//	         if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 2 4 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" d "<<d<<"\n";
	    //case5 0 1 0 0
	    en1=min(en1,branch_2(single2_w,w,i+1,j-2,k+1,l-1,N2,d,lowend,highend, iintercept, islope)+jdangle+gap+2*data->eparam[5]+2*data->eparam[10]
		     +data->eparam[6]);
//	       if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 2 5 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" d "<<d<<"\n";
	     //case6 0 1 0 1
	     en1=min(en1,branch_2(single2_w,w,i+1,j-2,k+1,l-2,N2,d,lowend,highend, iintercept, islope)+jdangle+ldangle+2*data->eparam[5]+2*data->eparam[10]
		     +2*data->eparam[6]);
//	              if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 2 6 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" d "<<d<<"\n";
	     //case7 0 1 1 0
	     en1=min(en1,branch_2(single2_w,w,i+1,j-2,k+2,l-1,N2,d,lowend,highend, iintercept, islope)+jdangle+kdangle+2*gap+2*data->eparam[5]+2*data->eparam[10]
		     +2*data->eparam[6]);
//	       if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 2 7 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" d "<<d<<"\n";
	      //case8 0 1 1 1
	     en1=min(en1,branch_2(single2_w,w,i+1,j-2,k+2,l-2,N2,d,lowend,highend, iintercept, islope)+jdangle+kdangle+ldangle+gap+2*data->eparam[5]+2*data->eparam[10]
		       +3*data->eparam[6]);
//	          if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 2 8 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" d "<<d<<"\n";
	      //case9 1 0 0 0
	     en1=min(en1,branch_2(single2_w,w,i+2,j-1,k+1,l-1,N2,d,lowend,highend, iintercept, islope)+idangle+gap+2*data->eparam[5]+2*data->eparam[10]
		      +data->eparam[6]);
//	              if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 2 9 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" d "<<d<<"\n";
	     //case10 1 0 0 1
	     en1=min(en1,branch_2(single2_w,w,i+2,j-1,k+1,l-2,N2,d,lowend,highend, iintercept, islope)+idangle+ldangle+2*gap+2*data->eparam[5]+2*data->eparam[10]
		      +2*data->eparam[6]);
//	        if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 2 10 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" d "<<d<<"\n";
	      //case11 1 0 1 0
	      en1=min(en1,branch_2(single2_w,w,i+2,j-1,k+2,l-1,N2,d,lowend,highend, iintercept, islope)+idangle+kdangle+2*data->eparam[5]+2*data->eparam[10]
		      +2*data->eparam[6]);
//	          if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 2 11 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" d "<<d<<"\n";
	      //case12 1 0 1 1
	      en1=min(en1,branch_2(single2_w,w,i+2,j-1,k+2,l-2,N2,d,lowend,highend, iintercept, islope)+idangle+kdangle+ldangle+gap+2*data->eparam[5]+2*data->eparam[10]
		      +3*data->eparam[6]);
//	         if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 2 12 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" d "<<d<<"\n";
	      //case13 1 1 0 0
	      en1=min(en1,branch_2(single2_w,w,i+2,j-2,k+1,l-1,N2,d,lowend,highend, iintercept, islope)+idangle+jdangle+2*gap+2*data->eparam[5]+2*data->eparam[10]
		      +2*data->eparam[6]);
//	          if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 2 13 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" d "<<d<<"\n";
	      //case14 1 1 0 1
	      en1=min(en1,branch_2(single2_w,w,i+2,j-2,k+1,l-2,N2,d,lowend,highend, iintercept, islope)+idangle+jdangle+ldangle+gap+2*data->eparam[5]+2*data->eparam[10]
		      +3*data->eparam[6]);
//	         if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 2 14 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" d "<<d<<"\n";
	      //case15 1 1 1 0
	      en1=min(en1,branch_2(single2_w,w,i+2,j-2,k+2,l-1,N2,d,lowend,highend, iintercept, islope)+idangle+jdangle+kdangle+gap+2*data->eparam[5]+2*data->eparam[10]
		      +3*data->eparam[6]);
//	         if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 2 15 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" d "<<d<<"\n";
	      //case16 1 1 1 1
	      en1=min(en1,branch_2(single2_w,w,i+2,j-2,k+2,l-2,N2,d,lowend,highend, iintercept, islope)+idangle+jdangle+kdangle+ldangle+2*data->eparam[5]+2*data->eparam[10]
		      +4*data->eparam[6]);
//	         if(i==91&&j==203&&k==109&&l==195)cerr<<"en1 2 16 "<<en1+ penalty(i,j,ct1,data) + penalty(k,l,ct2,data)<<" d "<<d<<"\n";
	 }

       }
#else
#endif
       en1 = en1 + penalty(i,j,ct1,data) + penalty(k,l,ct2,data);
    
       v->f(i,j,k,l)=min(en1,v->f(i,j,k,l));
       en3 = DYNALIGN_INFINITY;
       en3m = DYNALIGN_INFINITY; // Only used if force is on, but it's
        // simpler just to do the assign than to
        // check force.
       //     if(i==91&&j==203&&k==109&&l==195)cerr<<"v->f(i,j,k,l) "<<v->f(i,j,k,l)<<"\n";

        //first precalculate some internal loop free energies:
        for (c=i+1;c<=i+maxloop&&c<j/*-minloop*/&&c<=N;++c) {
          for (d=j-1;d>=j-maxloop&&d>c/*+minloop*/;--d) {
            if (c-i+j-d-2<maxinternal) {
              einternal[c-i][j-d] = erg2(i, j, c, d, ct1, data,
                                         force ? fce1[jref(i,c,N)][iref(i,c,N)] : 0,
                                         force ? fce1[jref(d,j,N)][iref(d,j,N)] : 0);
            }
          }
        }

        for (c=k+1;c<=k+maxloop&&c<l/*-minloop*/&&c<=N2;++c) {
          for (d=l-1;d>=l-maxloop&&d>c/*+minloop*/;--d) {
            if (c-k+l-d-2<maxinternal) {
              einternal2[c-k][l-d] = erg2(k, l, c, d, ct2, data,
                                          force ? fce2[jref(k,c,N2)][iref(k,c,N2)] : 0,
                                          force ? fce2[jref(d,l,N2)][iref(d,l,N2)] : 0);
            }
          }
        }
#ifdef DYNALIGN_II
        //Now consider internal loops/one stacked and other not
        for (c=i;c<=i+maxloop&&c<j/*-minloop*/&&c<=N;++c) {
          for (d=j;d>=j-maxloop&&d>c/*+minloop*/;--d) {//it seems there should be another constrain for d that when j>N, d>N report to Dave16
            //if (d>N) ap = Ndiff;
            //else ap=0;
            starte=max(k,lowend[c]);
            ende=min(min(min(k+maxloop,l),highend[c]),N2);
            for (e=starte/*max(k+1,c-maxsep)*/;e<=ende/*k+maxloop&&e<l&&(e-c)<=maxsep&&e<=N2*/;++e) {
              endf=min(l,highend[d]);
              startf=max(max((l-maxloop),e+1),lowend[d]);
              for (f=endf/*min(l-1,d+maxsep-ap)*/;f>=startf/*l-maxloop&&f>e&&f>=d-maxsep-ap*/;--f) {

                if (c-i+j-d-2<maxinternal&&e-k+l-f-2<maxinternal&&(((d>N)&&(f>N2))||((d<=N)&&(f<=N2)))) {
                  {
                    // Only used if force is on, but it's simpler just to do
                    // the assign than to check force.
                    helicalstack1=false;
                    helicalstack2=false;
		    stuck=false;
                  }

		  if(c==i&&d==j&&(e!=k&&f!=l)){
		    stuck=true;
		    en1=0;
		    if (force && modification) {
		      en1m = en1;
		     }
		  }

		  else if (c==i+1&&d==j-1&&(i!=N)&&(j-1!=N)) {//ask Dave is j-1!=N safe, will it assign an internal loop value to a stacking base pair? report to Dave13
		    en1 = erg1(i,j,i+1,j-1,ct1, data);//en1 is the energy of stacking base pair in seq1

		    if (force && modification) {
		      //calculate en1m
		      if (mod1[i]||mod1[j]) {
			//check for a GU exception
			if ((ct1->numseq[i]==3&&ct1->numseq[j]==4)||
			    (ct1->numseq[i]==4&&ct1->numseq[j]==3)||
			    (ct1->numseq[i+1]==3&&ct1->numseq[j-1]==4)||
			    (ct1->numseq[i+1]==4&&ct1->numseq[j-1]==3)||
			    (ct1->numseq[i-1]==3&&ct1->numseq[j+1]==4)||
			    (ct1->numseq[i-1]==4&&ct1->numseq[j+1]==3))//report to Dave10 what if i-1 and j+1 are unpaired nucleotides?
			  en1m=en1;
			else en1m=DYNALIGN_INFINITY;
		      }
		      else en1m=en1;
		      helicalstack1=true;
		    }
                  }
                  else if(c!=i&&d!=j){
		      //internal loop
		      en1 = einternal[c-i][j-d];
		      if (force && modification) {
			en1m=en1;
		      }
		    }
		  else continue;
		  
		 

		  if(e==k&&f==l&&(c!=i&&d!=j)){
		    stuck=true;
		    en2=0;
		    if (force && modification) {
		      en2m = en2;
		    }
		  }
		  
		  else if(e==k+1&&f==l-1&&(k!=N2)&&(l-1!=N2)){
		    en2 = erg1(k,l,k+1,l-1,ct2, data);
		    if (force && modification) {
		      if (mod2[k]||mod2[l]) {
			//check for a GU exception
			if ((ct2->numseq[k]==3&&ct2->numseq[l]==4)||
			    (ct2->numseq[k]==4&&ct2->numseq[l]==3)||
			    (ct2->numseq[k+1]==3&&ct2->numseq[l-1]==4)||
			    (ct2->numseq[k+1]==4&&ct2->numseq[l-1]==3)||
			    (ct2->numseq[k-1]==3&&ct2->numseq[l+1]==4)||
			    (ct2->numseq[k-1]==4&&ct2->numseq[l+1]==3))
			  en2m = en2;
			else en2m = DYNALIGN_INFINITY;
		      }
		      else en2m = en2;
		      helicalstack2=true;
		    }
		  }
		  
		  else if(e!=k&&f!=l){
		    en2 = einternal2[e-k][l-f];
		    if (force && modification) {
		      en2m = en2;
		    }
		  }
		  else continue;

                  if (force) {
                    if (modification) {
                      if (helicalstack1||helicalstack2) {
                        en3 = min(en3,en1+en2+vmod->f(c,d,e,f)/*vmod[d][iref(c,d,N)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l)));//report to Dave12 if it's modification, why assign any value to en3 at all? and what's the point of the meaning of this?
                        en3m = min(en3m,en1m+en2m+vmod->f(c,d,e,f)/*vmod[d][iref(c,d,N)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l)));
                      }
		      //report to Dave15 newly added stuff recommend Dave that we should write two separate vmod arrays.
		      else if(stuck){
			en3 = min(en3,en1+en2+v->f(c,d,e,f)/*v[d][iref(c,d,N)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l)));
                        en3m = min(en3m,en1+en2+vmod->f(c,d,e,f)/*v[d][iref(c,d,N)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l)));
		      }
		      
                      else{
                        en3 = min(en3,en1+en2+v->f(c,d,e,f)/*v[d][iref(c,d,N)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l)));
                        en3m = min(en3m,en1+en2+v->f(c,d,e,f)/*v[d][iref(c,d,N)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l)));
		
                      }

                    } else {
                      en3 = min(en3,en1+en2+v->f(c,d,e,f)/*v[d][iref(c,d,N)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l)));
                    }
                  } else {
                    en3 = min(en3,en1+en2+v->f(c,d,e,f)/*v[d][iref(c,d,N)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l)));
		    //		    if(i==20&&j==41&&k==23&&l==38&&en3==-1328){cerr<<c<<" "<<d<<" "<<e<<" "<<f<<" "<<en3<<"\n";}
                  }
                }
              }
            }
          }
        }
	//	if(i==62&&j==120&&k==61&&l==130)cerr<<"v->f(i,j,k,l) "<<v->f(i,j,k,l)<<"\n";
//	if(i==91&&j==203&&k==109&&l==195)cerr<<"v->f(i,j,k,l) "<<en3<<"\n";
#else
	        //Now consider internal loops/one stacked and other not
//	bool singleinsert=true;
	for (c=i+1;c<=i+maxloop&&c<j/*-minloop*/&&c<=N;++c) {
          for (d=j-1;d>=j-maxloop&&d>c/*+minloop*/;--d) {
            //if (d>N) ap = Ndiff;
            //else ap=0;
            starte=max(k+1,lowend[c]);
            ende=min(min(min(k+maxloop,l),highend[c]),N2);
            for (e=starte/*max(k+1,c-maxsep)*/;e<=ende/*k+maxloop&&e<l&&(e-c)<=maxsep&&e<=N2*/;++e) {
              endf=min(l-1,highend[d]);
              startf=max(max((l-maxloop),e+1),lowend[d]);
              for (f=endf/*min(l-1,d+maxsep-ap)*/;f>=startf/*l-maxloop&&f>e&&f>=d-maxsep-ap*/;--f) {

                if (c-i+j-d-2<maxinternal&&e-k+l-f-2<maxinternal&&(((d>N)&&(f>N2))||((d<=N)&&(f<=N2)))) {

                  {
                    // Only used if force is on, but it's simpler just to do
                    // the assign than to check force.
                    helicalstack1=false;
                    helicalstack2=false;
                  }

                  if (c==i+1&&d==j-1&&(i!=N)&&(j-1!=N)) {
                    //seq1 is helical stack
                    en1 = erg1(i,j,i+1,j-1,ct1, data);

                    if (force && modification) {
                      //calculate en1m
                      if (mod1[i]||mod1[j]) {
                        //check for a GU exception
                        if ((ct1->numseq[i]==3&&ct1->numseq[j]==4)||
                            (ct1->numseq[i]==4&&ct1->numseq[j]==3)||
                            (ct1->numseq[i+1]==3&&ct1->numseq[j-1]==4)||
                            (ct1->numseq[i+1]==4&&ct1->numseq[j-1]==3)||
                            (ct1->numseq[i-1]==3&&ct1->numseq[j+1]==4)||
                            (ct1->numseq[i-1]==4&&ct1->numseq[j+1]==3))
                          en1m=en1;
                        else en1m=DYNALIGN_INFINITY;
                      }
                      else en1m=en1;
                      helicalstack1=true;
                    }
                  }
                  else {
                    //internal loop
                    en1 = einternal[c-i][j-d];
                    //erg2(i,j,c,d,ct1,data,fce1[jref(i,c,N)][iref(i,c,N)],fce1[jref(d,j,N)][iref(d,j,N)]);
                    //erg2(i,j,c,d,ct1,data,0,0);

                    if (force && modification) {
                      en1m=en1;
                    }
                  }

                  if (e==k+1&&f==l-1&&k!=N2&&(l-1!=N2)) {
                    //seq2 is helical stack
                    en2 = erg1(k,l,k+1,l-1,ct2, data);
                    if (force && modification) {
                      if (mod2[k]||mod2[l]) {
                        //check for a GU exception
                        if ((ct2->numseq[k]==3&&ct2->numseq[l]==4)||
                            (ct2->numseq[k]==4&&ct2->numseq[l]==3)||
                            (ct2->numseq[k+1]==3&&ct2->numseq[l-1]==4)||
                            (ct2->numseq[k+1]==4&&ct2->numseq[l-1]==3)||
                            (ct2->numseq[k-1]==3&&ct2->numseq[l+1]==4)||
                            (ct2->numseq[k-1]==4&&ct2->numseq[l+1]==3))
                          en2m = en2;
                        else en2m = DYNALIGN_INFINITY;
                      }
                      else en2m = en2;
                      helicalstack2=true;
                    }

                    // also allow single base pair insertions into one sequence only
                    if (singleinsert&&c==i+2&&d==j-2&&j-2>i+2&&
                        inc[ct1->numseq[i+1]][ct1->numseq[j-1]]&&
                        inc[ct1->numseq[i+2]][ct1->numseq[j-2]]&&i+1!=N&&
                        (j-1!=N)&&(j-2!=N)

                        && (!force ||
                            (!(fce1[jref(i+1,j-1,N)][iref(i+1,j-1,N)] & SINGLE) &&
                             !(fce1[jref(i+1,j-1,N)][iref(i+1,j-1,N)] & NOPAIR) &&
                             !mod1[i+1] &&
                             !mod1[j-1]
                             )
                            )
                        ) {

                      en1 = min(en1,erg1(i,j,i+1,j-1,ct1, data)+
                                erg1(i+1,j-1,i+2,j-2,ct1,data));
                    }
                  } else {
                    //internal loop

                    en2 = einternal2[e-k][l-f];
                    //erg2(k,l,e,f,ct2,data,fce2[jref(k,e,N2)][iref(k,e,N2)],fce2[jref(f,l,N2)][iref(f,l,N2)]);
                    //erg2(k,l,e,f,ct2,data,0,0);
                    if (force && modification) {
                      en2m = en2;
                    }

                    //also allow single base pair insertions into one sequence only
                    if (singleinsert&&e==k+2&&f==l-2&&l-2>k+2&&c==i+1&&d==j-1&&
                        inc[ct2->numseq[k+1]][ct2->numseq[l-1]]&&
                        inc[ct2->numseq[k+2]][ct2->numseq[l-2]]&&
                        inc[ct1->numseq[i+1]][ct1->numseq[j-1]]&&
                        i!=N&&(j-1!=N)&&k!=N2&&(l-1!=N2)&&
                        (k+1!=N2)&&(l-2!=N2)

                        && (!force ||
                            (!(fce2[jref(k+1,l-1,N2)][iref(k+1,l-1,N2)] & SINGLE) &&
                             !(fce2[jref(k+1,l-1,N2)][iref(k+1,l-1,N2)] & NOPAIR) &&
                             !mod2[k+1] &&
                             !mod2[l-1]
                             )
                            )
                        ) {

                      en2 = min(en2,erg1(k,l,k+1,l-1,ct2, data)+
                                erg1(k+1,l-1,k+2,l-2,ct2,data));
                    }
                  }

                  if (force) {
                    if (modification) {
                      if (helicalstack1||helicalstack2) {
                        en3 = min(en3,en1+en2+vmod->f(c,d,e,f)/*vmod[d][iref(c,d,N)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l)));
                        en3m = min(en3m,en1m+en2m+vmod->f(c,d,e,f)/*vmod[d][iref(c,d,N)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l)));
                      }
                      else {
                        en3 = min(en3,en1+en2+v->f(c,d,e,f)/*v[d][iref(c,d,N)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l)));
                        en3m = min(en3m,en1+en2+v->f(c,d,e,f)/*v[d][iref(c,d,N)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l)));
                      }

                    } else {
                      en3 = min(en3,en1+en2+v->f(c,d,e,f)/*v[d][iref(c,d,N)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l)));
                    }
                  } else {
                    en3 = min(en3,en1+en2+v->f(c,d,e,f)/*v[d][iref(c,d,N)][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l)));
                  }
                }
              }
            }
          }
        }
#endif
#ifdef DYNALIGN_II

	int stacking_elongation_1=0;
	int stacking_elongation_2=0;
	int alignment_elongation=0;
	
	for(c=1;c<=max_elongation;c++){
	  
	  if((c+1>maxloop)||(i+c+1>=j-c-1)||(i+c+1>N)||(k+c+1>=l-c-1)||(k+c+1>N2))break;
	  if(j>N)if((j-c-1<=N+1)||(l-c-1<=N2+1))break;
	  
	  if(k+c<lowend[i+c]||k+c>highend[i+c]||l-c<lowend[j-c]||l-c>highend[j-c])break;
	  if(k+c+1<lowend[i+c+1]||k+c+1>highend[i+c+1]||l-c-1<lowend[j-c-1]||l-c-1>highend[j-c-1])break;
	  if (alignmentforced) {
	    if ((forcealign[0][i+c]>0&&forcealign[0][i+c]!=k+c)||
		(forcealign[0][j-c]>0&&forcealign[0][j-c]!=l-c)||
		(forcealign[1][k+c]>0&&forcealign[1][k+c]!=i+c)||
		(forcealign[1][l-c]>0&&forcealign[1][l-c]!=j-c)) {
	      break;
	    }
	    
	    else {
        //scan through to make sure that the base pair isn't precluding an alignment
	      for (g=i+c+1;g<j-c;g++) {
		if ((forcealign[0][g]>0)&&(forcealign[0][g]<k+c||forcealign[0][g]>l-c)) {
		  break;
		}
	      }
	      
	      for (g=k+c+1;g<l-c;g++) {
		if ((forcealign[1][g]>0)&&(forcealign[1][g]<i+c||forcealign[1][g]>j-c)) {
		  break;
		}
	      }
	    }

	    if ((forcealign[0][i+c+1]>0&&forcealign[0][i+c+1]!=k+c+1)||
		(forcealign[0][j-c-1]>0&&forcealign[0][j-c-1]!=l-c-1)||
		(forcealign[1][k+c+1]>0&&forcealign[1][k+c+1]!=i+c+1)||
		(forcealign[1][l-c-1]>0&&forcealign[1][l-c-1]!=j-c-1)) {
	      break;
	    }
	    
	    else {
        //scan through to make sure that the base pair isn't precluding an alignment
	      for (g=i+c+2;g<j-c-1;g++) {
		if ((forcealign[0][g]>0)&&(forcealign[0][g]<k+c+1||forcealign[0][g]>l-c-1)) {
		  break;
		}
	      }
	      
	      for (g=k+c+2;g<l-c-1;g++) {
		if ((forcealign[1][g]>0)&&(forcealign[1][g]<i+c+1||forcealign[1][g]>j-c-1)) {
		  break;
		}
	      }
	    }

	  }

	  alignment_elongation=c;
	  
	}

       	for(c=1;c<=alignment_elongation;c++){
	  if(!inc[ct1->numseq[i+c]][ct1->numseq[j-c]])break;
	  if (force)if(fce1[jref(i+c,j-c,N)][iref(i+c,j-c,N)]&SINGLE||fce1[jref(i+c,j+c,N)][iref(i+c,j+c,N)]&NOPAIR)break;
	  
	  if (force && modification) {
	    if (mod1[i+c]||mod1[j-c]) {
	      //check for a GU exception
	      if ((ct1->numseq[i+c]==3&&ct1->numseq[j-c]==4)||
		  (ct1->numseq[i+c]==4&&ct1->numseq[j-c]==3)||
		  (ct1->numseq[i+c+1]==3&&ct1->numseq[j-c-1]==4)||
		  (ct1->numseq[i+c+1]==4&&ct1->numseq[j-c-1]==3)||
		  (ct1->numseq[i+c-1]==3&&ct1->numseq[j-c+1]==4)||
		  (ct1->numseq[i+c-1]==4&&ct1->numseq[j-c+1]==3));
	      else break;
	    }
	  }
	  
	  if (ct1->templated) {
	    if (j<=N) {
	      d=j-c;
	      e=i+c;
	    }
	    else {
	      e=j-c-N;
	      d=i+c;
	    }
	    
	    if (!ct1->tem[d][e]) {
	      break;
	    }
	  }
	  stacking_elongation_1=c;	  
	  
	  

	}

	    
       	for(c=1;c<=alignment_elongation;c++){
	  if(!inc[ct2->numseq[k+c]][ct2->numseq[l-c]])break;
	  if (force)if(fce2[jref(k+c,l-c,N2)][iref(k+c,l-c,N2)]&SINGLE||fce2[jref(k+c,l-c,N2)][iref(k+c,l-c,N2)]&NOPAIR)break;
	  
	  if (force && modification) {
	    if (mod2[k+c]||mod2[l-c]) {
	      //check for a GU exception
	      if ((ct1->numseq[k+c]==3&&ct1->numseq[l-c]==4)||
		  (ct1->numseq[k+c]==4&&ct1->numseq[l-c]==3)||
		  (ct1->numseq[k+c+1]==3&&ct1->numseq[l-c-1]==4)||
		  (ct1->numseq[k+c+1]==4&&ct1->numseq[l-c-1]==3)||
		  (ct1->numseq[k+c-1]==3&&ct1->numseq[l-c+1]==4)||
		  (ct1->numseq[k+c-1]==4&&ct1->numseq[l-c+1]==3));
	      else break;
	    }
	  }

	  if (ct2->templated) {
	    if (j<=N) {
	      d=l-c;
	      e=k+c;
	    }
	    else {
	      e=l-c-N2;
	      d=k+c;
	    }
	    
	    if (!ct2->tem[d][e]) {
	      break;
	    }
	  }
	  stacking_elongation_2=c;	  
	  
	}	  

	//	cerr<<"i "<<i<<" j "<<j<<" k "<<k<<" l "<<l<<" stacking_elongation_1 "<<stacking_elongation_1<<" stacking_elongation_2 "<<stacking_elongation_2<<"\n";
	//	if(i==20&&j==41&&k==23&&l==38){cerr<<stacking_elongation_1<<" "<<stacking_elongation_2<<" "<<"\n";}
	en1_temp=0;
	for(c=1;c<=stacking_elongation_1;c++){
	  if((!inc[ct1->numseq[i+c+1]][ct1->numseq[j-c-1]])||(!inc[ct2->numseq[k+c+1]][ct2->numseq[l-c-1]]))continue;
	  en1_temp+=erg1(i+c-1,j-c+1,i+c,j-c,ct1,data);
	  en1=en1_temp+erg1(i+c,j-c,i+c+1,j-c-1,ct1,data);
	  en2=einternal2[c+1][c+1];
	  
	  if(force&&modification){
	    en3=min(en3,en1+en2+vmod->f(i+c+1,j-c-1,k+c+1,l-c-1));
	      if (mod1[i]||mod1[j]) {
			//check for a GU exception
		if ((ct1->numseq[i]==3&&ct1->numseq[j]==4)||
		    (ct1->numseq[i]==4&&ct1->numseq[j]==3)||
		    (ct1->numseq[i+1]==3&&ct1->numseq[j-1]==4)||
		    (ct1->numseq[i+1]==4&&ct1->numseq[j-1]==3)||
		    (ct1->numseq[i-1]==3&&ct1->numseq[j+1]==4)||
		    (ct1->numseq[i-1]==4&&ct1->numseq[j+1]==3))//report to Dave10 what if i-1 and j+1 are unpaired nucleotides?
		  en3m=min(en3m,en1+en2+vmod->f(i+c+1,j-c-1,k+c+1,l-c-1));
	      }
	      else en3m=min(en3m,en1+en2+vmod->f(i+c+1,j-c-1,k+c+1,l-c-1));
	  }

	 
	  else {
	    //  cerr<<"en1 "<<en1<<" en2 "<<en2<<" i "<<i<<" j "<<j<<" k "<<k<<" l "<<l<<" c "<<c<<"\n";
	    en3=min(en3,en1+en2+v->f(i+c+1,j-c-1,k+c+1,l-c-1));
	   }
	  // if(i==20&&j==41&&k==23&&l==38){cerr<<en3<<"\n";}
	}

	en1_temp=0;
	for(c=1;c<=stacking_elongation_2;c++){
	  if((!inc[ct1->numseq[i+c+1]][ct1->numseq[j-c-1]])||(!inc[ct2->numseq[k+c+1]][ct2->numseq[l-c-1]]))continue;
	  en1_temp+=erg1(k+c-1,l-c+1,k+c,l-c,ct2,data);
	  en2=en1_temp+erg1(k+c,l-c,k+c+1,l-c-1,ct2,data);
	  en1=einternal[c+1][c+1];
	  
	  if(force&&modification){
	    en3=min(en3,en1+en2+vmod->f(i+c+1,j-c-1,k+c+1,l-c-1));
	    if (mod1[k]||mod1[l]) {
			//check for a GU exception
	      if ((ct2->numseq[k]==3&&ct2->numseq[l]==4)||
		  (ct2->numseq[k]==4&&ct2->numseq[l]==3)||
		  (ct2->numseq[k+1]==3&&ct2->numseq[l-1]==4)||
		  (ct2->numseq[k+1]==4&&ct2->numseq[l-1]==3)||
		  (ct2->numseq[k-1]==3&&ct2->numseq[l+1]==4)||
		  (ct2->numseq[k-1]==4&&ct2->numseq[l+1]==3))//report to Dave10 what if i-1 and j+1 are unpaired nucleotides?
		en3m=min(en3m,en1+en2+vmod->f(i+c+1,j-c-1,k+c+1,l-c-1));
	    }
	    else en3m=min(en3m,en1+en2+vmod->f(i+c+1,j-c-1,k+c+1,l-c-1));
	  }  
	  else{
	    //  cerr<<"en1 "<<en1<<" en2 "<<en2<<" i "<<i<<" j "<<j<<" k "<<k<<" l "<<l<<" c "<<c<<"\n";
	     en3=min(en3,en1+en2+v->f(i+c+1,j-c-1,k+c+1,l-c-1));
	  }
	}
#else
#endif

        if (force && modification) {
          /*if (j<=N)*/ vmod->f(i,j,k,l)/*vmod[j][i][a][b]*/=min(en3m,v->f(i,j,k,l)/*v[j][i][a][b]*/);
          //else vmod[j][i-j+N][a][b]=min(en3m,v[j][i-j+N][a][b]);
        }

        /*if (j<=N)*/ v->f(i,j,k,l)/*v[j][i][a][b]*/=min(en3,v->f(i,j,k,l)/*v[j][i][a][b]*/);
        //else v[j][i-j+N][a][b]=min(en3,v[j][i-j+N][a][b]);
//	if(i==91&&j==203&&k==109&&l==195)cerr<<"v->f(i,j,k,l) "<<v->f(i,j,k,l)<<"\n";
	
	

    }
    else {
      v->f(i,j,k,l)/*v[j][iref(i,j,N)][a][b]*/ = DYNALIGN_INFINITY;
      //a bp wasn't allowed between either i and j or k and l
    }



  /////////////////////////////////////////////////////////
  //Fill w, the best energy for a fragment in a multiloop:

  //save a temporaray value for w in en1, a register short int


  //consider the possibilities for adding nucleotides to existing fragments
  //		i	j	k	l

  //2		0	0	0	1
  //3		0	0	1	0
  //4		0	0	1	1z
  //5		0	1	0	0
  //6		0	1	0	1z
  //7		0	1	1	0
  //8		0	1	1	1
  //9		1	0	0	0
  //10	1	0	0	1
  //11	1	0	1	0z
  //12	1	0	1	1
  //13	1	1	0	0z
  //14	1	1	0	1
  //15	1	1	1	0
  //16	1	1	1	1z

  //remember than j<=N is included fragment and j>N is excluded fragment
  en1=DYNALIGN_INFINITY;

  //case 6 0	1	0	1
  if (j-1 != N && l-1 != N2 && jldecrement &&
      (!force ||
       (!dbl1[j] && !dbl2[l])
       )
      ) {
    en1 = w->f(i,j-1,k,l-1)/*w[j-1][iref(i,j-1,N)][a][b]*/+2*data->eparam[6];

    //case 16 1	1	1	1
    if (j-1 > i+1 && i != N && k != N2 &&
        ikincrement && jldecrement &&
        (!force ||
         (!dbl1[i] && !dbl2[k])
         )
        ) {
      en1 = min(en1,w->f(i+1,j-1,k+1,l-1)/*w[j-1][iref(i+1,j-1,N)][a][b]*/+ 4*data->eparam[6]);
    }
  }

  //case 11	1	0	1	0
  if (i != N && k != N2 && ikincrement &&
      (!force ||
       (!dbl1[i] && !dbl2[k])
       )
      ) {
    en1 = min(en1,w->f(i+1,j,k+1,l)/*w[j][iref(i+1,j,N)][a][b]*/+2*data->eparam[6]);
  }
  if (k>=lowend[i+1]&&k<=highend[i+1]/*a>=1*/) {
    //case 9		1	0	0	0
    if (i != N &&
        (!force ||
         (!dbl1[i])
         )
        ) {
      en1 = min(en1,w->f(i+1,j,k,l)/*w[j][iref(i+1,j,N)][a-1][b]*/+data->eparam[6]+gap);

      //case 14	1	1	0	1
      if ((j-1>i+1)&&(j-1!=N)) {
        if (l-1 != N2 && jldecrement &&
            (!force ||
             (!dbl1[j]&&!dbl2[l])
             )
            ) {
          en1 = min(en1,w->f(i+1,j-1,k,l-1)/*w[j-1][iref(i+1,j-1,N)][a-1][b]*/+3*data->eparam[6]+gap);
        }
        if (l <= highend[j-1]&&l>=lowend[j-1]/*b+1<(2*maxseparation+2)*/ &&
            (!force ||
             (!dbl1[j])
             )
            ) {

          //case 13 1	1	0	0
          en1 = min(en1,w->f(i+1,j-1,k,l)/*w[j-1][iref(i+1,j-1,N)][a-1][b+1]*/+2*data->eparam[6]+2*gap);
        }
      }
    }

    if (l-1 >= lowend[j]/*b>=1*/ && i != N && l-1 != N2 &&
        (!force ||
         (!dbl1[i] && !dbl2[l])
         )
        ) {
      //case 10 1	0	0	1
      en1 = min(en1,w->f(i+1,j,k,l-1)/*w[j][iref(i+1,j,N)][a-1][b-1]*/+2*data->eparam[6]+2*gap);
    }
  }
  if (j-1 != N &&
      (l <= highend[j-1] && l>=lowend[j-1])) {
    //case 5		0	1	0	0
    if (!force ||
        (!dbl1[j])
        ) {
      en1 = min(en1,w->f(i,j-1,k,l)/*w[j-1][iref(i,j-1,N)][a][b+1]*/+data->eparam[6]+gap);
    }
    //case 15	1	1	1	0
    if (k!=N2) {
      if (j-1 > i+1 && i != N && ikincrement &&
          (!force ||
           (!dbl1[i] && !dbl1[j] && !dbl2[k])
           )
          ) {
        en1 = min(en1,w->f(i+1,j-1,k+1,l)/*w[j-1][iref(i+1,j-1,N)][a][b+1]*/+3*data->eparam[6]+gap);
      }
      if (k+1 <= highend[i]/*a+1<(2*maxseparation+2)*/ &&
          (!force ||
           (!dbl1[j] && !dbl2[k])
           )
          ) {
        //case 7		0	1	1	0
        en1=min(en1,w->f(i,j-1,k+1,l)/*w[j-1][iref(i,j-1,N)][a+1][b+1]*/+2*data->eparam[6]+2*gap);

      }
    }
  }
  if (k+1 <= highend[i]/*a+1<(2*maxseparation+2)*/ && k != N2 &&
      (!force ||
       (!dbl2[k])
       )
      ) {
    //case 3		0	0	1	0
    en1 = min(en1,w->f(i,j,k+1,l)/*w[j][iref(i,j,N)][a+1][b]*/+data->eparam[6]+gap);

    if (l-1 != N2 &&
        (!force ||
         (!dbl2[l])
         )
        ) {
      //case 8		0	1	1	1
      if (j-1 != N && jldecrement &&
          (!force ||
           (!dbl1[j])
           )
          ) {
        en1 = min(en1,w->f(i,j-1,k+1,l-1)/*w[j-1][iref(i,j-1,N)][a+1][b]*/+3*data->eparam[6]+gap);
      }
      if (l-1 >= lowend[j]/*b>=1*/) {
        //case 4		0	0	1	1
        en1 = min(en1,w->f(i,j,k+1,l-1)/*w[j][iref(i,j,N)][a+1][b-1]*/+2*data->eparam[6]+2*gap);
      }
    }
  }
  if (l-1 >= lowend[j]/*(b>=1)*/ && l-1 != N2 &&
      (!force ||
       (!dbl2[l])
       )
      ) {
    //case 2		0	0	0	1
    en1=min(en1,w->f(i,j,k,l-1)/*w[j][iref(i,j,N)][a][b-1]*/+data->eparam[6]+gap);

    //case12
    //12	1	0	1	1
    if (i != N && k != N2 && ikincrement &&
        (!force ||
         (!dbl1[i] && !dbl2[k])
         )
        ) {
      en1 = min(en1,w->f(i+1,j,k+1,l-1)/*w[j][iref(i+1,j,N)][a][b-1]*/+3*data->eparam[6]+gap);
    }
  }

  //Consider the case where none of the four nucs (i,j,k,l) are paired
  //Consider whether any of the 4 nucleotides are stacked on a helix
  //There are 16 cases:
  //consider them in this order (0 unstacked, 1 stacked)
  //		i	j	k	l
  //1		0	0	0	0
  //2		0	0	0	1
  //3		0	0	1	0
  //4		0	0	1	1
  //5		0	1	0	0
  //6		0	1	0	1
  //7		0	1	1	0
  //8		0	1	1	1
  //9		1	0	0	0
  //10	1	0	0	1
  //11	1	0	1	0
  //12	1	0	1	1
  //13	1	1	0	0
  //14	1	1	0	1
  //15	1	1	1	0
  //16	1	1	1	1


  //note:
  //	a = k-i+maxsep;
  //	b = l-j+maxsep;

  //	so when addressing i+1 => address a-1 to keep k unchanged
  //	and simlarly: j-1 => b+1 to keep l unchanged

  idangle=edangle5(i+1,j,i,ct1,data);
  jdangle=edangle3(j-1,i,j,ct1,data);
  kdangle=edangle5(k+1,l,k,ct2,data);
  ldangle=edangle3(l-1,k,l,ct2,data);
  ijdangle=edangle5(i+1,j-1,i,ct1,data)+edangle3(j-1,i+1,j,ct1,data);
  kldangle=edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data);

  //case 1 - nothing stacked:
  en1=min(en1,v->f(i,j,k,l)/*v[j][iref(i,j,N)][a][b]*/+2*data->eparam[10]+penalty(i,j,ct1,data)+penalty(k,l,ct2,data));

  //case 6 - j and l stacked
  if ((j-1!=N)&&(l!=N2+1)&&jldecrement) {
    en1 = min(en1,v->f(i,j-1,k,l-1)/*v[j-1][iref(i,j-1,N)][a][b]*/+2*data->eparam[6]+jdangle
              +ldangle+2*data->eparam[10]+penalty(i,j-1,ct1,data)+penalty(k,l-1,ct2,data));
  }
  //case 11 - i and k stacked
  if ((i!=N)&&(k!=N2)&&ikincrement) {
    en1 = min(en1,v->f(i+1,j,k+1,l)/*v[j][iref(i+1,j,N)][a][b]*/+2*data->eparam[6]+kdangle+
              idangle+2*data->eparam[10]+penalty(i+1,j,ct1,data)+penalty(k+1,l,ct2,data));

    //case 16 - i, j, k, and l stacked
    if (j-1>i+1&&(j-1!=N)&&(l!=N2+1)&&ikincrement&&jldecrement) {
      en1=min(en1,v->f(i+1,j-1,k+1,l-1)/*v[j-1][iref(i+1,j-1,N)][a][b]*/+4*data->eparam[6]
              +kldangle
              +ijdangle
              +2*data->eparam[10]
              +penalty(i+1,j-1,ct1,data)+penalty(k+1,l-1,ct2,data));
    }
  }

  if (l-1>=lowend[j]/*(b-1>=0)*/&&(l!=N2+1)) {
    //case 2 - l stacked
    en1=min(en1,v->f(i,j,k,l-1)/*v[j][iref(i,j,N)][a][b-1]*/+data->eparam[6]+ldangle+2*data->eparam[10]+gap+
            penalty(i,j,ct1,data)+penalty(k,l-1,ct2,data));

    if (k+1<=highend[i]/*(a+1<2*maxsep+2)*/&&(k!=N2)) {
      //case 4 - l and k stacked
      en1=min(en1,v->f(i,j,k+1,l-1)/*v[j][iref(i,j,N)][a+1][b-1]*/+2*data->eparam[6]
              +kldangle+2*data->eparam[10]+2*gap+
              penalty(i,j,ct1,data)+penalty(k+1,l-1,ct2,data));
    }
  }

  if (k+1<=highend[i]/*(a+1<2*maxsep+2)*/&&(k!=N2)) {



    //case 3 - k stacked
    en1 = min(en1,v->f(i,j,k+1,l)/*v[j][iref(i,j,N)][a+1][b]*/+data->eparam[6]+kdangle+2*data->eparam[10]+gap
              +penalty(i,j,ct1,data)+penalty(k+1,l,ct2,data));

    if (j!=N+1) {
      if (l<=highend[j-1]&&l>=lowend[j-1]) {
        //case 7 - j and k stacked
        en1 = min(en1,v->f(i,j-1,k+1,l)/*v[j-1][iref(i,j-1,N)][a+1][b+1]*/+2*data->eparam[6]+jdangle
                  +kdangle+2*data->eparam[10]+2*gap+
                  penalty(i,j-1,ct1,data)+penalty(k+1,l,ct2,data));

      }

      //case 8 - j, k, and l stacked:
      if (l!=N2+1&&jldecrement) {
        en1 = min(en1,v->f(i,j-1,k+1,l-1)/*v[j-1][iref(i,j-1,N)][a+1][b]*/+3*data->eparam[6]+jdangle
                  +kldangle+2*data->eparam[10]+gap
                  +penalty(i,j-1,ct1,data)+penalty(k+1,l-1,ct2,data));
      }

    }


  }

  if (l<=highend[j-1]&&l>=lowend[j-1]&&(j!=N+1)) {
    //case 5 - j stacked
    en1 = min(en1,v->f(i,j-1,k,l)/*v[j-1][iref(i,j-1,N)][a][b+1]*/+data->eparam[6]+jdangle+2*data->eparam[10]+gap
              +penalty(i,j-1,ct1,data)+penalty(k,l,ct2,data));

    //case 15 - i,j, and k stacked
    if ((j-1>i+1)&&(i!=N)) {
      if (k!=N2&&ikincrement) {
        en1 = min(en1,v->f(i+1,j-1,k+1,l)/*v[j-1][iref(i+1,j-1,N)][a][b+1]*/+3*data->eparam[6]+kdangle
                  +ijdangle
                  +2*data->eparam[10]+gap
                  +penalty(i+1,j-1,ct1,data)+penalty(k+1,l,ct2,data));
      }


      //case 13 - i and j stacked


      if (k>=lowend[i+1]&&k<=highend[i+1]) {
        en1 = min(en1,v->f(i+1,j-1,k,l)/*v[j-1][iref(i+1,j-1,N)][a-1][b+1]*/+2*data->eparam[6]+ijdangle
                  +2*data->eparam[10]+2*gap
                  +penalty(i+1,j-1,ct1,data)+penalty(k,l,ct2,data));
      }

    }

  }

  if (k>=lowend[i+1]&&k<=highend[i+1]&&(i!=N)) {
    //case 9 - i alone is stacked
    en1 = min(en1,v->f(i+1,j,k,l)/*v[j][iref(i+1,j,N)][a-1][b]*/+data->eparam[6]+idangle+2*data->eparam[10]+gap
              +penalty(i+1,j,ct1,data)+penalty(k,l,ct2,data));


    //case 14 - i, j, and l stacked
    if ((j-1>i+1)&&(j!=N+1)&&(l!=N2+1)&&jldecrement) {
      en1=min(en1,v->f(i+1,j-1,k,l-1)/*v[j-1][iref(i+1,j-1,N)][a-1][b]*/+3*data->eparam[6]+ldangle
              +ijdangle
              +2*data->eparam[10]+gap
              +penalty(i+1,j-1,ct1,data)+penalty(k,l-1,ct2,data));
    }
  }


  if (l-1>=lowend[j]/*b-1>=0*/&&(i!=N)&&(l!=N2+1)) {
    //case 12 - i, k, and l stacked:
    if (k!=N2&&ikincrement) {
      en1=min(en1,v->f(i+1,j,k+1,l-1)/*v[j][iref(i+1,j,N)][a][b-1]*/+3*data->eparam[6]+kldangle
              +idangle+2*data->eparam[10]+gap
              +penalty(i+1,j,ct1,data)+penalty(k+1,l-1,ct2,data));
    }
    if (k>=lowend[i+1]&&k<=highend[i+1]) {
      //case 10 - l and i stacked
      en1=min(en1,v->f(i+1,j,k,l-1)/*v[j][iref(i+1,j,N)][a-1][b-1]*/+2*data->eparam[6]+ldangle
              +idangle+2*data->eparam[10]+2*gap+
              penalty(i+1,j,ct1,data)+penalty(k,l-1,ct2,data));
    }
  }

  //calculate the free energy of 2 fragments merged:
  for (c=i+minloop;c<j-minloop;++c) {
    //if (c>N) {
    //	kp = Ndiff;
    //}
    //else kp=0;
    startd=max(k+minloop,lowend[c]);
    endd=min(l-minloop,highend[c]);
    for (d=startd/*max(k+minloop,c-maxsep-kp)*/;
         d<=endd/*d<l-minloop&&d<=c+maxsep-kp*/;++d) {

      //e = d-c+maxsep+kp;

      if ((c!=N)&&(d!=N2)&&d+1>=lowend[c+1]&&d+1<=highend[c+1]) {
        if (c<N&&d<N2) {
          en1 = min(en1,w->f(i,c,k,d)+w->f(c+1,j,d+1,l)/*w[c][iref(i,c,N)][a][e]+w[j][iref(c+1,j,N)][e][b]*/);
        }
        else {
          if (d>N2&&c>N) {
            en1 = min(en1,w->f(i,c,k,d)+w->f(c+1,j,d+1,l)/*w[c][iref(i,c,N)][a][e]+w[j-N][iref(c+1-N,j-N,N)][e][b]*/);
          }
        }
	//report to Dave6 what if c and d bigger than N and N2, the indexes are gonna be screwed up
	//report to Dave8 Dec.26 what about alignment envelope in w array?
      }

    }
  }
 
  //report to Dave5
#ifdef DYNALIGN_II

  for (c=i+minloop;c+1<j-minloop&&c!=N;++c) {
    //  if((c<N)&&k>=lowend[c+1]&&k<=highend[c+1])en1=min(en1,w->f(c+1,j,k,l)+single1_w.f(i,c)+(c-i+1)*gap);
      en1=min(en1,branch_1(single1_w,w,i,j,k,l,N,c,lowend,highend, iintercept, islope));
    
  }
  // if(i==28&&j==45&&k==28&&l==60)cerr<<"w->f(i,j,k,l) "<<en1<<"\n";
  for (d=k+minloop;d+1<l-minloop&&d!=N2;++d) {
  //  if((c<N)&&k>=lowend[c+1]&&k<=highend[c+1])en1=min(en1,w->f(c+1,j,k,l)+single1_w.f(i,c)+(c-i+1)*gap);
    en1=min(en1,branch_2(single2_w,w,i,j,k,l,N2,d,lowend,highend, iintercept, islope));

  }
#else
#endif






  w->f(i,j,k,l)=en1;

} //end of dynalignstep


//traceback a single conserved structure (either interior or exterior fragments)
//return an int that indicates whether an error occurred
#ifdef DYNALIGN_II
int  dyntrace(short i, short j, short a, short b, structure *ct1, structure *ct2,
              short structnum, short *alignment,
              dynalignarray *w, varray *v, wendarray *w3, wendarray *w5,
        short *lowend, short *highend,datatable *data, short islope, short iintercept, short gap, dynalignarray *vmod,DynProgArray<integersize> *single1_w,DynProgArray<integersize> *single2_w,DynProgArray<integersize> *single1_v,DynProgArray<integersize> *single2_v,DynProgArray<integersize> *single1_wmb,DynProgArray<integersize> *single2_wmb,DynProgArray<integersize> *single1_we,DynProgArray<integersize> *single2_we,bool *single1_lfce,bool *single2_lfce,bool *single1_mod,bool *single2_mod,forceclass *single1_fce,forceclass *single2_fce,
              bool local, int max_elongation,bool *mod1, bool *mod2, bool modification,
	      char **fce1, char **fce2, bool alignmentforced, short **forcealign,bool force,bool startopen)
#else
int  dyntrace(short i, short j, short a, short b, structure *ct1, structure *ct2,
              short structnum, short *alignment,
              dynalignarray *w, varray *v, wendarray *w3, wendarray *w5,
              short *lowend, short *highend,datatable *data, short gap, dynalignarray *vmod,
              bool local, bool startopen)
#endif

 {

  short int (*edangle5)(int i, int j, int ip, structure *ct, datatable *data) = edangle5noforce;
  short int (*edangle3)(int i, int j, int ip, structure *ct, datatable *data) = edangle3noforce;

  short k,l,c,d,en1,en2,e,f,g,kp,ip,jp,dp,lp;
  stackclass single1_stack;
  stackclass single2_stack;
  dynalignstackclass stack;
  bool open,closed,if_vmod;
  integersize en3;
#ifdef DYNALIGN_II
  integersize en1_temp;
#else
#endif
  register bool found;
  int constantclosure;
  bool ikincrement,jldecrement,ikincrement2,jldecrement2;
  int error = 0;


  //modification tracks whether chemical modification considerations apply
  bool helicalstack,stuck;

#ifndef DYNALIGN_II
  bool modification;
  modification = ( (ct1->GetNumberofModified() > 0) || (ct2->GetNumberofModified() > 0) );
#else
#endif
  vector< vector<bool> > inc = data->pairing;

  //put the whole fragment on the stack
  //a and b originally referred to positions in the large arrays -
  //now a and b refer to the nucleotide positions
  if (startopen) stack.push(j,b,1,0,w5->f(j,b),startopen);
  else stack.push(i,j,a,b,v->f(i,j,a,b),startopen);


#ifdef DYNALIGN_II
  while (stack.pull(&i,&j, &k, &l /*&a,&b*/, &en3, &open, &if_vmod))
#else
  while (stack.pull(&i,&j, &k, &l /*&a,&b*/, &en3, &open))
#endif
 {
    //cout << "pulled "<<i<< " "<<j<<" "<<a<<" "<<b<<"\n"<<flush;
//       cerr<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<en3<<" "<<open<<" "<<if_vmod<<" trace!\n";


    if (!open) {
      if (modification) {
        if (en3==v->f(i,j,k,l)/*v[j][iref(i,j,ct1->GetSequenceLength())][a][b]*/) closed = true;
        else if (en3==vmod->f(i,j,k,l)/*vmod[j][iref(i,j,ct1->GetSequenceLength())][a][b]*/) closed = true;
        else closed = false;
	
      }
      else {

        if (en3==v->f(i,j,k,l)) closed = true;
        else closed = false;

      }

    }

    if (open) {
      //This is from w3 or w5:
      if (k==ct1->GetSequenceLength()) {
        ip = i;
        kp = j;
        a = i;
        b = j;
	

	//ip and kp are now indexes of w3

        //ap = j;
        //we are dealing with a fragment ending in ct1->GetSequenceLength(), i.e. w3
        while (ip<ct1->GetSequenceLength()) {
			  found = false;
			//First check if there is no structure in the 3' ends.
			//This stopping rule is needed because the normal whittle away at the ends might
				//not work with the lowend/highend scheme for allowing aligned nucleotides.
			if (local) {
				  //local alignment calculation
				if(en3==0) {
					found = true;
					ip=ct1->GetSequenceLength();

				}

			}
			else {
				 //global alignment calculation
				if (en3==gap*(abs((ct1->GetSequenceLength()-ip)-(ct2->GetSequenceLength()-kp)))) {
					found =true;
					ip=ct1->GetSequenceLength();

				}

			}




          //check if it is safe to increment i and k
          ikincrement =
            kp+1 >= lowend[ip+1] &&
            kp+1 <= highend[ip+1];
	  //report to Dave6 why do kp+1 and ip+1 need to be aligned?
          //kp = ip + ap -maxsep;
          if (!found) {
            if (en3 == w3->f(ip+1,kp+1)) {//adding one nuc on both that is unpaired and unstacked
              found = true;
              en3 = w3->f(ip+1,kp+1);
              ip++;
              kp++;//

            }
          }

          if (!found){
            if (en3 == w3->f(ip+1,kp)+gap) {
              found = true;
              en3 = w3->f(ip+1,kp);
              ip++;


            }

          }


          if (!found) {
            if (en3==w3->f(ip,kp+1)+gap){
              found = true;
              en3 = w3->f(ip,kp+1);
              kp++;

            }
          }

          if (local) {

            if (!found){
              if (en3 == w3->f(ip+1,kp)) {
                found = true;
                en3 = w3->f(ip+1,kp);
                ip++;


              }

            }


            if (!found) {
              if (en3==w3->f(ip,kp+1)){
                found = true;
                en3 = w3->f(ip,kp+1);
                kp++;

              }
            }

          }



          for (jp=ct1->GetSequenceLength();jp>=ip+minloop;jp--) {
            dp = max(1/*lowlimit(jp,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/,kp+minloop);
            for (lp=ct2->GetSequenceLength();lp>=dp&&!found;lp--) {
              //bp = lp-jp+maxsep;

              //check whether mine[i][a] is split so that the lowest free energy is
              //an exterior fragment from 1 to j and a helix from j+1 to i;

              //check all possible alignments (in index a(k) and b(l))

              //must consider whether an unpaired nucleotide is stacked
              //is stacked onto each of the nucleotides in a helix
              //There are 16 cases:
              //consider them in this order (0 unstacked, 1 stacked)
              //		jp+1	ip	lp+1	k
              //1		0		0	0		0
              //2		0		0	0		1
              //3		0		0	1		0
              //4		0		0	1		1
              //5		0		1	0		0
              //6		0		1	0		1
              //7		0		1	1		0
              //8		0		1	1		1
              //9		1		0	0		0
              //10		1		0	0		1
              //11		1		0	1		0
              //12		1		0	1		1
              //13		1		1	0		0
              //14		1		1	0		1
              //15		1		1	1		0
              //16		1		1	1		1

              //note that for these exterior loops:
              //j<i
              //l<k

              //if (lp+1>=lowend[jp+1]&&lp+1<=highend[jp+1]/*lowlimit(jp+1,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/) {

                jldecrement =
                  lp-1 <= highend[jp-1]/*highlimit(jp-1, maxsep, ct1->GetSequenceLength(), ct2->GetSequenceLength())*/ &&
                  lp-1 >= lowend[jp-1]/*lowlimit(jp-1, maxsep, ct1->GetSequenceLength(), ct2->GetSequenceLength())*/;

				bool ipkp = kp>=lowend[ip]&&kp<=highend[ip];
				bool jplp = lp>=lowend[jp]&&lp<=highend[jp];

                //no stacking
                //case 1:
				if (ipkp&&jplp) {
					if (en3==w3->f(jp+1,lp+1)+v->f(ip,jp,kp,lp)+
						penalty(ip,jp,ct1,data)+penalty(kp,lp,ct2,data)) {

					  found = true;
					  stack.push(ip,jp,kp,lp,v->f(ip,jp,kp,lp));
					  en3=w3->f(jp+1,lp+1);
					  ip=jp+1;
					  kp=lp+1;;

					}
				}

                //case 6:
                if (!found&&ikincrement&&jplp) if (en3==w3->f(jp+1,lp+1)+v->f(ip+1,jp,kp+1,lp)+edangle5(ip+1,jp,ip,ct1,data)+
                                             edangle5(kp+1,lp,kp,ct2,data)+
                                             penalty(ip+1,jp,ct1,data)+penalty(kp+1,lp,ct2,data)){

                  found = true;
                  stack.push(ip+1,jp,kp+1,lp,v->f(ip+1,jp,kp+1,lp));
                  en3 = w3->f(jp+1,lp+1);
                  ip=jp+1;
                  kp=lp+1;


                }

                //case 11	1	0	1	0
                if (!found&&jldecrement&&ipkp) if(en3==w3->f(jp+1,lp+1)+v->f(ip,jp-1,kp,lp-1)+edangle3(jp-1,ip,jp,ct1,data)+
                                            edangle3(lp-1,kp,lp,ct2,data)+penalty(ip,jp-1,ct1,data)+penalty(lp-1,kp,ct2,data)) {

                  found = true;
                  stack.push(ip,jp-1,kp,lp-1,v->f(ip,jp-1,kp,lp-1));
                  en3 = w3->f(jp+1,lp+1);
                  ip = jp+1;
                  kp=lp+1;

                }


                //case 16	1	1	1	1
                if (!found&&ikincrement&&jldecrement) if(en3==w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/+v->f(ip+1,jp-1,kp+1,lp-1)/*v[jp-1][ip+1][ap][bp]*/+edangle3(jp-1,ip+1,jp,ct1,data)+
                                                         edangle5(ip+1,jp-1,ip,ct1,data)+edangle3(lp-1,kp+1,lp,ct2,data)+
                                                         edangle5(kp+1,lp-1,kp,ct2,data)+penalty(ip+1,jp-1,ct1,data)+penalty(lp-1,kp+1,ct2,data)) {

                  found = true;
                  stack.push(ip+1,jp-1,kp+1,lp-1,v->f(ip+1,jp-1,kp+1,lp-1));
                  en3 = w3->f(jp+1,lp+1);
                  ip = jp+1;
                  kp=lp+1;

                }


                if (kp+1<=highend[ip]&&kp+1>=lowend[ip]/*highlimit(ip,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/&&!found) {//report to Dave Nov.6 adding kp+1>=lowend[ip]
                  //case 2:
                  if (jplp) if (en3==w3->f(jp+1,lp+1)+v->f(ip,jp,kp+1,lp)+edangle5(kp+1,lp,kp,ct2,data)+
                      penalty(ip,jp,ct1,data)+penalty(kp+1,lp,ct2,data)+gap){

		      found = true;
		      stack.push(ip,jp,kp+1,lp,v->f(ip,jp,kp+1,lp));
		      en3=w3->f(jp+1,lp+1);
		      ip = jp+1;
		      kp=lp+1;
		      
		    }

                  //case 12	1	0	1	1
                  if (!found&&jldecrement) if(en3==w3->f(jp+1,lp+1)+v->f(ip,jp-1,kp+1,lp-1)+edangle3(jp-1,ip,jp,ct1,data)+
                                              edangle3(lp-1,kp+1,lp,ct2,data) + edangle5(kp+1,lp-1,kp,ct2,data)+//yinghan: report to Dave Nov.6 change kp->kp+1
                                              penalty(ip,jp-1,ct1,data)+penalty(kp+1,lp-1,ct2,data)+gap){

                    found = true;
                    stack.push(ip,jp-1,kp+1,lp-1,v->f(ip,jp-1,kp+1,lp-1));
                    en3 = w3->f(jp+1,lp+1);
                    ip = jp+1;
                    kp=lp+1;


                  }

                  if (lp-1>=lowend[jp]&&lp-1<=highend[jp]/*lowlimit(jp,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/&&!found) {
                    //case 4
                    if (en3==w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/+v->f(ip,jp,kp+1,lp-1)/*v[jp][ip][ap+1][bp-1]*/+edangle5(kp+1,lp-1,kp,ct2,data)+
                        edangle3(lp-1,kp+1,lp,ct2,data)+penalty(ip,jp,ct1,data)+penalty(kp+1,lp-1,ct2,data)+2*gap) {

                      found = true;
                      stack.push(ip,jp,kp+1/*ap+1*/,lp-1/*bp-1*/,v->f(ip,jp,kp+1,lp-1)/*v[jp][ip][ap+1][bp-1]*/);
                      en3 = w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/;
                      ip = jp+1;
                      kp=lp+1;//ap = bp;

                    }
                  }

                }


                if (lp-1>=lowend[jp]&&lp-1<=highend[jp]/*lowlimit(jp,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/&&!found) {
                  //case 3:
                  if (ipkp) if (en3==w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/+v->f(ip,jp,kp,lp-1)/*v[jp][ip][ap][bp-1]*/+edangle3(lp-1,kp,lp,ct2,data)+
                      penalty(ip,jp,ct1,data)+penalty(lp-1,kp,ct2,data)+gap) {

                    found = true;
                    stack.push(ip,jp,kp/*ap*/,lp-1/*bp-1*/,v->f(ip,jp,kp,lp-1)/*v[jp][ip][ap][bp-1]*/);
                    en3=w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/;
                    ip = jp+1;
                    kp=lp+1;//ap = bp;

                  }

                  //case 8:		0	1	1	1
                  if (!found&&ikincrement) if (en3==w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/+v->f(ip+1,jp,kp+1,lp-1)/*v[jp][ip+1][ap][bp-1]*/+edangle5(ip+1,jp,ip,ct1,data)+
                                               edangle3(lp-1,kp+1,lp,ct2,data)+edangle5(kp+1,lp-1,kp,ct2,data)+
                                               penalty(ip+1,jp,ct1,data)+penalty(lp-1,kp+1,ct2,data)+gap) {

                    found = true;
                    stack.push(ip+1,jp,kp+1/*ap*/,lp-1/*bp-1*/,v->f(ip+1,jp,kp+1,lp-1)/*v[jp][ip+1][ap][bp-1]*/);
                    en3=w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/;
                    ip = jp+1;
                    kp=lp+1;//ap = bp;


                  }

                  if (kp>=lowend[ip+1]&&kp<=highend[ip+1]&&!found) {
                    //case 7		0	1	1	0
                    if (en3==w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/+v->f(ip+1,jp,kp,lp-1)/*v[jp][ip+1][ap-1][bp-1]*/+edangle3(lp-1,kp,lp,ct2,data)+
                        edangle5(ip+1,jp,ip,ct1,data)+penalty(ip+1,jp,ct1,data)+penalty(lp-1,kp,ct2,data)+2*gap) {

                      found = true;
                      stack.push(ip+1,jp,kp/*ap-1*/,lp-1/*bp-1*/,v->f(ip+1,jp,kp,lp-1)/*v[jp][ip+1][ap-1][bp-1]*/);
                      en3=w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/;
                      ip = jp+1;
                      kp=lp+1;//ap = bp;

                    }
                  }


                }


                if (kp>=lowend[ip+1]&&kp<=highend[ip+1]&&!found) {
                  //case5: 0 1 0 0
                  if (jplp) if (en3==w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/+v->f(ip+1,jp,kp,lp)/*v[jp][ip+1][ap-1][bp]*/+edangle5(ip+1,jp,ip,ct1,data)+
                      penalty(ip+1,jp,ct1,data)+penalty(kp,lp,ct2,data)+gap) {

                    found = true;
                    stack.push(ip+1,jp,kp/*ap-1*/,lp/*bp*/,v->f(ip+1,jp,kp,lp)/*v[jp][ip+1][ap-1][bp]*/);
                    en3=w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/;
                    ip = jp+1;
                    kp=lp+1;//ap = bp;

                  }

                  //case 15	1	1	1	0
                  if (!found&&jldecrement) if (en3==w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/+v->f(ip+1,jp-1,kp,lp-1)/*v[jp-1][ip+1][ap-1][bp]*/ + edangle3(jp-1,ip+1,jp,ct1,data)+
                                               edangle5(ip+1,jp-1,ip,ct1,data)+edangle3(lp-1,kp,lp,ct2,data)+
                                               penalty(ip+1,jp-1,ct1,data)+penalty(lp-1,kp,ct2,data)+gap) {

                    found = true;
                    stack.push(ip+1,jp-1,kp/*ap-1*/,lp-1/*bp*/,v->f(ip+1,jp-1,kp,lp-1)/*v[jp-1][ip+1][ap-1][bp]*/);
                    en3=w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/;
                    ip = jp+1;
                    kp=lp+1;//ap = bp;

                  }

                  if (lp<=highend[jp-1]&&lp>=lowend[jp-1]/*highlimit(jp-1,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/&&!found) {
                    //case 13	1	1	0	0
                    if (en3==w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/+v->f(ip+1,jp-1,kp,lp)/*v[jp-1][ip+1][ap-1][bp+1]*/+edangle3(jp-1,ip+1,jp,ct1,data)+
                        edangle5(ip+1,jp-1,ip,ct1,data)+
                        penalty(ip+1,jp-1,ct1,data)+penalty(kp,lp,ct2,data)+2*gap) {


                      found = true;
                      stack.push(ip+1,jp-1,/*ap-1*/kp,lp/*bp+1*/,v->f(ip+1,jp-1,kp,lp)/*v[jp-1][ip+1][ap-1][bp+1]*/);
                      en3=w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/;
                      ip = jp+1;
                      kp=lp+1;//ap = bp;

                    }

                  }


                }
                if (lp<=highend[jp-1]&&lp>=lowend[jp-1]/*highlimit(jp-1,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/&&!found) {
                  //case 9		1	0	0	0
                  if(ipkp) if (en3==w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/+v->f(ip,jp-1,kp,lp)/*v[jp-1][ip][ap][bp+1]*/+edangle3(jp-1,ip,jp,ct1,data)+
                      penalty(ip,jp-1,ct1,data)+penalty(kp,lp,ct2,data)+gap) {

                    found = true;
                    stack.push(ip,jp-1,kp/*ap*/,lp/*bp+1*/,v->f(ip,jp-1,kp,lp)/*v[jp-1][ip][ap][bp+1]*/);
                    en3=w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/;
                    ip = jp+1;
                    kp=lp+1;//ap = bp;

                  }


                  //case 14	1	1	0	1
                  if (!found&&ikincrement) if (en3==w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/+v->f(ip+1,jp-1,kp+1,lp)/*v[jp-1][ip+1][ap][bp+1]*/+edangle3(jp-1,ip+1,jp,ct1,data)+
                                               edangle5(ip+1,jp-1,ip,ct1,data)+edangle5(kp+1,lp,kp,ct2,data)+
                                               penalty(ip+1,jp-1,ct1,data)+penalty(kp+1,lp,ct2,data)+gap) {

                    found = true;
                    stack.push(ip+1,jp-1,kp+1/*ap*/,lp/*bp+1*/,v->f(ip+1,jp-1,kp+1,lp)/*v[jp-1][ip+1][ap][bp+1]*/);
                    en3=w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/;
                    ip=jp+1;
                    kp=lp+1;//ap = bp;

                  }


                  if (kp+1<=highend[ip]&&kp+1>=lowend[ip]/*highlimit(ip,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/&&!found) {//report to Dave Nov. 6 add kp+1>=lowend[ip]
                    //case 10	1	0	0	1
                    if (en3==w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/+v->f(ip,jp-1,kp+1,lp)/*v[jp-1][ip][ap+1][bp+1]*/+edangle3(jp-1,ip,jp,ct1,data)+
                        edangle5(kp+1,lp,kp,ct2,data)+penalty(ip,jp-1,ct1,data)+penalty(kp+1,lp,ct2,data)+2*gap) {


                      found = true;
                      stack.push(ip,jp-1,kp+1/*ap+1*/,lp/*bp+1*/,v->f(ip,jp-1,kp+1,lp)/*v[jp-1][ip][ap+1][bp+1]*/);
                      en3=w3->f(jp+1,lp+1)/*w3[jp+1][bp]*/;
                      ip=jp+1;
                      kp=lp+1;//ap = bp;

                    }

                  }



                }


              }
            //}
          }
	
	//report to Dave4
#ifdef DYNALIGN_II

	    for (jp=ct1->GetSequenceLength();jp>=ip+minloop&&!found;jp--) {
	      if(!found&&en3==single1_we->f(ip,jp)+w3->f(jp+1,kp)+gap_branch(jp-ip+1)){
//		cerr<<"w31 ip="<<ip<<" jp="<<jp<<" kp="<<kp<<"\n";
		trace(ct1,data,ip,jp,single1_v,single1_w,single1_wmb,NULL,NULL,single1_lfce,single1_fce,NULL,NULL,single1_mod,single1_we,single1_we->f(ip,jp),1,0);
		found=true;
		en3=w3->f(jp+1,kp);
		ip=jp+1;
	      }
	      
	    }
	    
	    for (lp=ct2->GetSequenceLength();lp>=kp+minloop&&!found;lp--) {
	      if(!found&&en3==single2_we->f(kp,lp)+w3->f(ip,lp+1)+gap_branch(lp-kp+1)){
		
//		cerr<<"w32 ip="<<ip<<" lp="<<lp<<" kp="<<kp<<"\n";
		trace(ct2,data,kp,lp,single2_v,single2_w,single2_wmb,NULL,NULL,single2_lfce,single2_fce,NULL,NULL,single2_mod,single2_we,single2_we->f(kp,lp),1,0);
		found=true;
		en3=w3->f(ip,lp+1);
		kp=lp+1;
	      }
	    }	

#else
#endif	  
	
	if (!found) {
		  //A traceback error occurred
//		  cerr << "Traceback error at w3\n";
		  error = 14;
                  //	  return error;
		  break;

		}
 

		}
	    
	
	
      }


      
      else {
        //dealing with w5

        //a = j;
        k = j;
        while(i>0) {
          found = false;
          //k = a+i-maxsep;
//	  cerr<<"w5->(i,k) "<<i<<" "<<k<<"\n";
	  //i and k are indexes of w5 array


          //although this is really a decrement...
	  if(i>1){
	    ikincrement =
	      k-1 <= highend[i-1]/*highlimit(i-1, maxsep, ct1->GetSequenceLength(), ct2->GetSequenceLength())*/ &&
	      k-1 >= lowend[i-1]/*lowlimit(i-1, maxsep, ct1->GetSequenceLength(), ct2->GetSequenceLength())*/;
	  }
		  //check if there is no structure in the 5' direction.
		  //This is important because the usual whittle away might not work with the highlimit/lowlimit tracking
			//of what can be aligned...
		  if (local) {

			  if (en3==0) {

				  //If energy is zero, no need to traceback farther
				  found = true;
				  i=0;
			  }

		  }
		  else {



			  if (en3 == gap* (abs(i - k))) {

				  //If energy is accounted with gaps, then no need for more traceback
				  found = true;
				  i=0;

			  }

		  }


          if (!found) {
            if (en3==w5->f(i-1,k-1)/*w5[i-1][a]*/) {
              found = true;
              en3 = w5->f(i-1,k-1)/*w5[i-1][a]*/;
              i--;
              k--;//

            }
          }

		  if (!found) {
			  //if (k<=highend[i-1]&&k>=lowend[i-1]) {
				if (en3==w5->f(i-1,k)/*w5[i-1][a+1]*/+gap){
				  found = true;
				  en3 = w5->f(i-1,k);/*w5[i-1][a+1]*/
				  i--;
				  //a++;


				}

			  //}
		  }



          //if (k-1>=lowend[i]&&k-1<=highend[i]&&!found) {
		  if(!found){
		    if(en3==w5->f(i,k-1)/*w5[i][a-1]*/+gap) {
		      found = true;
		      en3 = w5->f(i,k-1)/*w5[i][a-1]*/;
		      k--;//a--;
		      
		    }
		  }
          //}

          else if (local) {
            if (!found&&i>0) {
              if (en3==w5->f(i-1,k)){
                found = true;
                en3 = w5->f(i-1,k);
                i--;



              }

            }



            //if (k-1>=lowend[i]&&k-1<=highend[i]&&!found) {

              if(en3==w5->f(i,k-1)) {
                found = true;
                en3 = w5->f(i,k-1);
                k--;

              }

            //}

          }

          for (j=0;j+minloop<i&&!found;j++) {
            for (l=0;l+minloop<k&&!found;l++) {
              //b = l-j+maxsep;

              //check whether w5[i][a] is split so that the lowest free energy is
              //an exterior fragment from 1 to j and a helix from j+1 to i;

              //check all possible alignments (in index a(k) and b(l))



              //no stacking
              //if(en3== w5[j][b]+v[i][j+1][b][a]) {
              //	found = true;

              //	stack.push(j+1,i,b,a,v[i][j+1][b][a]);
              //	en3 = w5[j][b];
              //	i = j;
              //	a = b;


              //}

              //check whether w5[i][a] is split so that the lowest free energy is
              //an exterior fragment from 1 to j and a helix from j+1 to i;

              //check all possible alignments (in index a(k) and b(l))

              //must consider whether an unpaired nucleotide is stacked
              //is stacked onto each of the nucleotides in a helix
              //There are 16 cases:
              //consider them in this order (0 unstacked, 1 stacked)
              //		j+1	i	l+1	k
              //1		0	0	0	0
              //2		0	0	0	1
              //3		0	0	1	0
              //4		0	0	1	1
              //5		0	1	0	0
              //6		0	1	0	1
              //7		0	1	1	0
              //8		0	1	1	1
              //9		1	0	0	0
              //10		1	0	0	1
              //11		1	0	1	0
              //12		1	0	1	1
              //13		1	1	0	0
              //14		1	1	0	1
              //15		1	1	1	0
              //16		1	1	1	1

              //note that for these exterior loops:
              //j<i
              //l<k

              //although this is really an increment...
              jldecrement =
                l+1 >= lowend[j+1]/*lowlimit(j+1, maxsep, ct1->GetSequenceLength(), ct2->GetSequenceLength())*/ &&
                l+1 <= highend[j+1]/*highlimit(j+1, maxsep, ct1->GetSequenceLength(), ct2->GetSequenceLength())*/;
              //although this is really an increment...
              jldecrement2 =
                l+2 >= lowend[j+2]/*lowlimit(j+2, maxsep, ct1->GetSequenceLength(), ct2->GetSequenceLength())*/ &&
                l+2 <= highend[j+2]/*highlimit(j+2, maxsep, ct1->GetSequenceLength(), ct2->GetSequenceLength())*/;

			  bool ik = k >= lowend[i] && k <=highend[i];
			  //bool jl = l >= lowend[j] && l <=highend[j];

              //no stacking
              //case 1:
              if (jldecrement&&ik) {
                if (en3 == w5->f(j,l)/*w5[j][b]*/+v->f(j+1,i,l+1,k)/*v[i][j+1][b][a]*/+penalty(i,j+1,ct1,data)+penalty(k,l+1,ct2,data)) {
                  found = true;
                  en3 = w5->f(j,l)/*w5[j][b]*/;
                  stack.push(j+1,i,l+1/*b*/,k/*a*/,v->f(j+1,i,l+1,k)/*v[i][j+1][b][a]*/);
                  i=j;
                  k=l;
                  //a=b;

                }
              }
              //case 6:

              if (!found&&i>1&&ikincrement&&jldecrement) {
                if(en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+1,i-1,l+1,k-1)/*v[i-1][j+1][b][a]*/+edangle3(i-1,j+1,i,ct1,data)+
                   edangle3(k-1,l+1,k,ct2,data)+penalty(i-1,j+1,ct1,data)+penalty(k-1,l+1,ct2,data)){

                  found = true;
                  en3 = w5->f(j,l)/*w5[j][b]*/;
                  stack.push(j+1,i-1,l+1/*b*/,k-1/*a*/,v->f(j+1,i-1,l+1,k-1)/*v[i-1][j+1][b][a]*/);
                  i = j;
                  k=l;
                  //a=b;
                }
              }

              //case 11	1	0	1	0
              if (!found&&jldecrement2&&ik) {
                if (en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+2,i,l+2,k)/*v[i][j+2][b][a]*/+edangle5(j+2,i,j+1,ct1,data)+
                    edangle5(l+2,k,l+1,ct2,data)+penalty(i,j+2,ct1,data)+penalty(l+2,k,ct2,data)) {

                  stack.push(j+2,i,l+2/*b*/,k/*a*/,v->f(j+2,i,l+2,k)/*v[i][j+2][b][a]*/);
                  found = true;
                  en3=w5->f(j,l)/*w5[j][b]*/;
                  i = j;
                  k=l;
                  //a=b;


                }
              }


              //case 16	1	1	1	1
              if (!found&&i>1&&ikincrement&&jldecrement2) {
                if(en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+2,i-1,l+2,k-1)/*v[i-1][j+2][b][a]*/+edangle5(j+2,i-1,j+1,ct1,data)+
                   edangle3(i-1,j+2,i,ct1,data)+edangle5(l+2,k-1,l+1,ct2,data)+
                   edangle3(k-1,l+2,k,ct2,data)+penalty(i-1,j+2,ct1,data)+penalty(k-1,l+2,ct2,data)) {

                  stack.push(j+2,i-1,l+2/*b*/,k-1/*a*/,v->f(j+2,i-1,l+2,k-1)/*v[i-1][j+2][b][a]*/);
                  found = true;
                  en3=w5->f(j,l)/*w5[j][b]*/;
                  i = j;
                  k=l;
                  //a=b;



                }
              }



              if (!found&&k-1>=lowend[i]&&k-1<=highend[i]) {
                //case 2: //gap added
                if (jldecrement) {
                  if (en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+1,i,l+1,k-1)/*v[i][j+1][b][a-1]*/+edangle3(k-1,l+1,k,ct2,data)
                      +penalty(i,j+1,ct1,data)+penalty(k-1,l+1,ct2,data)+gap) {

                    stack.push(j+1,i,l+1/*b*/,k-1/*a-1*/,v->f(j+1,i,l+1,k-1)/*v[i][j+1][b][a-1]*/);
                    found = true;
                    en3=w5->f(j,l)/*w5[j][b]*/;
                    i = j;
                    k=l;
                    //a=b;



                  }
                }

                //case 12	1	0	1	1
                if(!found&&jldecrement2) {

                  if (en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+2,i,l+2,k-1)/*v[i][j+2][b][a-1]*/+edangle5(j+2,i,j+1,ct1,data)+
                      edangle5(l+2,k-1,l+1,ct2,data) + edangle3(k-1,l+2,k,ct2,data)+
                      penalty(i,j+2,ct1,data)+penalty(k-1,l+2,ct2,data)+gap) {

                    stack.push(j+2,i,l+2/*b*/,k-1/*a-1*/,v->f(j+2,i,l+2,k-1)/*v[i][j+2][b][a-1]*/);
                    found = true;
                    en3=w5->f(j,l)/*w5[j][b]*/;
                    i = j;
                    k=l;
                    //a=b;



                  }
                }

                if (!found&&l+2<=highend[j+1]/*highlimit(j+1,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/&&(l+2>=lowend[j+1]/*lowlimit(j+1,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/)) {
                  //case 4
                  if	(en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+1,i,l+2,k-1)/*v[i][j+1][b+1][a-1]*/+edangle3(k-1,l+2,k,ct2,data)+
                       edangle5(l+2,k-1,l+1,ct2,data)+penalty(i,j+1,ct1,data)+penalty(k-1,l+2,ct2,data)+2*gap) {

                    stack.push(j+1,i,l+2/*b+1*/,k-1/*a-1*/,v->f(j+1,i,l+2,k-1)/*v[i][j+1][b+1][a-1]*/);
                    found = true;
                    en3=w5->f(j,l)/*w5[j][b]*/;
                    i = j;
                    k=l;
                    //a=b;



                  }

                }

              }


              if (!found&&l+2<=highend[j+1]/*highlimit(j+1,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/&&l+2>=lowend[j+1]/*lowlimit(j+1,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/) {
                //case 3:
                if(ik) if(en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+1,i,l+2,k)/*v[i][j+1][b+1][a]*/+edangle5(l+2,k,l+1,ct2,data)+
                   penalty(i,j+1,ct1,data)+penalty(l+2,k,ct2,data)+gap) {

                  stack.push(j+1,i,l+2/*b+1*/,k/*a*/,v->f(j+1,i,l+2,k)/*v[i][j+1][b+1][a]*/);
                  found = true;
                  en3=w5->f(j,l)/*w5[j][b]*/;
                  i = j;
                  k=l;
                  //a=b;



                }

                //case 8:		0	1	1	1
                if (!found&&i>1&&ikincrement) {
                  if(en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+1,i-1,l+2,k-1)/*v[i-1][j+1][b+1][a]*/+edangle3(i-1,j+1,i,ct1,data)+
                     edangle5(l+2,k-1,l+1,ct2,data)+edangle3(k-1,l+2,k,ct2,data)+

                     penalty(i-1,j+1,ct1,data)+penalty(k-1,l+2,ct2,data)+gap) {
                    stack.push(j+1,i-1,l+2/*b+1*/,k-1/*a*/,v->f(j+1,i-1,l+2,k-1)/*v[i-1][j+1][b+1][a]*/);
                    found = true;
                    en3=w5->f(j,l)/*w5[j][b]*/;
                    i = j;
                    k=l;
                    //a=b;



                  }
                }

                if (!found&&i>1&&k<=highend[i-1]&&k>=lowend[i-1]) {
                  //case 7		0	1	1	0
                  if(en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+1,i-1,l+2,k)/*v[i-1][j+1][b+1][a+1]*/+edangle5(l+2,k,l+1,ct2,data)+
                     edangle3(i-1,j+1,i,ct1,data)+penalty(i-1,j+1,ct1,data)+penalty(l+2,k,ct2,data)+2*gap) {

                    stack.push(j+1,i-1,l+2/*b+1*/,k/*a+1*/,v->f(j+1,i-1,l+2,k)/*v[i-1][j+1][b+1][a+1]*/);
                    found = true;
                    en3=w5->f(j,l)/*w5[j][b]*/;
                    i = j;
                    k=l;
                    //a=b;



                  }
                }
              }


			if (!found&&i>1&&k<=highend[i-1]&&k>=lowend[i-1]) {
                //case5:
                if (jldecrement) {
                  if(en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+1,i-1,l+1,k)/*v[i-1][j+1][b][a+1]*/+edangle3(i-1,j+1,i,ct1,data)+
                     penalty(i-1,j+1,ct1,data)+penalty(k,l+1,ct2,data)+gap) {

                    stack.push(j+1,i-1,l+1/*b*/,k/*a+1*/,v->f(j+1,i-1,l+1,k)/*v[i-1][j+1][b][a+1]*/);
                    found = true;
                    en3=w5->f(j,l)/*w5[j][b]*/;
                    i = j;
                    k=l;
                    //a=b;

                  }

                }

                //case 15	1	1	1	0
                if (!found&&jldecrement2) {
                  if(en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+2,i-1,l+2,k)/*v[i-1][j+2][b][a+1]*/ + edangle5(j+2,i-1,j+1,ct1,data)+
                     edangle3(i-1,j+2,i,ct1,data)+edangle5(l+2,k,l+1,ct2,data)+
                     penalty(i-1,j+2,ct1,data)+penalty(l+2,k,ct2,data)+gap) {

                    stack.push(j+2,i-1,l+2/*b*/,k/*a+1*/,v->f(j+2,i-1,l+2,k)/*v[i-1][j+2][b][a+1]*/);
                    found = true;
                    en3=w5->f(j,l)/*w5[j][b]*/;
                    i = j;
                    k=l;
                    //a=b;



                  }
                }

                if (!found&&l+1>=lowend[j+2]&&l+1<=highend[j+2]/*lowlimit(j+2,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/) {
                  //case 13	1	1	0	0
                  if(en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+2,i-1,l+1,k)/*v[i-1][j+2][b-1][a+1]*/+edangle5(j+2,i-1,j+1,ct1,data)+
                     edangle3(i-1,j+2,i,ct1,data)+penalty(i-1,j+2,ct1,data)+penalty(k,l+1,ct2,data)+2*gap) {

                    stack.push(j+2,i-1,l+1/*b-1*/,k/*a+1*/,v->f(j+2,i-1,l+1,k)/*v[i-1][j+2][b-1][a+1]*/);
                    found = true;
                    en3=w5->f(j,l)/*w5[j][b]*/;
                    i = j;
                    k=l;
                    //a=b;




                  }

                }

              }
              if (!found&&l+1>=lowend[j+2]&&l+1<=highend[j+2]/*lowlimit(j+2,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/) {
                //case 9		1	0	0	0
                if(ik&&en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+2,i,l+1,k)/*v[i][j+2][b-1][a]*/+edangle5(j+2,i,j+1,ct1,data)+
                   penalty(i,j+2,ct1,data)+penalty(k,l+1,ct2,data)+gap) {

                  stack.push(j+2,i,l+1/*b-1*/,k/*a*/,v->f(j+2,i,l+1,k)/*v[i][j+2][b-1][a]*/);
                  found = true;
                  en3=w5->f(j,l)/*w5[j][b]*/;
                  i = j;
                  k=l;
                  //a=b;




                }

                //case 14	1	1	0	1
                
                  if(!found&&i>1&&ikincrement&&en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+2,i-1,l+1,k-1)/*v[i-1][j+2][b-1][a]*/+edangle5(j+2,i-1,j+1,ct1,data)+
                     edangle3(i-1,j+2,i,ct1,data)+edangle3(k-1,l+1,k,ct2,data)+
                     penalty(i-1,j+2,ct1,data)+penalty(k-1,l+1,ct2,data)+gap) {

                    stack.push(j+2,i-1,l+1/*b-1*/,k-1/*a*/,v->f(j+2,i-1,l+1,k-1)/*v[i-1][j+2][b-1][a]*/);
                    found = true;
                    en3=w5->f(j,l)/*w5[j][b]*/;
                    i = j;
                    k=l;
                    //a=b;



                  }
                

                if (!found&&k-1>=lowend[i]&&k-1<=highend[i]/*lowlimit(i,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/) {
                  //case 10	1	0	0	1
                  if(en3==w5->f(j,l)/*w5[j][b]*/+v->f(j+2,i,l+1,k-1)/*v[i][j+2][b-1][a-1]*/+edangle5(j+2,i,j+1,ct1,data)+
                     edangle3(k-1,l+1,k,ct2,data)+penalty(i,j+2,ct1,data)+penalty(k-1,l+1,ct2,data)+2*gap) {

                    stack.push(j+2,i,l+1/*b-1*/,k-1/*a-1*/,v->f(j+2,i,l+1,k-1)/*v[i][j+2][b-1][a-1]*/);
                    found = true;
                    en3=w5->f(j,l)/*w5[j][b]*/;
                    i = j;
                    k=l;
                    //a=b;



                  }
                }

              }


            }

		  }

	  //report to Dave2
#ifdef DYNALIGN_II

	  for (j=0;j+minloop<i&&!found;++j) {
	    if(en3==w5->f(j,k)+single1_we->f(j+1,i)+gap_branch(i-j)){
//	      cerr<<"w51 i="<<i<<" j="<<j<<" k="<<k<<"\n";
	      trace(ct1,data,j+1,i,single1_v,single1_w,single1_wmb,NULL,NULL,single1_lfce,single1_fce,NULL,NULL,single1_mod,single1_we,single1_we->f(j+1,i),1,0);
	      found=true;
	      en3=w5->f(j,k);
	      i=j;
	    }
	  }
	
	  for (l=0;l+minloop<k&&!found;++l) {
	    if(en3==w5->f(i,l)+single2_we->f(l+1,k)+gap_branch(k-l)){
//	      cerr<<"w52 i="<<i<<" l="<<l<<" k="<<k<<"\n";
	      trace(ct2,data,l+1,k,single2_v,single2_w,single2_wmb,NULL,NULL,single2_lfce,single2_fce,NULL,NULL,single2_mod,single2_we,single2_we->f(l+1,k),1,0);
	      found=true;
	      en3=w5->f(i,l);
	      k=l;
	    }
	  }

#else
#endif

          if (!found) {
				//a traceback error occurred
            //     cerr << "Traceback error at w5!\n";
			  error = 14;
                          //		  return error;
            break;

          }

        }
      }

    }
    else if (closed) {
      //      if(i==54&&j==61&&k==66&&l==74)cerr<<"closed!\n";
      //check if it is safe to increment i and k
      ikincrement =
        k+1 >= lowend[i+1]/*lowlimit(i+1, maxsep, ct1->GetSequenceLength(), ct2->GetSequenceLength())*/ &&
        k+1 <= highend[i+1]/*highlimit(i+1, maxsep, ct1->GetSequenceLength(), ct2->GetSequenceLength())*/;

      ikincrement2 =
        k+2 >= lowend[i+2]/*lowlimit(i+2, maxsep, ct1->GetSequenceLength(), ct2->GetSequenceLength())*/ &&
        k+2 <= highend[i+2]/*highlimit(i+2, maxsep, ct1->GetSequenceLength(), ct2->GetSequenceLength())*/;

      //check to see if it is safe to decrement j and l
      jldecrement =
        l-1 <= highend[j-1]/*highlimit(j-1, maxsep, ct1->GetSequenceLength(), ct2->GetSequenceLength())*/ &&
        l-1 >= lowend[j-1]/*lowlimit(j-1, maxsep, ct1->GetSequenceLength(), ct2->GetSequenceLength())*/;

      jldecrement2 =
        l-2 <= highend[j-2]/*highlimit(j-2, maxsep, ct1->GetSequenceLength(), ct2->GetSequenceLength())*/ &&
        l-2 >= lowend[j-2]/*lowlimit(j-2, maxsep, ct1->GetSequenceLength(), ct2->GetSequenceLength())*/;


      //i and j are paired and aligned to a and b

      if (j<=ct1->GetSequenceLength()) {
       ct1->SetPair(i,j,ct1->GetNumberofStructures());
		
		  ct2->SetPair(k,l,ct2->GetNumberofStructures());
	if(alignment[i]==-1&&alignment[j]==-1){
	  alignment[i]=0;
	  alignment[j]=0;
	}
	else{
	  alignment[i]=k/*a+i-maxsep*/;
	  alignment[j]=l/*b+j-maxsep*/;
	}
      }
      else {
        //j (and l) are > N (and N2)
      
		  ct1->SetPair(i,j-ct1->GetSequenceLength(),ct1->GetNumberofStructures());
		
		  ct2->SetPair(k,l-ct2->GetSequenceLength(),ct2->GetNumberofStructures());
      	if(alignment[i]==-1&&alignment[j-ct1->GetSequenceLength()]==-1){
	  alignment[i]=0;
	  alignment[j-ct1->GetSequenceLength()]=0;
	}
	else{
	  alignment[i]=k/*a+i-maxsep*/;
	  alignment[j-ct1->GetSequenceLength()]=l-ct2->GetSequenceLength()/*b+j-maxsep*/;
	}

      }

      //now find the next pair:

      //      if(i==54&&j==61&&k==66&&l==74)if(en3==erg3(i,j,ct1,data,0)+erg3(k,l,ct2,data,0)+gap*abs(j-i-l+k)){cerr<<"lala!what's wrong?\n";}


      if (en3!=erg3(i,j,ct1,data,0)+erg3(k,l,ct2,data,0)+gap*abs(j-i-l+k)){
        //the fragment does not close hairpins, therefore, the internal
        // fragment does need to be characterized
	
        //first check for internal loop:
        found= false;
#ifdef DYNALIGN_II
        //Now consider internal loops/one stacked and other not
        for (c=i;c<=i+maxloop&&c<j&&!found&&c<=ct1->GetSequenceLength();c++) {
          for (d=j;d>=j-maxloop&&d>c;d--) {
            //if (d>ct1->GetSequenceLength()) ap = ct1->GetSequenceLength()-ct2->GetSequenceLength();
            //else ap =0;
            for (e=max(k,lowend[c]/*lowlimit(c,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/);
                 e<=k+maxloop&&e<l&&e<=highend[c]/*highlimit(c,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/&&e<=ct2->GetSequenceLength()&&!found;e++) {
	      
              for (f=min(l,highend[d]/*highlimit(d,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/);f>=l-maxloop&&f>e&&f>=lowend[d]/*lowlimit(d,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/
                     /*d-maxsep-ap*/&&!found;f--) {
	
                if (c-i+j-d-2<maxinternal&&e-k+l-f-2<maxinternal&&(((d>ct1->GetSequenceLength())&&(f>ct2->GetSequenceLength()))||((d<=ct1->GetSequenceLength())&&(f<=ct2->GetSequenceLength())))) {//newest//report to Dave 16 change of parantheses 
	
		  helicalstack=false;
		  stuck=false;

		  if(c==i&&d==j&&(e!=k&&f!=l)){
		
		    en1=0;
		    if (modification) stuck=true;
		  }
		  
		  //		  	  /*we should remove this part? report to Dave14
		  else if (c==i+1&&d==j-1&&(i!=ct1->GetSequenceLength())&&(j-1!=ct1->GetSequenceLength())
		      //&&report to Dave 17 unnecessary?		      ((d>ct1->GetSequenceLength()&&f>ct2->GetSequenceLength())||(d<=ct1->GetSequenceLength()&&f<=ct2->GetSequenceLength()))
		      ) {
                    //seq1 is helical stack
                    en1 = erg1(i,j,i+1,j-1,ct1, data);
                    if (modification) helicalstack=true;
		    
		    
                  }
                  else  if(c!=i&&d!=j){
                    //internal loop

                    en1 = erg2(i,j,c,d,ct1,data,0,0);

                  }
		  else continue;

		  if(e==k&&f==l&&(c!=i&&d!=j)){
		    en2=0;
		    if (modification) stuck=true;
		  }
		  
                  else if (e==k+1&&f==l-1&&k!=ct2->GetSequenceLength()&&(l-1!=ct2->GetSequenceLength())) {
                    //seq2 is helical stack
                    en2 = erg1(k,l,k+1,l-1,ct2, data);
                    if (modification) helicalstack=true;
                  }
                  else if(e!=k&&f!=l){
                    //internal loop
                    en2 = erg2(k,l,e,f,ct2,data,0,0);

                  }
		  else continue;

                  if (helicalstack) {
                    if(en3== en1+en2+vmod->f(c,d,e,f)+gap*(abs(c-i+k-e)+abs(j-d+f-l))){
		      found = true;
                      stack.push(c,d,e,f,vmod->f(c,d,e,f),false,true);
		      if(e==k&&f==l&&(c!=i&&d!=j)){
			if(j<=ct1->GetSequenceLength()){alignment[c]=-1;alignment[d]=-1;}
			else{alignment[c]=-1;alignment[d-ct1->GetSequenceLength()]=-1;}
		      }
		    }
		 
		    

                  }
		  else if(stuck&&if_vmod){
		     if(en3== en1+en2+vmod->f(c,d,e,f)+gap*(abs(c-i+k-e)+abs(j-d+f-l))){
                      found = true;
                      stack.push(c,d,e,f,vmod->f(c,d,e,f),false,true);
		      if(e==k&&f==l&&(c!=i&&d!=j)){
			if(j<=ct1->GetSequenceLength()){alignment[c]=-1;alignment[d]=-1;}
			else{alignment[c]=-1;alignment[d-ct1->GetSequenceLength()]=-1;}
		      }
		    }
		  
		  
		  }
                  else if((c!=i&&d!=j)||(e!=k&&f!=l)){
                    if(en3== en1+en2+v->f(c,d,e,f)+gap*(abs(c-i+k-e)+abs(j-d+f-l))){
		      //   cerr<<c<<" "<<d<<" "<<e<<" "<<f<<" closed!\n";
                      found = true;
                      stack.push(c,d,e,f,v->f(c,d,e,f));
		      if(e==k&&f==l&&(c!=i&&d!=j)){
			if(j<=ct1->GetSequenceLength()){alignment[c]=-1;alignment[d]=-1;}
			else{alignment[c]=-1;alignment[d-ct1->GetSequenceLength()]=-1;}
		      }
		    }
		  
		    /*    report to Dave14 remove this part*/
                }


              }
            }
          }
        }
	}
#else
   //Now consider internal loops/one stacked and other not
        for (c=i+1;c<=i+maxloop&&c<j&&!found&&c<=ct1->GetSequenceLength();c++) {
          for (d=j-1;d>=j-maxloop&&d>c;d--) {
            //if (d>ct1->GetSequenceLength()) ap = ct1->GetSequenceLength()-ct2->GetSequenceLength();
            //else ap =0;
            for (e=max(k+1,lowend[c]/*lowlimit(c,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/);
                 e<=k+maxloop&&e<l&&e<=highend[c]/*highlimit(c,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/&&e<=ct2->GetSequenceLength()&&!found;e++) {

              for (f=min(l-1,highend[d]/*highlimit(d,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/);f>=l-maxloop&&f>e&&f>=lowend[d]/*lowlimit(d,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/
                     /*d-maxsep-ap*/&&!found;f--) {

                if (c-i+j-d-2<maxinternal&&e-k+l-f-2<maxinternal&&((d>ct1->GetSequenceLength())&&(f>ct2->GetSequenceLength()))||((d<=ct1->GetSequenceLength())&&(f<=ct2->GetSequenceLength()))) {
                  helicalstack=false;
                  if (c==i+1&&d==j-1&&(i!=ct1->GetSequenceLength())&&(j-1!=ct1->GetSequenceLength())&&
                      ((d>ct1->GetSequenceLength()&&f>ct2->GetSequenceLength())||(d<=ct1->GetSequenceLength()&&f<=ct2->GetSequenceLength()))) {
                    //seq1 is helical stack
                    en1 = erg1(i,j,i+1,j-1,ct1, data);
                    if (modification) helicalstack=true;


                  }
                  else {
                    //internal loop

                    en1 = erg2(i,j,c,d,ct1,data,0,0);

                  }

                  if (e==k+1&&f==l-1&&k!=ct2->GetSequenceLength()&&(l-1!=ct2->GetSequenceLength())) {
                    //seq2 is helical stack
                    en2 = erg1(k,l,k+1,l-1,ct2, data);
                    if (modification) helicalstack=true;
                  }
                  else {
                    //internal loop
                    en2 = erg2(k,l,e,f,ct2,data,0,0);

                  }

                  if (helicalstack) {
                    if(en3== en1+en2+vmod->f(c,d,e,f)/*vmod[d][iref(c,d,ct1->GetSequenceLength())][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l))){
                      found = true;
                      stack.push(c,d,e/*-c+maxsep*/,f/*-d+maxsep+ap*/,vmod->f(c,d,e,f)/*vmod[d][iref(c,d,ct1->GetSequenceLength())][e-c+maxsep][ f-d+maxsep+ap]*/);

                    }


                    else if (c==i+2&&d==j-2&&j-2>i+2&&
                             inc[ct1->numseq[i+1]][ct1->numseq[j-1]]&&
                             inc[ct1->numseq[i+2]][ct1->numseq[j-2]]&&
                             e==k+1&&f==l-1&&
                             inc[ct2->numseq[k+1]][ct2->numseq[l-1]]
                             &&i!=ct1->GetSequenceLength()&&i+1!=ct1->GetSequenceLength()&&
                             (j-1!=ct1->GetSequenceLength())&&(j-2!=ct1->GetSequenceLength())
                             &&k!=ct2->GetSequenceLength()&&(l-1!=ct2->GetSequenceLength())) {

                      en1 = erg1(i,j,i+1,j-1,ct1,data)+
                        erg1(i+1,j-1,i+2,j-2,ct1,data);

                      if (en3==en1+en2+vmod->f(c,d,e,f)/*vmod[d][iref(c,d,ct1->GetSequenceLength())][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l))) {
                        //base pairs with single base pair insertion in ct1
                        //ct1->basepr[ct1->GetNumberofStructures()][i+1]=j-1;
                        //ct1->basepr[ct1->GetNumberofStructures()][j-1]=i+1;
                        if (j<=ct1->GetSequenceLength()) {
                        	ct1->SetPair(i+1,j-1,ct1->GetNumberofStructures());

                        }
                        else {
                          //j (and l) are > N (and N2)
                       	ct1->SetPair(i+1,j-1-ct1->GetSequenceLength(),ct1->GetNumberofStructures());



                        }
                        stack.push(c,d,e/*-c+maxsep*/,f/*-d+maxsep+ap*/,v->f(c,d,e,f)/*v[d][iref(c,d,ct1->GetSequenceLength())][e-c+maxsep][f-d+maxsep+ap]*/);
                        found = true;


                      }




                    }
                    else if (e==k+2&&f==l-2&&l-2>k+2&&c==i+1&&d==j-1&&
                             inc[ct2->numseq[k+1]][ct2->numseq[l-1]]&&
                             inc[ct2->numseq[k+2]][ct2->numseq[l-2]]&&
                             inc[ct1->numseq[i+1]][ct1->numseq[j-1]]&&
                             i!=ct1->GetSequenceLength()&&(j-1!=ct1->GetSequenceLength())&&k!=ct2->GetSequenceLength()&&(l-1!=ct2->GetSequenceLength())&&
                             (k+1!=ct2->GetSequenceLength())&&(l-2!=ct2->GetSequenceLength())) {

                      en2= erg1(k,l,k+1,l-1,ct2,data)+
                        erg1(k+1,l-1,k+2,l-2,ct2,data);

                      if (en3==en1+en2+vmod->f(c,d,e,f)/*vmod[d][iref(c,d,ct1->GetSequenceLength())][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l))) {
                        //base pairs with single base pair insertion in ct1
                        //ct2->basepr[ct2->GetNumberofStructures()][k+1]=l-1;
                        //ct2->basepr[ct2->GetNumberofStructures()][l-1]=k+1;
                        if (l<=ct2->GetSequenceLength()) {
                            ct2->SetPair(k+1,l-1,ct2->GetNumberofStructures());

                         
                        }
                        else {
                          //j (and l) are > N (and N2)

                        	ct2->SetPair(k+1,l-1-ct2->GetSequenceLength(),ct2->GetNumberofStructures());


                        }
                        stack.push(c,d,e/*-c+maxsep*/,f/*-d+maxsep+ap*/,v->f(c,d,e,f)/*v[d][iref(c,d,ct1->GetSequenceLength())][e-c+maxsep][f-d+maxsep+ap]*/);

                        found = true;

                      }


                    }


                  }
                  else {
                    if(en3== en1+en2+v->f(c,d,e,f)/*v[d][iref(c,d,ct1->GetSequenceLength())][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l))){
                      found = true;
                      stack.push(c,d,e/*-c+maxsep*/,f/*-d+maxsep+ap*/,v->f(c,d,e,f)/*v[d][iref(c,d,ct1->GetSequenceLength())][e-c+maxsep][ f-d+maxsep+ap]*/);



                    }

                    else if (c==i+2&&d==j-2&&j-2>i+2&&
                             inc[ct1->numseq[i+1]][ct1->numseq[j-1]]&&
                             inc[ct1->numseq[i+2]][ct1->numseq[j-2]]&&
                             e==k+1&&f==l-1&&
                             inc[ct2->numseq[k+1]][ct2->numseq[l-1]]
                             &&i!=ct1->GetSequenceLength()&&i+1!=ct1->GetSequenceLength()&&
                             (j-1!=ct1->GetSequenceLength())&&(j-2!=ct1->GetSequenceLength())
                             &&k!=ct2->GetSequenceLength()&&(l-1!=ct2->GetSequenceLength())) {

                      en1 = erg1(i,j,i+1,j-1,ct1,data)+
                        erg1(i+1,j-1,i+2,j-2,ct1,data);

                      if (en3==en1+en2+v->f(c,d,e,f)/*v[d][iref(c,d,ct1->GetSequenceLength())][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l))) {
                        //base pairs with single base pair insertion in ct1
                        //ct1->basepr[ct1->GetNumberofStructures()][i+1]=j-1;
                        //ct1->basepr[ct1->GetNumberofStructures()][j-1]=i+1;
                        if (j<=ct1->GetSequenceLength()) {
                     	ct1->SetPair(i+1,j-1,ct1->GetNumberofStructures());
                          

                        }
                        else {
                          //j (and l) are > N (and N2)
                        	ct1->SetPair(i+1,j-1-ct1->GetSequenceLength(),ct1->GetNumberofStructures());



                        }
                        stack.push(c,d,e/*-c+maxsep*/,f/*-d+maxsep+ap*/,v->f(c,d,e,f)/*v[d][iref(c,d,ct1->GetSequenceLength())][e-c+maxsep][f-d+maxsep+ap]*/);
                        found = true;


                      }




                    }
                    else if (e==k+2&&f==l-2&&l-2>k+2&&c==i+1&&d==j-1&&
                             inc[ct2->numseq[k+1]][ct2->numseq[l-1]]&&
                             inc[ct2->numseq[k+2]][ct2->numseq[l-2]]&&
                             inc[ct1->numseq[i+1]][ct1->numseq[j-1]]&&
                             i!=ct1->GetSequenceLength()&&(j-1!=ct1->GetSequenceLength())&&k!=ct2->GetSequenceLength()&&(l-1!=ct2->GetSequenceLength())&&
                             (k+1!=ct2->GetSequenceLength())&&(l-2!=ct2->GetSequenceLength())) {

                      en2= erg1(k,l,k+1,l-1,ct2,data)+
                        erg1(k+1,l-1,k+2,l-2,ct2,data);

                      if (en3==en1+en2+v->f(c,d,e,f)/*v[d][iref(c,d,ct1->GetSequenceLength())][e-c+maxsep][f-d+maxsep+ap]*/+gap*(abs(c-i+k-e)+abs(j-d+f-l))) {
                        //base pairs with single base pair insertion in ct1
                        //ct2->basepr[ct2->GetNumberofStructures()][k+1]=l-1;
                        //ct2->basepr[ct2->GetNumberofStructures()][l-1]=k+1;
                        if (l<=ct2->GetSequenceLength()) {

                  	ct2->SetPair(k+1,l-1,ct2->GetNumberofStructures());

                        }
                        else {
                          //j (and l) are > N (and N2)

                        	ct2->SetPair(k+1,l-1-ct2->GetSequenceLength(),ct2->GetNumberofStructures());


                        }
                        stack.push(c,d,e/*-c+maxsep*/,f/*-d+maxsep+ap*/,v->f(c,d,e,f)/*v[d][iref(c,d,ct1->GetSequenceLength())][e-c+maxsep][f-d+maxsep+ap]*/);

                        found = true;

                      }


                    }
                  }
                }

              }
            }
          }
        }
#endif
#ifdef DYNALIGN_II

	int stacking_elongation_1=0;
	int stacking_elongation_2=0;
	int alignment_elongation=0;
	
	for(c=1;!found&&c<=max_elongation;c++){
	  
	  if((c+1>maxloop)||(i+c+1>=j-c-1)||(i+c+1>ct1->GetSequenceLength())||(k+c+1>=l-c-1)||(k+c+1>ct2->GetSequenceLength()))break;
	  if(j>ct1->GetSequenceLength())if((j-c-1<=ct1->GetSequenceLength()+1)||(l-c-1<=ct2->GetSequenceLength()+1))break;
	  
	  if(k+c<lowend[i+c]||k+c>highend[i+c]||l-c<lowend[j-c]||l-c>highend[j-c])break;
	  if(k+c+1<lowend[i+c+1]||k+c+1>highend[i+c+1]||l-c-1<lowend[j-c-1]||l-c-1>highend[j-c-1])break;

	  if (alignmentforced) {
	    if ((forcealign[0][i+c]>0&&forcealign[0][i+c]!=k+c)||
		(forcealign[0][j-c]>0&&forcealign[0][j-c]!=l-c)||
		(forcealign[1][k+c]>0&&forcealign[1][k+c]!=i+c)||
		(forcealign[1][l-c]>0&&forcealign[1][l-c]!=j-c)) {
	      break;
	    }
	    
	    else {
        //scan through to make sure that the base pair isn't precluding an alignment
	      for (g=i+c+1;g<j-c;g++) {
		if ((forcealign[0][g]>0)&&(forcealign[0][g]<k+c||forcealign[0][g]>l-c)) {
		  break;
		}
	      }
	      
	      for (g=k+c+1;g<l-c;g++) {
		if ((forcealign[1][g]>0)&&(forcealign[1][g]<i+c||forcealign[1][g]>j-c)) {
		  break;
		}
	      }
	    }

	    if ((forcealign[0][i+c+1]>0&&forcealign[0][i+c+1]!=k+c+1)||
		(forcealign[0][j-c-1]>0&&forcealign[0][j-c-1]!=l-c-1)||
		(forcealign[1][k+c+1]>0&&forcealign[1][k+c+1]!=i+c+1)||
		(forcealign[1][l-c-1]>0&&forcealign[1][l-c-1]!=j-c-1)) {
	      break;
	    }
	    
	    else {
        //scan through to make sure that the base pair isn't precluding an alignment
	      for (g=i+c+2;g<j-c-1;g++) {
		if ((forcealign[0][g]>0)&&(forcealign[0][g]<k+c+1||forcealign[0][g]>l-c-1)) {
		  break;
		}
	      }
	      
	      for (g=k+c+2;g<l-c-1;g++) {
		if ((forcealign[1][g]>0)&&(forcealign[1][g]<i+c+1||forcealign[1][g]>j-c-1)) {
		  break;
		}
	      }
	    }

	  }



	  alignment_elongation=c;
	  
	}

       	for(c=1;!found&&c<=alignment_elongation;c++){
	  if(!inc[ct1->numseq[i+c]][ct1->numseq[j-c]])break;
	  if (force)if(fce1[jref(i+c,j-c,ct1->GetSequenceLength())][iref(i+c,j-c,ct1->GetSequenceLength())]&SINGLE||fce1[jref(i+c,j+c,ct1->GetSequenceLength())][iref(i+c,j+c,ct1->GetSequenceLength())]&NOPAIR)break;
	  
	  if (force && modification) {
	    if (mod1[i+c]||mod1[j-c]) {
	      //check for a GU exception
	      if ((ct1->numseq[i+c]==3&&ct1->numseq[j-c]==4)||
		  (ct1->numseq[i+c]==4&&ct1->numseq[j-c]==3)||
		  (ct1->numseq[i+c+1]==3&&ct1->numseq[j-c-1]==4)||
		  (ct1->numseq[i+c+1]==4&&ct1->numseq[j-c-1]==3)||
		  (ct1->numseq[i+c-1]==3&&ct1->numseq[j-c+1]==4)||
		  (ct1->numseq[i+c-1]==4&&ct1->numseq[j-c+1]==3));
	      else break;
	    }
	  }
	  
	  if (ct1->templated) {
	    if (j<=ct1->GetSequenceLength()) {
	      d=j-c;
	      e=i+c;
	    }
	    else {
	      e=j-c-ct1->GetSequenceLength();
	      d=i+c;
	    }
	    
	    if (!ct1->tem[d][e]) {
	      break;
	    }
	  }
	  stacking_elongation_1=c;	  
	  
	  

	}

	    
       	for(c=1;!found&&c<=alignment_elongation;c++){
	  if(!inc[ct2->numseq[k+c]][ct2->numseq[l-c]])break;
	  if (force)if(fce2[jref(k+c,l-c,ct2->GetSequenceLength())][iref(k+c,l-c,ct2->GetSequenceLength())]&SINGLE||fce2[jref(k+c,l-c,ct2->GetSequenceLength())][iref(k+c,l-c,ct2->GetSequenceLength())]&NOPAIR)break;
	  
	  if (force && modification) {
	    if (mod2[k+c]||mod2[l-c]) {
	      //check for a GU exception
	      if ((ct2->numseq[k+c]==3&&ct2->numseq[l-c]==4)||
		  (ct2->numseq[k+c]==4&&ct2->numseq[l-c]==3)||
		  (ct2->numseq[k+c+1]==3&&ct2->numseq[l-c-1]==4)||
		  (ct2->numseq[k+c+1]==4&&ct2->numseq[l-c-1]==3)||
		  (ct2->numseq[k+c-1]==3&&ct2->numseq[l-c+1]==4)||
		  (ct2->numseq[k+c-1]==4&&ct2->numseq[l-c+1]==3));
	      else break;
	    }
	  }

	  if (ct2->templated) {
	    if (j<=ct1->GetSequenceLength()) {
	      d=l-c;
	      e=k+c;
	    }
	    else {
	      e=l-c-ct2->GetSequenceLength();
	      d=k+c;
	    }
	    
	    if (!ct2->tem[d][e]) {
	      break;
	    }
	  }
	  stacking_elongation_2=c;	  
	  
	}	  

	en1_temp=0;
	for(c=1;!found&&c<=stacking_elongation_1;c++){
  if((!inc[ct1->numseq[i+c+1]][ct1->numseq[j-c-1]])||(!inc[ct2->numseq[k+c+1]][ct2->numseq[l-c-1]]))continue;
	  en1_temp+=erg1(i+c-1,j-c+1,i+c,j-c,ct1,data);
	  en1=en1_temp+erg1(i+c,j-c,i+c+1,j-c-1,ct1,data);
	  en2=erg2(k,l,k+c+1,l-c-1,ct2,data,0,0);
	  
	  if(modification){
	    if(en3==en1+en2+vmod->f(i+c+1,j-c-1,k+c+1,l-c-1)){
	      found=true;
	      for(d=1;d<=c;d++){
		alignment[i+d]=k+d;
		if(j<=ct1->GetSequenceLength()){
		  alignment[j-d]=l-d;
//		  ct1->basepr[ct1->GetNumberofStructures()][i+d]=j-d;
//		  ct1->basepr[ct1->GetNumberofStructures()][j-d]=i+d;
                  ct1->SetPair(i+d, j-d, ct1->GetNumberofStructures());
		}
		else{
		  alignment[j-d-ct1->GetSequenceLength()]=l-d-ct2->GetSequenceLength();
//		  ct1->basepr[ct1->GetNumberofStructures()][i+d]=j-d-ct1->GetSequenceLength();
//		  ct1->basepr[ct1->GetNumberofStructures()][j-d-ct1->GetSequenceLength()]=i+d;
                  ct1->SetPair(i+d, j-d-ct1->GetSequenceLength(), ct1->GetNumberofStructures());
		}
	      }
	      stack.push(i+c+1,j-c-1,k+c+1,l-c-1,vmod->f(i+c+1,j-c-1,k+c+1,l-c-1),false,true);
	    }
	  }
	  else {
	    if(en3==en1+en2+v->f(i+c+1,j-c-1,k+c+1,l-c-1)){
	      found=true;
	      for(d=1;d<=c;d++){
		alignment[i+d]=k+d;
		if(j<=ct1->GetSequenceLength()){
		  alignment[j-d]=l-d;
//		  ct1->basepr[ct1->GetNumberofStructures()][i+d]=j-d;
//		  ct1->basepr[ct1->GetNumberofStructures()][j-d]=i+d;
                  ct1->SetPair(i+d, j-d, ct1->GetNumberofStructures());
		}
		else{
		  alignment[j-d-ct1->GetSequenceLength()]=l-d-ct2->GetSequenceLength();
//		  ct1->basepr[ct1->GetNumberofStructures()][i+d]=j-d-ct1->GetSequenceLength();
//		  ct1->basepr[ct1->GetNumberofStructures()][j-d-ct1->GetSequenceLength()]=i+d;
                  ct1->SetPair(i+d, j-d-ct1->GetSequenceLength(), ct1->GetNumberofStructures());
		}
	      }
	      stack.push(i+c+1,j-c-1,k+c+1,l-c-1,v->f(i+c+1,j-c-1,k+c+1,l-c-1));
	    }
	  }
	 
	}
	      
	en1_temp=0;
	for(c=1;!found&&c<=stacking_elongation_2;c++){
	  if((!inc[ct1->numseq[i+c+1]][ct1->numseq[j-c-1]])||(!inc[ct2->numseq[k+c+1]][ct2->numseq[l-c-1]]))continue;
	  
	  en1_temp+=erg1(k+c-1,l-c+1,k+c,l-c,ct2,data);
	  en2=en1_temp+erg1(k+c,l-c,k+c+1,l-c-1,ct2,data);
	  en1=erg2(i,j,i+c+1,j-c-1,ct1,data,0,0);
	  
	  if(modification){
	    if(en3==en1+en2+vmod->f(i+c+1,j-c-1,k+c+1,l-c-1)){
	      found=true;
	      for(d=1;d<=c;d++){
		alignment[i+d]=k+d;
		if(j<=ct1->GetSequenceLength()){
		  alignment[j-d]=l-d;
//		  ct2->basepr[ct2->GetNumberofStructures()][k+d]=l-d;
//		  ct2->basepr[ct2->GetNumberofStructures()][l-d]=k+d;
                  ct2->SetPair(k+d, l-d, ct2->GetNumberofStructures());
		}
		else{
		  alignment[j-d-ct1->GetSequenceLength()]=l-d-ct2->GetSequenceLength();
//		  ct2->basepr[ct2->GetNumberofStructures()][k+d]=l-d-ct2->GetSequenceLength();
//		  ct2->basepr[ct2->GetNumberofStructures()][l-d-ct2->GetSequenceLength()]=k+d;
                  ct2->SetPair(k+d, l-d-ct2->GetSequenceLength(), ct2->GetNumberofStructures());
		}
	      }
	      stack.push(i+c+1,j-c-1,k+c+1,l-c-1,vmod->f(i+c+1,j-c-1,k+c+1,l-c-1),false,true);
	    }
	  }
	  else {
	    if(en3==en1+en2+v->f(i+c+1,j-c-1,k+c+1,l-c-1)){
	      
	      found=true;
	      for(d=1;d<=c;d++){
		alignment[i+d]=k+d;
		if(j<=ct1->GetSequenceLength()){
		  alignment[j-d]=l-d;
//		  ct2->basepr[ct2->GetNumberofStructures()][k+d]=l-d;
//		  ct2->basepr[ct2->GetNumberofStructures()][l-d]=k+d;
                  ct2->SetPair(k+d, l-d, ct2->GetNumberofStructures());
		}
		else{
		  alignment[j-d-ct1->GetSequenceLength()]=l-d-ct2->GetSequenceLength();
//		  ct2->basepr[ct2->GetNumberofStructures()][k+d]=l-d-ct2->GetSequenceLength();
//		  ct2->basepr[ct2->GetNumberofStructures()][l-d-ct2->GetSequenceLength()]=k+d;
                  ct2->SetPair(k+d, l-d-ct2->GetSequenceLength(), ct2->GetNumberofStructures());
		}
	      }
	      stack.push(i+c+1,j-c-1,k+c+1,l-c-1,v->f(i+c+1,j-c-1,k+c+1,l-c-1));
	    }
	  }
	 
	}
#else
#endif
        constantclosure = penalty(i,j,ct1,data)+penalty(k,l,ct2,data);
        if (j>ct1->GetSequenceLength()&&!found) {
          //consider an exterior loop closed by i-j and k-l
          //There are 16 cases:
          //consider them in this order (0 unstacked, 1 stacked)
          //		i	j	k	l
          //1		0	0	0	0
          //2		0	0	0	1
          //3		0	0	1	0
          //4		0	0	1	1
          //5		0	1	0	0
          //6		0	1	0	1
          //7		0	1	1	0
          //8		0	1	1	1
          //9		1	0	0	0
          //10		1	0	0	1
          //11		1	0	1	0
          //12		1	0	1	1
          //13		1	1	0	0
          //14		1	1	0	1
          //15		1	1	1	0
          //16		1	1	1	1

		  //allow an exception for jldecrement2 if this reaches the end of the sequence
		  if (j-2==ct1->GetSequenceLength()) jldecrement2=true;

          //case 1, no stacks
          //if (jldecrement) {
			  //if (ikincrement) {
		  //  if (k<ct2->GetSequenceLength()) {
				  if (en3==w5->f(j-1-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b]*/+
						w3->f(i+1,k+1)/*w3[i+1][a]*/+constantclosure) {
						stack.push(j-1-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength()/*b*/,1,0,w5->f(j-1-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b]*/,true);
						stack.push(i+1,k+1/*a*/,ct1->GetSequenceLength(),0,w3->f(i+1,k+1)/*w3[i+1][a]*/,true);
						found = true;
					}

				  //	  }
				  //	  else if (k==ct2->GetSequenceLength()&&local) {
				  //		if (en3==w5->f(j-1-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b]*/+
				  //			+constantclosure) {
				  //			stack.push(j-1-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength()/*b*/,1,0,w5->f(j-1-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b]*/,true);

				  //			found = true;
				  //		}

				  //  }
				  //		  else if (k==ct2->GetSequenceLength()) {
				  //		  if (en3==w5->f(j-1-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b]*/+
				  //			+constantclosure+gap*(ct1->GetSequenceLength()-i)) {
				  //				stack.push(j-1-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength()/*b*/,1,0,w5->f(j-1-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b]*/,true);

				  //			found = true;
				  //		}
	//
				  //}

				  //

          //case 16, all stack
          if ((i<ct1->GetSequenceLength())&&(j>ct1->GetSequenceLength()+1)&&!found&&k+2<=ct2->GetSequenceLength()+1&&j-2-ct1->GetSequenceLength()>=0&&l-2-ct2->GetSequenceLength()>=0) {
            if (en3==w5->f(j-2-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b]*/+w3->f(i+2,k+2)/*w3[i+2][a]*/+edangle5(j,i,j-1,ct1,data)+edangle5(l,k,l-1,ct2,data)
                +edangle3(i,j,i+1,ct1,data)+edangle3(k,l,k+1,ct2,data)+constantclosure) {

              stack.push(j-2-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength()/*b*/,1,0,w5->f(j-2-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b]*/,true);
              stack.push(i+2,k+2/*a*/,ct1->GetSequenceLength(),0,w3->f(i+2,k+2)/*w3[i+2][a]*/,true);
              found = true;

            }

          }
          //case 6, j and l stacked
          if (j>ct1->GetSequenceLength()+1&&!found&&j-2-ct1->GetSequenceLength()>=0&&l-2-ct2->GetSequenceLength()>=0) {
	    //  if (ikincrement) {
				  if (en3==w5->f(j-2-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b]*/+w3->f(i+1,k+1)/*w3[i+1][a]*/+edangle5(j,i,j-1,ct1,data)+edangle5(l,k,l-1,ct2,data)+constantclosure) {

					stack.push(j-2-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength()/*b*/,1,0,w5->f(j-2-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b]*/,true);
					stack.push(i+1,k+1/*a*/,ct1->GetSequenceLength(),0,w3->f(i+1,k+1)/*w3[i+1][a]*/,true);
					found = true;
				  }

				  //    }
				  //  else if (k==ct2->GetSequenceLength()&&local) {
				  //	  if (en3==w5->f(j-2-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b]*/+edangle5(j,i,j-1,ct1,data)+edangle5(l,k,l-1,ct2,data)+constantclosure) {
				  
				  //		stack.push(j-2-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength()/*b*/,1,0,w5->f(j-2-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b]*/,true);

				  //		found = true;
				  //	  }

				  //  }
				  //	  else if (k==ct2->GetSequenceLength()) {
				  //	  if (en3==w5->f(j-2-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b]*/+edangle5(j,i,j-1,ct1,data)+edangle5(l,k,l-1,ct2,data)+constantclosure+gap*(ct1->GetSequenceLength()-i)) {
				  //
				  //		stack.push(j-2-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength()/*b*/,1,0,w5->f(j-2-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b]*/,true);
				  //
				  //		found = true;
				  //	  }

				  //  }
				    	  }

          //case 11, i and k stacked
          if ((i<ct1->GetSequenceLength())&&!found&&k+2<=ct2->GetSequenceLength()+1) {
            if (en3==w5->f(j-1-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b]*/
                +w3->f(i+2,k+2)/*w3[i+2][a]*/+edangle3(i,j,i+1,ct1,data)+edangle3(k,l,k+1,ct2,data)+constantclosure) {

              stack.push(j-1-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength()/*b*/,1,0,w5->f(j-1-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b]*/,true);
              stack.push(i+2,k+2/*a*/,ct1->GetSequenceLength(),0,w3->f(i+2,k+2)/*w3[i+2][a]*/,true);
              found = true;


            }
          }

		  if (l-2-ct2->GetSequenceLength()>=0) if (!found) {
            //case 2, l stacked
		      //	  if (k<ct2->GetSequenceLength()) {
				  if (en3==w5->f(j-1-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b-1]*/+
					w3->f(i+1,k+1)/*w3[i+1][a]*/+edangle5(l,k,l-1,ct2,data)+gap+constantclosure) {
					stack.push(j-1-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength()/*b-1*/,1,0,w5->f(j-1-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b-1]*/,true);
					stack.push(i+1,k+1/*a*/,ct1->GetSequenceLength(),0,w3->f(i+1,k+1)/*w3[i+1][a]*/,true);
					found = true;

				}
				  //}
				  //else if (k==ct2->GetSequenceLength()&&local) {
				  //if (en3==w5->f(j-1-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b-1]*/+
				  //	+edangle5(l,k,l-1,ct2,data)+gap+constantclosure) {
				  //	stack.push(j-1-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength()/*b-1*/,1,0,w5->f(j-1-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b-1]*/,true);

				  //	found = true;
				  //			  }

				  //		  }
				  //		  else if (k==ct2->GetSequenceLength()){
				  //			if (en3==w5->f(j-1-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b-1]*/+
				  //				+edangle5(l,k,l-1,ct2,data)+gap+constantclosure+gap*(ct1->GetSequenceLength()-i)) {
				  //			stack.push(j-1-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength()/*b-1*/,1,0,w5->f(j-1-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b-1]*/,true);
				  //
				  //			found = true;
					//		}

				  //	  }

            //case 12, i, k, and l stacked
            if ((i<ct1->GetSequenceLength())&&!found&&k+2<=ct2->GetSequenceLength()+1) {

              if (en3==w5->f(j-1-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b-1]*/+
                  w3->f(i+2,k+2)/*w3[i+2][a]*/+edangle5(l,k,l-1,ct2,data)+edangle3(i,j,i+1,ct1,data)
                  +edangle3(k,l,k+1,ct2,data)+gap+constantclosure) {

                stack.push(j-1-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength()/*b-1*/,1,0,w5->f(j-1-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b-1]*/,true);
                stack.push(i+2,k+2/*a*/,ct1->GetSequenceLength(),0,w3->f(i+2,k+2)/*w3[i+2][a]*/,true);
                found = true;

              }
            }
            if (k+2<=ct2->GetSequenceLength()+1/*highlimit(i+1,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/&&!found) {
              //case 4, l and k stacked
              if (en3==w5->f(j-1-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b-1]*/+
                  w3->f(i+1,k+2)/*w3[i+1][a+1]*/+edangle5(l,k,l-1,ct2,data)+
                  edangle3(k,l,k+1,ct2,data)+2*gap+constantclosure){

                stack.push(j-1-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength()/*b-1*/,1,0,w5->f(j-1-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b-1]*/,true);
                stack.push(i+1,k+2/*a+1*/,ct1->GetSequenceLength(),0,w3->f(i+1,k+2)/*w3[i+1][a+1]*/,true);
                found = true;

              }
            }

            if ((i<ct1->GetSequenceLength())/*lowlimit(i+2,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())/*a>0*/&&!found) {
              //case 10, i and l stacked
	      //		if (k<ct2->GetSequenceLength()) {
					if (en3==w5->f(j-1-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b-1]*/+
						w3->f(i+2,k+1)/*w3[i+2][a-1]*/+edangle5(l,k,l-1,ct2,data)+
						edangle3(i,j,i+1,ct1,data)+2*gap+constantclosure){

					  stack.push(j-1-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength()/*b-1*/,1,0,w5->f(j-1-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b-1]*/,true);
					  stack.push(i+2,k+1/*a-1*/,ct1->GetSequenceLength(),0,w3->f(i+2,k+1)/*w3[i+2][a-1]*/,true);
					  found = true;


					}
					//	}
					//		else if (k==ct2->GetSequenceLength()&&local) {
					//			if (en3==w5->f(j-1-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b-1]*/+
					//			+edangle5(l,k,l-1,ct2,data)+
					//			edangle3(i,j,i+1,ct1,data)+2*gap+constantclosure){

					//		  stack.push(j-1-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength()/*b-1*/,1,0,w5->f(j-1-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b-1]*/,true);

					//		  found = true;


					//		}
					//	}
					//	else if (k==ct2->GetSequenceLength()) {
					//		if (en3==w5->f(j-1-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b-1]*/+
					//			+edangle5(l,k,l-1,ct2,data)+
					//			edangle3(i,j,i+1,ct1,data)+2*gap+constantclosure+gap*(ct1->GetSequenceLength()-i-1)){

					//		  stack.push(j-1-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength()/*b-1*/,1,0,w5->f(j-1-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b-1]*/,true);

					//		  found = true;


					//		}
					//	}

            }


          }//

          if (k+2<=ct2->GetSequenceLength()+1/*highlimit(i+1,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())/*a<2*maxsep+2*/&&!found) {
            //case 3, k alone stacked
			//if (jldecrement) {
				if (en3==w5->f(j-1-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b]*/+
					w3->f(i+1,k+2)/*w3[i+1][a+1]*/+edangle3(k,l,k+1,ct2,data)+gap+constantclosure){

				  stack.push(j-1-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength()/*b*/,1,0,w5->f(j-1-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b]*/,true);
				  stack.push(i+1,k+2/*a+1*/,ct1->GetSequenceLength(),0,w3->f(i+1,k+2)/*w3[i+1][a+1]*/,true);
				  found = true;

				}
		    //}

            //case 8, j, k, and l stacked
            if ((j>ct1->GetSequenceLength()+1)&&!found&&k+2<=ct2->GetSequenceLength()+1&&l-2-ct2->GetSequenceLength()>=0) {
              if (en3==w5->f(j-2-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b]*/+
                  w3->f(i+1,k+2)/*w3[i+1][a+1]*/+edangle3(k,l,k+1,ct2,data)
                  +edangle5(j,i,j-1,ct1,data)+edangle5(l,k,l-1,ct2,data)+gap+constantclosure) {

                stack.push(j-2-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength()/*b*/,1,0,w5->f(j-2-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b]*/,true);
                stack.push(i+1,k+2/*a+1*/,ct1->GetSequenceLength(),0,w3->f(i+1,k+2)/*w3[i+1][a+1]*/,true);
                found = true;

              }
            }

            if (j-2-ct1->GetSequenceLength()>=0) if (!found) {
              //case 7, j and k stacked
              if ((j>ct1->GetSequenceLength()+1)&&k+2<=ct2->GetSequenceLength()+1) {

                if (en3==w5->f(j-2-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b+1]*/+
                    w3->f(i+1,k+2)/*w3[i+1][a+1]*/+edangle3(k,l,k+1,ct2,data)+
                    edangle5(j,i,j-1,ct1,data)+2*gap+constantclosure){

                  stack.push(j-2-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength()/*b+1*/,1,0,w5->f(j-2-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b+1]*/,true);
                  stack.push(i+1,k+2/*a+1*/,ct1->GetSequenceLength(),0,w3->f(i+1,k+2)/*w3[i+1][a+1]*/,true);
                  found = true;

                }
              }
            }
          }

		  if (j-2-ct1->GetSequenceLength()>=0) if (!found) {
            //case 5, j stacked
		      //f ((k<ct2->GetSequenceLength())&&(j>ct1->GetSequenceLength()+1)) {
              if (en3==w5->f(j-2-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b+1]*/+
                  w3->f(i+1,k+1)/*w3[i+1][a]*/+edangle5(j,i,j-1,ct1,data)+gap+constantclosure) {

                stack.push(j-2-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength()/*b+1*/,1,0,w5->f(j-2-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b+1]*/,true);
                stack.push(i+1,k+1/*a*/,ct1->GetSequenceLength(),0,w3->f(i+1,k+1)/*w3[i+1][a]*/,true);
                found = true;

              }
	      // }
	      //		else if ((k==ct2->GetSequenceLength()&&local)&&(j>ct1->GetSequenceLength()+1)) {
	      //    if (en3==w5->f(j-2-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b+1]*/+
	      //     +edangle5(j,i,j-1,ct1,data)+gap+constantclosure) {

	      //    stack.push(j-2-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength()/*b+1*/,1,0,w5->f(j-2-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b+1]*/,true);

	      //    found = true;

	      //   }
	      //   }
	      //		else if ((k==ct2->GetSequenceLength())&&(j>ct1->GetSequenceLength()+1)) {
	      //   if (en3==w5->f(j-2-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b+1]*/+
	      //    +edangle5(j,i,j-1,ct1,data)+gap+constantclosure+gap*(ct1->GetSequenceLength()-i)) {

	      //   stack.push(j-2-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength()/*b+1*/,1,0,w5->f(j-2-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b+1]*/,true);

	      //  found = true;

	      //}
	      // }

            //case 15, i, j, and k stacked
            if ((i<ct1->GetSequenceLength())&&(j>ct1->GetSequenceLength()+1)&&!found&&k+2<=ct2->GetSequenceLength()+1) {

              if (en3==w5->f(j-2-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b+1]*/+
                  w3->f(i+2,k+2)/*w3[i+2][a]*/+edangle5(j,i,j-1,ct1,data)
                  +edangle3(i,j,i+1,ct1,data)+edangle3(k,l,k+1,ct2,data)+gap+constantclosure) {

                stack.push(j-2-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength()/*b+1*/,1,0,w5->f(j-2-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b+1]*/,true);
                stack.push(i+2,k+2/*a*/,ct1->GetSequenceLength(),0,w3->f(i+2,k+2)/*w3[i+2][a]*/,true);
                found = true;


              }

            }
            if (!found) {
              //case 13, j and i stacked
              if ((i<ct1->GetSequenceLength())&&(j>ct1->GetSequenceLength()+1)) {
		//	  if (k<ct2->GetSequenceLength()) {
					  if (en3==w5->f(j-2-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b+1]*/+
						w3->f(i+2,k+1)/*w3[i+2][a-1]*/+edangle5(j,i,j-1,ct1,data)+
						edangle3(i,j,i+1,ct1,data)+2*gap+constantclosure){

						stack.push(j-2-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength()/*b+1*/,1,0,w5->f(j-2-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b+1]*/,true);
						stack.push(i+2,k+1/*a-1*/,ct1->GetSequenceLength(),0,w3->f(i+2,k+1)/*w3[i+2][a-1]*/,true);
						found = true;
					  }

					  //   }
					  //	  else if (k==ct2->GetSequenceLength()&&local) {
					  //		  if (en3==w5->f(j-2-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b+1]*/
		  //			+edangle5(j,i,j-1,ct1,data)+
					  //			edangle3(i,j,i+1,ct1,data)+2*gap+constantclosure){

					  //			stack.push(j-2-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength()/*b+1*/,1,0,w5->f(j-2-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b+1]*/,true);

					  //			found = true;
					  //		  }

					  //        }
					  //	  else if (k==ct2->GetSequenceLength()) {
					  //		  if (en3==w5->f(j-2-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b+1]*/
					  //			+edangle5(j,i,j-1,ct1,data)+
					  //			edangle3(i,j,i+1,ct1,data)+2*gap+constantclosure+gap*(ct1->GetSequenceLength()-i-1)){

					  //			stack.push(j-2-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength()/*b+1*/,1,0,w5->f(j-2-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b+1]*/,true);

					  //			found = true;
					  //		  }

					  //         }

              }
            }
          }

          if (/*i+2<ct1->GetSequenceLength()&&*/!found) {
            //case 9, i stacked
            if ((i<ct1->GetSequenceLength())) {

	      //	if (k+1<=ct2->GetSequenceLength()) {
					if (en3==w5->f(j-1-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b]*/+
						w3->f(i+2,k+1)/*w3[i+2][a-1]*/+edangle3(i,j,i+1,ct1,data)+gap+constantclosure) {

						stack.push(j-1-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength()/*b*/,1,0,w5->f(j-1-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b]*/,true);
						stack.push(i+2,k+1/*a-1*/,ct1->GetSequenceLength(),0,w3->f(i+2,k+1)/*w3[i+2][a-1]*/,true);
						found = true;

					}
					//	}
					//		else if (k==ct2->GetSequenceLength()&&local) {
					//			if (en3==w5->f(j-1-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b]*/
		//				+edangle3(i,j,i+1,ct1,data)+gap+constantclosure) {

		//				stack.push(j-1-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength()/*b*/,1,0,w5->f(j-1-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b]*/,true);

					//			found = true;
					//		}

					//	}
					///	else if (k==ct2->GetSequenceLength()) {
					//		if (en3==w5->f(j-1-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b]*/
		//			+edangle3(i,j,i+1,ct1,data)+gap+constantclosure+gap*(ct1->GetSequenceLength()-i-1)) {

					//			stack.push(j-1-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength()/*b*/,1,0,w5->f(j-1-ct1->GetSequenceLength(),l-1-ct2->GetSequenceLength())/*w5[j-1-ct1->GetSequenceLength()][b]*/,true);

					//		found = true;
					//	}

					//		}

            }
            //case 14, i, j, and l stacked
			if ((i<ct1->GetSequenceLength())&&(j>ct1->GetSequenceLength()+1)&&!found&&j-2-ct1->GetSequenceLength()>=0&&l-2-ct2->GetSequenceLength()>=0) {
			  //	if (k+1<=ct2->GetSequenceLength()) {
					if (en3==w5->f(j-2-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b]*/+
						w3->f(i+2,k+1)/*w3[i+2][a-1]*/+edangle3(i,j,i+1,ct1,data)
						+edangle5(j,i,j-1,ct1,data)+edangle5(l,k,l-1,ct2,data)+gap+constantclosure) {

					stack.push(j-2-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength()/*b*/,1,0,w5->f(j-2-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b]*/,true);
					stack.push(i+2,k+1/*a-1*/,ct1->GetSequenceLength(),0,w3->f(i+2,k+1)/*w3[i+2][a-1]*/,true);
					found = true;

					}
					//	}
					//		else if (k==ct2->GetSequenceLength()&&local) {
					//			if (en3==w5->f(j-2-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b]*/+
					//			edangle3(i,j,i+1,ct1,data)
					//			+edangle5(j,i,j-1,ct1,data)+edangle5(l,k,l-1,ct2,data)+gap+constantclosure) {

					//			stack.push(j-2-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength()/*b*/,1,0,w5->f(j-2-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b]*/,true);

					//			found = true;

					//		}
					// 	}
					//		else if (k==ct2->GetSequenceLength()) {
					//		if (en3==w5->f(j-2-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b]*/+
					//			edangle3(i,j,i+1,ct1,data)
					//			 +edangle5(j,i,j-1,ct1,data)+edangle5(l,k,l-1,ct2,data)+gap+constantclosure+gap*(ct1->GetSequenceLength()-i-1)) {

					//			stack.push(j-2-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength()/*b*/,1,0,w5->f(j-2-ct1->GetSequenceLength(),l-2-ct2->GetSequenceLength())/*w5[j-2-ct1->GetSequenceLength()][b]*/,true);

					//			found = true;
					//		}

					//		}
            }
          }

		  //reset jl decrement2 for remainder of code... i.e. undo the exception above
		  if (j-2==ct1->GetSequenceLength()) jldecrement2 = l-2 <= highend[j-2] && l-2 >= lowend[j-2];


		}

        //Now consider multiloop:
        //junction closed by i-j pair aligned with a and b
        //calculate the free energy of 2 fragments merged:
        if (i!=ct1->GetSequenceLength()&&j!=ct1->GetSequenceLength()+1&&k!=ct2->GetSequenceLength()&&l!=ct2->GetSequenceLength()+1) {
          for (c=i+minloop+1;c<j-minloop&&!found;c++) {
            //if (c>ct1->GetSequenceLength()) kp=ct1->GetSequenceLength()-ct2->GetSequenceLength();
            //else kp = 0;
            for (d=max(k+minloop+1,lowend[c]/*lowlimit(c,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())/*c-maxsep-kp*/);
                 d<l-minloop&&d<=highend[c]/*highlimit(c,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())/*(c+maxsep-kp*/&&!found;d++) {
              //e = d-c+maxsep+kp;

              if ((((c<ct1->GetSequenceLength())&&d<ct2->GetSequenceLength())||((c>ct1->GetSequenceLength())&&(d>ct2->GetSequenceLength())))
				  &&(c!=ct1->GetSequenceLength())&&(d!=ct2->GetSequenceLength())&&d+1>=lowend[c+1]/*lowlimit(c+1,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/
				  &&d+1<=highend[c+1]&&d+1>=lowend[c+1]) {


                //must consider whether an unpaired nucleotide is stacked
                //is stacked onto each of the nucleotides in a helix
                //There are 16 cases:
                //consider rthem in this order (0 unstacked, 1 stacked)
                //		i	j	k	l
                //1		0	0	0	0
                //2		0	0	0	1
                //3		0	0	1	0
                //4		0	0	1	1
                //5		0	1	0	0
                //6		0	1	0	1
                //7		0	1	1	0
                //8		0	1	1	1
                //9		1	0	0	0
                //10	1	0	0	1
                //11	1	0	1	0
                //12	1	0	1	1
                //13	1	1	0	0
                //14	1	1	0	1
                //15	1	1	1	0
                //16	1	1	1	1

                //case 1 - no stacks
                if (ikincrement&&jldecrement) {
                  if(en3== w->f(i+1,c,k+1,d)/*w[c][iref(i+1,c,ct1->GetSequenceLength())][a][e]*/
                     +w->f(c+1,j-1,d+1,l-1)/*w[jref(c+1,j-1,ct1->GetSequenceLength())][iref(c+1,j-1,ct1->GetSequenceLength())][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                     +constantclosure) {

                    stack.push(i+1,c,k+1,d,w->f(i+1,c,k+1,d)/*w[c][iref(i+1,c,ct1->GetSequenceLength())][a][e]*/);
                    stack.push(ideref(c+1,j-1,ct1->GetSequenceLength())/*ideref(c+1,j-1,ct1->GetSequenceLength())*/,jderef(c+1,j-1,ct1->GetSequenceLength())/*jderef(c+1,j-1,ct1->GetSequenceLength())*/,ideref(d+1,l-1,ct2->GetSequenceLength())/*e*/,jderef(d+1,l-1,ct2->GetSequenceLength())/*b*/,w->f(c+1,j-1,d+1,l-1)/*w[jref(c+1,j-1,ct1->GetSequenceLength())][iref(c+1,j-1,ct1->GetSequenceLength())][e][b]*/);
                    found = true;
                  }

                }

                //case 16 - all four stacked
                if (!found&&ikincrement2&&jldecrement2) {
                  if((en3 ==w->f(i+2,c,k+2,d)/*w[c][iref(i+2,c,ct1->GetSequenceLength())][a][e]*/
                      +w->f(c+1,j-2,d+1,l-2)/*w[jref(c+1,j-2,ct1->GetSequenceLength())][iref(c+1,j-2,ct1->GetSequenceLength())][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                      +4*data->eparam[6]+edangle5(l,k,l-1,ct2,data)
                      +edangle3(k,l,k+1,ct2,data)
                      +edangle5(j,i,j-1,ct1,data)
                      +edangle3(i,j,i+1,ct1,data)+constantclosure)&&
                     ((i+1)!=ct1->GetSequenceLength())&&((j-1)!=ct1->GetSequenceLength())&&((k+1)!=ct2->GetSequenceLength())&&((l-1)!=ct2->GetSequenceLength())) {

                    found = true;
                    stack.push(i+2,c,k+2/*a*/,d/*e*/,w->f(i+2,c,k+2,d)/*w[c][iref(i+2,c,ct1->GetSequenceLength())][a][e]*/);
                    stack.push(ideref(c+1,j-2,ct1->GetSequenceLength()),jderef(c+1,j-2,ct1->GetSequenceLength()),ideref(d+1,l-2,ct2->GetSequenceLength())/*e*/,jderef(d+1,l-2,ct2->GetSequenceLength())/*b*/,w->f(c+1,j-2,d+1,l-2)/*w[jref(c+1,j-2,ct1->GetSequenceLength())][iref(c+1,j-2,ct1->GetSequenceLength())][e][b]*/);

                  }
                }

                //case 6 - j and l stacked:
                if (ikincrement&&jldecrement2&&!found) {
                  if((en3 == w->f(i+1,c,k+1,d)/*w[c][iref(i+1,c,ct1->GetSequenceLength())][a][e]*/
                      +w->f(c+1,j-2,d+1,l-2)/*w[jref(c+1,j-2,ct1->GetSequenceLength())][iref(c+1,j-2,ct1->GetSequenceLength())][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                      +2*data->eparam[6]+edangle5(j,i,j-1,ct1,data)
                      +edangle5(l,k,l-1,ct2,data)+constantclosure)&&((j-1)!=ct1->GetSequenceLength())&&((l-1)!=ct2->GetSequenceLength())){

                    found = true;
                    stack.push(i+1,c,k+1/*a*/,d/*e*/,w->f(i+1,c,k+1,d)/*w[c][iref(i+1,c,ct1->GetSequenceLength())][a][e]*/);
                    stack.push(ideref(c+1,j-2,ct1->GetSequenceLength()),jderef(c+1,j-2,ct1->GetSequenceLength()),ideref(d+1,l-2,ct2->GetSequenceLength())/*e*/,jderef(d+1,l-2,ct2->GetSequenceLength())/*b*/,w->f(c+1,j-2,d+1,l-2)/*w[jref(c+1,j-2,ct1->GetSequenceLength())][iref(c+1,j-2,ct1->GetSequenceLength())][e][b]*/);
                  }
                }

                //case 11 - i and k stacked
                if (ikincrement2&&jldecrement&&!found) {
                  if((en3 ==w->f(i+2,c,k+2,d)/*w[c][iref(i+2,c,ct1->GetSequenceLength())][a][e]*/
                      +w->f(c+1,j-1,d+1,l-1)/*w[jref(c+1,j-1,ct1->GetSequenceLength())][iref(c+1,j-1,ct1->GetSequenceLength())][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                      +2*data->eparam[6]+edangle3(k,l,k+1,ct2,data)+
                      edangle3(i,j,i+1,ct1,data)+constantclosure)&&((i+1)!=ct1->GetSequenceLength())&&((k+1)!=ct2->GetSequenceLength())){

                    found = true;
                    stack.push(i+2,c,k+2/*a*/,d/*e*/,w->f(i+2,c,k+2,d)/*w[c][iref(i+2,c,ct1->GetSequenceLength())][a][e]*/);
                    stack.push(ideref(c+1,j-1,ct1->GetSequenceLength()),jderef(c+1,j-1,ct1->GetSequenceLength()),ideref(d+1,l-1,ct2->GetSequenceLength())/*e*/,jderef(d+1,l-1,ct2->GetSequenceLength())/*b*/,w->f(c+1,j-1,d+1,l-1)/*w[jref(c+1,j-1,ct1->GetSequenceLength())][iref(c+1,j-1,ct1->GetSequenceLength())][e][b]*/);
                  }
                }

                if (l-2>=lowend[j-1]&&l-2<=highend[j-1]/*lowlimit(j-1,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())/*(b-1>=0*/&&!found&&((l-1)!=ct2->GetSequenceLength())) {
                  //case 2 - stack on l
                  if (ikincrement) {
                    if(en3== w->f(i+1,c,k+1,d)/*w[c][iref(i+1,c,ct1->GetSequenceLength())][a][e]*/
                       +w->f(c+1,j-1,d+1,l-2)/*w[jref(c+1,j-1,ct1->GetSequenceLength())][iref(c+1,j-1,ct1->GetSequenceLength())][e][b-1]*/+2*data->eparam[5]+2*data->eparam[10]
                       +data->eparam[6]+edangle5(l,k,l-1,ct2,data)+gap+constantclosure){

                      found = true;
                      stack.push(i+1,c,k+1/*a*/,d/*e*/,w->f(i+1,c,k+1,d)/*w[c][iref(i+1,c,ct1->GetSequenceLength())][a][e]*/);
                      stack.push(ideref(c+1,j-1,ct1->GetSequenceLength()),jderef(c+1,j-1,ct1->GetSequenceLength()),ideref(d+1,l-2,ct2->GetSequenceLength())/*e*/,jderef(d+1,l-2,ct2->GetSequenceLength())/*b-1*/,w->f(c+1,j-1,d+1,l-2)/*w[jref(c+1,j-1,ct1->GetSequenceLength())][iref(c+1,j-1,ct1->GetSequenceLength())][e][b-1]*/);

                    }
                  }

                  //case 12 - i, k, and l stacked
                  if (!found&&ikincrement2) {
                    if((en3 ==w->f(i+2,c,k+2,d)/*w[c][iref(i+2,c,ct1->GetSequenceLength())][a][e]*/
                        +w->f(c+1,j-1,d+1,l-2)/*w[jref(c+1,j-1,ct1->GetSequenceLength())][iref(c+1,j-1,ct1->GetSequenceLength())][e][b-1]*/+2*data->eparam[5]+2*data->eparam[10]
                        +3*data->eparam[6]+edangle5(l,k,l-1,ct2,data)
                        +edangle3(i,j,i+1,ct1,data)+edangle3(k,l,k+1,ct2,data)+gap
                        +constantclosure)&&((i+1)!=ct1->GetSequenceLength())){


                      found = true;
                      stack.push(i+2,c,k+2/*a*/,d/*e*/,w->f(i+2,c,k+2,d)/*w[c][iref(i+2,c,ct1->GetSequenceLength())][a][e]*/);
                      stack.push(ideref(c+1,j-1,ct1->GetSequenceLength()),jderef(c+1,j-1,ct1->GetSequenceLength()),ideref(d+1,l-2,ct2->GetSequenceLength())/*e*/,jderef(d+1,l-2,ct2->GetSequenceLength())/*b-1*/,w->f(c+1,j-1,d+1,l-2)/*w[jref(c+1,j-1,ct1->GetSequenceLength())][iref(c+1,j-1,ct1->GetSequenceLength())][e][b-1]*/);


                    }
                  }

                  if (k+1>=lowend[i+2]&&k+1<=highend[i+2]/*lowlimit(i+2,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())/*a-1>0*/&&!found&&(i+1)!=ct1->GetSequenceLength()) {
                    //case 10 - l and i stacked
                    if(en3 ==w->f(i+2,c,k+1,d)/*w[c][iref(i+2,c,ct1->GetSequenceLength())][a-1][e]*/
                       +w->f(c+1,j-1,d+1,l-2)/*w[jref(c+1,j-1,ct1->GetSequenceLength())][iref(c+1,j-1,ct1->GetSequenceLength())][e][b-1]*/+2*data->eparam[5]+2*data->eparam[10]
                       +2*data->eparam[6]+edangle5(l,k,l-1,ct2,data)
                       +edangle3(i,j,i+1,ct1,data)+2*gap+constantclosure){


                      found = true;
                      stack.push(i+2,c,k+1/*a-1*/,d/*e*/,w->f(i+2,c,k+1,d)/*w[c][iref(i+2,c,ct1->GetSequenceLength())][a-1][e]*/);
                      stack.push(ideref(c+1,j-1,ct1->GetSequenceLength()),jderef(c+1,j-1,ct1->GetSequenceLength()),ideref(d+1,l-2,ct2->GetSequenceLength())/*e*/,jderef(d+1,l-2,ct2->GetSequenceLength())/*b-1*/,w->f(c+1,j-1,d+1,l-2)/*w[jref(c+1,j-1,ct1->GetSequenceLength())][iref(c+1,j-1,ct1->GetSequenceLength())][e][b-1]*/);


                    }

                  }

                  if (k+2<=highend[i+1]&&k+2>=lowend[i+1]/*highlimit(i+1,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())/*a+1<2*maxsep+2*/&&!found&&((k+1)!=ct2->GetSequenceLength())) {
                    //case 4 - k and l stacked
                    if(en3 ==w->f(i+1,c,k+2,d)/*w[c][iref(i+1,c,ct1->GetSequenceLength())][a+1][e]*/
                       +w->f(c+1,j-1,d+1,l-2)/*w[jref(c+1,j-1,ct1->GetSequenceLength())][iref(c+1,j-1,ct1->GetSequenceLength())][e][b-1]*/+2*data->eparam[5]+2*data->eparam[10]
                       +2*data->eparam[6]+edangle5(l,k,l-1,ct2,data)
                       +edangle3(k,l,k+1,ct2,data)+2*gap+constantclosure){

                      found = true;
                      stack.push(i+1,c,k+2/*a+1*/,d/*e*/,w->f(i+1,c,k+2,d)/*w[c][iref(i+1,c,ct1->GetSequenceLength())][a+1][e]*/);
                      stack.push(ideref(c+1,j-1,ct1->GetSequenceLength()),jderef(c+1,j-1,ct1->GetSequenceLength()),ideref(d+1,l-2,ct2->GetSequenceLength())/*e*/,jderef(d+1,l-2,ct2->GetSequenceLength())/*b-1*/,w->f(c+1,j-1,d+1,l-2)/*w[jref(c+1,j-1,ct1->GetSequenceLength())][iref(c+1,j-1,ct1->GetSequenceLength())][e][b-1]*/);


                    }




                  }





                }
                if (k+1>=lowend[i+2]&&k+1<=highend[i+2]/*lowlimit(i+2,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())/*a-1>0*/&&!found&&(i+1!=ct1->GetSequenceLength())) {
                  //case 14 - i, j, and l stacked
                  if (jldecrement2) {
                    if((en3 ==w->f(i+2,c,k+1,d)/*w[c][iref(i+2,c,ct1->GetSequenceLength())][a-1][e]*/
                        +w->f(c+1,j-2,d+1,l-2)/*w[jref(c+1,j-2,ct1->GetSequenceLength())][iref(c+1,j-2,ct1->GetSequenceLength())][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                        +3*data->eparam[6]+edangle5(l,k,l-1,ct2,data)
                        +edangle5(j,i,j-1,ct1,data)
                        +edangle3(i,j,i+1,ct1,data)+gap+constantclosure)&&(j-1!=ct1->GetSequenceLength())){


                      found = true;
                      stack.push(i+2,c,k+1/*a-1*/,d/*e*/,w->f(i+2,c,k+1,d)/*w[c][iref(i+2,c,ct1->GetSequenceLength())][a-1][e]*/);
                      stack.push(ideref(c+1,j-2,ct1->GetSequenceLength()),jderef(c+1,j-2,ct1->GetSequenceLength()),ideref(d+1,l-2,ct2->GetSequenceLength())/*e*/,jderef(d+1,l-2,ct2->GetSequenceLength())/*b*/,w->f(c+1,j-2,d+1,l-2)/*w[jref(c+1,j-2,ct1->GetSequenceLength())][iref(c+1,j-2,ct1->GetSequenceLength())][e][b]*/);

                    }
                  }

                  //case 9 - i stacked
                  if (jldecrement&&!found) {
                    if(en3 ==w->f(i+2,c,k+1,d)/*w[c][iref(i+2,c,ct1->GetSequenceLength())][a-1][e]*/
                       +w->f(c+1,j-1,d+1,l-1)/*w[jref(c+1,j-1,ct1->GetSequenceLength())][iref(c+1,j-1,ct1->GetSequenceLength())][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                       +data->eparam[6]+edangle3(i,j,i+1,ct1,data)+gap+constantclosure){

                      found = true;
                      stack.push(i+2,c,k+1/*a-1*/,d/*e*/,w->f(i+2,c,k+1,d)/*w[c][iref(i+2,c,ct1->GetSequenceLength())][a-1][e]*/);
                      stack.push(ideref(c+1,j-1,ct1->GetSequenceLength()),jderef(c+1,j-1,ct1->GetSequenceLength()),ideref(d+1,l-1,ct2->GetSequenceLength())/*e*/,jderef(d+1,l-1,ct2->GetSequenceLength())/*b*/,w->f(c+1,j-1,d+1,l-1)/*w[jref(c+1,j-1,ct1->GetSequenceLength())][iref(c+1,j-1,ct1->GetSequenceLength())][e][b]*/);


                    }
                  }

                  if (!found&&l-1<=highend[j-2]&&l-1>=lowend[j-2]/*highlimit(j-2,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())/*(b+1<2*maxsep+2)*/&&(j-1!=ct1->GetSequenceLength())) {
                    //case 13 - i and j stacked
                    if(en3 ==w->f(i+2,c,k+1,d)/*w[c][iref(i+2,c,ct1->GetSequenceLength())][a-1][e]*/
                       +w->f(c+1,j-2,d+1,l-1)/*w[jref(c+1,j-2,ct1->GetSequenceLength())][iref(c+1,j-2,ct1->GetSequenceLength())][e][b+1]*/+2*data->eparam[5]+2*data->eparam[10]
                       +2*data->eparam[6]+edangle5(j,i,j-1,ct1,data)
                       +edangle3(i,j,i+1,ct1,data)+2*gap+constantclosure){


                      found = true;
                      stack.push(i+2,c,k+1/*a-1*/,d/*e*/,w->f(i+2,c,k+1,d)/*w[c][iref(i+2,c,ct1->GetSequenceLength())][a-1][e]*/);
                      stack.push(ideref(c+1,j-2,ct1->GetSequenceLength()),jderef(c+1,j-2,ct1->GetSequenceLength()),ideref(d+1,l-1,ct2->GetSequenceLength())/*e*/,jderef(d+1,l-1,ct2->GetSequenceLength())/*b+1*/,w->f(c+1,j-2,d+1,l-1)/*w[jref(c+1,j-2,ct1->GetSequenceLength())][iref(c+1,j-2,ct1->GetSequenceLength())][e][b+1]*/);

                    }

                  }

                }

                if (k+2<=highend[i+1]&&k+2>=lowend[i+1]/*highlimit(i+1,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())/*a+1<2*maxsep+2*/&&!found&&(k+1!=ct2->GetSequenceLength())) {
                  //case 8 - j, k, and l stacked
                  if (jldecrement2) {
                    if((en3==w->f(i+1,c,k+2,d)/*w[c][iref(i+1,c,ct1->GetSequenceLength())][a+1][e]*/
                        +w->f(c+1,j-2,d+1,l-2)/*w[jref(c+1,j-2,ct1->GetSequenceLength())][iref(c+1,j-2,ct1->GetSequenceLength())][e][b]*/+2*data->eparam[5]+2*data->eparam[10]
                        +3*data->eparam[6]+edangle5(j,i,j-1,ct1,data)
                        +edangle5(l,k,l-1,ct2,data)+edangle3(k,l,k+1,ct2,data)+gap
                        +constantclosure)&&(j-1!=ct1->GetSequenceLength())){

                      found = true;
                      stack.push(i+1,c,k+2/*a+1*/,d/*e*/,w->f(i+1,c,k+2,d)/*w[c][iref(i+1,c,ct1->GetSequenceLength())][a+1][e]*/);
                      stack.push(ideref(c+1,j-2,ct1->GetSequenceLength()),jderef(c+1,j-2,ct1->GetSequenceLength()),ideref(d+1,l-2,ct2->GetSequenceLength())/*e*/,jderef(d+1,l-2,ct2->GetSequenceLength())/*b*/,w->f(c+1,j-2,d+1,l-2)/*w[jref(c+1,j-2,ct1->GetSequenceLength())][iref(c+1,j-2,ct1->GetSequenceLength())][e][b]*/);

                    }
                  }

                  //case 3 - stack on k
                  if (jldecrement&&!found) {
                    if(en3 ==w->f(i+1,c,k+2,d)/*w[c][iref(i+1,c,ct1->GetSequenceLength())][a+1][e]*/+2*data->eparam[10]
                       +w->f(c+1,j-1,d+1,l-1)/*w[jref(c+1,j-1,ct1->GetSequenceLength())][iref(c+1,j-1,ct1->GetSequenceLength())][e][b]*/+2*data->eparam[5]
                       +data->eparam[6]+edangle3(k,l,k+1,ct2,data)+gap+constantclosure){

                      found = true;
                      stack.push(i+1,c,k+2/*a+1*/,d/*e*/,w->f(i+1,c,k+2,d)/*w[c][iref(i+1,c,ct1->GetSequenceLength())][a+1][e]*/);
                      stack.push(ideref(c+1,j-1,ct1->GetSequenceLength()),jderef(c+1,j-1,ct1->GetSequenceLength()),ideref(d+1,l-1,ct2->GetSequenceLength())/*e*/,jderef(d+1,l-1,ct2->GetSequenceLength())/*b*/,w->f(c+1,j-1,d+1,l-1)/*w[jref(c+1,j-1,ct1->GetSequenceLength())][iref(c+1,j-1,ct1->GetSequenceLength())][e][b]*/);

                    }
                  }

                  if (!found&&l-1<=highend[j-2]&&l-1>=lowend[j-2]/*highlimit(j-2,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())/*b+1<2*maxsep+2*/) {
                    //case 7 - j and k stacked
                    if ((en3 ==w->f(i+1,c,k+2,d)/*w[c][iref(i+1,c,ct1->GetSequenceLength())][a+1][e]*/
                         +w->f(c+1,j-2,d+1,l-1)/*w[jref(c+1,j-2,ct1->GetSequenceLength())][iref(c+1,j-2,ct1->GetSequenceLength())][e][b+1]*/+2*data->eparam[5]+2*data->eparam[10]
                         +2*data->eparam[6]+edangle5(j,i,j-1,ct1,data)
                         +edangle3(k,l,k+1,ct2,data)+2*gap+constantclosure)&&(j-1!=ct1->GetSequenceLength())){

                      found = true;
                      stack.push(i+1,c,k+2/*a+1*/,d/*e*/,w->f(i+1,c,k+2,d)/*w[c][iref(i+1,c,ct1->GetSequenceLength())][a+1][e]*/);
                      stack.push(ideref(c+1,j-2,ct1->GetSequenceLength()),jderef(c+1,j-2,ct1->GetSequenceLength()),ideref(d+1,l-1,ct2->GetSequenceLength())/*e*/,jderef(d+1,l-1,ct2->GetSequenceLength())/*b+1*/,w->f(c+1,j-2,d+1,l-1)/*w[jref(c+1,j-2,ct1->GetSequenceLength())][iref(c+1,j-2,ct1->GetSequenceLength())][e][b+1]*/);

                    }

                  }


                }
                if (l-1<=highend[j-2]&&l-1>=lowend[j-2]/*highlimit(j-2,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())/*b+1<2*maxsep+2*/&&!found&&(j-1!=ct1->GetSequenceLength())) {
                  //case 15 - i, j, and k stacked
                  if (ikincrement2) {
                    if((en3 ==w->f(i+2,c,k+2,d)/*w[c][iref(i+2,c,ct1->GetSequenceLength())][a][e]*/
                        +w->f(c+1,j-2,d+1,l-1)/*w[jref(c+1,j-2,ct1->GetSequenceLength())][iref(c+1,j-2,ct1->GetSequenceLength())][e][b+1]*/+2*data->eparam[5]+2*data->eparam[10]
                        +3*data->eparam[6]+edangle3(k,l,k+1,ct2,data)
                        +edangle5(j,i,j-1,ct1,data)
                        +edangle3(i,j,i+1,ct1,data)+gap+constantclosure)&&(i+1!=ct1->GetSequenceLength())&&(k+1!=ct2->GetSequenceLength())){


                      found = true;
                      stack.push(i+2,c,k+2/*a*/,d/*e*/,w->f(i+2,c,k+2,d)/*w[c][iref(i+2,c,ct1->GetSequenceLength())][a][e]*/);
                      stack.push(ideref(c+1,j-2,ct1->GetSequenceLength()),jderef(c+1,j-2,ct1->GetSequenceLength()),ideref(d+1,l-1,ct2->GetSequenceLength())/*e*/,jderef(d+1,l-1,ct2->GetSequenceLength())/*b+1*/,w->f(c+1,j-2,d+1,l-1)/*w[jref(c+1,j-2,ct1->GetSequenceLength())][iref(c+1,j-2,ct1->GetSequenceLength())][e][b+1]*/);

                    }
                  }

                  //case 5 - j stacked
                  if (!found&&ikincrement) {
                    if(en3 ==w->f(i+1,c,k+1,d)/*w[c][iref(i+1,c,ct1->GetSequenceLength())][a][e]*/
                       +w->f(c+1,j-2,d+1,l-1)/*w[jref(c+1,j-2,ct1->GetSequenceLength())][iref(c+1,j-2,ct1->GetSequenceLength())][e][b+1]*/+2*data->eparam[5]+2*data->eparam[10]
                       +data->eparam[6]+edangle5(j,i,j-1,ct1,data)+gap+constantclosure){

                      found = true;
                      stack.push(i+1,c,k+1/*a*/,d/*e*/,w->f(i+1,c,k+1,d)/*w[c][iref(i+1,c,ct1->GetSequenceLength())][a][e]*/);
                      stack.push(ideref(c+1,j-2,ct1->GetSequenceLength()),jderef(c+1,j-2,ct1->GetSequenceLength()),ideref(d+1,l-1,ct2->GetSequenceLength())/*e*/,jderef(d+1,l-1,ct2->GetSequenceLength())/*b+1*/,w->f(c+1,j-2,d+1,l-1)/*w[jref(c+1,j-2,ct1->GetSequenceLength())][iref(c+1,j-2,ct1->GetSequenceLength())][e][b+1]*/);

                    }
                  }
                }

              }




            }
          }
        }

#ifdef DYNALIGN_II


	//	 if(i==45&&j==93&&k==45&&l==105)cerr<<"constantclosure "<<constantclosure<<"\n";
	if (i!=ct1->GetSequenceLength()&&j!=ct1->GetSequenceLength()+1&&k!=ct2->GetSequenceLength()&&l!=ct2->GetSequenceLength()+1){
	  
	  //newest
	  short int idangle = edangle3(i,j,i+1,ct1,data);
	  short int jdangle = edangle5(j,i,j-1,ct1,data);
	  short int kdangle = edangle3(k,l,k+1,ct2,data);
	  short int ldangle = edangle5(l,k,l-1,ct2,data);

	  for (c=i+1+minloop;!found&&c+1<j-1-minloop&&c!=ct1->GetSequenceLength();++c) {

	    //case1 0 0 0 0
	    trace_branch_1(i+1,j-1,k+1,l-1,c,found,ct1,data,single1_v, single1_w, single1_wmb, 
			   single1_lfce, single1_fce, single1_mod,w,stack,lowend,highend,en3,2*data->eparam[5]+2*data->eparam[10]+constantclosure, iintercept, islope);

	    //case2 0 0 0 1
	    trace_branch_1(i+1,j-1,k+1,l-2,c,found,ct1,data,single1_v, single1_w, single1_wmb, 
			   single1_lfce, single1_fce, single1_mod,w,stack,lowend,highend,en3,ldangle+gap+2*data->eparam[5]+2*data->eparam[10]
		    +data->eparam[6]+constantclosure, iintercept, islope);

	    //case3 0 0 1 0
	    trace_branch_1(i+1,j-1,k+2,l-1,c,found,ct1,data,single1_v, single1_w, single1_wmb, 
			   single1_lfce, single1_fce, single1_mod,w,stack,lowend,highend,en3,kdangle+gap+2*data->eparam[5]+2*data->eparam[10]
			   +data->eparam[6]+constantclosure, iintercept, islope);
	    
	    //case4 0 0 1 1
	    trace_branch_1(i+1,j-1,k+2,l-2,c,found,ct1,data,single1_v, single1_w, single1_wmb, 
			   single1_lfce, single1_fce, single1_mod,w,stack,lowend,highend,en3,kdangle+ldangle+2*gap+2*data->eparam[5]+2*data->eparam[10]
			   +2*data->eparam[6]+constantclosure, iintercept, islope);
	    //case5 0 1 0 0
	    trace_branch_1(i+1,j-2,k+1,l-1,c,found,ct1,data,single1_v, single1_w, single1_wmb, 
			   single1_lfce, single1_fce, single1_mod,w,stack,lowend,highend,en3,jdangle+gap+2*data->eparam[5]+2*data->eparam[10]
			   +data->eparam[6]+constantclosure, iintercept, islope);	    

	    //case6 0 1 0 1
	    trace_branch_1(i+1,j-2,k+1,l-2,c,found,ct1,data,single1_v, single1_w, single1_wmb, 
			   single1_lfce, single1_fce, single1_mod,w,stack,lowend,highend,en3,jdangle+ldangle+2*data->eparam[5]+2*data->eparam[10]
		     +2*data->eparam[6]+constantclosure, iintercept, islope);

	    //case7 0 1 1 0
	    trace_branch_1(i+1,j-2,k+2,l-1,c,found,ct1,data,single1_v, single1_w, single1_wmb, 
			   single1_lfce, single1_fce, single1_mod,w,stack,lowend,highend,en3,jdangle+kdangle+2*gap+2*data->eparam[5]+2*data->eparam[10]
		      +2*data->eparam[6]+constantclosure, iintercept, islope);

	    //case8 0 1 1 1
	    trace_branch_1(i+1,j-2,k+2,l-2,c,found,ct1,data,single1_v, single1_w, single1_wmb, 
			   single1_lfce, single1_fce, single1_mod,w,stack,lowend,highend,en3,jdangle+kdangle+ldangle+gap+2*data->eparam[5]+2*data->eparam[10]
		       +3*data->eparam[6]+constantclosure, iintercept, islope);

	    //case9 1 0 0 0
	    trace_branch_1(i+2,j-1,k+1,l-1,c,found,ct1,data,single1_v, single1_w, single1_wmb, 
			   single1_lfce, single1_fce, single1_mod,w,stack,lowend,highend,en3,idangle+gap+2*data->eparam[5]+2*data->eparam[10]
			   +data->eparam[6]+constantclosure, iintercept, islope);

	    //case10 1 0 0 1
	    trace_branch_1(i+2,j-1,k+1,l-2,c,found,ct1,data,single1_v, single1_w, single1_wmb, 
			   single1_lfce, single1_fce, single1_mod,w,stack,lowend,highend,en3,idangle+ldangle+2*gap+2*data->eparam[5]+2*data->eparam[10]
		      +2*data->eparam[6]+constantclosure, iintercept, islope);
	    
	    //case11 1 0 1 0
	    trace_branch_1(i+2,j-1,k+2,l-1,c,found,ct1,data,single1_v, single1_w, single1_wmb, 
			   single1_lfce, single1_fce, single1_mod,w,stack,lowend,highend,en3,idangle+kdangle+2*data->eparam[5]+2*data->eparam[10]
		      +2*data->eparam[6]+constantclosure, iintercept, islope);

	    //case12 1 0 1 1
	    trace_branch_1(i+2,j-1,k+2,l-2,c,found,ct1,data,single1_v, single1_w, single1_wmb, 
			   single1_lfce, single1_fce, single1_mod,w,stack,lowend,highend,en3,idangle+kdangle+ldangle+gap+2*data->eparam[5]+2*data->eparam[10]
		      +3*data->eparam[6]+constantclosure, iintercept, islope);

	    //case13 1 1 0 0
	    trace_branch_1(i+2,j-2,k+1,l-1,c,found,ct1,data,single1_v, single1_w, single1_wmb, 
			   single1_lfce, single1_fce, single1_mod,w,stack,lowend,highend,en3,idangle+jdangle+2*gap+2*data->eparam[5]+2*data->eparam[10]
		      +2*data->eparam[6]+constantclosure, iintercept, islope);
	    
	    //case14 1 1 0 1
	    trace_branch_1(i+2,j-2,k+1,l-2,c,found,ct1,data,single1_v, single1_w, single1_wmb, 
			   single1_lfce, single1_fce, single1_mod,w,stack,lowend,highend,en3,idangle+jdangle+ldangle+gap+2*data->eparam[5]+2*data->eparam[10]
			   +3*data->eparam[6]+constantclosure, iintercept, islope);

	    //case15 1 1 1 0
	    trace_branch_1(i+2,j-2,k+2,l-1,c,found,ct1,data,single1_v, single1_w, single1_wmb, 
			   single1_lfce, single1_fce, single1_mod,w,stack,lowend,highend,en3,idangle+jdangle+kdangle+gap+2*data->eparam[5]+2*data->eparam[10]
			   +3*data->eparam[6]+constantclosure, iintercept, islope);
	    
	    //case16 1 1 1 1
	    trace_branch_1(i+2,j-2,k+2,l-2,c,found,ct1,data,single1_v, single1_w, single1_wmb, 
			   single1_lfce, single1_fce, single1_mod,w,stack,lowend,highend,en3,idangle+jdangle+kdangle+ldangle+2*data->eparam[5]+2*data->eparam[10]
			   +4*data->eparam[6]+constantclosure, iintercept, islope);
	  }

	  
	  for (d=k+1+minloop;!found&&d+1<l-1-minloop&&d!=ct2->GetSequenceLength();++d) {

	    //case1 0 0 0 0
	    trace_branch_2(i+1,j-1,k+1,l-1,d,found,ct2,data,single2_v, single2_w, single2_wmb, 
			   single2_lfce, single2_fce, single2_mod,w,stack,lowend,highend,en3,2*data->eparam[5]+2*data->eparam[10]+constantclosure, iintercept, islope);

	    //case2 0 0 0 1
	    trace_branch_2(i+1,j-1,k+1,l-2,d,found,ct2,data,single2_v, single2_w, single2_wmb, 
			   single2_lfce, single2_fce, single2_mod,w,stack,lowend,highend,en3,ldangle+gap+2*data->eparam[5]+2*data->eparam[10]
		    +data->eparam[6]+constantclosure, iintercept, islope);

	    //case3 0 0 1 0
	    trace_branch_2(i+1,j-1,k+2,l-1,d,found,ct2,data,single2_v, single2_w, single2_wmb, 
			   single2_lfce, single2_fce, single2_mod,w,stack,lowend,highend,en3,kdangle+gap+2*data->eparam[5]+2*data->eparam[10]
			   +data->eparam[6]+constantclosure, iintercept, islope);
	    
	    //case4 0 0 1 1
	    trace_branch_2(i+1,j-1,k+2,l-2,d,found,ct2,data,single2_v, single2_w, single2_wmb, 
			   single2_lfce, single2_fce, single2_mod,w,stack,lowend,highend,en3,kdangle+ldangle+2*gap+2*data->eparam[5]+2*data->eparam[10]
			   +2*data->eparam[6]+constantclosure, iintercept, islope);
	    //case5 0 1 0 0
	    trace_branch_2(i+1,j-2,k+1,l-1,d,found,ct2,data,single2_v, single2_w, single2_wmb, 
			   single2_lfce, single2_fce, single2_mod,w,stack,lowend,highend,en3,jdangle+gap+2*data->eparam[5]+2*data->eparam[10]
			   +data->eparam[6]+constantclosure, iintercept, islope);
	    
	    //case6 0 1 0 1
	    trace_branch_2(i+1,j-2,k+1,l-2,d,found,ct2,data,single2_v, single2_w, single2_wmb, 
			   single2_lfce, single2_fce, single2_mod,w,stack,lowend,highend,en3,jdangle+ldangle+2*data->eparam[5]+2*data->eparam[10]
		     +2*data->eparam[6]+constantclosure, iintercept, islope);

	    //case7 0 1 1 0
	    trace_branch_2(i+1,j-2,k+2,l-1,d,found,ct1,data,single2_v, single2_w, single2_wmb, 
			   single2_lfce, single2_fce, single2_mod,w,stack,lowend,highend,en3,jdangle+kdangle+2*gap+2*data->eparam[5]+2*data->eparam[10]
		      +2*data->eparam[6]+constantclosure, iintercept, islope);

	    //case8 0 1 1 1
	    trace_branch_2(i+1,j-2,k+2,l-2,d,found,ct2,data,single2_v, single2_w, single2_wmb, 
			   single2_lfce, single2_fce, single2_mod,w,stack,lowend,highend,en3,jdangle+kdangle+ldangle+gap+2*data->eparam[5]+2*data->eparam[10]
		       +3*data->eparam[6]+constantclosure, iintercept, islope);

	    //case9 1 0 0 0
	    trace_branch_2(i+2,j-1,k+1,l-1,d,found,ct2,data,single2_v, single2_w, single2_wmb, 
			   single2_lfce, single2_fce, single2_mod,w,stack,lowend,highend,en3,idangle+gap+2*data->eparam[5]+2*data->eparam[10]
			   +data->eparam[6]+constantclosure, iintercept, islope);

	    //case10 1 0 0 1
	    trace_branch_2(i+2,j-1,k+1,l-2,d,found,ct2,data,single2_v, single2_w, single2_wmb, 
			   single2_lfce, single2_fce, single2_mod,w,stack,lowend,highend,en3,idangle+ldangle+2*gap+2*data->eparam[5]+2*data->eparam[10]
		      +2*data->eparam[6]+constantclosure, iintercept, islope);
	    
	    //case11 1 0 1 0
	    trace_branch_2(i+2,j-1,k+2,l-1,d,found,ct2,data,single2_v, single2_w, single2_wmb, 
			   single2_lfce, single2_fce, single2_mod,w,stack,lowend,highend,en3,idangle+kdangle+2*data->eparam[5]+2*data->eparam[10]
		      +2*data->eparam[6]+constantclosure, iintercept, islope);

	    //case12 1 0 1 1
	    trace_branch_2(i+2,j-1,k+2,l-2,d,found,ct2,data,single2_v, single2_w, single2_wmb, 
			   single2_lfce, single2_fce, single2_mod,w,stack,lowend,highend,en3,idangle+kdangle+ldangle+gap+2*data->eparam[5]+2*data->eparam[10]
		      +3*data->eparam[6]+constantclosure, iintercept, islope);

	    //case13 1 1 0 0
	    trace_branch_2(i+2,j-2,k+1,l-1,d,found,ct2,data,single2_v, single2_w, single2_wmb, 
			   single2_lfce, single2_fce, single2_mod,w,stack,lowend,highend,en3,idangle+jdangle+2*gap+2*data->eparam[5]+2*data->eparam[10]
		      +2*data->eparam[6]+constantclosure, iintercept, islope);
	    
	    //case14 1 1 0 1
	    trace_branch_2(i+2,j-2,k+1,l-2,d,found,ct2,data,single2_v, single2_w, single2_wmb, 
			   single2_lfce, single2_fce, single2_mod,w,stack,lowend,highend,en3,idangle+jdangle+ldangle+gap+2*data->eparam[5]+2*data->eparam[10]
			   +3*data->eparam[6]+constantclosure, iintercept, islope);

	    //case15 1 1 1 0
	    trace_branch_2(i+2,j-2,k+2,l-1,d,found,ct2,data,single2_v, single2_w, single2_wmb, 
			   single2_lfce, single2_fce, single2_mod,w,stack,lowend,highend,en3,idangle+jdangle+kdangle+gap+2*data->eparam[5]+2*data->eparam[10]
			   +3*data->eparam[6]+constantclosure, iintercept, islope);
	    
	    //case16 1 1 1 1
	    trace_branch_2(i+2,j-2,k+2,l-2,d,found,ct2,data,single2_v, single2_w, single2_wmb, 
			   single2_lfce, single2_fce, single2_mod,w,stack,lowend,highend,en3,idangle+jdangle+kdangle+ldangle+2*data->eparam[5]+2*data->eparam[10]
			   +4*data->eparam[6]+constantclosure, iintercept, islope);
	  }

	}
#else
#endif
        if (!found) {
			//a traceback error occurred
            //     cerr << "Traceback error in v!!\n";
	  //	cerr<<"closed!\n";
			error = 14;
          return error;
        }


      }




    }	
    else {//v[i][j][a][b]!=en3 so we need to search for w solutions
      found = false;
     
      //check if it is safe to increment i and k
      ikincrement =
        k+1 >= lowend[i+1]/*lowlimit(i+1, maxsep, ct1->GetSequenceLength(), ct2->GetSequenceLength())*/ &&
        k+1 <= highend[i+1]/*highlimit(i+1, maxsep, ct1->GetSequenceLength(), ct2->GetSequenceLength())*/;


      //check to see if it is safe to decrement j and l
      jldecrement =
        l-1 <= highend[j-1]/*highlimit(j-1, maxsep, ct1->GetSequenceLength(), ct2->GetSequenceLength())*/ &&
        l-1 >= lowend[j-1]/*lowlimit(j-1, maxsep, ct1->GetSequenceLength(), ct2->GetSequenceLength())*/;


      //case 6
      if (jldecrement) {
        if (en3==w->f(i,j-1,k,l-1)/*w[j-1][iref(i,j-1,ct1->GetSequenceLength())][a][b]*/+2*data->eparam[6]&&(j-1!=ct1->GetSequenceLength())&&(l-1!=ct2->GetSequenceLength())){
          found = true;
          stack.push(i,j-1,k/*a*/,l-1,w->f(i,j-1,k,l-1)/*w[j-1][iref(i,j-1,ct1->GetSequenceLength())][a][b]*/);
        }
      }
      //case 11
      if (!found&&ikincrement){
        if (en3==w->f(i+1,j,k+1,l)/*w[j][iref(i+1,j,ct1->GetSequenceLength())][a][b]*/+2*data->eparam[6]&&i!=ct1->GetSequenceLength()&&k!=ct2->GetSequenceLength()) {

          found = true;
          stack.push(i+1,j,k+1,l,w->f(i+1,j,k+1,l));
        }
      }
      //case 16
      if (!found&&ikincrement&&jldecrement) {
        if (j-1>i+1&&i!=ct1->GetSequenceLength()&&k!=ct2->GetSequenceLength()) {

          if (en3 ==w->f(i+1,j-1,k+1,l-1)/*w[j-1][iref(i+1,j-1,ct1->GetSequenceLength())][a][b]*/+ 4*data->eparam[6]) {
            found = true;
            stack.push(i+1,j-1,k+1,l-1,w->f(i+1,j-1,k+1,l-1));
          }
        }
      }



      if (/*a>=1*/k>=lowend[i+1]&&k<=highend[i+1]/*lowlimit(i+1,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/&&!found) {
        //case 9
        if(en3==w->f(i+1,j,k,l)/*w[j][iref(i+1,j,ct1->GetSequenceLength())][a-1][b]*/+data->eparam[6]+gap&&(i!=ct1->GetSequenceLength())){
          found = true;
          stack.push(i+1,j,k,l,w->f(i+1,j,k,l));

        }

        //case 14
        else if((j-1>i+1)&&(j-1!=ct1->GetSequenceLength())&&jldecrement) {

          if(en3==w->f(i+1,j-1,k,l-1)/*w[j-1][iref(i+1,j-1,ct1->GetSequenceLength())][a-1][b]*/+3*data->eparam[6]+gap){
            found = true;
            stack.push(i+1,j-1,k,l-1,w->f(i+1,j-1,k,l-1));
          }

        }

        if (/*b+1<(2*maxsep+2)*/l<=highend[j-1]&&l>=lowend[j-1]/*highlimit(j-1,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/&&!found&&(j-1>i+1)&&(j-1!=ct1->GetSequenceLength())) {
          //case 13
          if(en3==w->f(i+1,j-1,k,l)/*w[j-1][iref(i+1,j-1,ct1->GetSequenceLength())][a-1][b+1]*/+2*data->eparam[6]+2*gap){
            found = true;
            stack.push(i+1,j-1,k,l,w->f(i+1,j-1,k,l));

          }

        }
        if (/*b>=1*/l-1>=lowend[j]/*lowlimit(j,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/&&!found&&(i!=ct1->GetSequenceLength())&&(l-1!=ct2->GetSequenceLength())) {
          //case 10
          if(en3==w->f(i+1,j,k,l-1)/*w[j][iref(i+1,j,ct1->GetSequenceLength())][a-1][b-1]*/+2*data->eparam[6]+2*gap){
            found = true;
            stack.push(i+1,j,k,l-1,w->f(i+1,j,k,l-1));

          }

        }
      }
      if (/*b+1<(2*maxsep+2)*/l<=highend[j-1]&&l>=lowend[j-1]/*highlimit(j-1,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/&&!found&&(j-1!=ct1->GetSequenceLength())) {

        //case 5
        if(en3==w->f(i,j-1,k,l)/*w[j-1][iref(i,j-1,ct1->GetSequenceLength())][a][b+1]*/+data->eparam[6]+gap){
          found = true;
          stack.push(i,j-1,k,l,w->f(i,j-1,k,l));

        }



        else if (k+1<=highend[i]/*highlimit(i,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())/*a+1<(2*maxsep+2)*/&&(k!=ct2->GetSequenceLength())&&(j-1!=ct1->GetSequenceLength())) {
          //case 7
          if(en3==w->f(i,j-1,k+1,l)/*w[j-1][iref(i,j-1,ct1->GetSequenceLength())][a+1][b+1]*/+2*data->eparam[6]+2*gap){
            found = true;
            stack.push(i,j-1,k+1,l,w->f(i,j-1,k+1,l));

          }

        }

        //case 15
        else if((j-1>i+1)&&(k!=ct2->GetSequenceLength())&&(j-1!=ct1->GetSequenceLength())&&ikincrement) {
          if(en3==w->f(i+1,j-1,k+1,l)/*w[j-1][iref(i+1,j-1,ct1->GetSequenceLength())][a][b+1]*/+3*data->eparam[6]+gap) {
            found = true;
            stack.push(i+1,j-1,k+1,l,w->f(i+1,j-1,k+1,l));
          }

        }

      }
      if (/*a+1<(2*maxsep+2)*/k+1<=highend[i]/*highlimit(i,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/&&!found&&(k!=ct2->GetSequenceLength())) {
        //case 3
        if(en3==w->f(i,j,k+1,l)/*w[j][iref(i,j,ct1->GetSequenceLength())][a+1][b]*/+data->eparam[6]+gap){
          found = true;
          stack.push(i,j,k+1,l,w->f(i,j,k+1,l));

        }

        //case 8
        else if (jldecrement){
          if (en3==w->f(i,j-1,k+1,l-1)/*w[j-1][iref(i,j-1,ct1->GetSequenceLength())][a+1][b]*/+3*data->eparam[6]+gap&&(l-1!=ct2->GetSequenceLength())) {
            found = true;
            stack.push(i,j-1,k+1,l-1,w->f(i,j-1,k+1,l-1));
          }
        }

        if (!found&&/*b>=1*/l-1>=lowend[j]/*lowlimit(j,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/&&(l-1!=ct2->GetSequenceLength())) {
          //case 4
          if(en3==w->f(i,j,k+1,l-1)/*w[j][iref(i,j,ct1->GetSequenceLength())][a+1][b-1]*/+2*data->eparam[6]+2*gap){
            found = true;
            stack.push(i,j,k+1,l-1,w->f(i,j,k+1,l-1));

          }
        }
      }
      if (/*b>=1*/l-1>=lowend[j]/*lowlimit(j,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/&&!found&&(l-1!=ct2->GetSequenceLength())) {
        //case 2
        if(en3==w->f(i,j,k,l-1)/*w[j][iref(i,j,ct1->GetSequenceLength())][a][b-1]*/+data->eparam[6]+gap) {
          found = true;
          stack.push(i,j,k,l-1,w->f(i,j,k,l-1));
        }

        //case12
        else if (ikincrement){
          if ((en3==w->f(i+1,j,k+1,l-1)/*w[j][iref(i+1,j,ct1->GetSequenceLength())][a][b-1]*/+3*data->eparam[6]+gap&&i!=ct1->GetSequenceLength())&&(k!=ct2->GetSequenceLength())){
            found = true;
            stack.push(i+1,j,k+1,l-1,w->f(i+1,j,k+1,l-1));
          }
        }
      }
      //Consider the case where none of the four nucs (i,j,k,l) are paired
      //Consider whether any of the 4 nucleotides are stacked on a helix
      //There are 16 cases:
      //consider them in this order (0 unstacked, 1 stacked)
      //		i	j	k	l
      //1		0	0	0	0
      //2		0	0	0	1
      //3		0	0	1	0
      //4		0	0	1	1
      //5		0	1	0	0
      //6		0	1	0	1
      //7		0	1	1	0
      //8		0	1	1	1
      //9		1	0	0	0
      //10	1	0	0	1
      //11	1	0	1	0
      //12	1	0	1	1
      //13	1	1	0	0
      //14	1	1	0	1
      //15	1	1	1	0
      //16	1	1	1	1


      //case 1 - nothing stacked:
      if(en3==v->f(i,j,k,l)/*v[j][iref(i,j,ct1->GetSequenceLength())][a][b]*/+2*data->eparam[10]+penalty(i,j,ct1,data)+penalty(k,l,ct2,data)&&!found) {
        found = true;
        stack.push(i,j,k,l,v->f(i,j,k,l));

      }



      //case 6 - j and l stacked
      if (!found&&jldecrement) {
        if(en3==v->f(i,j-1,k,l-1)/*v[j-1][iref(i,j-1,ct1->GetSequenceLength())][a][b]*/+2*data->eparam[6]+edangle3(j-1,i,j,ct1,data)
           +edangle3(l-1,k,l,ct2,data)+2*data->eparam[10]+penalty(i,j-1,ct1,data)+penalty(l-1,k,ct2,data)
           &&(j-1!=ct1->GetSequenceLength())&&(l!=ct2->GetSequenceLength()+1)) {

          found = true;
          stack.push(i,j-1,k,l-1,v->f(i,j-1,k,l-1));
        }
      }

      //case 11 - i and k stacked
      if(!found&&ikincrement) {
        if ((i!=ct1->GetSequenceLength())&&(k!=ct2->GetSequenceLength())&&en3==v->f(i+1,j,k+1,l)/*v[j][iref(i+1,j,ct1->GetSequenceLength())][a][b]*/+2*data->eparam[6]+edangle5(k+1,l,k,ct2,data)+
            edangle5(i+1,j,i,ct1,data)+2*data->eparam[10]+penalty(i+1,j,ct1,data)+penalty(k+1,l,ct2,data)) {

          found = true;
          stack.push(i+1,j,k+1,l,v->f(i+1,j,k+1,l));
        }

      }

      //case 16 - i, j, k, and l stacked
      if(!found&&(j-1>i+1)&&(i!=ct1->GetSequenceLength())&&(k!=ct2->GetSequenceLength())&&(j-1!=ct1->GetSequenceLength())&&(l!=ct2->GetSequenceLength()+1)&&ikincrement&&jldecrement) {
        if (en3==v->f(i+1,j-1,k+1,l-1)/*v[j-1][iref(i+1,j-1,ct1->GetSequenceLength())][a][b]*/+4*data->eparam[6]+edangle3(l-1,k+1,l,ct2,data)
            +edangle5(k+1,l-1,k,ct2,data)
            +edangle3(j-1,i+1,j,ct1,data)
            +edangle5(i+1,j-1,i,ct1,data)+2*data->eparam[10]+
            penalty(i+1,j-1,ct1,data)+penalty(k+1,l-1,ct2,data)) {


          found = true;
          stack.push(i+1,j-1,k+1,l-1,v->f(i+1,j-1,k+1,l-1));

        }


      }



      if (!found&&/*b-1>=0*/l-1>=lowend[j]/*lowlimit(j,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/&&(l!=ct2->GetSequenceLength()+1)) {
        //case 2 - l stacked
        if(en3==v->f(i,j,k,l-1)/*v[j][iref(i,j,ct1->GetSequenceLength())][a][b-1]*/+data->eparam[6]+edangle3(l-1,k,l,ct2,data)+2*data->eparam[10]+gap+
           penalty(i,j,ct1,data)+penalty(l-1,k,ct2,data)) {


          found = true;
          stack.push(i,j,k,l-1,v->f(i,j,k,l-1));

        }

        else if (/*(a+1<2*maxsep+2)*/k+1<=highend[i]/*highlimit(i,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/&&(k!=ct2->GetSequenceLength())) {
          //case 4 - l and k stacked
          if(en3==v->f(i,j,k+1,l-1)/*v[j][iref(i,j,ct1->GetSequenceLength())][a+1][b-1]*/+2*data->eparam[6]+edangle3(l-1,k+1,l,ct2,data)
             +edangle5(k+1,l-1,k,ct2,data)+2*data->eparam[10]+2*gap
             +penalty(i,j,ct1,data)+penalty(l-1,k+1,ct2,data)) {

            found = true;
            stack.push(i,j,k+1,l-1,v->f(i,j,k+1,l-1));

          }
        }




      }

      if (!found&&k+1<=highend[i]/*highlimit(i,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())/*a+1<2*maxsep+2*/&&(k!=ct2->GetSequenceLength())) {

        //case 8 - j, k, and l stacked:
        if (jldecrement) {
          if(en3 ==v->f(i,j-1,k+1,l-1)/*v[j-1][iref(i,j-1,ct1->GetSequenceLength())][a+1][b]*/+3*data->eparam[6]+edangle3(j-1,i,j,ct1,data)
             +edangle3(l-1,k+1,l,ct2,data)+edangle5(k+1,l-1,k,ct2,data)+2*data->eparam[10]+gap
             +penalty(i,j-1,ct1,data)+penalty(k+1,l-1,ct2,data)&&(l!=ct2->GetSequenceLength()+1)&&(j-1!=ct1->GetSequenceLength())) {


            found = true;
            stack.push(i,j-1,k+1,l-1,v->f(i,j-1,k+1,l-1));
          }
        }

        //case 3 - k stacked
        if(!found&&en3==v->f(i,j,k+1,l)/*v[j][iref(i,j,ct1->GetSequenceLength())][a+1][b]*/+data->eparam[6]+edangle5(k+1,l,k,ct2,data)+2*data->eparam[10]+gap+
           penalty(i,j,ct1,data)+penalty(k+1,l,ct2,data)) {


          found = true;
          stack.push(i,j,k+1,l,v->f(i,j,k+1,l));

        }

        if (!found&&l<=highend[j-1]&&l>=lowend[j-1]/*highlimit(j-1,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())/*b+1<2*maxsep+2*/&&(j!=ct1->GetSequenceLength()+1)) {
          //case 7 - j and k stacked
          if(en3==v->f(i,j-1,k+1,l)/*v[j-1][iref(i,j-1,ct1->GetSequenceLength())][a+1][b+1]*/+2*data->eparam[6]+edangle3(j-1,i,j,ct1,data)
             +edangle5(k+1,l,k,ct2,data)+2*data->eparam[10]+2*gap+
             penalty(i,j-1,ct1,data)+penalty(k+1,l,ct2,data)) {


            found = true;
            stack.push(i,j-1,k+1,l,v->f(i,j-1,k+1,l));

          }

        }


      }

      if (l<=highend[j-1]&&l>=lowend[j-1]/*highlimit(j-1,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())/*b+1<2*maxsep+2*/&&!found&&(j!=ct1->GetSequenceLength()+1)) {

        //case 5 - j stacked
        if(en3==v->f(i,j-1,k,l)/*v[j-1][iref(i,j-1,ct1->GetSequenceLength())][a][b+1]*/+data->eparam[6]+edangle3(j-1,i,j,ct1,data)+2*data->eparam[10]+gap+
           penalty(i,j-1,ct1,data)+penalty(k,l,ct2,data)) {


          found = true;
          stack.push(i,j-1,k,l,v->f(i,j-1,k,l));

        }

        //case 15 - i,j, and k stacked
        else if (ikincrement) {
          if(en3==v->f(i+1,j-1,k+1,l)/*v[j-1][iref(i+1,j-1,ct1->GetSequenceLength())][a][b+1]*/+3*data->eparam[6]+edangle5(k+1,l,k,ct2,data)
             +edangle3(j-1,i+1,j,ct1,data)
             +edangle5(i+1,j-1,i,ct1,data)+2*data->eparam[10]+gap+
             penalty(i+1,j-1,ct1,data)+penalty(k+1,l,ct2,data)&&(i!=ct1->GetSequenceLength())) {

            found = true;
            stack.push(i+1,j-1,k+1,l,v->f(i+1,j-1,k+1,l));
          }

        }

        if (!found&&/*a>=1*/k>=lowend[i+1]&&k<=highend[i+1]/*lowlimit(i+1,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/&&(j-1>i+1)&&(i!=ct1->GetSequenceLength())) {
          //case 13 - i and j stacked

          if(en3==v->f(i+1,j-1,k,l)/*v[j-1][iref(i+1,j-1,ct1->GetSequenceLength())][a-1][b+1]*/+2*data->eparam[6]+edangle3(j-1,i+1,j,ct1,data)
             +edangle5(i+1,j-1,i,ct1,data)+2*data->eparam[10]+2*gap+
             penalty(i+1,j-1,ct1,data)+penalty(k,l,ct2,data)) {


            found = true;
            stack.push(i+1,j-1,k,l,v->f(i+1,j-1,k,l));

          }
        }

      }

      if (!found&&k>=lowend[i+1]&&k<=highend[i+1]/*lowlimit(i+1,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())/*a>=1*/&&(i!=ct1->GetSequenceLength())) {
        //case 9 - i alone is stacked
        if(en3==v->f(i+1,j,k,l)/*v[j][iref(i+1,j,ct1->GetSequenceLength())][a-1][b]*/+data->eparam[6]+edangle5(i+1,j,i,ct1,data)+2*data->eparam[10]+gap+
           penalty(i+1,j,ct1,data)+penalty(k,l,ct2,data)) {
          found = true;
          stack.push(i+1,j,k,l,v->f(i+1,j,k,l));


        }


        //case 14 - i, j, and l stacked
        else if(j-1>i+1&&(j!=ct1->GetSequenceLength()+1)&&(l!=ct2->GetSequenceLength()+1)&&jldecrement) {
          if(en3==v->f(i+1,j-1,k,l-1)/*v[j-1][iref(i+1,j-1,ct1->GetSequenceLength())][a-1][b]*/+3*data->eparam[6]+edangle3(l-1,k,l,ct2,data)
             +edangle3(j-1,i+1,j,ct1,data)
             +edangle5(i+1,j-1,i,ct1,data)+2*data->eparam[10]+gap+
             penalty(i+1,j-1,ct1,data)+penalty(k,l-1,ct2,data)) {

            found = true;
            stack.push(i+1,j-1,k,l-1,v->f(i+1,j-1,k,l-1));

          }


        }
      }


      if (!found&&l-1>=lowend[j]/*lowlimit(j,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())/*b-1>=0*/&&(i!=ct1->GetSequenceLength())&&(l!=ct2->GetSequenceLength()+1)) {
        //case 12 - i, k, and l stacked:
        if (ikincrement) {
          if(en3==v->f(i+1,j,k+1,l-1)/*v[j][iref(i+1,j,ct1->GetSequenceLength())][a][b-1]*/+3*data->eparam[6]+edangle3(l-1,k+1,l,ct2,data)
             +edangle5(i+1,j,i,ct1,data)+edangle5(k+1,l-1,k,ct2,data)+2*data->eparam[10]+gap+
             penalty(i+1,j,ct1,data)+penalty(l-1,k+1,ct2,data)&&(k!=ct2->GetSequenceLength())) {


            found = true;

			stack.push(i+1,j,k+1,l-1,v->f(i+1,j,k+1,l-1));

          }
        }
        if (!found&&k>=lowend[i+1]&&k<=highend[i+1]/*lowlimit(i+1,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())/*a>=1*/) {
          //case 10 - l and i stacked
          if(en3==v->f(i+1,j,k,l-1)/*v[j][iref(i+1,j,ct1->GetSequenceLength())][a-1][b-1]*/+2*data->eparam[6]+edangle3(l-1,k,l,ct2,data)
             +edangle5(i+1,j,i,ct1,data)+2*data->eparam[10]+2*gap+
             penalty(i+1,j,ct1,data)+penalty(k,l-1,ct2,data)) {

            found = true;
            stack.push(i+1,j,k,l-1,v->f(i+1,j,k,l-1));

          }
        }
      }

      //calculate the free energy of 2 fragments merged:
      for (c=i+minloop;c<j-minloop&&!found;c++) {
        //if (c>ct1->GetSequenceLength()) kp = ct1->GetSequenceLength()-ct2->GetSequenceLength();
        //else kp = 0;
        for (d=max(k+minloop,/*c-maxsep-kp*/lowend[c]/*lowlimit(c,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/);
             d<l-minloop&&d<=/*c+maxsep-kp*/highend[c]/*highlimit(c,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/&&!found;d++) {

          //e = d-c+maxsep+kp;
          if ((c!=ct1->GetSequenceLength())&&(d!=ct2->GetSequenceLength())&&d+1>=lowend[c+1]&&d+1<=highend[c+1]) {
            if (c<ct1->GetSequenceLength()) {
              if (d<ct2->GetSequenceLength()) {
				  if (d+1>=lowend[c+1]&&d+1<=highend[c+1]) {
					if(en3 ==w->f(i,c,k,d)/*w[c][i][a][e]*/+w->f(c+1,j,d+1,l)/*w[j][iref(c+1,j,ct1->GetSequenceLength())][e][b]*/) {
						found = true;
						stack.push(i,c,k,d,w->f(i,c,k,d));
						stack.push(c+1,j,d+1,l,w->f(c+1,j,d+1,l));

					}
				  }

              }
            }
            else {
              if (d>ct2->GetSequenceLength()) {
                if(en3 ==w->f(i,c,k,d)/*w[c][iref(i,c,ct1->GetSequenceLength())][a][e]*/+
                   w->f(c+1-ct1->GetSequenceLength(),j-ct1->GetSequenceLength(),d+1-ct2->GetSequenceLength(),l-ct2->GetSequenceLength())
                   /*w[j-ct1->GetSequenceLength()][iref(c+1-ct1->GetSequenceLength(),j-ct1->GetSequenceLength(),ct1->GetSequenceLength())][e][b]*/) {
                  found = true;
                  stack.push(i,c,k,d,w->f(i,c,k,d));
                  stack.push(c+1-ct1->GetSequenceLength(),j-ct1->GetSequenceLength(),d+1-ct2->GetSequenceLength(),l-ct2->GetSequenceLength(),w->f(c+1-ct1->GetSequenceLength(),j-ct1->GetSequenceLength(),d+1-ct2->GetSequenceLength(),l-ct2->GetSequenceLength()));
		  //report to Dave6
                }
              }
            }
          }
        }
      }
#ifdef DYNALIGN_II

      //report to Dave5
      //      if(i==28&&j==45&&k==28&&l==60)cerr<<"w->f(i,j,k,l) in traceback"<<en3<<"\n";
      for (c=i+minloop;!found&&c+1<j-minloop&&c!=ct1->GetSequenceLength();++c) {
    //  if((c<N)&&k>=lowend[c+1]&&k<=highend[c+1])en1=min(en1,w->f(c+1,j,k,l)+single1_w.f(i,c)+(c-i+1)*gap);
	trace_branch_1(i,j,k,l,c,found,ct1,data,single1_v, single1_w, single1_wmb, 
		       single1_lfce, single1_fce, single1_mod,w,stack,lowend,highend,en3,0, iintercept, islope);
 
      }

      for (d=k+minloop;!found&&d+1<l-minloop&&d!=ct2->GetSequenceLength();++d) {
	trace_branch_2(i,j,k,l,d,found,ct2,data,single2_v, single2_w, single2_wmb, 
		       single2_lfce, single2_fce, single2_mod,w,stack,lowend,highend,en3,0, iintercept, islope);
 
      }

#else
#endif


      if (!found) {
		  //a traceback error occurred
//	   cerr << "Traceback error in w!\n";
		  error = 14;
        return error;
      }
    }
  }
  //end of while (stack.pull(&i,&j, &k, &l /*&a,&b*/, &en3, &open))
  return error;
}


//Perform multiple trackbacks.
//Return an int that indicates whether an error occurred (0=no error).
#ifdef DYNALIGN_II
int  dyntraceback(short maxtracebacks, short window, short awindow,
                  short percentsort,
                  varray *v, dynalignarray *w, wendarray *w3, wendarray *w5,
                  /*bool **pair,*/DynProgArray<integersize> *single1_w,DynProgArray<integersize> *single2_w,DynProgArray<integersize> *single1_v,DynProgArray<integersize> *single2_v,DynProgArray<integersize> *single1_wmb,DynProgArray<integersize> *single2_wmb,DynProgArray<integersize> *single1_we,DynProgArray<integersize> *single2_we,bool *single1_lfce,bool *single2_lfce,bool *single1_mod,bool *single2_mod,forceclass *single1_fce,forceclass *single2_fce,
                  structure *ct1, structure *ct2, short **alignment,
                  short *lowend, short *highend, short int islope, short int iintercept, short int gapincrease,
                  datatable *data,
                  integersize lowest, dynalignarray *vmod,
                  bool local,int max_elongation,bool *mod1, bool *mod2, bool modification,
                  char **fce1, char **fce2, bool alignmentforced, short **forcealign,bool force) 
#else
int  dyntraceback(short maxtracebacks, short window, short awindow,
                  short percentsort,
                  varray *v, dynalignarray *w, wendarray *w3, wendarray *w5,
                  /*bool **pair,*/
                  structure *ct1, structure *ct2, short **alignment,
                  short *lowend, short *highend, short int gapincrease,
                  datatable *data, bool singleinsert,
                  integersize lowest, dynalignarray *vmod,
                  bool local)
#endif
{

  short int i,j,k,l,gap,N,N2,index;
  int ip,jp,kp,cp;
  integersize en1;
  integersize crit;
  dynalignheap heap;
  int point;
  bool **mark1, **mark2, **marka;
  //  bool modification;
  int error = 0;
  int localerror;

#ifndef DYNALIGN_II
  bool modification;
   if (ct1->GetNumberofModified()>0||ct2->GetNumberofModified()>0) modification = true;
   else modification = false;
#else
#endif
  vector< vector<bool> > inc = data->pairing;

  //declarations for traceback portion:


  gap = gapincrease;
  N = ct1->GetSequenceLength();
  N2 = ct2->GetSequenceLength();

  if ((local&&lowest==0)||(!local&&(lowest==gap*abs(ct1->GetSequenceLength()-ct2->GetSequenceLength())))) {

	//lowest is zero because there is no structure, the ends are single-stranded
	//In this case, revise lowest using the varrays to determine optimal folding (>0 DG)
	lowest = DYNALIGN_INFINITY;
	for (i=1;i<=N;i++) {
		for (j=i+minloop;j<=N;j++) {
		  for (k=max(/*i-maxsep*/lowend[i]/*lowlimit(i,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/,1);k<=(min(ct2->GetSequenceLength(),/*i+maxsep*/highend[i]/*highlimit(i,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/));k++) {
			for (l=max(/*j-maxsep*/lowend[j]/*lowlimit(j,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/,k);l<=(min(ct2->GetSequenceLength(),/*j+maxsep*/highend[j]/*highlimit(j,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/));l++) {



			  if (modification) {
				if (vmod->f(i,j,k,l)/*vmod[j][i][a][b]*/+vmod->f(j,i+N,l,k+N2)/*vmod[i+N][j-i][b][a]*/<=lowest)
				  lowest = vmod->f(i,j,k,l)+vmod->f(j,i+N,l,k+N2);
			  }
			  else {
				  if (v->f(i,j,k,l)+v->f(j,i+N,l,k+N2)<=lowest)
						lowest =v->f(i,j,k,l)+v->f(j,i+N,l,k+N2);
			  }

			}//loop over l
		  }//loop over k
		}//loop over j
	  }//loop over i


  }

  if (lowest>=DYNALIGN_INFINITY) {
	//There is no way to get a common structure at all

    //Set the energies to DYNALIGN_INFINITY
	  ct1->SetEnergy(1,DYNALIGN_INFINITY);
	  ct2->SetEnergy(1,DYNALIGN_INFINITY);
	  
    //return with empty CT files
   ct1->AddStructure();
    ct2->AddStructure();


    //initialize the alignment:
    for (index=0;index<=ct1->GetSequenceLength();index++) alignment[0][index]=0;
    return 0;

  }



  mark1 = new bool *[N+1];
  for (i=0;i<=N;i++) mark1[i] = new bool [i];
  mark2 = new bool *[N2+1];
  for (i=0;i<=N2;i++) mark2[i] = new bool [i];

  marka = new bool *[N+1];
  for (i=0;i<=N;i++) marka[i] = new bool [highend[i]-lowend[i]+1/*2*maxsep+2*/];

  for (i=0;i<=N;i++) {
    for (j=0;j<i;j++) {
      mark1[i][j] = false;
    }
  }
  for (i=0;i<=N2;i++) {
    for (j=0;j<i;j++) {
      mark2[i][j] = false;
    }
  }
  for (i=0;i<=N;i++) {
    for (j=0;j<highend[i]-lowend[i]+1/*2*maxsep+2*/;j++) {
      marka[i][j] = false;
    }
	//change the pointer for marka for faster indexing:
	marka[i]-=lowend[i];
  }





  //now traceback:

  //For suboptimal structure prediction, build a heap of possible traceback starts

  crit = lowest;
  if (percentsort>0) crit = crit + abs((integersize)((float) crit*(((float) percentsort)/100.0 )));
  else crit = crit - percentsort;


  for (i=1;i<=N;i++) {
    for (j=i+minloop;j<=N;j++) {
      for (k=max(/*i-maxsep*/lowend[i]/*lowlimit(i,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/,1);k<=(min(ct2->GetSequenceLength(),/*i+maxsep*/highend[i]/*highlimit(i,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/));k++) {
        for (l=max(/*j-maxsep*/lowend[j]/*lowlimit(j,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/,k);l<=(min(ct2->GetSequenceLength(),/*j+maxsep*/highend[j]/*highlimit(j,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/));l++) {

          //use a and b when refering to the energy arrays
          //a = k-i+maxsep;
          //b = l-j+maxsep;

          //if (/*(a+N2-N>=0)&&(a+N2-N<=2*maxsep)*/) {
          if (modification) {
            if (vmod->f(i,j,k,l)/*vmod[j][i][a][b]*/+vmod->f(j,i+N,l,k+N2)/*vmod[i+N][j-i][b][a]*/<=crit)
              heap.push(i,j,k,l,vmod->f(i,j,k,l)+vmod->f(j,i+N,l,k+N2));
          }
          else if (v->f(i,j,k,l)/*v[j][i][a][b]*/+v->f(j,i+N,l,k+N2)/*v[i+N][j-i][b][a]*/<=crit)
            heap.push(i,j,k,l,v->f(i,j,k,l)+v->f(j,i+N,l,k+N2));
          //}

        }
      }
    }

  }

  //make a heap out of dynalign heap:
  for (ip=1;ip<heap.size;ip++) {
    jp = ip;
    kp = jp/2;
    while ((heap.peak(jp)<heap.peak(kp))&&kp>=0){
      heap.swap(jp,kp);
      jp = jp/2;
      kp = kp/2;
    }

  }
  //sort the dynalign heap with heap sort
  for (ip = heap.size-2;ip>=0;ip--) {
    heap.swap(ip+1,0);

    jp = 0;
    cp = 1;
    while (cp<=ip) {
      if (cp!=ip) {
        if (heap.peak(cp+1)<heap.peak(cp)) cp++;

      }
      if (heap.peak(cp)<heap.peak(jp)) {
        heap.swap(cp,jp);
        jp = cp;
        cp = 2*cp;
      }
      else cp = ip + 1;
    }

  }


  //call dyntrace to determine suboptimal structures:

  point = heap.size;
 
  while (ct1->GetNumberofStructures()<maxtracebacks&&point>0) {

    point--;


    heap.read(point,&i,&j,&k,&l,&en1);
    //k = a+i-maxsep;
    //l = b+j-maxsep;
    // cerr<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<en1<<" heap trace!\n";
    if (!mark1[j][i]||!mark2[l][k]||!marka[i][k/*reference(i, k, N,N2,maxsep)/*a*/]||!marka[j][l/*reference(j, l, N,N2,maxsep)/*b*/]) {
    	//Need to add the next structure for traceback
		ct1->AddStructure();
		ct2->AddStructure();

		ct1->SetCtLabel(ct1->GetSequenceLabel(),ct1->GetNumberofStructures());
		ct2->SetCtLabel(ct2->GetSequenceLabel(),ct2->GetNumberofStructures());
		
	
      //initialize the alignment:
      for (index=0;index<=ct1->GetSequenceLength();index++) alignment[ct1->GetNumberofStructures()-1][index]=0;

      //dyntrace(57,126,61,131,ct1,ct2,ct1->GetNumberofStructures(),alignment[ct1->GetNumberofStructures()-1],w,v,w3,w5,maxsep,data,gap,vmod);
//       cerr<<"start!\n";

#ifdef DYNALIGN_II
      localerror = dyntrace(i,j,k,l,
			    ct1,ct2,ct1->GetNumberofStructures(),
			    alignment[ct1->GetNumberofStructures()-1],
			    w,v,w3,w5,lowend,highend,data,islope,iintercept,gap,vmod,single1_w,single2_w,single1_v,single2_v,single1_wmb,single2_wmb,single1_we,single2_we,single1_lfce,single2_lfce,single1_mod,single2_mod,single1_fce,single2_fce,
			    local,max_elongation,mod1, mod2, modification,
                         fce1, fce2,
                         alignmentforced, forcealign,force);
#else
         localerror = dyntrace(i,j,k,l,
               ct1,ct2,ct1->GetNumberofStructures(),
               alignment[ct1->GetNumberofStructures()-1],
               w,v,w3,w5,lowend,highend,data,gap,vmod,
               local);
#endif
//      cerr<<"over!\n";
      if (localerror!=0)error = localerror;
#ifdef DYNALIGN_II
      localerror = dyntrace(j,i+N,l,k+N2,
			    ct1,ct2,ct1->GetNumberofStructures(),
			    alignment[ct1->GetNumberofStructures()-1],
			    w,v,w3,w5,lowend,highend,data,islope,iintercept,gap,vmod,single1_w,single2_w,single1_v,single2_v,single1_wmb,single2_wmb,single1_we,single2_we,single1_lfce,single2_lfce,single1_mod,single2_mod,single1_fce,single2_fce,
			    local,max_elongation,mod1, mod2, modification,
                         fce1, fce2,
                         alignmentforced, forcealign,force);
#else
      localerror = dyntrace(j,i+N,l,k+N2,
                            ct1,ct2,ct1->GetNumberofStructures(),
                            alignment[ct1->GetNumberofStructures()-1],
                            w,v,w3,w5,lowend,highend,data,gap,vmod,
                            local);
#endif
//      cerr<<"really over!\n";
      if (localerror!=0) error = localerror;
       ct1->SetEnergy(ct1->GetNumberofStructures(),en1);
	  ct2->SetEnergy(ct2->GetNumberofStructures(),en1);
      




      for (i=1;i<N;i++) {
        if (ct1->GetPair(i,ct1->GetNumberofStructures())>i) {
          if (!mark1[ct1->GetPair(i,ct1->GetNumberofStructures())][i]) {

            for (j=i-window;j<=i+window;j++) {
              for (k=ct1->GetPair(i,ct1->GetNumberofStructures())-window;
                   k<=ct1->GetPair(i,ct1->GetNumberofStructures())+window;k++) {

                if (j<k&&j>0&&k<=ct1->GetSequenceLength()) {
                  mark1[k][j] = true;
                }
              }
            }
          }
        }
      }

      for (i=1;i<N2;i++) {
        if (ct2->GetPair(i,ct2->GetNumberofStructures())>i) {
          if (!mark2[ct2->GetPair(i,ct2->GetNumberofStructures())][i]) {

            for (j=i-window;j<=i+window;j++) {
              for (k=ct2->GetPair(i,ct2->GetNumberofStructures())-window;
                   k<=ct2->GetPair(i,ct2->GetNumberofStructures())+window;k++) {

                if (j<k&&j>0&&k<=ct2->GetSequenceLength()) {
                  mark2[k][j] = true;
                }
              }
            }
          }
        }
      }

      //cerr<<"2 really really over!\n";
      for (i=1;i<=N;i++) {
        if (alignment[ct1->GetNumberofStructures()-1][i]) {
          if (!marka[i][alignment[ct1->GetNumberofStructures()-1][i]]) {/*[reference(i, alignment[ct1->GetNumberofStructures()-1][i], N,N2,maxsep)]/*alignment[ct1->GetNumberofStructures()-1][i]-i+maxsep]*/

            for (j=i-awindow;j<=i+awindow;j++) {
			  for (k=alignment[ct1->GetNumberofStructures()-1][i]-awindow;k<=alignment[ct1->GetNumberofStructures()-1][i]+awindow;k++) {
			  //for (k=reference(i, alignment[ct1->GetNumberofStructures()-1][i], N,N2,maxsep)/*alignment[ct1->GetNumberofStructures()-1][i]-i+maxsep*/-awindow;
              //     k<=reference(i, alignment[ct1->GetNumberofStructures()-1][i], N,N2,maxsep)/*alignment[ct1->GetNumberofStructures()-1][i]-i+maxsep*/+awindow;
              //     k++) {

                if (j>0&&j<=ct1->GetSequenceLength()&&k>=lowend[j]&&k<=highend[j]) {
                  marka[j][k] = true;
                }
              }
            }
          }
        }
      }
      // cerr<<"1 really really over!\n";
    }
   
  }
  //cerr<<"3 really really over!\n";
  for (i=0;i<=N;i++) delete[] mark1[i];
  delete[] mark1;
  for (i=0;i<=N2;i++) delete[] mark2[i];;
  delete[] mark2;
  for (i=0;i<=N;i++) {
	  //move back the pointer:
	  marka[i]+=lowend[i];
	  delete[] marka[i];
  }
  delete[] marka;

  //return the error code
  return error;
}

void alignout(short** align, const char *aout, structure *ct1, structure *ct2) {
  ofstream out;
  short i,j,last,next,index;
  char *line1,*line2,*line3;

  line1 = new char [ct1->GetSequenceLength()+ct2->GetSequenceLength()+100];
  line2 = new char [ct1->GetSequenceLength()+ct2->GetSequenceLength()+100];
  line3 = new char [ct1->GetSequenceLength()+ct2->GetSequenceLength()+100];

  out.open(aout);
  //   cerr<<"GetNumberofStructures() "<<ct1->GetNumberofStructures()<<"\n";
  for (index=0;index<ct1->GetNumberofStructures();index++) {
    strcpy(line1,"");
    strcpy(line2,"");
    strcpy(line3,"");
    //   cerr<<"in the GetNumberofStructures() cycle! "<<index<<"\n";
    //   for (i=1;i<=ct1->GetSequenceLength();i++) {
    //     cerr<<i<<" "<<align[index][i]<<"\n";
    //   }
    //    cerr<<"\n";
    last = 0;
    for (i=1;i<=ct1->GetSequenceLength();i++) {
      //  cerr<<"in the nucleotide cycle! "<<i<<"\n";
      if (last==ct2->GetSequenceLength()) {
	
        //nothing more can be put down in line 2 for the second sequence
	
        line1[strlen(line1)+1]='\0';
        line1[strlen(line1)]=ct1->nucs[i];
        strcat(line2,"-");
        strcat(line3," ");

      }
      else if (align[index][i]>0) {
        //this nucleotide is aligned to something
        while (align[index][i]!=last+1) {
          //need to lay nucleotides down in second sequence
          strcat(line1,"-");
          last++;
          line2[strlen(line2)+1]='\0';
          line2[strlen(line2)]=ct2->nucs[last];
          strcat(line3," ");

        }
        //now put i in alignment with last+1
        line1[strlen(line1)+1]='\0';
        line1[strlen(line1)]=ct1->nucs[i];
        last++;
        line2[strlen(line2)+1]='\0';
        line2[strlen(line2)]=ct2->nucs[last];
        strcat(line3,"^");

      }
      else {//align[i]=0
        next=0;
        for (j=i+1;j<=ct1->GetSequenceLength()&&next==0;j++) next = align[index][j];
        if (next==last+1) {
          //no nucs from seq2 can be put down next to seq1
          line1[strlen(line1)+1]='\0';
          line1[strlen(line1)]=ct1->nucs[i];
          strcat(line2,"-");
          strcat(line3," ");

        }
        else {
          line1[strlen(line1)+1]='\0';
          line1[strlen(line1)]=ct1->nucs[i];
          last++;
          line2[strlen(line2)+1]='\0';
          line2[strlen(line2)]=ct2->nucs[last];
          strcat(line3," ");

        }


      }


    }

    //put down any nucs from seq2 that did not go down:
    for (i=last+1;i<=ct2->GetSequenceLength();i++) {
      strcat(line1,"-");
      line2[strlen(line2)+1]='\0';
      line2[strlen(line2)]=ct2->nucs[i];
      strcat(line3," ");
    }

    out<<"Alignment #"<<index+1<<" Score= "<<ct1->GetEnergy(index+1)<<"\n";
    out <<line1<<"\n"<<line2<<"\n"<<line3<<"\n\n\n";
  }

  out.close();
  delete[] line1;
  delete[] line2;
  delete[] line3;
}

void parse(structure *ct, char *seq1, char *seq2, datatable* data) {
  short int *seqone, *seqtwo;
  vector< vector<bool> > inc = data->pairing;
  short int i,j;
  //char temp;

  seqone = new short int [ct->GetSequenceLength()+1];
  seqtwo = new short int [ct->GetSequenceLength()+1];

  j = 1;
  for (i=0;i<(short) (strlen(seq1));i++) {
   
    if (seq1[i]!='-') {
      seqone[j] = ct->GetThermodynamicDataTable()->basetonum(seq1[i]);
      
      seqtwo[j] = ct->GetThermodynamicDataTable()->basetonum(seq2[i]);
      j++;
    }
  }

  for (i=1;i<=ct->GetSequenceLength();i++) {
    for (j=i+1;j<=ct->GetSequenceLength();j++) {
      if (inc[seqone[i]][seqone[j]]&&inc[seqtwo[i]][seqtwo[j]]) {
        ct->tem[j][i] = true;
      }
else {
        ct->tem[j][i] = false;
      }
    }
  }
  delete[] seqone;
  delete[] seqtwo;
}


#ifdef DYNALIGN_II
void opendynalignsavefile(const char *filename, structure *ct1, structure *ct2,
                        varray *v, dynalignarray *w, dynalignarray *vmod,
                          wendarray *w3, wendarray *w5, DynProgArray<integersize> *single1_w,DynProgArray<integersize> *single2_w,DynProgArray<integersize> *single1_v,DynProgArray<integersize> *single2_v,DynProgArray<integersize> *single1_wmb,DynProgArray<integersize> *single2_wmb,DynProgArray<integersize> *single1_we,DynProgArray<integersize> *single2_we,bool *single1_lfce,bool *single2_lfce,bool *single1_mod,bool *single2_mod,forceclass *single1_fce,forceclass *single2_fce,int *single1_vmin,int *single2_vmin,datatable *data,int *max_elongation,
			  short *maxsep, short *islope, short *iintercept, short *gap,
                          short *lowest, bool *local, bool **allowed_alignments,
						  short *lowend, short *highend)
#else
void opendynalignsavefile(const char *filename, structure *ct1, structure *ct2,
                        varray *v, dynalignarray *w, dynalignarray *vmod,
                          wendarray *w3, wendarray *w5, datatable *data,
                          bool *singleinsert, short *maxsep, short *gap,
                          short *lowest, bool *local, bool **allowed_alignments,
						  short *lowend, short *highend)
#endif
 {
  short i,j,k,l,a,b,c,d,I;
  int temp;
  bool optimalonly,modification;

  ifstream sav(filename,ios::binary);
  int versioncheck;
  read(&sav, &versioncheck);
 

  //start with structure information
  read(&sav,&temp);//read the modification/optimalonly key - this has been noted previously
  //because it is needed for array allocation
  //The flag returns:
  //0 - suboptimal, no modifications
  //1 - suboptimal, with modifications
  //2 - optimal only, no modifications
  //3 - optimal only, with modifications
  if (temp==1||temp==3) modification = true;
  else modification = false;

  int localint;
  read(&sav,&(localint));
  ct1->allocate(localint);
  read(&sav,&(localint));
  ct2->allocate(localint);
  read(&sav,maxsep);
  


  //Read the structure-class specific items:
  openstructuresave(&sav,ct1);
  openstructuresave(&sav,ct2);

  //read datatable
  read(&sav, data);

  //now read the folding data:
  optimalonly = temp==2 || temp==3;


  read(&sav,gap);
  read(&sav,lowest);
#ifdef DYNALIGN_II
  read(&sav,max_elongation);
#else
  read(&sav,singleinsert);
#endif
 
#ifdef DYNALIGN_II
    read(&sav,islope);
    read(&sav,iintercept);
#else
#endif

	//read allowed_alignments if maxsep < 0
	if (*maxsep < 0) {

		for (i=0;i<=ct1->GetSequenceLength();++i) {
			for (j=0;j<=ct2->GetSequenceLength();++j) {
				read(&sav,&(allowed_alignments[i][j]));
			}
		}
	}

	//now calculate the lowend and highend arrays that indicate the first and last index to which a nucleotide
		//from sequence 1 can be aligned.

	if (*maxsep>0) {
		//Using the traditional M parameter to constrain the alignment space
		for (i=0;i<2*ct1->GetSequenceLength();++i) {
			lowend[i]=lowlimit(i,*maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength());
			highend[i]=highlimit(i,*maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength());
		}
	}
	else {
		//allowed_alignments must be defined to constrain the alignment space
		//need to determine lowend and highend
		for (i=0;i<2*ct1->GetSequenceLength();++i) {
			lowend[i]=lowlimit(i,allowed_alignments,ct1->GetSequenceLength(),ct2->GetSequenceLength());
			highend[i]=highlimit(i,allowed_alignments,ct1->GetSequenceLength(),ct2->GetSequenceLength());
		}


	}

	//Allocate v, w, vmod, w3, w5
	v->allocate(ct1->GetSequenceLength(),ct2->GetSequenceLength(),lowend,highend,ct1->tem,optimalonly);
	w->allocate(ct1->GetSequenceLength(),ct2->GetSequenceLength(),lowend,highend,optimalonly);
	if (modification) vmod->allocate(ct1->GetSequenceLength(),ct2->GetSequenceLength(),lowend,highend,optimalonly);

	if (!optimalonly) w3->allocate(ct1->GetSequenceLength(), ct2->GetSequenceLength(), lowend,highend);
	w5->allocate(ct1->GetSequenceLength(), ct2->GetSequenceLength(), lowend,highend);


  for (i=0;i<=ct1->GetSequenceLength();++i) {

    if (temp==2||temp==3) I = ct1->GetSequenceLength();
    else I = i+ct1->GetSequenceLength()-1;

    for (j=i;j<=I;++j) {


      for (k=lowend[i];k<=highend[i];k++) {


        for (l=lowend[j];l<=highend[j];l++) {
          if (j>ct1->GetSequenceLength()) {
            b = i;
            a = j-ct1->GetSequenceLength();
          }
          else {
            b = j;
            a = i;
          }

          if (ct1->tem[b][a]) read(&sav,&(v->f(i,j,k,l)));
          read(&sav,&(w->f(i,j,k,l)));

          if (modification) read(&sav,&(vmod->f(i,j,k,l)));


        }
      }

    }
  }





  for (i=0;i<=ct1->GetSequenceLength()+1;++i) {
    for (j=0;j<=ct2->GetSequenceLength()+1;++j) {
      if (temp < 2) read(&sav,&(w3->f(i,j))); //do not attemptto read w3 if this is an optimal-only file
      read(&sav,&(w5->f(i,j)));

    }

  }
  if (temp < 2) for (j=0/*lowlimit(ct1->GetSequenceLength()+1,*maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/;j<=ct2->GetSequenceLength()+1/*highlimit(ct1->GetSequenceLength()+1,*maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/;j++) read(&sav,&(w3->f(ct1->GetSequenceLength()+1,j)));


  read(&sav,&temp);
  if (temp==0) *local=false;
  else *local=true;

	//allocate everything
 
 
	
#ifdef DYNALIGN_II
  short vers;//save file version
	//start with structure information
	int sequencelength1;
	read(&sav,&(sequencelength1));
	read(&sav,&(ct1->intermolecular));

	
	int pairnumber;
	read(&sav,&(pairnumber));

	int pair5,pair3;
	for (i=0;i<pairnumber;i++) {
		read(&sav,&(pair5));
		
		read(&sav,&(pair3));
		ct1->AddPair(pair5,pair3);
	}
	read(&sav,&(pairnumber));
	for (i=0;i<pairnumber;i++) {
		read(&sav,&(pair5));
		read(&sav,&(pair3));

		ct1->AddForbiddenPair(pair5,pair3);
	}
	for (i=0;i<=ct1->GetSequenceLength();i++) {
		
		read(&sav,&(ct1->hnumber[i]));
		sav.read(&(ct1->nucs[i]),1);

	}

	for (i=0;i<=2*ct1->GetSequenceLength();i++) read(&sav,&(ct1->numseq[i]));
	
	//Read the number of nucleotides forced double stranded
	read(&sav,&(pairnumber));

	//Read and set double stranded nucleotides
	for (i=0;i<pairnumber;i++) {
		int doubles;
		read(&sav,&(doubles));

		ct1->AddDouble(doubles);

	}

	

	//Read the number of nucleotides that are not allowed to pair
	read(&sav,&(pairnumber));

	//Read the nucleotides that cannot pair and store the information.
	for (i=0;i<pairnumber;i++) {
		int singles;
		read(&sav,&(singles));

		ct1->AddSingle(singles);

	}


	//Read the number of chemically modified nucleotides
	read(&sav,&(pairnumber));

	//Read and store the chemically modified nucleotides
	for (i=0;i<pairnumber;i++) {

		int modified;
		read(&sav,&(modified));

		ct1->AddModified(modified);

	}
	
	//Read the number of Us in GU pairs
	read(&sav,&(pairnumber));

	//Read the Us in GU pairs, and store information
	for (i=0;i<pairnumber;i++) {
		
		int GU;
		read(&sav,&(GU));  

		ct1->AddGUPair(GU);

	}
	
	string label;
	read(&sav,&(label));
	ct1->SetSequenceLabel(label);

	read(&sav,&(ct1->templated));
	if (ct1->templated) {

		ct1->allocatetem();
		for (i=0;i<=ct1->GetSequenceLength();i++) {
			for (j=0;j<=i;j++) read(&sav,&(ct1->tem[i][j]));	

		}

	}

	read(&sav,&ct1->shaped);
	if (ct1->shaped) {
		ct1->SHAPE = new double [2*ct1->GetSequenceLength()+1];
		for (i=0;i<=2*ct1->GetSequenceLength();i++) read(&sav,&(ct1->SHAPE[i]));

	}
	

	//now read the array class data for v, w, and wmb:
	//now write the array class data for v, w, and wmb:
	for (i=0;i<=ct1->GetSequenceLength();i++) {
            for (j=0;j<=ct1->GetSequenceLength();j++) {
			read(&sav,&(single1_v->dg[i][j+i]));
			read(&sav,&(single1_w->dg[i][j+i]));
			read(&sav,&(single1_wmb->dg[i][j+i]));
                        read(&sav,&(single1_we->dg[i][j+i]));
			readsinglechar(&sav,&(single1_fce->dg[i][j]));
		}
	}

	for (i=0;i<=2*ct1->GetSequenceLength();i++) {
		read(&sav,&(single1_lfce[i]));
		read(&sav,&(single1_mod[i]));
	}

	read(&sav, single1_vmin);


	int sequencelength2;
	read(&sav,&(sequencelength2));
	read(&sav,&(ct1->intermolecular));

	
	read(&sav,&(pairnumber));

	for (i=0;i<pairnumber;i++) {
		read(&sav,&(pair5));
		
		read(&sav,&(pair3));
		ct2->AddPair(pair5,pair3);
	}
	read(&sav,&(pairnumber));
	for (i=0;i<pairnumber;i++) {
		read(&sav,&(pair5));
		read(&sav,&(pair3));

		ct2->AddForbiddenPair(pair5,pair3);
	}
	for (i=0;i<=ct2->GetSequenceLength();i++) {
		
		read(&sav,&(ct2->hnumber[i]));
		sav.read(&(ct2->nucs[i]),1);

	}

	for (i=0;i<=2*ct2->GetSequenceLength();i++) read(&sav,&(ct2->numseq[i]));
	
	//Read the number of nucleotides forced double stranded
	read(&sav,&(pairnumber));

	//Read and set double stranded nucleotides
	for (i=0;i<pairnumber;i++) {
		int doubles;
		read(&sav,&(doubles));

		ct2->AddDouble(doubles);

	}

	

	//Read the number of nucleotides that are not allowed to pair
	read(&sav,&(pairnumber));

	//Read the nucleotides that cannot pair and store the information.
	for (i=0;i<pairnumber;i++) {
		int singles;
		read(&sav,&(singles));

		ct2->AddSingle(singles);

	}


	//Read the number of chemically modified nucleotides
	read(&sav,&(pairnumber));

	//Read and store the chemically modified nucleotides
	for (i=0;i<pairnumber;i++) {

		int modified;
		read(&sav,&(modified));

		ct2->AddModified(modified);

	}
	
	//Read the number of Us in GU pairs
	read(&sav,&(pairnumber));

	//Read the Us in GU pairs, and store information
	for (i=0;i<pairnumber;i++) {
		
		int GU;
		read(&sav,&(GU));  

		ct2->AddGUPair(GU);

	}
	
	read(&sav,&(label));
	ct2->SetSequenceLabel(label);

	read(&sav,&(ct2->templated));
	if (ct2->templated) {

		ct2->allocatetem();
		for (i=0;i<=ct2->GetSequenceLength();i++) {
			for (j=0;j<=i;j++) read(&sav,&(ct2->tem[i][j]));	

		}

	}

	read(&sav,&ct2->shaped);
	if (ct2->shaped) {
		ct2->SHAPE = new double [2*ct2->GetSequenceLength()+1];
		for (i=0;i<=2*ct2->GetSequenceLength();i++) read(&sav,&(ct2->SHAPE[i]));

	}
	

	//now read the array class data for v, w, and wmb:
	//now write the array class data for v, w, and wmb:
	for (i=0;i<=ct2->GetSequenceLength();i++) {
            for (j=0;j<=ct2->GetSequenceLength();j++) {
			read(&sav,&(single2_v->dg[i][j+i]));
			read(&sav,&(single2_w->dg[i][j+i]));
			read(&sav,&(single2_wmb->dg[i][j+i]));
                        read(&sav,&(single2_we->dg[i][j+i]));
			readsinglechar(&sav,&(single2_fce->dg[i][j]));
		}
	}

	for (i=0;i<=2*ct2->GetSequenceLength();i++) {
		read(&sav,&(single2_lfce[i]));
		read(&sav,&(single2_mod[i]));
	}

	read(&sav, single2_vmin);


#else
#endif

  sav.close();
}

 
                                                            


//Refold dynalign using save file information
//return an int that indicates whther an error occurred (0 = no error)
#ifdef DYNALIGN_II
int refolddynalign(const char* filename, structure *ct1, structure *ct2,
                    short **alignment, short maxtracebacks,
		   short window, short awindow, short percentsort,bool forced,short** forcealign)
#else
int refolddynalign(const char* filename, structure *ct1, structure *ct2,
                    short **alignment, short maxtracebacks,
		   short window, short awindow, short percentsort)
#endif
 {


  short i,j,ip,kp,ipp,kpp,cp;
#ifdef DYNALIGN_II
  int max_elongation;
  short islope,iintercept;
#else
  bool singleinsert;
#endif

  short maxsep,gap,lowest,index,*lowend,*highend;
  //integersize ****v,****w,****vmod,**w3,**w5;
  dynalignarray *w,*vmod;
  varray *v;
  wendarray *w3,*w5;
  datatable *data;
  int modificationflag;
  bool optimalonly;
  bool local,found;
  bool **allowed_alignments;
  int error;

#ifdef DYNALIGN_II
  char **fce1;
  char **fce2;

  //dbl1 and dbl2 flag nucleotides that must be double-stranded in seq1 and seq2, respectively
  //false means a nucleotide is not forced double
  //true indicates a nucleotide is forced double
  bool *dbl1,*dbl2;

  //mod indicates whether a nucleotide is chemically modified
  //mod1 for seq1 and mod2 or seq2
  bool *mod1,*mod2;
  bool modification;
  bool alignmentforced;
  int single1_vmin=0;
  int single2_vmin=0;
  integersize *single1_w5,*single1_w3,*single2_w5,*single2_w3;
  bool *single1_lfce,*single2_lfce,*single1_mod,*single2_mod;
  // DynProgArray<integersize> *single1_w2,*single1_wmb2,*single2_w2,*single2_wmb2;
  
  //allocate everything
  DynProgArray<integersize> single1_w(ct1->GetSequenceLength()),single2_w(ct2->GetSequenceLength());
  DynProgArray<integersize> single1_v(ct1->GetSequenceLength()),single2_v(ct2->GetSequenceLength());
  DynProgArray<integersize> single1_wmb(ct1->GetSequenceLength()),single2_wmb(ct2->GetSequenceLength());
  DynProgArray<integersize> single1_we(ct1->GetSequenceLength(),0),single2_we(ct2->GetSequenceLength(),0);
  forceclass single1_fce(ct1->GetSequenceLength()),single2_fce(ct2->GetSequenceLength());


  single1_lfce = new bool [2*ct1->GetSequenceLength()+1];
  single2_lfce = new bool [2*ct2->GetSequenceLength()+1];
  single1_mod = new bool [2*ct1->GetSequenceLength()+1];
  single2_mod = new bool [2*ct2->GetSequenceLength()+1];
 single1_w5 = new integersize [ct1->GetSequenceLength()+1];
  single1_w3 = new integersize [ct1->GetSequenceLength()+2];
  single2_w5 = new integersize [ct2->GetSequenceLength()+1];
  single2_w3 = new integersize [ct2->GetSequenceLength()+2];
  //the following are used for chemical modfication cases
#else
#endif

  data = new datatable();

  //open the save file to peek at the sizes needed to allocate arrays:
  ifstream sav(filename,ios::binary);

  read(&sav, &modificationflag); //note that this includes modification and
  //optimal only data

  if (modificationflag==2||modificationflag==3) optimalonly=true;
  else optimalonly=false;

  //start with structure information
   int length1,length2;
  read(&sav,&(length1));
  read(&sav,&(length2));
  read(&sav,&maxsep);
  sav.close();



  if (maxsep<0) {
	//This means that the allowed alignments were constrained by the HMM alignment
	  //probs and not by the traditional M parameter.
	  //Therefore, space must be allocated for the allowed_alignment arrays.

	  allowed_alignments=new bool *[length1+1];
	  for (i=0;i<=length1;i++) allowed_alignments[i]=new bool [length2+1];

  }
  else allowed_alignments=NULL;
  //fill the low and highend arrays:
  //allocate space:
  lowend = new short [2*length1];
  highend = new short [2*length1];

  //DHM (11/30/06): move these allocations to where the file is read
  if (modificationflag==1||modificationflag==3) {

    vmod = new dynalignarray();


  }
  else vmod=NULL;

  v = new varray();
  w = new dynalignarray();




  if (!optimalonly) w3 = new wendarray();
  else w3= NULL;
  w5 = new wendarray();


#ifdef DYNALIGN_II
  opendynalignsavefile(filename, ct1, ct2, v, w, vmod, w3, w5, &single1_w,&single2_w,&single1_v,&single2_v,&single1_wmb,&single2_wmb,&single1_we,&single2_we,
		       single1_lfce,single2_lfce,single1_mod,single2_mod,&single1_fce,&single2_fce,&single1_vmin,&single2_vmin,data, &max_elongation,  &maxsep, &islope, &iintercept, &gap, &lowest,&local,allowed_alignments,lowend,highend);
#else
    opendynalignsavefile(filename, ct1, ct2, v, w, vmod, w3, w5, data, &singleinsert,  &maxsep, &gap, &lowest,&local,allowed_alignments,lowend,highend);
#endif

#ifdef DYNALIGN_II
 if (forced) {

    //allocate the fce array for constraints
    //The following definitions are bitwise applied to the fce array:
    //SINGLE applies to any fce(i,j) s.t. i or j should be single stranded
    //PAIR applies to any fce(i,j) where i is paired to j
    //NOPAIR applies to any fce(i,j) where either i or j is paired to
    //  another nucleotide or i and j are forbidden to pair
    //DUBLE applies to any fce(i.j) where an nuc, k, i<k<j is double stranded
    //INTER applies to any fce(i,j) such that some k, i<k<j, is the virtual linker
    //  used for intermolecular folding
    //INTER applies to any fce(i,i) s.t. i is the center of the virtual linker
    //The above terms are defined in define.h

    //note that INTER is not currently supported in Dynalign
    //fce = new dynalignarray(N,N2,maxseparation);


    //allocate space for single sequence configuration in fce1 and fce2
   int N=length1;
   int N2=length2;
   fce1 = new char *[2*N+1];

    for (i=0;i<=2*N;++i) {
      if (i<=N) {
        fce1[i] = new char [i+1];
        for (j=0;j<=i;++j) fce1[i][j]=0;
      }
      else {
        fce1[i] = new char [2*N-i+1];
        for (j=0;j<=2*N-i;++j) fce1[i][j]=0;
      }
    }

    fce2 = new char *[2*N2+1];

    for (i=0;i<=2*N2;++i) {
      if (i<=N2) {
        fce2[i] = new char [i+1];
        for (j=0;j<=i;++j) fce2[i][j]=0;
      }
      else {
        fce2[i] = new char [2*N2-i+1];
        for (j=0;j<=2*N2-i;++j) fce2[i][j]=0;
      }
    }

    dbl1 = new bool [2*N];
    mod1=new bool [2*N];
    for (i=0;i<2*N;++i) {
      dbl1[i]=false;
      mod1[i]=false;
    }
    dbl2 = new bool [2*N2];
    mod2 = new bool [2*N2];
    for (i=0;i<2*N2;++i) {
      dbl2[i]=false;
      mod2[i]=false;
    }

    //store the locations of dbl1 and dbl2 in ct1 and ct2, respectively
    //this facilitates checks done in functions edangle5 and edangle3
    ct1->fcedbl=dbl1;
    ct2->fcedbl=dbl2;

    //now assign values to fce according to the above key using the function dynalignforce
    dynalignforce(ct1,ct2/*,fce*/,fce1,fce2,dbl1,dbl2,mod1,mod2);

    if (ct1->GetNumberofModified()>0||ct2->GetNumberofModified()>0) {
      modification = true;
      //For chemical modification, a second v array, vmod is required
      vmod = new dynalignarray(N,N2,lowend,highend);

    }
    else {
      modification =false;
      vmod = NULL;
    }

    alignmentforced = (forcealign != NULL);

  } else {
    vmod = NULL;
	dbl1=NULL;
	dbl2=NULL;
	alignmentforced=false;
	fce1=NULL;fce2=NULL;
	modification=false;
	mod1=NULL;
	mod2=NULL;
  }
#else
#endif

  //double up the sequences as part of suboptimal traceback support:
  for (i=1;i<=ct1->GetSequenceLength();i++) {
    ct1->numseq[(ct1->GetSequenceLength())+i] = ct1->numseq[i];
  }
  for (i=1;i<ct2->GetSequenceLength();i++) ct2->numseq[ct2->GetSequenceLength()+i]=ct2->numseq[i];




 
  if (!optimalonly) {
#ifdef DYNALIGN_II
    error = dyntraceback(maxtracebacks,window,awindow,percentsort,
                 v,w,w3,w5,&single1_w,&single2_w,&single1_v,&single2_v,&single1_wmb,&single2_wmb,&single1_we,&single2_we,
			 single1_lfce,single2_lfce,single1_mod,single2_mod,&single1_fce,&single2_fce,
			 ct1,ct2,alignment,lowend,highend,islope, iintercept,gap,data,lowest,vmod,local,max_elongation,mod1, mod2, modification,
                         fce1, fce2,
                         alignmentforced, forcealign,forced);
#else
     error = dyntraceback(maxtracebacks,window,awindow,percentsort,
                 v,w,w3,w5,
                 ct1,ct2,alignment,lowend,highend,gap,data,singleinsert,lowest,vmod,local);
#endif
  } else {

    //need to preset the structure and alignment info to zeros

    for (index=0;index<=ct1->GetSequenceLength();index++) alignment[0][index]=0;

    //ct1->numofstructures=1;
    //ct2->numofstructures=1;
    ct1->AddStructure();
    ct2->AddStructure();

    //calculate the total energy
    ct1->SetEnergy(1,lowest);
    ct2->SetEnergy(1,lowest);
    

    found = false;
    for (ip=1;ip<=ct1->GetSequenceLength()&&!found;++ip) {



      cp = min(ct2->GetSequenceLength(),highlimit(ip,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength()));
      for (kp=max(1,lowlimit(ip,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength()));kp<=cp&&!found;++kp) {
        if (local) {
          if (lowest==w5->f(ip,kp)) {
            ipp=ip;
            kpp=kp;
            found=true;
          }

        }
        else if (abs((ct1->GetSequenceLength()-ip)-(ct2->GetSequenceLength()-kp))*gap+w5->f(ip,kp)==lowest) {
          ipp=ip;
          kpp=kp;
          found=true;
        }
      }
    }


#ifdef DYNALIGN_II
    error = dyntrace(1, ipp, 1, kpp, ct1, ct2, 1, alignment[0],
                     w, v, w3, w5, lowend, highend,data, islope, iintercept, gap, vmod, &single1_w,&single2_w,&single1_v,&single2_v,&single1_wmb,&single2_wmb,&single1_we,&single2_we,single1_lfce,single2_lfce,single1_mod,single2_mod,&single1_fce,&single2_fce,local,max_elongation,mod1, mod2, modification,
                         fce1, fce2,
			    alignmentforced, forcealign,forced,true);
#else
      error = dyntrace(1, ipp, 1, kpp, ct1, ct2, 1, alignment[0],
             w, v, w3, w5, lowend, highend,data, gap, vmod, local, true);
#endif
  }



  delete v;
  delete w;
  delete w3;
  delete w5;

  delete lowend;
  delete highend;

  if (maxsep<0) {

		for (i=0;i<=length1;i++) delete[] allowed_alignments[i];
		delete[] allowed_alignments;
  }

  if (modificationflag==1||modificationflag==3) delete vmod;



  //delete[] pair[0];
  //delete[] pair[1];
  //delete[] pair;
  delete data;

  return error;
}



//Fold a single sequence to decide what pairs should be allowed in a subsequent dynalign calculation.
void templatefromfold(structure *ct, datatable *data, int singlefold_subopt_percent) {
  integersize *w5,*w3;
  bool *lfce,*mod;
  short crit,i,j;
  int vmin;
  DynProgArray<integersize> *w2,*wmb2;

  //allocate everything
  DynProgArray<integersize> w(ct->GetSequenceLength());
  DynProgArray<integersize> v(ct->GetSequenceLength());
  DynProgArray<integersize> wmb(ct->GetSequenceLength());
#ifdef DYNALIGN_II
  DynProgArray<integersize> we(ct->GetSequenceLength(),0);
#else
#endif
  forceclass fce(ct->GetSequenceLength());


  lfce = new bool [2*ct->GetSequenceLength()+1];
  mod = new bool [2*ct->GetSequenceLength()+1];

  for (i=0;i<=2*ct->GetSequenceLength();i++) {
    lfce[i] = false;
    mod[i] = false;
  }

  w5 = new integersize [ct->GetSequenceLength()+1];
  w3 = new integersize [ct->GetSequenceLength()+2];


  for (i=0;i<=ct->GetSequenceLength();i++) {
    w5[i] = 0;
    w3[i] = 0;

  }
  w3[ct->GetSequenceLength()+1] = 0;

  if (ct->intermolecular) {
    w2 = new DynProgArray<integersize>(ct->GetSequenceLength());
    wmb2 = new DynProgArray<integersize>(ct->GetSequenceLength());



  }
  else {
    w2 = NULL;
    wmb2 = NULL;

  }

  force(ct,&fce,lfce);
  vmin=DYNALIGN_INFINITY;

  //perform the fill steps:(i.e. fill arrays v and w.)
#ifdef DYNALIGN_II
  fill(ct, v, w, wmb, fce, vmin,lfce, mod,w5, w3, false, data, w2, wmb2, &we);
#else
    fill(ct, v, w, wmb, fce, vmin,lfce, mod,w5, w3, false, data, w2, wmb2);
#endif
  //find the dots here:


	crit= ((short) (abs(vmin)*(float (singlefold_subopt_percent)/100.0)));

  crit = crit + vmin;


  i = 1;
  j = 2;
  while (i<(ct->GetSequenceLength())) {

    if ((v.f(i,j)+v.f(j,i+ct->GetSequenceLength()))>crit) {
      //don't allow this pair
      ct->tem[j][i]=false;

    }
    j++;
    if (j>ct->GetSequenceLength()) {
      i++;
      j=i+1;
    }
  }

  delete[] lfce;
  delete[] mod;

  delete[] w5;
  delete[] w3;

  if (ct->intermolecular) {
    delete w2;
    delete wmb2;


  }
}

//load a dynalign saved file and determine all possible pairs
//in structures within percentdots of the lowest free energy structure.  These pairs are then
//the only allowed pairs for dyalign structure prediction for the passed structure.
int templatefromdsv(structure *cttemplate, const char *savefilename, float maxdsvchange, int maxpairs)
{

	datatable *data;

	short i,j,k,l,crit;
	bool local,modification;
	short maxsep,gap,lowest,*lowend,*highend;

	dynalignarray *w,*vmod;
	varray *v;
	wendarray *w3,*w5;

	bool optimalonly;

	  int modificationflag;
#ifdef DYNALIGN_II
        int max_elongation;
        short islope, iintercept;
#else
        bool singleinsert;
#endif

	dynalignheap heap;

	bool **allowed_alignments;


	structure *ct1, *ct2;
	ct1 = new structure();
	ct2 = new structure();

#ifdef DYNALIGN_II
	int single1_vmin=0;
	int single2_vmin=0;
	integersize *single1_w5,*single1_w3,*single2_w5,*single2_w3;
	bool *single1_lfce,*single2_lfce,*single1_mod,*single2_mod;
  // DynProgArray<integersize> *single1_w2,*single1_wmb2,*single2_w2,*single2_wmb2;
  
  //allocate everything
	DynProgArray<integersize> single1_w(ct1->GetSequenceLength()),single2_w(ct2->GetSequenceLength());
	DynProgArray<integersize> single1_v(ct1->GetSequenceLength()),single2_v(ct2->GetSequenceLength());
	DynProgArray<integersize> single1_wmb(ct1->GetSequenceLength()),single2_wmb(ct2->GetSequenceLength());
	DynProgArray<integersize> single1_we(ct1->GetSequenceLength(),0),single2_we(ct2->GetSequenceLength(),0);
	forceclass single1_fce(ct1->GetSequenceLength()),single2_fce(ct2->GetSequenceLength());


	single1_lfce = new bool [2*ct1->GetSequenceLength()+1];
	single2_lfce = new bool [2*ct2->GetSequenceLength()+1];
	single1_mod = new bool [2*ct1->GetSequenceLength()+1];
	single2_mod = new bool [2*ct2->GetSequenceLength()+1];
	single1_w5 = new integersize [ct1->GetSequenceLength()+1];
	single1_w3 = new integersize [ct1->GetSequenceLength()+2];
	single2_w5 = new integersize [ct2->GetSequenceLength()+1];
	single2_w3 = new integersize [ct2->GetSequenceLength()+2];
#else
#endif
	data = new datatable();

	//open the save file to peek at the sizes needed to allocate arrays:
	ifstream sav(savefilename,ios::binary);

	int versioncheck;
	read(&sav, &versioncheck);
	


	//read a flag that indicates whether chemical modification constraints were used
	read(&sav, &modificationflag);

	if (modificationflag==2||modificationflag==3) optimalonly=true;
	else optimalonly=false;

	//start with structure information needed for array allocation
	int length1,length2;
	read(&sav,&(length1));
	read(&sav,&(length2));
	read(&sav,&maxsep);
	sav.close();

	if (maxsep<0) {
	//This means that the allowed alignments were constrained by the HMM alignment
	  //probs and not by the traditional M parameter.
	  //Therefore, space must be allocated for the allowed_alignment arrays.

	  allowed_alignments=new bool *[length1+1];
	  for (i=0;i<=length1;i++) allowed_alignments[i]=new bool [length2+1];

	}
	else allowed_alignments=NULL;

	//allocate the arrays
	if (modificationflag==1||modificationflag==3) {
		vmod = new dynalignarray();
	} else {
		vmod = NULL;
	}

	v = new varray();
	w = new dynalignarray();

	if (!optimalonly) w3 = new wendarray();
	else w3= NULL;
	w5 = new wendarray();


	lowend = new short [2*length1];
	highend = new short [2*length1];

	//open the save file
	// this function does not work on iced for some reason

#ifdef DYNALIGN_II
	opendynalignsavefile(savefilename, ct1, ct2, v, w, vmod, w3, w5, &single1_w,&single2_w,&single1_v,&single2_v,&single1_wmb,&single2_wmb,&single1_we,&single2_we,
                             single1_lfce,single2_lfce,single1_mod,single2_mod,&single1_fce,&single2_fce,&single1_vmin,&single2_vmin,data, &max_elongation,  &maxsep, &islope, &iintercept, &gap, &lowest,&local,allowed_alignments,lowend,highend);
#else
	opendynalignsavefile(savefilename, ct1, ct2, v, w, vmod, w3, w5, data, &singleinsert,  &maxsep, &gap, &lowest,&local,allowed_alignments,lowend,highend);
#endif
	//double up the sequences as part of suboptimal traceback support:
	for (i=1;i<=ct1->GetSequenceLength();i++) {
		ct1->numseq[(ct1->GetSequenceLength())+i] = ct1->numseq[i];
	}
	for (i=1;i<ct2->GetSequenceLength();i++) ct2->numseq[ct2->GetSequenceLength()+i]=ct2->numseq[i];

	if (ct1->GetNumberofModified()>0||ct2->GetNumberofModified()>0) modification = true;
	else modification = false;

	crit = lowest;
	/*for (i=1;i<=ct1->GetSequenceLength();i++) {
		for (k=max(lowend[i],1);k<=min(highend[i],ct2->GetSequenceLength());k++) {



				if (local) {
					if ((w5->f(i,k))<crit) {
						crit = w5->f(i,k);
					}
				}
				else {
					if ((w5->f(i,k)+gap*abs(ct1->GetSequenceLength()-i-(ct2->GetSequenceLength()-k))  )<crit) {
						crit = w5->f(i,k)+gap*abs(ct1->GetSequenceLength()-i-(ct2->GetSequenceLength()-k));
					}
				}

		}
	}*/

	//Check for case when energy is zero because an empty structure is better than one with pairs.
		//In this case, calculate a different crit value:
	if ((local&&lowest==0)||(!local&&(lowest==gap*abs(ct1->GetSequenceLength()-ct2->GetSequenceLength())))) {
		lowest = DYNALIGN_INFINITY;
		for (i=1;i<=ct1->GetSequenceLength();i++) {
			for (j=i+minloop;j<=ct1->GetSequenceLength();j++) {
			  for (k=max(/*i-maxsep*/lowend[i]/*lowlimit(i,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/,1);k<=(min(ct2->GetSequenceLength(),/*i+maxsep*/highend[i]/*highlimit(i,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/));k++) {
				for (l=max(/*j-maxsep*/lowend[j]/*lowlimit(j,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/,k);l<=(min(ct2->GetSequenceLength(),/*j+maxsep*/highend[j]/*highlimit(j,maxsep,ct1->GetSequenceLength(),ct2->GetSequenceLength())*/));l++) {



				  if (modification) {
					if (vmod->f(i,j,k,l)/*vmod[j][i][a][b]*/+vmod->f(j,i+ct1->GetSequenceLength(),l,k+ct2->GetSequenceLength())/*vmod[i+N][j-i][b][a]*/<=lowest)
					  lowest = vmod->f(i,j,k,l)+vmod->f(j,i+ct1->GetSequenceLength(),l,k+ct2->GetSequenceLength());
				  }
				  else {
					  if (v->f(i,j,k,l)+v->f(j,i+ct1->GetSequenceLength(),l,k+ct2->GetSequenceLength())<=lowest)
							lowest =v->f(i,j,k,l)+v->f(j,i+ct1->GetSequenceLength(),l,k+ct2->GetSequenceLength());
				  }

				}//loop over l
			  }//loop over k
			}//loop over j
		  }//loop over i
		crit = lowest;
	}


//zane 07/20/2009:
	int minimum = crit; // lowest free energy * 10, usually a negative number
        float stepsize = 1 * conversionfactor;
	int dsv_energy = minimum + stepsize; //It's gonna be the energy threshold defined by maxpairs
	unsigned int PairsCount = 0; // count the number of pairs, with each of which the lowest energy structure has energy lower than dsv_energy
        //cout << maxpairs << endl;
        while ( PairsCount < maxpairs ){
          PairsCount=0;
          for (i=1;i<=ct1->GetSequenceLength();i++)
            for (j=i+minloop;j<=ct1->GetSequenceLength();j++) {
              bool lower = false;

              for (k=max(lowend[i],1);k<=(min(ct2->GetSequenceLength(),highend[i]));k++) {
                for (l=max(lowend[j],k);l<=(min(ct2->GetSequenceLength(),highend[j]));l++)
                  if (v->f(i,j,k,l)+v->f(j,i+ct1->GetSequenceLength(),l,k+ct2->GetSequenceLength())<dsv_energy){
                    lower = true;//the lowest energy structure( with i paired with j, k paired with l, i and j aligned to k and l respectively) is found to be lower than dsv_energy
                    break;//terminate l loop when this lowest energy structure is lower than dsv_energy
                  }
                if(lower)
                  break; //terminate k loop when this lowest energy structure is lower than dsv_energy
              }

              if(lower) // count this i-j pair in PairsCount
                ++PairsCount;
            }
          dsv_energy += stepsize;
          if (dsv_energy >= DYNALIGN_INFINITY) break;
        }

        dsv_energy -= stepsize;

        if (minimum < 0)
            maxdsvchange = minimum - minimum * maxdsvchange / 100;
        else
            maxdsvchange = minimum + minimum * maxdsvchange / 100;

        crit = maxdsvchange < dsv_energy ? dsv_energy : maxdsvchange;
        //crit = (short) ((float) crit -  ((float) crit*(((float) maxdsvchange)/100.0 )));
        //cout << maxdsvchange << ' ' << dsv_energy << ' ' << crit << endl;
	////
	for (i=1;i<=ct1->GetSequenceLength();i++) {
		for (j=i+minloop;j<=ct1->GetSequenceLength();j++) {

			cttemplate->tem[j][i]=false; //the default is to not allow the pair...
			for (k=max(lowend[i],1);k<=(min(ct2->GetSequenceLength(),highend[i]));k++) {
				for (l=max(lowend[j],k);l<=(min(ct2->GetSequenceLength(),highend[j]));l++) {
						if (v->f(i,j,k,l)+v->f(j,i+ct1->GetSequenceLength(),l,k+ct2->GetSequenceLength())<=crit)
							cttemplate->tem[j][i]=true;
				}
			}
		}
	}

	delete w3;
	delete w5;
	delete w;
	delete v;

	if (modificationflag==1||modificationflag==3) delete vmod;

	if (maxsep<0) {

		for (i=0;i<=length1;i++) delete[] allowed_alignments[i];
		delete[] allowed_alignments;
	}

	delete data;
	delete ct1;
	delete ct2;
	return 0;
}

//Using a set of known pairs for ct that have been previously loaded, set those as the only allowed pairs for a Dynalign calculation.
void templatefromct(structure *ct)
{
	short i,j;

	for (i=1;i<=ct->GetSequenceLength();i++) {
		for (j=i+minloop;j<=ct->GetSequenceLength();j++) {

			if (ct->GetPair(i)==j) ct->tem[j][i]=true; //allow this pair because it occurs in the structure
			else ct->tem[j][i]=false; //the default is to not allow the pair...

		}
	}


}


//Read alignment constraints from disk
void readalignmentconstraints(const char *filename, short **forcealign, structure *ct1, structure *ct2) {
	ifstream in;
	int i,j;

	in.open(filename);


	in >> i;
	in >> j;
	while (i!=-1) {
		forcealign[0][i]=j;
		forcealign[1][j]=i;

	        in >> i;
	        in >> j;
	}

	in.close();


}

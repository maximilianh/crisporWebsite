/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005,2006
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 */

#ifndef DYNALIGN_H
#define DYNALIGN_H

#include "TProgressDialog.h"
#include "dynalignarray.h"
#include "rna_library.h"
#include "structure.h"
#include "dynalignstackclass.h"
#ifdef _WINDOWS_GUI
	#include "../Windows_interface_deprecated/platform.h" 
#else
	#include "platform.h"
#endif //_WINDOWS_GUI
#include "varray.h"
#include "wendarray.h"
#include "DynProgArray.h"
#include "forceclass.h"

short int edangle5noforce(int i,int j,int ip,structure* ct,datatable* data);
short int edangle3noforce(int i,int j,int ip,structure* ct,datatable* data);
short int edangle5force(int i,int j,int ip,structure* ct,datatable* data);
short int edangle3force(int i,int j,int ip,structure* ct,datatable* data);
short int ebp(int i,int j,int ip,int jp,structure *ct, datatable *data);

#ifdef DYNALIGN_II
short int branch_1(DynProgArray<integersize> *single1_w, dynalignarray *w,int i, int j, int k, int l,int N,int c,short *lowend, short *highend, short iintercept, short islope);
short int branch_2(DynProgArray<integersize> *single2_w, dynalignarray *w,int i, int j, int k, int l,int N2,int d,short *lowend, short *highend, short iintercept, short islope);
void trace_branch_1(int i,int j,int k,int l,int c,bool& found,structure *ct1,datatable *data,DynProgArray<integersize> *single1_v, DynProgArray<integersize> *single1_w, DynProgArray<integersize> *single1_wmb, 
		    bool *single1_lfce, forceclass *single1_fce, bool *single1_mod,dynalignarray *w,dynalignstackclass& stack,short *lowend,short *highend,short int energy,short int addition, short iintercept, short islope);
void trace_branch_2(int i,int j,int k,int l,int d,bool& found,structure *ct2,datatable *data,DynProgArray<integersize> *single2_v, DynProgArray<integersize> *single2_w, DynProgArray<integersize> *single2_wmb, 
		    bool *single2_lfce, forceclass *single2_fce, bool *single2_mod,dynalignarray *w,dynalignstackclass& stack,short *lowend,short *highend,short int energy,short int addition, short iintercept, short islope);
#else
#endif
// The following functions, iref, jref, jderef, and ideref are used
// when addressing the large arrays
short iref(short i, short j, short N);
short jref(short i, short j, short N);
short ideref(short i, short j, short N);
short jderef(short i, short j, short N);

// Return the most 5' nucleotide in sequence 2 that can be aligned to
// nucleotide i in sequence 1, using the M constraint:
short lowlimit(short i, short M, short N1, short N2);
//overload to use a bool array of allowed alignments:
short lowlimit(short i, bool **allowed_alignments, short N1, short N2);

// Return the most 3' nucleotide in sequence 2 that can be aligned to
// nucleotide i in sequence 1, using the M constraint:
short highlimit(short i, short M, short N1, short N2);
//overload to use a bool array of allowed alignments:
short highlimit(short i, bool **allowed_alignments, short N1, short N2);

// Return a reference to a position in array based on i, k, N, N2, and M
/* short reference(short i, short k, short N, short N2, short M); */


// Encode the fact that nucleotide #nopair must be single stranded
void dynalignfceunpaired(structure *ct,char **fce,int nopair);

void dynforcedbl(int dbl,structure* ct,char **fce,bool *lineardbl);

// Register in array v that x and y must pair
void dynforcepair(int x,int y,structure *ct,char **v);

// Register in fce array that x must be paired to a G, as in FMN
// cleavage
void dynforcepairg(int x,structure *ct,char **fce);

// Register constraints into arrays
void dynalignforce(structure *ct1, structure *ct2,
                   /*short ****fce,*/ char **fce1, char **fce2, bool *dbl1, 
                   bool *dbl2, bool *mod1, bool *mod2);

/* This function takes 2 sequences and creates an alignment
 *
 * It requires maxseparation (the largest difference in position that
 * is permissible to align), maxgap (the maximum length of a gap),
 * comp (a bonus for compensating changes), gapstart (an affine gap
 * penalty for starting a gap), and gapincrease (an affine gap penalty
 * for each gapped position) datatable is the thermodynamic data
 * files.
 *
 * maxtracebacks, window, awindow, and percent sort apply to
 * suboptimal tracebacks
 *
 * window is the window parameter adapted from Zuker, Science, 1989,
 * and applied to structure dimension
 *
 * awindow is similar to window, but applied to the alignment dimension
 *
 * maxtracebacks is the maximum number of suboptimal structures and alignments
 *
 * percent sort is the maximum percent difference in score from the
 * lowest score (score = energy of stucture for seq 1 + energy of
 * structure for seq 2 + gap penalties)
 *
 * energyonly is a bool that, when true, indicates that dynalign is
 * just being used to calculate the lowest free energy possible for a
 * structure common to two sequences.  It halves the calculation.
 * ct1->energy[1] holds the value alignment, maxtracebacks, window,
 * awindow, percentsort are then not used in calculation.  This is
 * prototyped to false.
 *
 * local is a bool that, when true, indicates that a local alignment
 * calculation is being performed.  This is prototyped to false.
 *
 * force is a bool that will cause the algorithm to use experimental
 * constraints
 *
 * allowed_alignments is a bool array that summarizes the allowed nucleotide alignments between the two sequences.
 *
 *Return an error code that indicates whether an error occured.  (0=no error, 14=traceback error). 
 *
 */
#ifdef DYNALIGN_II
int dynalign(structure *ct1, structure *ct2, short **alignment,
             short int maxseparation, short int islope, short int iintercept, short int gapincrease, datatable *data, 
              short maxtracebacks, short window, short awindow, short percentsort, short **forcealign, int max_elongation,bool **allowed_alignments=NULL, 
              ProgressHandler *progress=NULL, const char *Savefile=NULL,
              bool energyonly = false, bool local = false, bool force = false, short numProcessors = 1);
#else
int dynalign(structure *ct1, structure *ct2, short **alignment,
              short int maxseparation, short int gapincrease, datatable *data, bool singleinsert, 
              short maxtracebacks, short window, short awindow, short percentsort, short **forcealign, bool **allowed_alignments=NULL, 
              ProgressHandler *progress=NULL, const char *Savefile=NULL,
              bool energyonly = false, bool local = false, bool force = false, short numProcessors = 1);
#endif
//calculate a single point in the v and w arrays -- allowing constraints if force is true

#ifdef DYNALIGN_II
void dynalignstep(structure *ct1, structure *ct2, datatable *data, varray *v, dynalignarray *w, wendarray *w5, wendarray *w3,  
                  dynalignarray *vmod, DynProgArray<integersize> *single1_w, DynProgArray<integersize> *single2_w,bool *mod1, bool *mod2, bool modification,
                  char **fce1, char **fce2, bool alignmentforced, short **forcealign, bool *dbl1, bool *dbl2,int max_elongation,
                  int i, int j, int k, int l, int N, int N2,
                  short *lowend, short *highend, int islope, int iintercept, int gap,
                  bool force,bool local,
                  short int (*edangle5)(int i, int j, int ip, structure *ct, datatable *data),
                  short int (*edangle3)(int i, int j, int ip, structure *ct, datatable *data));
#else

void dynalignstep(structure *ct1, structure *ct2, datatable *data, varray *v, dynalignarray *w, wendarray *w5, wendarray *w3,  
                  dynalignarray *vmod, bool *mod1, bool *mod2, bool modification,
                  char **fce1, char **fce2, bool alignmentforced, short **forcealign, bool *dbl1, bool *dbl2,  
                  int i, int j, int k, int l, int N, int N2,
                  short *lowend, short *highend, int gap, bool singleinsert,
                  bool force,bool local,
                  short int (*edangle5)(int i, int j, int ip, structure *ct, datatable *data),
                  short int (*edangle3)(int i, int j, int ip, structure *ct, datatable *data));
#endif
/* dyntrace traces back secondary structures using the arrays filled
 * in dynalign().
 *
 * Traceback a structure with the i-j pair is startopen=false
 * Otherwise, do an optimal traceback by providing i = 1, j = N1, a =
 * 1, b = N2
 *
 * Tracebacks are performed from nucleotides i to j in sequence 1 and
 * k to l in sequence 2.
 * Basepairs are stored in structure #structnum.  Alignemnt data is
 * stored in alignemnt.
 * Startopen can be used to do optimal traceback by settting as true
 * and passing i = 1, j = N, a = 1, b = N2.
 * local = true indicates local alignment.
 *
 *Return an error code that indicates whether an error occurred.  (0=no error, 14 = tracaeback error)
 *
 */
#ifdef DYNALIGN_II
int dyntrace(short i, short j, short a, short b, structure *ct1, structure *ct2,
	     short structnum, short *alignment,
	     dynalignarray *w, varray *v, wendarray *w3, wendarray *w5, 
	     short *lowend, short *highend,datatable *data, short islope, short iintercept, short gap, dynalignarray *vmod,DynProgArray<integersize> *single1_w,DynProgArray<integersize> *single2_w,DynProgArray<integersize> *single1_v,DynProgArray<integersize> *single2_v,DynProgArray<integersize> *single1_wmb,DynProgArray<integersize> *single2_wmb,DynProgArray<integersize> *single1_we,DynProgArray<integersize> *single2_we,bool *single1_lfce,bool *single2_lfce,bool *single1_mod,bool *single2_mod,forceclass *single1_fce,forceclass *single2_fce,
	     bool local, int max_elongation,bool *mod1, bool *mod2, bool modification,
	     char **fce1, char **fce2, bool alignmentforced, short **forcealign,bool force,bool startopen=false);

/* Percentsort is the maximum % difference in energy in suboptimal
 * structures if > 0 otherwise, percentsort is the maximum kcal/mol*10
 * difference in energy.
  * Return an int that indicates an errorcode.  (0= no error, 14 = traceback error)
 */
int dyntraceback(short maxtracebacks, short window, short awindow,
                  short percentsort,
                  varray *v, dynalignarray *w, wendarray *w3, wendarray *w5,DynProgArray<integersize> *single1_w,DynProgArray<integersize> *single2_w,DynProgArray<integersize> *single1_v,DynProgArray<integersize> *single2_v,DynProgArray<integersize> *single1_wmb,DynProgArray<integersize> *single2_wmb,DynProgArray<integersize> *single1_we,DynProgArray<integersize> *single2_we,bool *single1_lfce,bool *single2_lfce,bool *single1_mod,bool *single2_mod,forceclass *single1_fce,forceclass *single2_fce,
                  /*bool **pair,*/
                  structure *ct1, structure *ct2, short **alignment,
                 short *lowend, short *highend, short int islope, short int iintercept, short int gapincrease,
                  datatable *data,
                  integersize lowest, dynalignarray *vmod,
		 bool local,int max_elongation,bool *mod1, bool *mod2, bool modification,
                  char **fce1, char **fce2, bool alignmentforced, short **forcealign,bool force);
#else

int dyntrace(short i, short j, short a, short b, structure *ct1, structure *ct2,
              short structnum, short *alignment,
              dynalignarray *w, varray *v, wendarray *w3, wendarray *w5, 
              short *lowend, short *highend,datatable *data, short gap, dynalignarray *vmod,
              bool local, bool startopen=false);

/* Percentsort is the maximum % difference in energy in suboptimal
 * structures if > 0 otherwise, percentsort is the maximum kcal/mol*10
 * difference in energy.
  * Return an int that indicates an errorcode.  (0= no error, 14 = traceback error)
 */
int dyntraceback(short maxtracebacks, short window, short awindow,
                  short percentsort,
                  varray *v, dynalignarray *w, wendarray *w3, wendarray *w5,
                  /*bool **pair,*/
                  structure *ct1, structure *ct2, short **alignment,
                  short *lowend, short *highend, short int gapincrease,
                  datatable *data, bool singleinsert,
                  integersize lowest, dynalignarray *vmod,
                  bool local);
#endif

void alignout(short** align, const char *aout, structure *ct1, structure *ct2);

/* This function parses out alignement data in seq1 and seq2 and saves
 * the info in ct.tem[][] -- initialize ct.tem before calling this
 * function.
 */
void parse(structure *ct, char *seq1, char *seq2, datatable *data);

// Open a sav file for generating a dot plot or re-doing suboptimal
// structure prediction.
//This function takes care of calling allocate for both structure classes (ct1 and ct2).
#ifdef DYNALIGN_II
void opendynalignsavefile(const char *filename, structure *ct1, structure *ct2,
                          varray *v, dynalignarray *w, dynalignarray *vmod,
                          wendarray *w3, wendarray *w5, DynProgArray<integersize> *single1_w,DynProgArray<integersize> *single2_w,DynProgArray<integersize> *single1_v,DynProgArray<integersize> *single2_v,DynProgArray<integersize> *single1_wmb,DynProgArray<integersize> *single2_wmb,DynProgArray<integersize> *single1_we,DynProgArray<integersize> *single2_we,bool *single1_lfce,bool *single2_lfce,bool *single1_mod,bool *single2_mod,forceclass *single1_fce,forceclass *single2_fce,int *single1_vmin,int *single2_vmin,datatable *data,
                          int *max_elongation, short *maxsep, short *islope, short *iintercept, short *gap,
                          short *lowest, bool *local, bool **allowed_alignments,
						  short *lowend, short *highend);
#else

void opendynalignsavefile(const char *filename, structure *ct1, structure *ct2,
                          varray *v, dynalignarray *w, dynalignarray *vmod,
                          wendarray *w3, wendarray *w5, datatable *data,
                          bool *singleinsert, short *maxsep, short *gap,
                          short *lowest, bool *local, bool **allowed_alignments,
						  short *lowend, short *highend);
#endif


int traceback_debugger(const char* filename, structure *ct1, structure *ct2,int i,int j,int k,int l, short **alignment, short maxtracebacks,
		       short window, short awindow, short percentsort,bool forced,short** forcealign,bool outer);
//Redo a dynalign calculation using saved data from filename.
//Return an int that indicates errors.  (0=no error, 14 = traceback error).
#ifdef DYNALIGN_II
int refolddynalign(const char* filename, structure *ct1, structure *ct2,
                    short **alignment, short maxtracebacks,
                       short window, short awindow, short percentsort,bool forced, short** forcealign);
#else
int refolddynalign(const char* filename, structure *ct1, structure *ct2,
                    short **alignment, short maxtracebacks,
		   short window, short awindow, short percentsort);
#endif
/* Load a saved file from single sequence free energy minimization and
 * determine all possible pairs in structures within percentdots of
 * the lowest free energy structure.  These pairs are then the only
 * allowed pairs for dyalign structure prediction for the passed
 * structure.
 *
 * Use a single sequence save file to restrict pairs that will be
 * allowed in a Dynalign structure prediction.
 *
 * Deprecated.  Josh Keegan 2006-06-15
 */
//void templatefromsave(structure *cttemplate, char *savename);

/* Use single sequence free energy minimization to determine all
 * possible pairs in structures within percentdots of the lowest free
 * energy structure.  These pairs are then the only allowed pairs for
 * dyalign structure prediction for the passed structure.
 *
 * Use single sequence secondary structure prediction to restrict
 * pairs that will be allowed in a Dynalign structure prediction.
 *
 * singlefold_subopt_percent = maximum % difference in free energy change for allowed pairs
 */
void templatefromfold(structure *ct, datatable *data, int singlefold_subopt_percent=PERCENTDOTS);

//Below are two function that provide alternatives to templating the allowed pairs in the first sequence:

//This function determines which pairs should be allowed for ct by Dynalign, based on a previous Dynalign calculation.
int templatefromdsv(structure *ct, const char *savename, float maxdsvchange, int maxpairs);
//This function deteremined which pairs should be allowed for ct by Dynalign, based on a known structure.
void templatefromct(structure *ct);

//read alignment contraints from disk
void readalignmentconstraints(const char *filename, short **forcealign, structure *ct1, structure *ct2);


//Determine an alignment envelope of nucleotide pairs that are allowed to be aligned, uses HMM
void calculate_coinc_probs_env(structure *ct1, structure *ct2, bool **allowed_alignments, short **forcealign);

// inline implementations of select prototypes above

inline short int ebp(int i,int j,int ip,int jp,structure *ct, datatable *data) {
  return data->stack[(ct->numseq[i])][(ct->numseq[j])]
    [(ct->numseq[ip])][(ct->numseq[jp])];
}

inline short iref(short i, short j, short N) {
  if (j<=N) return i;
  else {
    if (i>N) return i-N;
    else return (i-j+N);
  }
}

inline short jref(short i, short j, short N) {
  if (i>N) return (j-N);
  else return j;
}

inline short ideref(short i, short j, short N) {
  if (i>N&&j>N) return i-N;
  else return i;
}

inline short jderef(short i, short j, short N) {
  if (i>N&&j>N) return j-N;
  else return j;
}

#endif

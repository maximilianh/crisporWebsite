/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Authors: Chris Connett and Andrew Yohn, 2006
 */

#include "rankconsumer.h"

#include "dynalign.h"
#include "dynalignarray.h"
#include "rankmanager.h"
#include "rna_library.h"
#include "structure.h"
#include "varray.h"
#include "wendarray.h"
#include "workslice.h"

void *rankconsumer(void *arg) {
  rankconsumerargs *rca = static_cast<rankconsumerargs*>(arg);

  RankManager *manager = rca->manager;
  structure *ct1 = rca->ct1;
  structure *ct2 = rca->ct2;


  datatable *data = rca->data;
  varray *v = rca->v;
  dynalignarray *w = rca->w;
  wendarray *w5 = rca->w5;
  wendarray *w3 = rca->w3;
  dynalignarray *vmod = rca->vmod;
  bool *mod1 = rca->mod1;
  bool *mod2 = rca->mod2;
  bool modification = rca->modification;
  char **fce1 = rca->fce1;
  char **fce2 = rca->fce2;
  bool alignmentforced = rca->alignmentforced;
  short **forcealign = rca->forcealign;
  bool *dbl1 = rca->dbl1;
  bool *dbl2 = rca->dbl2;
  short int *lowend = rca->lowend;
  short int *highend = rca->highend;
  int gap = rca->gap;
#ifdef DYNALIGN_II
  DynProgArray<integersize> *single1_w = rca->single1_w;
  DynProgArray<integersize> *single2_w = rca->single2_w;
  int max_elongation = rca->max_elongation;
  short islope = rca->islope;
  short iintercept = rca->iintercept;
#else
  bool singleinsert = rca->singleinsert;
#endif
  bool force = rca->force;
  bool local = rca->local;
  short int (*edangle5)(int i, int j, int ip, structure *ct, datatable *data) = rca->edangle5;
  short int (*edangle3)(int i, int j, int ip, structure *ct, datatable *data) = rca->edangle3;

  int N = ct1->GetSequenceLength();
  int N2 = ct2->GetSequenceLength();
  
  workslice slice;
  while (slice = manager->getWorkslice(), slice.begin != slice.end) {
    while (slice.begin != slice.end) {
#ifdef DYNALIGN_II
        dynalignstep(ct1, ct2, data,
                     v, w, w5, w3,
                     vmod, single1_w, single2_w, mod1, mod2, modification,
                     fce1, fce2,
                     alignmentforced, forcealign,
                     dbl1, dbl2, max_elongation,
                     slice.begin->i, slice.begin->j,
                     slice.begin->k, slice.begin->l,
                     N, N2, lowend, highend, islope, iintercept, gap, 
                     force,local,
                     edangle5, edangle3);
#else
      dynalignstep(ct1, ct2, data,
                   v, w, w5, w3,
                   vmod, mod1, mod2, modification,
                   fce1, fce2,
                   alignmentforced, forcealign,
                   dbl1, dbl2,
                   slice.begin->i, slice.begin->j,
                   slice.begin->k, slice.begin->l,
                   N, N2, lowend, highend, gap, singleinsert,
                   force,local,
                   edangle5, edangle3);
#endif
 slice.begin++;
    }
  }

  return NULL;
}

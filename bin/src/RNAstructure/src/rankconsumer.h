/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Authors: Chris Connett and Andrew Yohn, 2006
 */

#ifndef RANKCONSUMER_H
#define RANKCONSUMER_H

#include "DynProgArray.h"
#include "dynalignarray.h"
#include "rankmanager.h"
#include "rna_library.h"
#include "structure.h"
#include "varray.h"
#include "wendarray.h"
#include "workslice.h"

typedef struct {
  RankManager *manager;

#ifdef DYNALIGN_II
  DynProgArray<integersize> *single1_w;
  DynProgArray<integersize> *single2_w;
  int max_elongation;
  short islope, iintercept;
#else
#endif

  datatable *data;
  dynalignarray *vmod;
  dynalignarray *w;
  structure *ct1, *ct2;
  varray *v;
  wendarray *w5, *w3;

  bool *mod1, *mod2, *dbl1, *dbl2;
  short int *lowend, *highend;

  char **fce1, **fce2;
  short **forcealign;

  short int (*edangle5)(int i, int j, int ip, structure *ct, datatable *data);
  short int (*edangle3)(int i, int j, int ip, structure *ct, datatable *data);

  int gap;
  bool modification, alignmentforced;
#ifndef DYNALIGN_II
    bool singleinsert;
#else
#endif
  bool force, local;  
} rankconsumerargs;

void *rankconsumer(void *arg);

#endif

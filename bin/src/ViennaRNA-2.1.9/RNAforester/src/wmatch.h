#ifndef _WMATCH_H_
#define _WMATCH_H_

#ifdef __cplusplus
extern "C"
{
#endif
  
  int *Weighted_Match (void *gptr,int type,int maximize);
  
  // pairs.c
  void PAIR (int *outcome);
  void MERGE_PAIRS (int v);
  void LINK_PATH (int e);
  void INSERT_PAIR ();

  // pointer.c
  void POINTER (int u, int v, int e);
  void SCAN (int x, int del);

  // readgraph.c
  void SetUp (void* gptr, int type);
  void SetStandard(Graph graph);
  void SetEuclid(EuclidGraph graph);
  void SetMatrix(MatrixGraph graph);

  // termc.c
  void SET_BOUNDS ();
  void UNPAIR_ALL ();

  // unpairs.c
  void UNPAIR (int oldbase, int oldmate);
  void REMATCH (int firstmate, int e) ;
  void UNLINK (int oldbase);
  
#ifdef __cplusplus
}
#endif

#endif

#include <math.h>
//#include "stdafx.h"
#include "rna_library.h"
#include "pclass.h"
//#include "platform.h"





#define MaxStructure 10000000 // the max number of structures explored


//This version of the alltrace header is used by intermolecular.* for OligoWalk calculations.  

void alltrace(bool writect,int *ctenergy,structure* ct,datatable* data, short percentdelta, short absolutedelta, ProgressHandler* update, char* save);
//void readalltrace(char *filename, structure *ct, 
//			 short *w5,  
//			 atDynProgArray *v, atDynProgArray *w, atDynProgArray *wmb, atDynProgArray *wmbl, atDynProgArray *wl, atDynProgArray *wcoax,
//			 atDynProgArray *w2, atDynProgArray *wmb2, forceclass *fce, bool *lfce, bool *mod, datatable *data);

void realltrace(bool writect,int *ctenergy,char *savefilename, structure *ct, short percentdelta, short absolutedelta, char *ctname = NULL);


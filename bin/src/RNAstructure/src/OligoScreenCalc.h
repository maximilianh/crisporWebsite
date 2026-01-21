#if !defined(OLIGOSCREENCALC_H)
#define OLIGOSCREENCALC_H


#include "intermolecular.h"
#include "rna_library.h"

//Perform the OligoScreen calculation.
//infilename is the input list of oligos.
//outfilename is the output file to be written with a report.
//data is a a datatable that is pre-filled with DNA or RNA parameters as needed.
//hybriddata needs to contain DNA-RNA stack parameters if oligos are DNA.  It should = NULL for RNA oligos.
int OligoScreenCalc(const char *infilename, const char *outfilename, datatable *table, rddata *hybriddata);




#endif
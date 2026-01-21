
//Use the precompiler to make sure the class definition is not included more than once.
#if !defined(BIMOL_H)
#define BIMOL_H


#include "structure.h"
#include "rna_library.h"

//prototype for bimol, a bimolecular folding routine.
void bimol(structure *ct1, structure *ct2, structure *ct3, int maxloop, int maxtracebacks, int percent, int windowsize, 
           datatable *data);

//prototype for accessfold, bimolecular folding without intramolecular pairs, using a partition function heuristic to model accessibility to binding.
void accessfold(structure *ct1, structure *ct2, structure *ct3, int maxloop, int maxtracebacks, int percent, int windowsize, 
           datatable *data, double gamma, bool IsRNA, double temperature );


#endif

#if !defined (findp)

#define findp


#include "structure.h"
#define maxpseudo 100


void findpseudo(structure *ct, int structnum, int *npseudo, int *pseudo);

//return true if the structure #StructureNumber has a pseudoknot or false otherwise
bool HasPseudo(structure *ct, int StructureNumber);

#endif
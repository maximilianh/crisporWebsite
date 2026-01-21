

#include "findpseudo.h"


//find pseudoknots, break them, and return the number of pseudoknots and their location

//NOTE: This function has been replaced by a much more robust function, RNA::BreakPseudoknot(const bool minimum_energy=true, const int structurenumber = 0);
void findpseudo(structure *ct, int structnum, int *npseudo, int *pseudo) {
	register int i,j,k,L;
	int cons,incons;


   *npseudo = 0;
  for (i=1;i<=ct->GetSequenceLength();i++) {
   	if (ct->GetPair(i,structnum)>i) {//basepair found
		for (j=ct->GetSequenceLength();j>ct->GetPair(i,structnum);j--) {

			if ((ct->GetPair(j,structnum)<ct->GetPair(i,structnum))&&(ct->GetPair(j,structnum)>i)) {
				//a pseudoknot has been found
				cons=1;
				incons=1;
				for (k=i+1;k<j;k++) {
					if (ct->GetPair(k,structnum)>k) {
						
						if(ct->GetPair(k,structnum)>ct->GetPair(i,structnum)) incons++;
						else cons++;

					

					}


				}
				if (cons>incons) {
					L = ct->GetPair(i,structnum);
					for (k=i;k<L;k++) {

						if (ct->GetPair(k,structnum)>ct->GetPair(i,structnum)) {
							(*npseudo)++;
							if ((*npseudo)<maxpseudo) pseudo[*npseudo]=k;
							(*npseudo)++;
							if ((*npseudo)<maxpseudo) pseudo[*npseudo]=ct->GetPair(k,structnum);
							ct->SetPair(k,0,structnum);
							ct->SetPair(ct->GetPair(k,structnum),0,structnum);
							
				

						}

					}


				}
				else {
					L = ct->GetPair(i,structnum);
					for (k=i;k<L;k++) {

						if (ct->GetPair(k,structnum)<=ct->GetPair(i,structnum)&&ct->GetPair(k,structnum)>0) {
							(*npseudo)++;
							if ((*npseudo)<maxpseudo) pseudo[*npseudo]=k;
							(*npseudo)++;
							if ((*npseudo)<maxpseudo) pseudo[*npseudo]=ct->GetPair(k,structnum);
							ct->SetPair(k,0,structnum);
							ct->SetPair(ct->GetPair(k,structnum),0,structnum);
							
				

						}

					}


				}


			}
		}
   }
}

}

//return true if the structure #StructureNumber has a pseudoknot or false otherwise
bool HasPseudo(structure *ct, int StructureNumber) {
	int i,j;

	//use the standard definition of a pseudoknot to find pseudoknots
	for (i=1;i<=ct->GetSequenceLength();i++) {
		if (ct->GetPair(i,StructureNumber)>i) {
			//found a pair
			for (j=i+1;j<ct->GetPair(i,StructureNumber);j++) {
				if (ct->GetPair(j,StructureNumber)>0&&ct->GetPair(j,StructureNumber)<i) return true;
				else if (ct->GetPair(j,StructureNumber)>ct->GetPair(i,StructureNumber)) return true;

			}
		}
	}


	return false;//no pseudoknot found

}

#include "forceclass.h"

// forceclass encapsulates a large 2-d arrays of char used by the
// dynamic algorithm to enforce folding constraints

forceclass::forceclass(int size) {
  Size = size;
  register int i,j;
  dg = new char *[size+1];

	for (i=0;i<=(size);i++)  {
    dg[i] = new char [size+1];
  }
  for (i=0;i<=size;i++) {
    for (j=0;j<size+1;j++) {
      dg[i][j] = 0;
             
    }
  }
}

forceclass::~forceclass() {
  for (int i = 0; i <= Size; i++) {
    delete[] dg[i];
  }

  delete[] dg;
}

//this function is used to set up the fce array for a base, x, that should be single-stranded
void forcesingle(int x,structure* ct,forceclass *v) {
	int i;

	for (i=x;i<x+(ct->GetSequenceLength());i++) {
		v->f(x,i)=v->f(x,i)|SINGLE;
	}
	for (i=1;i<=x;i++) {
		v->f(i,x)=v->f(i,x)|SINGLE;
	}
	for (i=x+1;i<=ct->GetSequenceLength();i++) {
		v->f(i,x+ct->GetSequenceLength()) = v->f(i,x+ct->GetSequenceLength())|SINGLE;
	}
}

void forcepair(int x,int y,structure *ct,forceclass *v) {
	v->f(x,y) = v->f(x,y)|PAIR;
	v->f(y,x+ct->GetSequenceLength())=v->f(y,x+ct->GetSequenceLength())|PAIR;
	forcedomain(x, y, ct, v);
}

void forcedomain(int x,int y,structure *ct,forceclass *v) {
	int i,j;
	for (i=y+1;i<=x-1+ct->GetSequenceLength();i++) {
		v->f(x,i) = v->f(x,i)|NOPAIR;
	}
	for (i=x;i<=y-1;i++) {
		v->f(x,i) = v->f(x,i)|NOPAIR;
	}
	for (i=1;i<=x-1;i++) {
		v->f(i,y) = v->f(i,y)|NOPAIR;
	}
	for (i=x+1;i<=y;i++) {
		v->f(i,y) = v->f(i,y)|NOPAIR;
	}
	for (i=1;i<=x-1;i++) {
		v->f(i,x) = v->f(i,x)|NOPAIR;
	}
	for (i=y+1;i<=ct->GetSequenceLength();i++) {
		v->f(i,y+ct->GetSequenceLength())=v->f(i,y+ct->GetSequenceLength())|NOPAIR;
	}
	for (i=y;i<=x-1+(ct->GetSequenceLength());i++) {
		v->f(y,i) = v->f(y,i)|NOPAIR;
	}
	for (i=(ct->GetSequenceLength())+x+1;i<=(ct->GetSequenceLength())+y-1;i++) {
		v->f(y,i) = v->f(y,i)|NOPAIR;
	}
	for (i=x+1;i<=y-1;i++) {
		v->f(i,x+ct->GetSequenceLength()) = v->f(i,x+ct->GetSequenceLength())|NOPAIR;
	}
	for (i=y+1;i<=ct->GetSequenceLength();i++) {
		v->f(i,x+ct->GetSequenceLength()) = v->f(i,x+ct->GetSequenceLength())|NOPAIR;
	}
	for (i=1;i<=x-1;i++) {
		for (j = x+1;j<=y-1;j++){
			v->f(i,j) = v->f(i,j)|NOPAIR;
		}
	}
	for (i=x+1;i<=y-1;i++) {
		for (j=y+1;j<=(ct->GetSequenceLength())+x-1;j++) {
			v->f(i,j) = v->f(i,j)|NOPAIR;
		}
	}
	for (i=y+1;i<=ct->GetSequenceLength();i++) {
		for (j=(ct->GetSequenceLength())+x+1;j<=(ct->GetSequenceLength())+y-1;j++) {
			v->f(i,j) = v->f(i,j)|NOPAIR;
		}
	}
}

void forcedbl(int dbl,structure* ct,forceclass *w,bool *v) {
	int i,j;

	v[dbl] = true;
	v[dbl+ct->GetSequenceLength()] = true;

	for(i=dbl+1;i<=ct->GetSequenceLength();i++) {
		for (j=1;j<dbl;j++) {
			w->f(j,i) = w->f(j,i)|DUBLE;
		}
	}
	for(j=(dbl+(ct->GetSequenceLength())-1);j>ct->GetSequenceLength();j--) {
		for (i=dbl+1;i<=ct->GetSequenceLength();i++) {
			w->f(i,j) = w->f(i,j)|DUBLE;
		}
	}

}

void forceinter(int dbl,structure* ct,forceclass *w) {
	int i,j;

	for(i=dbl+1;i<=ct->GetSequenceLength();i++) {
		for (j=1;j<dbl;j++) {
			w->f(j,i) = w->f(j,i)|INTER;
		}
	}
	for(j=(dbl+(ct->GetSequenceLength())-1);j>ct->GetSequenceLength();j--) {
		for (i=dbl+1;i<=ct->GetSequenceLength();i++) {
			w->f(i,j) = w->f(i,j)|INTER;
		}
	}
	for(i=dbl+1+ct->GetSequenceLength();i<=2*ct->GetSequenceLength();i++) {
		for (j=ct->GetSequenceLength();j<dbl+ct->GetSequenceLength();j++) {
			w->f(j,i) = w->f(j,i)|INTER;
		}
	}
}

void forceinterefn(int dbl,structure* ct,forceclass *w) {
	int i,j;
	for(i=dbl+1;i<=ct->GetSequenceLength();i++) {
		for (j=1;j<dbl;j++) {
			w->f(j,i) = w->f(j,i)|INTER;
		}
	}
}

//force is used to enforce the folding constraints specified in a structure
//The following definitions are bitwise applied to the fce array:
//SINGLE applies to any fce(i,j) s.t. i or j should be single stranded
//PAIR applies to any fce(i,j) where i is paired to j
//NOPAIR applies to any fce(i,j) where either i or j is paired to
//  another nucleotide or i and j are forbidden to pair
//DOUBLE applies to any fce(i.j) where an nuc, k, i<k<j is double stranded
//INTER applies to any fce(i,j) such that some k, i<k<j, is the virtual linker
//  used for intermolecular folding
//INTER applies to any fce(i,i) s.t. i is the center of the virtual linker
//The above terms are defined in define.h


//mod[i] is a bool array that is set to true if i is a nuc accessible to chemical
//  modification
//lfce[i] is a bool array that is set to true if i is double stranded


//force also double the sequence s.t. ct->numseq[i+N] = ct->numseq[i] where
//  N is the number of nucleotides

void force(structure *ct,forceclass *fce, bool *lfce){
	int i,j;
	register int number;

	number = ct->GetSequenceLength();

	for (i=0;i<ct->GetNumberofSingles();i++) {
		
		if(ct->GetSingle(i)<=ct->GetSequenceLength()) forcesingle(ct->GetSingle(i),ct,fce);
	}

	for (i=0;i<ct->GetNumberofPairs();i++) {
		if(ct->GetPair5(i)<=ct->GetSequenceLength()&&ct->GetPair3(i)<=ct->GetSequenceLength()) {
			forcepair(ct->GetPair5(i),ct->GetPair3(i),ct,fce);
			forcedbl(ct->GetPair5(i),ct,fce,lfce);
			forcedbl(ct->GetPair3(i),ct,fce,lfce);
		}
	}

	for (i=0;i<ct->GetNumberofDoubles();i++) {
		if (ct->GetDouble(i)<=ct->GetSequenceLength()) forcedbl(ct->GetDouble(i),ct,fce,lfce);
	}

	for (i=0;i<ct->GetNumberofDomains();i++) {
		if(ct->GetDomain5(i)<=ct->GetSequenceLength()&&ct->GetDomain3(i)<=ct->GetSequenceLength()) {
			forcedomain(ct->GetDomain5(i),ct->GetDomain3(i),ct,fce);
		}
	}

	//u's in gu pairs must be double stranded
	for (i=0;i<ct->GetNumberofGU();i++) {
		if (ct->GetGUpair(i)<=ct->GetSequenceLength()) forcedbl(ct->GetGUpair(i),ct,fce,lfce);
	}



	if (ct->intermolecular) {//this indicates an intermolecular folding
		//for (i=0;i<=3;i++) {
		//  forcesingle(ct->inter[i])//don't allow the intermolecular indicators to pair
		//}
		for (i=0;i<3;i++) {
			forceinter(ct->inter[i],ct,fce);
		}

		fce->f(ct->inter[1],ct->inter[1]) = fce->f(ct->inter[1],ct->inter[1])|INTER;
	}

	for (i=0;i<ct->GetNumberofForbiddenPairs();i++) {
		if(ct->GetForbiddenPair5(i)<=ct->GetSequenceLength()&&ct->GetForbiddenPair3(i)<=ct->GetSequenceLength()) fce->f(ct->GetForbiddenPair5(i),ct->GetForbiddenPair3(i)) = fce->f(ct->GetForbiddenPair5(i),ct->GetForbiddenPair3(i))|NOPAIR;
		if(ct->GetForbiddenPair5(i)<=ct->GetSequenceLength()&&ct->GetForbiddenPair3(i)<=ct->GetSequenceLength()) fce->f(ct->GetForbiddenPair3(i),ct->GetForbiddenPair5(i)+ct->GetSequenceLength())=fce->f(ct->GetForbiddenPair3(i),ct->GetForbiddenPair5(i)+ct->GetSequenceLength())|NOPAIR;
	}

	//Double up the sequence
	for (i = 1; i <= number; i++) 
	{
		ct->numseq[(number)+i] = ct->numseq[i];
	}
	
	// modify numseq in msa stored in ct 
	// might create problems for fold; will have to keep a check
	if (ct->number_of_sequences > 1) 
	{
		for (int j = 0; j < ct->number_of_sequences; j++)
		{
			for (int i = 1; i <= number; i++)
			{
				ct->get_individual_sequence(j)->numseq[number + i] = ct->get_individual_sequence(j)->numseq[i];
			}
		}
	}

	//Double up the sequence
	for (i = 1; i <= number; i++)
	{
		ct->numseq[(number)+i] = ct->numseq[i];
	}

	//The next section handles the case where base pairs are not
	//not allowed to form between nucs more distant
	//than ct->GetPairingDistanceLimit()
	if (ct->DistanceLimited()) {

		if (!ct->templated) ct->allocatetem();

		for (j=minloop+2;j<=ct->GetSequenceLength();j++) {
			for (i=1;i<j;i++) {
				if (j-i>=ct->GetPairingDistanceLimit()) ct->tem[j][i]=false;
			}
		}
	}
}

//returns true if pair i and j is not a GU pair and the most adjacent pairs are not GU
//	used with chemical modification data to make sure that some contributions are not over-counted
bool notgu(int i, int j, structure *ct) {
	if ((ct->IsNuc(i,'G')||ct->IsNuc(i,'g'))&&(ct->IsNuc(j,'U')||ct->IsNuc(j,'u'))) return false;
	else if ((ct->IsNuc(i,'U')||ct->IsNuc(i,'u'))&&(ct->IsNuc(j,'G')||ct->IsNuc(j,'g'))) return false;
	else if ((ct->IsNuc(i+1,'G')||ct->IsNuc(i+1,'g'))&&(ct->IsNuc(j-1,'U')||ct->IsNuc(j-1,'u'))) return false;
	else if ((ct->IsNuc(i+1,'U')||ct->IsNuc(i+1,'u'))&&(ct->IsNuc(j-1,'G')||ct->IsNuc(j-1,'g'))) return false;
	else if (i>1) {//make sure i isn't 1
		if ((ct->IsNuc(i-1,'G')||ct->IsNuc(i-1,'g'))&&(ct->IsNuc(j+1,'U')||ct->IsNuc(j+1,'u'))) return false;
		else if ((ct->IsNuc(i-1,'U')||ct->IsNuc(i-1,'u'))&&(ct->IsNuc(j+1,'G')||ct->IsNuc(j+1,'g'))) return false;
	}
	return true;
}

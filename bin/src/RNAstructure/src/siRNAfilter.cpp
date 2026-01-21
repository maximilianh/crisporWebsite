#include "siRNAfilter.h"

//=======================================================================
//siPREFILTER is a class to calculate the pre-filter score of functional siRNA
//Those siRNA with certain score will be chosed as functional candidate for folding target


/*
datatable DATA stores the RNA folding free energy parameters.
datatable dhDATA stores the RNA folding enthalpy change parameters.
usefilter: =0 -> filter not used by OligoWalk.  >0-> filter used.
scoreit: set to true to use empirical parameters for siRNA filtering.
size: "ct->numofbases - length + 2": length of target sequence - length of oligo + 2


*/

siPREFILTER::siPREFILTER(datatable& DATA,datatable& dhDATA,int usefilter,bool scoreit,int size, bool isdna) {
	data=&DATA;
	dhdata=&dhDATA;
	melt= new float[size+1];
	useit=usefilter;
	usescore=scoreit;
	if (useit != 0)	score = new int[size+1];
	if (useit != 0)	enddiff = new double[size+1];
	//initialize the free energy arrays:
	stack[0][1]=0;
	stack[0][2]=0;
	stack[0][3]=0;
	stack[0][4]=0;
	stack[1][0]=0;
	stack[2][0]=0;
	stack[3][0]=0;
	stack[4][0]=0;
	if (isdna) {
	stack[1][1]=-1.0;
	stack[1][2]=-2.1;
	stack[1][3]=-1.8;
	stack[1][4]=-0.9;
	stack[2][1]=-0.9;
	stack[2][2]=-2.1;
	stack[2][3]=-1.7;
	stack[2][4]=-0.9;
	stack[3][1]=-1.3;
	stack[3][2]=-2.7;
	stack[3][3]=-2.9;
	stack[3][4]=-1.1;
	stack[4][1]=-0.6;
	stack[4][2]=-1.5;
	stack[4][3]=-1.6;
	stack[4][4]=-0.2;
	end[0]=0;
	end[1]=0;
	end[2]=0;
	end[3]=0;
	end[4]=0;
	}
	else {
	stack[1][1]=-.93;
	stack[1][2]=-2.24;
	stack[1][3]=-2.08;
	stack[1][4]=-1.1;
	stack[2][1]=-2.11;
	stack[2][2]=-3.26;
	stack[2][3]=-2.36;
	stack[2][4]=-2.08;
	stack[3][1]=-2.35;
	stack[3][2]=-3.42;
	stack[3][3]=-3.26;
	stack[3][4]=-2.24;
	stack[4][1]=-1.33;
	stack[4][2]=-2.35;
	stack[4][3]=-2.11;
	stack[4][4]=-.93;
	end[0]=0;
	end[1]=0.45;
	end[2]=0;
	end[3]=0;
	end[4]=0.45;
	}
}

siPREFILTER::~siPREFILTER() {

	delete[] melt;
	if (useit != 0) delete[] score;
	if (useit != 0)	delete[] enddiff;

}

void siPREFILTER::count(structure *ct,int i,int test) {
	
	int length=ct->GetSequenceLength();
	score[i]=0; 

	//Criteria 0: unstable 5'  end of antisense strand(AS)
	//find if the sequence with a unstable 5' AS end, with a widows having 2 base pairs
	//5' end of the siRNA strand
	DG5=0;
	for (j=1;j<2;j++)	DG5+=stack[ct->numseq[j]][ct->numseq[j+1]];
	DG5+=end[ct->numseq[1]];
	//substract the AU penalty if the 5th begin with AU
	DG5-=end[ct->numseq[2]];
	//3' end of the siRNA strand
	DG3=0;
	for (j=length-1;j<length;j++) DG3+=stack[ct->numseq[j]][ct->numseq[j+1]];
	DG3+=end[ct->numseq[length]];
	//substract the AU penalty if the 5th begin with AU
	DG3-=end[ct->numseq[length-1]];
	//functional siRNA AS  would have a unfavorable(more positive) DG5	
	
	enddiff[i]=DG5-DG3;

	
if (usescore) {
	//Criteria 0: unstable 5'  end of antisense strand(AS)
	//-2 point
	if ( enddiff[i] <= 0) score[i]-=3; //unstrict value(0) is used  to prefilting siRNA

	//Criteria I: 30-52% G/C content 
	//+1 point
	int numofGC=0;
	for (j=1; j<=ct->GetSequenceLength(); j++) {
		if (ct->numseq[j] ==2 || ct->numseq[j]==3) {
			numofGC++;
		}
	}
	float contentofGC= (float)numofGC/(float)ct->GetSequenceLength();
	if (contentofGC >=0.3 && contentofGC <=0.52)	score[i]++;


	//Criteria II: at least 3 'A/U' bases at positions 15-19(sense strand)( 1-5 for antisense)
	//+ 1xN points (1 point for each 'A/U', at most 5 points in total)
	for (j=1; j<=5; j++) {
		if (ct->numseq[j] ==1 || ct->numseq[j]==4) {
			score[i]++;
		}
	}

	//Criteria III: Absence of internal repeats (Tm of potential internal hairpin is< 57 degree)
	//+1 point
	dynamic(ct,data,1,10,0,0);//predict the self structure of oligo
	int dg= (ct->GetEnergy(1));
	efn2(dhdata,ct);
	int dh= (ct->GetEnergy(1));
	float ds= (float)(dh-dg)/(float)310.15;
	float Tm;
	if (dg <0 && dh <0 && ds <0) {
		Tm= (float)dh/ds;
		if (Tm < 57+273.15)			score[i]++;
	}
	else if(dg>=0) {//for those having no self-structure, dg=0
		score[i]++;
		Tm=0;
	}
	else {
		Tm=0;
	}
	if (Tm==0)	melt[i]=0;
	else melt[i]=Tm-(float)273.15;// this will be used to report melting temperature
    
	//Criteria IV: An 'A' base at position 19 (sense strand) (1 for antisense)
	//+1 point
	if (ct->numseq[1] == 1)		score[i]++;

	//Criteria V: An 'A' base at position 3(sense strand) (17 for antisense)
	//+1 point
	if (ct->numseq[17] == 1)		score[i]++;

	//Criteria VI : A 'U' base at position 10 (sense strand) (10 for antisense)
	//+1 point
	if (ct->numseq[10] ==4)			score[i]++;

	//Criteria VII: A base other than 'G' or 'C' at 19 (sense strand) (1 for antisense)
	//-1 point  
	//deprecated
//	if (ct->numseq[1] == 3 || ct->numseq[1]==2)		score[i]--;

	//Criteria VIII: A base other than 'G' at position 13 (sense strand) (7 for antisense)
	//-1 point
	if (ct->numseq[7] ==3)		score[i]--;
}
//sites specified in oligo() by -test option to test the correlation
else if (test)	 	score[i]=99; 

}

#include "stdafx.h"
#include "dyndotplot.h"
#include "../src/defines.h"
#include "../src/structure.h"
#include "../src/dynalign.h"

using namespace std;

IMPLEMENT_DYNCREATE(CDynDotPlot, CDocument)
CDynDotPlot::CDynDotPlot(void)
{
}

CDynDotPlot::CDynDotPlot(CString Filename, int seq, CRNAstructureApp *app):PlotDoc(app)
{
	//Place all the energies in the array[j][i] array. (Inherited from PlotDoc.)
	//(Note that array is a double.)

	//In the CDynDotPlot, seq is either 1 or 2 to indicate the first or second sequence in the Dynalign calculation
	//A seperate CDynaDotPlot is needed for each sequence.  At this point, the save file is read twice for each plot.

	short i,j,k,l/*,a,b*/;
	bool singleinsert,local/***pair*/;
	short maxsep,gap,lowest;
	//integersize ****v,****w,**w3,**w5,****vmod;
	dynalignarray *w,*vmod;
	varray *v;
	wendarray *w3,*w5;
	datatable *data;
	short int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
	{0,1,0,1,0,0},{0,0,0,0,0,0}};
	int modificationflag;
	short *lowend,*highend;
	bool **allowed_alignments;

	data = new datatable();

	structure *ct1, *ct2;
	ct1 = new structure();
	ct2 = new structure();

	//open the save file to peek at the sizes needed to allocate arrays:
	ifstream sav(Filename.GetBuffer(0),ios::binary);


	read(&sav, &modificationflag);


	//start with structure information
	int length1;
	int length2;
	read(&sav,&(length1));
	read(&sav,&(length2));
	read(&sav,&maxsep);
	sav.close();


	low = INFINITE_ENERGY;
	high = 0;
	
	if (maxsep<0) {
	//This means that the allowed alignments were constrained by the HMM alignment
	  //probs and not by the traditional M parameter.
	  //Therefore, space must be allocated for the allowed_alignment arrays.

	  allowed_alignments=new bool *[length1+1];
	  for (i=0;i<=length1;i++) allowed_alignments[i]=new bool [length2+1];

	}
	else allowed_alignments=NULL;
	//fill the low and highend arrays:
	//allocate space:
	lowend = new short [2*length1];
	highend = new short [2*length1];


	if (modificationflag==1||modificationflag==3) {

		vmod = new dynalignarray();
		

	}
	else vmod=NULL;




	v = new varray();
	w = new dynalignarray();

	

	w5 = new wendarray();
	w3 = new wendarray();


	
    opendynalignsavefile(Filename.GetBuffer(0),ct1,ct2,v,w,vmod,w3,w5,data,&singleinsert,&maxsep,&gap,&lowest,&local,allowed_alignments,lowend,highend);

	//double up the sequences as part of suboptimal traceback support:
	for (i=1;i<=ct1->GetSequenceLength();i++) {
		ct1->numseq[(ct1->GetSequenceLength())+i] = ct1->numseq[i];
	}
	for (i=1;i<ct2->GetSequenceLength();i++) ct2->numseq[ct2->GetSequenceLength()+i]=ct2->numseq[i];
	
	//if (seq==1)ct.GetSequenceLength()=ct1->GetSequenceLength();
	//else ct.GetSequenceLength() = ct2->GetSequenceLength();


	//allocate everything for the plot doc
	xstart = 1;
	ystart = 1;
	
	string seqlabel;
	if (seq==1) {
		xstop = ct1->GetSequenceLength();
		ystop = ct1->GetSequenceLength();
		ct.allocate(ct1->GetSequenceLength());
		for (i=1;i<=ct1->GetSequenceLength();i++) {
			ct.numseq[i]=ct1->numseq[i];
			ct.nucs[i]=ct1->nucs[i];
	
		}
		seqlabel=ct1->GetSequenceLabel();
		//strcpy(ct.ctlabel[1],ct1->ctlabel[1]);


	}
	else {//seq==2 
		xstop = ct2->GetSequenceLength();
		ystop = ct2->GetSequenceLength();
		ct.allocate(ct2->GetSequenceLength());
		for (i=1;i<=ct2->GetSequenceLength();i++) {
			ct.numseq[i]=ct2->numseq[i];
			ct.nucs[i]=ct2->nucs[i];
	
		}
		
		seqlabel=ct2->GetSequenceLabel();
		//strcpy(ct.ctlabel[1],ct2->ctlabel[1]);

	}
	//There is a newline at the end of ctlabel, this needs to be removed
	seqlabel.erase(seqlabel.size()-1,1);
	//ct.ctlabel[1][strlen(ct.ctlabel[1])-1]='\0';
	ct.SetSequenceLabel(seqlabel);

	arrayvalues = new double *[ct.GetSequenceLength()+1];
	for (i=0;i<=ct.GetSequenceLength();i++) {
		arrayvalues[i] = new double [i+1];
	}
	for (i=0;i<=ct.GetSequenceLength();i++) {
		for (j=0;j<=i;j++) {
			arrayvalues[i][j]=(double) INFINITE_ENERGY;
		}
	}


	////
	for (i=1;i<=ct1->GetSequenceLength();i++) {
		for (j=i+minloop;j<=ct1->GetSequenceLength();j++) {
			for (k=max(lowend[i],1);k<=(min(ct2->GetSequenceLength(),highend[i]));k++) {
				for (l=max(lowend[j],k);l<=(min(ct2->GetSequenceLength(),highend[j]));l++) {
					
					//use a and b when refering to the energy arrays
					//a = k-i+maxsep;
					//b = l-j+maxsep;	

					//if ((a+ct2->numofbases-ct1->numofbases>=0)&&(a+ct2->numofbases-ct1->numofbases<=2*maxsep)) {

						if (seq==1) {
						
							arrayvalues[j][i]=min(arrayvalues[j][i],(double) (v->f(i,j,k,l)+v->f(j,i+ct1->GetSequenceLength(),l,k+ct2->GetSequenceLength()))/((PFPRECISION)conversionfactor));
							
							//allow single BP inserts
							if (i>1&&j<ct1->GetSequenceLength()&&k>1&&l<ct1->GetSequenceLength()&&k+1>=lowend[i+1]&&k+1<=highend[i+1]&&k<=highend[i-1]&&k>=lowend[i-1]
								&&l>=lowend[j+1]&&l<=highend[j+1]&&l-1<=highend[j-1]&&l-1>=lowend[j-1]) {
								arrayvalues[j][i]=min(arrayvalues[j][i],(double) (v->f(i+1,j-1,k+1,l-1)+v->f(j+1,(i-1)+ct1->GetSequenceLength(),l,k+ct2->GetSequenceLength())+
									erg1(i-1,j+1,i,j,ct1,data)+erg1(i,j,i+1,j-1,ct1,data)+erg1(k,l,k+1,l-1,ct2,data)+2*gap)/((double)conversionfactor));	
							}
							
							if (arrayvalues[j][i]<low) low = arrayvalues[j][i];



						}
						else { //seq==2
							arrayvalues[l][k]=min(arrayvalues[l][k],(double) (v->f(i,j,k,l)+v->f(j,i+ct1->GetSequenceLength(),l,k+ct2->GetSequenceLength()))/((double)conversionfactor));

							//allow single BP inserts
							if (i>1&&j<ct1->GetSequenceLength()&&k>1&&l<ct2->GetSequenceLength()&&k+1>=lowend[i+1]&&k+1<=highend[i+1]&&k-1<=highend[i]&&k-1>=lowend[i]
								&&l+1>=lowend[j]&&l+1<=highend[j]&&l-1<=highend[j-1]&&l-1>=lowend[j-1]) {
								arrayvalues[l][k]=min(arrayvalues[l][k],(double) (v->f(i+1,j-1,k+1,l-1)+v->f(j,i+ct1->GetSequenceLength(),l+1,k-1+ct2->GetSequenceLength())+
									erg1(k-1,l+1,k,l,ct2,data)+erg1(k,l,k+1,l-1,ct2,data)+erg1(i,j,i+1,j-1,ct1,data)+2*gap)/((double)conversionfactor));	
							}

							if (arrayvalues[l][k]<low) low = arrayvalues[l][k];
						}
						
					//}

				}
			}
		}

	}



	


	delete w3;
	delete w5;
	delete w;
	delete v;

	if (modificationflag) delete vmod;


	
	
	delete lowend;
	delete highend;

	if (maxsep<0) {
	  
		for (i=0;i<=ct1->GetSequenceLength();i++) delete[] allowed_alignments[i];
		delete[] allowed_alignments;
	}

	delete data;
	delete ct1;
	delete ct2;

	originallow = low;
	originalhigh = high;
   
	//select the ranges for each color in display
	colorranges();

	Filename+="  Sequence: ";
	Filename+=ct.GetSequenceLabel().c_str();
	message = " Free Energy (kcal/mol) ";
	SetTitle(Filename);
	outside = "infinity";
}

CDynDotPlot::~CDynDotPlot(void)
{
	int i;
	for (i=0;i<=ct.GetSequenceLength();i++) {
		delete[] arrayvalues[i];
	}
	delete[] arrayvalues;
}

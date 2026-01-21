#include "stdafx.h"
#include "dpdoc.h"
#include "../src/algorithm.h"
#include "../src/defines.h"

using namespace std;

IMPLEMENT_DYNCREATE(DPDoc, CDocument)

DPDoc::DPDoc(void)
{
}

DPDoc::DPDoc(CString Filename, CRNAstructureApp *app, bool istextfile):PlotDoc(app)
{
	//Place all the energies in the array[j][i] array. (Inherited from PlotDoc.)
	//(Note that array is a double.)
	short vers;//short for saving the save file version	
	int i,j;
	double ene;

	low = INFINITE_ENERGY;
	high = 0;

	if (!istextfile) {
		register int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
		{0,1,0,1,0,0},{0,0,0,0,0,0}};
		
		DynProgArray<integersize> *w2,*wmb2;
		integersize *w5,*w3;
		int vmin;
		bool *lfce,*mod;
		datatable *data;

		data = new datatable;

		//peek at the length of the sequence and whether the folding is intermolecular to allocate arrays:
		ifstream sav(Filename.GetBuffer(10),ios::binary);
	
		read(&sav,&vers);
		if (vers!=safiversion) {
			//Add some error handling here at a future time.
			//sav.close();

			//MessageBox("Incorrect save file version. Please refold the sequence first with this version of RNAstructure.");

			//return;
		}
		int sequencelength;
		read(&sav,&(sequencelength));
		read(&sav,&(ct.intermolecular));
		sav.close();

		//allocate everything
		ct.allocate(sequencelength);
		
		xstart = 1;
		xstop = ct.GetSequenceLength();
		ystart = 1;
		ystop = ct.GetSequenceLength();
		
	
		DynProgArray<integersize> w(ct.GetSequenceLength());
		DynProgArray<integersize> v(ct.GetSequenceLength());
		DynProgArray<integersize> wmb(ct.GetSequenceLength());	
		forceclass fce(ct.GetSequenceLength());


		lfce = new bool [2*ct.GetSequenceLength()+1];
		mod = new bool [2*ct.GetSequenceLength()+1];
	
		w5 = new integersize [ct.GetSequenceLength()+1];
		w3 = new integersize [ct.GetSequenceLength()+2];

		if (ct.intermolecular) {
			w2 = new DynProgArray<integersize>(ct.GetSequenceLength());
			wmb2 = new DynProgArray<integersize>(ct.GetSequenceLength());
		
			for (i=0;i<3;i++) read(&sav,&(ct.inter[i]));
		
		}
		else {
			w2 = NULL;
			wmb2 = NULL;
		}

		arrayvalues = new double *[ct.GetSequenceLength()+1];
		for (i=1;i<=ct.GetSequenceLength();i++) {
			arrayvalues[i] = new double [i+1];
		}
		for (i=1;i<=ct.GetSequenceLength();i++) {
			for (j=1;j<=i;j++) {
				arrayvalues[i][j]= (double) 0;
			}
		}
	
		readsav(Filename.GetBuffer(10), &ct, w2, wmb2, w5, w3, lfce, mod, data,
			 &v, &w, &wmb, &fce, &vmin);

		for (i=1;i<=ct.GetSequenceLength();i++) {
			for (j=i+1;j<=ct.GetSequenceLength();j++) {

				arrayvalues[j][i] = ((double) (v.f(i,j)+v.f(j,i+ct.GetSequenceLength())))/( (double) conversionfactor);
				if (arrayvalues[j][i]<low) low = arrayvalues[j][i];
				//if (array[j][i]>high) high = array[j][i];
			}
		}

		delete[] lfce;
		delete[] mod;
	



		delete[] w5;
		delete[] w3;


		if (ct.intermolecular) {
			delete w2;
			delete wmb2;
		}

		delete data;

	}
	else {//This is a plain text file
	
		ifstream *in;
		in = new ifstream();
		in->open(Filename.GetBuffer(0));
		
		int length;
		*in >> length;

		//foo is a string used to temporarily store some info 
		char foo[maxfil];
		
		*in >>foo;//read the "i" header
		*in >>foo;//read the "j" header
		*in >>foo;//read the type of data

		//set message to the type of data
		message = foo;
		
		ct.allocate(length);
		
		for (i=1;i<=ct.GetSequenceLength();i++) ct.nucs[i]='X';
		
		
		//allocate everything
		xstart = 1;
		xstop = ct.GetSequenceLength();
		ystart = 1;
		ystop = ct.GetSequenceLength();

		arrayvalues = new double *[ct.GetSequenceLength()+1];
		for (i=1;i<=ct.GetSequenceLength();i++) {
			arrayvalues[i] = new double [i+1];
		}
		for (i=1;i<=ct.GetSequenceLength();i++) {
			for (j=1;j<=i;j++) {
				arrayvalues[i][j]= INFINITE_ENERGY;
			}
		}

		

		//read the file once to determine maximum i and j;
		*in >> i;
		*in >> j;
		*in >> ene;
		
		
		while(!in->eof()) {
			if (i<j) arrayvalues[j][i]=ene;
			else arrayvalues[i][j]=ene;
			if (arrayvalues[j][i]<low) low = arrayvalues[j][i];
			if (arrayvalues[j][i]>high) high=arrayvalues[j][i];
			*in >> i;
			*in >> j;
			*in >> ene;

			
			

		}

		in->close();
		delete in;

		
	}


	originallow = low;
	originalhigh = high;
   
	//select the ranges for each color in display
	colorranges();

	
	SetTitle(Filename);
	outside = "infinity";
}

DPDoc::~DPDoc(void)
{
	int i;
	for (i=1;i<=ct.GetSequenceLength();i++) {
		delete[] arrayvalues[i];
	}
	delete[] arrayvalues;
}

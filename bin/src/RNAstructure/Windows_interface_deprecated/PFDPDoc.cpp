#include "stdafx.h"
#include "pfdpdoc.h"
#include <math.h>

using namespace std;

IMPLEMENT_DYNCREATE(PFDPDoc, CDocument)

PFDPDoc::PFDPDoc():PlotDoc() {

}

PFDPDoc::PFDPDoc(CString Filename, CRNAstructureApp *app):PlotDoc(app)
{
	filename = Filename;
	PFPRECISION logint;
	int i,j;
	low = 1;
	high = -INFINITE_ENERGY;
	short vers;

	//allocate the ct file by reading the save file:
	ifstream sav(filename.GetBuffer(10),ios::binary);

	
	read(&sav,&(vers));//read the version of the save file
		
	//If the version is not correct, bail out and set v to NULL to indicate the bailout
	if (vers!=pfsaveversion) {
		arrayvalues = NULL;	
		return;
	}
	int slen;
	read(&sav,&(slen));
	
	sav.close();
	//allocate everything
	ct.allocate(slen);
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
			arrayvalues[i][j]= 0;
		}
	}
	
	

	w = new pfunctionclass(ct.GetSequenceLength());
	v = new pfunctionclass(ct.GetSequenceLength());
	wmb = new pfunctionclass(ct.GetSequenceLength());
	wmbl = new pfunctionclass(ct.GetSequenceLength());
	wcoax = new pfunctionclass(ct.GetSequenceLength());
	wl = new pfunctionclass(ct.GetSequenceLength());
	wlc = new pfunctionclass(ct.GetSequenceLength());
	fce = new forceclass(ct.GetSequenceLength());

	w5 = new PFPRECISION [ct.GetSequenceLength()+1];
	w3 = new PFPRECISION [ct.GetSequenceLength()+2];

	lfce = new bool [2*ct.GetSequenceLength()+1];
    mod = new bool [2*ct.GetSequenceLength()+1];

	data = new pfdatatable();

	//load all the data from the pfsavefile:
	readpfsave(filename.GetBuffer(10), &ct, w5, w3,v, w, wmb,wl,wlc, wmbl, wcoax, fce,&scaling,mod,lfce,data);

	//fill array with the values for the plot:
	for (i=1;i<ct.GetSequenceLength();i++) {
		for (j=i+1;j<=ct.GetSequenceLength();j++) {

			if ((v->f(i,j)*v->f(j,i+ct.GetSequenceLength()))>0) {
				logint = calculateprobability(i,j,v,w5,&ct,data,lfce,mod,scaling,fce);
				arrayvalues[j][i] = -log10(logint);

				if (arrayvalues[j][i]<low) low = arrayvalues[j][i];
				if (arrayvalues[j][i]>high) high = arrayvalues[j][i];
			}
			else arrayvalues[j][i] = 1e20;
		}
	}

	originallow = low;
	originalhigh = high;
   
	//select the ranges for each color in display
	colorranges();

	
	message = " -log10(BP Probability) ";
	SetTitle(Filename);

	outside = "infinity";

}

PFDPDoc::~PFDPDoc(void)
{
	int i;

	if (arrayvalues!=NULL) {

		delete w;
		delete v;
		delete wmb;
		delete fce;
		delete[] w5;
		delete[] w3;
		for (i=1;i<=ct.GetSequenceLength();i++) {
			delete[] arrayvalues[i];
		}
		delete[] arrayvalues;
		delete[] lfce;
		delete[] mod;
		delete data;
		delete wmbl;
		delete wl;
		delete wcoax;
	}


}

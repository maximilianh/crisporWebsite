// PlotDoc.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "PlotDoc.h"
#include <math.h>


#define colorbrilliance 220 //describes how bright the dot colors will be

// PlotDoc

IMPLEMENT_DYNCREATE(PlotDoc, CDocument)

PlotDoc::PlotDoc()
{
	

	


}

PlotDoc::PlotDoc(CRNAstructureApp *App) {

	app = App;
	colors = 5;
	rgbcolors=NULL;
	showlegend = true;

}

BOOL PlotDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;
	return TRUE;
}

PlotDoc::~PlotDoc()
{
	delete[] rgbcolors;	


}


BEGIN_MESSAGE_MAP(PlotDoc, CDocument)
END_MESSAGE_MAP()


// PlotDoc diagnostics

#ifdef _DEBUG
void PlotDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void PlotDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG


// PlotDoc serialization

void PlotDoc::Serialize(CArchive& ar)
{
	if (ar.IsStoring())
	{
		// TODO: add storing code here
	}
	else
	{
		// TODO: add loading code here
	}
}

void PlotDoc::colorranges() {
	int i,divisions,per;

	//define the range of values each color denotes
	for (i=1;i<colors;i++) {
		colorrange[i] = low + (high-low)*fabs(((double) i)/((double)colors));
	}
	colorrange[colors]=high+0.000000000001;
	colorrange[0] = low-0.0000000001;
	

	///Define the colors themselves

	if (rgbcolors!=NULL) delete[] rgbcolors;


	rgbcolors = new COLORREF[colors];
	divisions = (colors-3)/2;

	per = colorbrilliance/(divisions+1);

	rgbcolors[0]=RGB(colorbrilliance,0,0);
	for (i=1;i<=divisions;i++) {
		rgbcolors[i] = RGB(colorbrilliance-per*(i),per*(i),0);
	}
	rgbcolors[divisions+1]=RGB(0,colorbrilliance,0);
	for (i=1;i<=divisions;i++) {
		rgbcolors[divisions+1+i] = RGB(0,colorbrilliance-per*(i),per*(i));
	}
	rgbcolors[colors-1]=RGB(0,0,colorbrilliance);
}



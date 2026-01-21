#pragma once
#include "../src/structure.h"
#define maxcolors 17
#include "RNAstructure.h"

// PlotDoc document

class PlotDoc : public CDocument
{
	DECLARE_DYNCREATE(PlotDoc)

public:
	PlotDoc();
	PlotDoc(CRNAstructureApp *App);
	virtual ~PlotDoc();
	virtual void Serialize(CArchive& ar);   // overridden for document i/o
	int xstart,xstop,ystart,ystop;
	double **arrayvalues;
	structure ct;
	int colors;
	COLORREF *rgbcolors;
	double colorrange[maxcolors+1];
	double low, high,originallow,originalhigh;
	CRNAstructureApp *app;
	CString message;
	bool showlegend;
	CString outside;

	void colorranges();

#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:
	virtual BOOL OnNewDocument();

	DECLARE_MESSAGE_MAP()
};

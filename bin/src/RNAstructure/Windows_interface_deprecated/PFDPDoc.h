#pragma once
#include "plotdoc.h"
#include "../src/pfunction.h"
#include "rnastructure.h"
#include "../src/phmm/utils/xmath/log/xlog_math.h"

class PFDPDoc :
	public PlotDoc
{
	DECLARE_DYNCREATE(PFDPDoc)
public:
	PFDPDoc(CString Filename, CRNAstructureApp *app);
	PFDPDoc();
	~PFDPDoc(void);

	CString filename;
	PFPRECISION *w5, *w3, scaling;
	pfunctionclass *v, *w,*wmb,*wmbl,*wl,*wcoax;
	forceclass *fce;
	bool *mod,*lfce;
	pfdatatable *data;

};

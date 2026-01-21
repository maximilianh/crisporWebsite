#pragma once
#include "plotdoc.h"

class DPDoc :
	public PlotDoc
{
	DECLARE_DYNCREATE(DPDoc)
public:
	DPDoc(void);

	//if istextfile is true, then the file isn't as savefile, but a plain text file
	DPDoc(CString Filename, CRNAstructureApp *app, bool istextfile=false);
	~DPDoc(void);
};

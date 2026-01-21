#pragma once
#include "plotdoc.h"
#include "../src/phmm/utils/xmath/log/xlog_math.h"


class CDynDotPlot :
	public PlotDoc
{
	DECLARE_DYNCREATE(CDynDotPlot)
public:
	CDynDotPlot(void);
	CDynDotPlot(CString Filename, int seq, CRNAstructureApp *app);
	~CDynDotPlot(void);
};

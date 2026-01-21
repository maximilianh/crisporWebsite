#pragma once
#include "plot.h"

class CBoxPlotView :
	public CPlot
{
	DECLARE_DYNCREATE(CBoxPlotView)
public:
	CBoxPlotView(void);
	~CBoxPlotView(void);
	DECLARE_MESSAGE_MAP()
	afx_msg void OnOutputOutputprobablestructure();
	afx_msg void OnOutputOutputtextfile();
	afx_msg void OnOutputStochastictraceback();
};

#pragma once

#include "plotdoc.h"

// CPlotRange dialog

class CPlotRange : public CDialog
{
	DECLARE_DYNAMIC(CPlotRange)

public:
	CPlotRange(PlotDoc *Plotdoc =NULL, CWnd* pParent = NULL);   // standard constructor
	virtual ~CPlotRange();

// Dialog Data
	enum { IDD = IDD_PLOTRANGE };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	
	afx_msg void OnBnClickedResetrange();
	afx_msg void OnBnClickedOk();
	PlotDoc *plotdoc;
	CWnd *parent;

	
	double m_min;
	double m_max;
};

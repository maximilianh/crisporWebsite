// PlotRange.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "PlotRange.h"


// CPlotRange dialog

IMPLEMENT_DYNAMIC(CPlotRange, CDialog)
CPlotRange::CPlotRange(PlotDoc *Plotdoc, CWnd* pParent /*=NULL*/)
	: CDialog(CPlotRange::IDD, pParent)
		
	
{
	if(Plotdoc!=NULL) {
		plotdoc = Plotdoc;
		m_min=plotdoc->low;
		m_max=plotdoc->high;
		parent = pParent;
	}
	
}

CPlotRange::~CPlotRange()
{
}

void CPlotRange::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);


	DDX_Text(pDX, IDC_EDITMIN, m_min);
	DDX_Text(pDX, IDC_EDITMAX, m_max);
}


BEGIN_MESSAGE_MAP(CPlotRange, CDialog)
	ON_BN_CLICKED(IDC_RESETRANGE, OnBnClickedResetrange)
	ON_BN_CLICKED(IDOK, OnBnClickedOk)
END_MESSAGE_MAP()


// CPlotRange message handlers

void CPlotRange::OnBnClickedResetrange()
{
	m_min = plotdoc->originallow;
	m_max = plotdoc->originalhigh;
	UpdateData(FALSE);
}

void CPlotRange::OnBnClickedOk()
{
	UpdateData(TRUE);
	plotdoc->low = m_min;
	plotdoc->high = m_max;
	plotdoc->colorranges();
	parent->Invalidate();
	OnOK();
}

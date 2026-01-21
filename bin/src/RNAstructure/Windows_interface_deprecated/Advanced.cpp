// Advanced.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "Advanced.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CAdvanced dialog


CAdvanced::CAdvanced(short *Dynnumber, short *Dynpercent, short *Dynwindow, CWnd* pParent /*=NULL*/)
	: CDialog(CAdvanced::IDD, pParent)
{
	//{{AFX_DATA_INIT(CAdvanced)
	m_dnumber = Dynnumber;
	m_dwindow = Dynwindow;
	m_dpercent = Dynpercent;
	//}}AFX_DATA_INIT
	dynnumber = Dynnumber;
	dynpercent = Dynpercent;
	dynwindow = Dynwindow;
	savenumber = *Dynnumber;
	savewindow = *Dynwindow;
	savepercent = *Dynpercent;

}


void CAdvanced::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CAdvanced)
	DDX_Text(pDX, IDC_ADVANCED_NUMBER, *m_dnumber);
	DDV_MinMaxInt(pDX, *m_dnumber, 1, 1000);
	DDX_Text(pDX, IDC_ADVANCED_WINDOW, *m_dwindow);
	DDV_MinMaxInt(pDX, *m_dwindow, 0, 32000);
	DDX_Text(pDX, IDC_ADVANCED_PERCENT, *m_dpercent);
	DDV_MinMaxInt(pDX, *m_dpercent, 0, 99);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CAdvanced, CDialog)
	//{{AFX_MSG_MAP(CAdvanced)
	ON_EN_CHANGE(IDC_ADVANCED_PERCENT, OnChangeParameter)
	ON_EN_CHANGE(IDC_ADVANCED_NUMBER, OnChangeParameter)
	ON_EN_CHANGE(IDC_ADVANCED_WINDOW, OnChangeParameter)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CAdvanced message handlers

void CAdvanced::OnOK() 
{
	

	CDialog::OnOK();
}

void CAdvanced::OnChangeParameter() {

	UpdateData(TRUE);

}

void CAdvanced::OnCancel() 
{
	//reset the values
	*dynnumber = savenumber;
	*dynwindow = savewindow;
	*dynpercent = savepercent;

	
	CDialog::OnCancel();
}


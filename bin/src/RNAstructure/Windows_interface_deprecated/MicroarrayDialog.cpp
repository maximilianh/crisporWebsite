// MicroarrayDialog.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "MicroarrayDialog.h"


// CMicroarrayDialog dialog

IMPLEMENT_DYNAMIC(CMicroarrayDialog, CDialog)
CMicroarrayDialog::CMicroarrayDialog(CWnd* pParent /*=NULL*/)
	: CDialog(CMicroarrayDialog::IDD, pParent)
	, m_start(0)
	, m_stop(0)
	, m_min(0)
{
}

CMicroarrayDialog::CMicroarrayDialog(structure * CT, CWnd* pParent /*=NULL*/)
	: CDialog(CMicroarrayDialog::IDD, pParent)
	, m_start(0)
	, m_stop(0)
	, m_min(0)
{
	ct = CT;
}

CMicroarrayDialog::~CMicroarrayDialog()
{
}

void CMicroarrayDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_EDIT1, m_start);
	DDX_Text(pDX, IDC_EDIT2, m_stop);
	DDX_Text(pDX, IDC_EDIT3, m_min);
}


BEGIN_MESSAGE_MAP(CMicroarrayDialog, CDialog)
	ON_BN_CLICKED(IDOK, OnBnClickedOk)
END_MESSAGE_MAP()


// CMicroarrayDialog message handlers

void CMicroarrayDialog::OnBnClickedOk()
{
	UpdateData(TRUE);
	ct->microstart[ct->nmicroarray]=m_start;
	ct->microstop[ct->nmicroarray]=m_stop;
	ct->microunpair[ct->nmicroarray]=m_min;

	ct->nmicroarray++;	


	OnOK();
}

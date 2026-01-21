// InternalDialog.cpp : implementation file
//

//Dialog box to allow user to change the maximum size of bulge/internal loops

#include "stdafx.h"
#include "RNAstructure.h"
#include "InternalDialog.h"


// CInternalDialog dialog

IMPLEMENT_DYNAMIC(CInternalDialog, CDialog)
CInternalDialog::CInternalDialog(int *Maxinternal, CWnd* pParent /*=NULL*/)
	: CDialog(CInternalDialog::IDD, pParent)
{
	maxinternal = Maxinternal;
	m_max=*maxinternal;
}

CInternalDialog::~CInternalDialog()
{
}

void CInternalDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_MAXSIZEEDIT, m_max);
}


BEGIN_MESSAGE_MAP(CInternalDialog, CDialog)
	ON_BN_CLICKED(IDOK, OnButtonOK)
	ON_BN_CLICKED(IDCANCEL, OnButtonCancel)
END_MESSAGE_MAP()


// CInternalDialog message handlers
void CInternalDialog::OnButtonOK()
{
	UpdateData(TRUE);
	*maxinternal = m_max;
	CDialog::OnOK();
	
}

void CInternalDialog::OnButtonCancel()
{
	CDialog::OnCancel();
	
}
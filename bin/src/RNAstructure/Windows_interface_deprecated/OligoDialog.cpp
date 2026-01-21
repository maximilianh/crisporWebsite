// OligoDialog.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "OligoDialog.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// OligoDialog dialog


COligoDialog::COligoDialog(CWnd* pParent /*=NULL*/)
	: CDialog(COligoDialog::IDD, pParent)
{
	//{{AFX_DATA_INIT(OligoDialog)
		// NOTE: the ClassWizard will add member initialization here
	//}}AFX_DATA_INIT
}


void COligoDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(OligoDialog)
		// NOTE: the ClassWizard will add DDX and DDV calls here
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(COligoDialog, CDialog)
	//{{AFX_MSG_MAP(OligoDialog)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// OligoDialog message handlers

void COligoDialog::OnCancel() 
{
	// TODO: Add extra cleanup here
	
	CDialog::OnCancel();
}

void COligoDialog::OnOK() 
{
	// TODO: Add extra validation here
	
	CDialog::OnOK();
}



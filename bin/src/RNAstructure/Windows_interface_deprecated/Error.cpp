// Error.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "Error.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CError dialog

//This dialog now only reports traceback errors.  And should not be used for anything else.


CError::CError(CWnd* pParent /*=NULL*/)
	: CDialog(CError::IDD, pParent)
{
	//{{AFX_DATA_INIT(CError)
		// NOTE: the ClassWizard will add member initialization here
	//}}AFX_DATA_INIT
}


void CError::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CError)
		// NOTE: the ClassWizard will add DDX and DDV calls here
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CError, CDialog)
	//{{AFX_MSG_MAP(CError)
		// NOTE: the ClassWizard will add message map macros here
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CError message handlers

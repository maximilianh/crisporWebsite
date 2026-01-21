
// TProgressDialog.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "TProgressDialog.h"


#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif



/////////////////////////////////////////////////////////////////////////////
// TProgressDialog




//caption is a pointer to the window's caption
TProgressDialog::TProgressDialog()
{
	CDialog::CDialog(IDD_PROGRESSDIALOG);
	

}



TProgressDialog::~TProgressDialog()
{

	//delete progressbar;
}


//BEGIN_MESSAGE_MAP(TProgressDialog, CDialog)
	//{{AFX_MSG_MAP(TProgressDialog)
	
	//}}AFX_MSG_MAP
//END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// TProgressDialog message handlers

//call this function with the percent completed
void TProgressDialog::update(int frac) {
	CProgressCtrl *progressbar;
	progressbar = (CProgressCtrl*) GetDlgItem( IDC_PROGRESSBAR );	
	progressbar->SetPos(frac);
	UpdateWindow();
}


//int TProgressDialog::OnCreate(LPCREATESTRUCT lpCreateStruct) 
//{
	//if (CMDIChildWnd::OnCreate(lpCreateStruct) == -1) return -1;
	//progressbar = new CProgressCtrl();
	//progressbar->Create(WS_CHILD|WS_VISIBLE|WS_BORDER|PBS_SMOOTH,CRect(20,60,200,110),this,IDC_PROGRESSBAR);
	//progressbar->Create(WS_CHILD|WS_VISIBLE|WS_BORDER|PBS_SMOOTH,CRect(5,15,180,45),this,IDC_PROGRESSBAR);
	//progressbar->SetRange(0,100);--default
	//progressbar->SetPos(0);--default
	//progressbar->SetStep(1);
//	return 0;
//}
	



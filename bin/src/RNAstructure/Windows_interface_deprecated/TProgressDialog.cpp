
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

IMPLEMENT_DYNCREATE(TProgressDialog, CMDIChildWnd)


//caption is a pointer to the window's caption
TProgressDialog::TProgressDialog()
{
	CMDIChildWnd::CMDIChildWnd();
	

}



TProgressDialog::~TProgressDialog()
{

	delete progressbar;
}


BEGIN_MESSAGE_MAP(TProgressDialog, CMDIChildWnd)
	//{{AFX_MSG_MAP(TProgressDialog)
	ON_WM_CREATE()
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// TProgressDialog message handlers

//call this function with the percent completed
void TProgressDialog::update(int frac) {
		
	progressbar->SetPos(frac);
	UpdateWindow();
}


int TProgressDialog::OnCreate(LPCREATESTRUCT lpCreateStruct) 
{
	if (CMDIChildWnd::OnCreate(lpCreateStruct) == -1) return -1;
	progressbar = new CProgressCtrl();
	
	progressbar->Create(WS_CHILD|WS_VISIBLE|WS_BORDER|PBS_SMOOTH,CRect(5,15,180,45),this,IDC_PROGRESSBAR);
	
	progressbar->SetStep(1);
	return 0;
}
	



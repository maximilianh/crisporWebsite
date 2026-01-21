// Splash.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "Splash.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif








//This is the backend functions, set up to run in a seperate thread
UINT TimerProc( LPVOID pParam ){    
	CTimeObject* pObject = (CTimeObject*) pParam;
	CSplash* psplash = (CSplash*) pObject->splash; 
    
	
	Sleep(5000);
	
	if (psplash->m_hWnd!=0)
	::PostMessage(psplash->m_hWnd,ID_DONE,0,0);

	//::PostMessage(psplash->m_hWnd,WM_CLOSE,0,0);
	
	

	return 0;   // thread completed successfully
}








/////////////////////////////////////////////////////////////////////////////
// CSplash dialog
IMPLEMENT_DYNAMIC(CSplash, CDialog)
CSplash::CSplash(CWnd* pParent /*=NULL*/)
	: CDialog(CSplash::IDD, pParent)
{
	parent = pParent;
	//{{AFX_DATA_INIT(CSplash)
		// NOTE: the ClassWizard will add member initialization here
	//}}AFX_DATA_INIT
	
}

CSplash::~CSplash(){
	//CDialog::~CDialog();//Removed 3/31/11
}


void CSplash::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CSplash)
		// NOTE: the ClassWizard will add DDX and DDV calls here
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CSplash, CDialog)
	ON_MESSAGE(ID_DONE, Done)
	//{{AFX_MSG_MAP(CSplash)
	ON_WM_LBUTTONDBLCLK()
	ON_WM_LBUTTONDOWN()
	
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CSplash message handlers
BOOL CSplash::OnInitDialog( ) {
	
	//BOOL rtn = CDialog::OnInitDialog();

	


	time.splash = this;
	AfxBeginThread(TimerProc,&time);
	///return rtn;
	return CDialog::OnInitDialog();
}

LRESULT CSplash::Done(WPARAM wParam, LPARAM lParam) {

	//DestroyWindow();
	//SendMessage(WM_CLOSE);
	
	//delete this;//with Erol


	OnOK();
	return 0;

}

void CSplash::OnLButtonDblClk(UINT nFlags, CPoint point) 
{
	
	CDialog::OnLButtonDblClk(nFlags, point);
	// Hide the dialog window:
	//DestroyWindow();
	//SendMessage(WM_CLOSE);
	
	//OnOK();
	ShowWindow(SW_HIDE);
	

		
	

	
}

void CSplash::OnLButtonDown(UINT nFlags, CPoint point) 
{

	CDialog::OnLButtonDown(nFlags, point);
	// Hide the dialog window:
	//DestroyWindow();
	//SendMessage(WM_CLOSE);
	//OnOK();
	ShowWindow(SW_HIDE);
	
	
}

void CSplash::OnOK() {
	
	//CDialog::OnOK();
	//CWnd::DestroyWindow();
	//this->SendMessage(WM_CLOSE);
	
	//parent->SendMessage(ID_KILLSPLASH); -This doesn't work either


	DestroyWindow();
	//delete this;//deleted 3/20/2011 by DHM - NOT needed in current Windows??
}


void CSplash::PostNcDestroy()
{
	delete this;
}
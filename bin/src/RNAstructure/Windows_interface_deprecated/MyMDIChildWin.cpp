// MyMDIChildWin.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "MyMDIChildWin.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CMyMDIChildWin

IMPLEMENT_DYNCREATE(CMyMDIChildWin, CMDIChildWnd)

CMyMDIChildWin::CMyMDIChildWin()
{
}

CMyMDIChildWin::~CMyMDIChildWin()
{
}


/*afx_msg int CMyMDIChildWin::OnCreate(LPCREATESTRUCT lpCreateStruct) {
	

	lpCreateStruct->style=
		WS_CHILD;// | WS_VISIBLE | WS_OVERLAPPED|WS_CAPTION|WS_SYSMENU|WS_MINIMIZEBOX;

	return CMDIChildWnd::OnCreate(lpCreateStruct);

}*/

BOOL CMyMDIChildWin::PreCreateWindow(CREATESTRUCT& cs)
{
	cs.style &= ~WS_THICKFRAME;
	cs.style |= WS_BORDER;

	//Use some large dimensions so that the window can snap down to the right size when the view requests ResizetoFit()
	cs.cx=768;
	cs.cy=1024;

// remove the minimize and maximize buttons
// so that the MDI child frame "snaps" to the property pages.

	cs.style &= ~(WS_MAXIMIZEBOX);

	return CMDIChildWnd::PreCreateWindow(cs);
}


BEGIN_MESSAGE_MAP(CMyMDIChildWin, CMDIChildWnd)
	//{{AFX_MSG_MAP(CMyMDIChildWin)
		// NOTE - the ClassWizard will add and remove mapping macros here.
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CMyMDIChildWin message handlers

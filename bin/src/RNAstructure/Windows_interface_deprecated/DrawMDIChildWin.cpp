// DrawMDIChildWin.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "DrawMDIChildWin.h"
#include <Afxwin.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CDrawMDIChildWin

IMPLEMENT_DYNCREATE(CDrawMDIChildWin, CMDIChildWnd)

CDrawMDIChildWin::CDrawMDIChildWin()
{
}

CDrawMDIChildWin::~CDrawMDIChildWin()
{
}

BOOL CDrawMDIChildWin::PreCreateWindow(CREATESTRUCT& cs)
{
	//cs.style &= ~WS_THICKFRAME;
	//cs.style |= WS_BORDER;

// remove the minimize and maximize buttons
// so that the MDI child frame "snaps" to the property pages.

	//cs.style &= ~(WS_MAXIMIZEBOX);

	return CMDIChildWnd::PreCreateWindow(cs);
}
BOOL CDrawMDIChildWin::Create( LPCTSTR lpszClassName, LPCTSTR lpszWindowName, DWORD dwStyle = WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW, const RECT& rect = rectDefault, CMDIFrameWnd* pParentWnd = NULL, CCreateContext* pContext = NULL ) {



	
	lpszWindowName="Draw Window Check";
	return CMDIChildWnd::Create( lpszClassName, lpszWindowName,dwStyle , rect , pParentWnd, pContext);
}




BEGIN_MESSAGE_MAP(CDrawMDIChildWin, CMDIChildWnd)
	//{{AFX_MSG_MAP(CDrawMDIChildWin)
		// NOTE - the ClassWizard will add and remove mapping macros here.
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDrawMDIChildWin message handlers

// DrawFrame.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "DrawFrame.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CDrawFrame

IMPLEMENT_DYNCREATE(CDrawFrame, CFrameWnd)

CDrawFrame::CDrawFrame()
{
}


/*BOOL Create( LPCTSTR lpszClassName, LPCTSTR lpszWindowName, DWORD dwStyle, 
			RECT& rect, CWnd* pParentWnd, UINT nID, CCreateContext* pContext) {



	return CFrameWnd::Create(lpszClassName,lpszWindowName,dwStyle, 
			rect, pParentWnd,nID, pContext) {

}*/

CDrawFrame::~CDrawFrame()
{
}


BEGIN_MESSAGE_MAP(CDrawFrame, CFrameWnd)
	//{{AFX_MSG_MAP(CDrawFrame)
		// NOTE - the ClassWizard will add and remove mapping macros here.
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDrawFrame message handlers

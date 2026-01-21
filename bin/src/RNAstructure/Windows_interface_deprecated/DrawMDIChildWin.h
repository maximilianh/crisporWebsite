#if !defined(AFX_DRAWMDICHILDWIN_H__9C9028E2_9A0E_11D4_9F32_00C0F02A5F5D__INCLUDED_)
#define AFX_DRAWMDICHILDWIN_H__9C9028E2_9A0E_11D4_9F32_00C0F02A5F5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// DrawMDIChildWin.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CDrawMDIChildWin frame

class CDrawMDIChildWin : public CMDIChildWnd
{
	DECLARE_DYNCREATE(CDrawMDIChildWin)
protected:
	CDrawMDIChildWin();           // protected constructor used by dynamic creation

// Attributes
public:
	
// Operations
public:
	BOOL PreCreateWindow(CREATESTRUCT& cs);
	BOOL Create(  LPCTSTR lpszClassName, LPCTSTR lpszWindowName, DWORD dwStyle , const RECT& rect, CMDIFrameWnd* pParentWnd , CCreateContext* pContext );



// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CDrawMDIChildWin)
	//}}AFX_VIRTUAL

// Implementation
protected:
	virtual ~CDrawMDIChildWin();

	// Generated message map functions
	//{{AFX_MSG(CDrawMDIChildWin)
		// NOTE - the ClassWizard will add and remove member functions here.
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_DRAWMDICHILDWIN_H__9C9028E2_9A0E_11D4_9F32_00C0F02A5F5D__INCLUDED_)

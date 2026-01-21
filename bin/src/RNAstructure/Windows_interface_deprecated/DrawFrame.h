#if !defined(AFX_DRAWFRAME_H__221D28C2_B68C_11D4_9F32_00C0F02A5F5D__INCLUDED_)
#define AFX_DRAWFRAME_H__221D28C2_B68C_11D4_9F32_00C0F02A5F5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// DrawFrame.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CDrawFrame frame

class CDrawFrame : public CFrameWnd
{
	DECLARE_DYNCREATE(CDrawFrame)
protected:
	CDrawFrame();           // protected constructor used by dynamic creation

// Attributes
public:

	//BOOL Create(LPCTSTR lpszClassName, LPCTSTR lpszWindowName, DWORD dwStyle, 
	//	RECT& rect, CWnd* pParentWnd, UINT nID, CCreateContext* pContext = NULL); 

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CDrawFrame)
	//}}AFX_VIRTUAL

// Implementation
protected:
	virtual ~CDrawFrame();

	// Generated message map functions
	//{{AFX_MSG(CDrawFrame)
		// NOTE - the ClassWizard will add and remove member functions here.
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_DRAWFRAME_H__221D28C2_B68C_11D4_9F32_00C0F02A5F5D__INCLUDED_)

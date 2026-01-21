#if !defined(AFX_MYMDICHILDWIN_H__CD74BB02_8F97_11D4_9F32_00C0F02A5F5D__INCLUDED_)
#define AFX_MYMDICHILDWIN_H__CD74BB02_8F97_11D4_9F32_00C0F02A5F5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// MyMDIChildWin.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CMyMDIChildWin frame

class CMyMDIChildWin : public CMDIChildWnd
{
	DECLARE_DYNCREATE(CMyMDIChildWin)
protected:
	CMyMDIChildWin();           // protected constructor used by dynamic creation

// Attributes
public:

// Operations
public:
	BOOL PreCreateWindow(CREATESTRUCT& cs);
	//afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CMyMDIChildWin)
	//}}AFX_VIRTUAL

// Implementation
protected:
	virtual ~CMyMDIChildWin();

	// Generated message map functions
	//{{AFX_MSG(CMyMDIChildWin)
		// NOTE - the ClassWizard will add and remove member functions here.
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_MYMDICHILDWIN_H__CD74BB02_8F97_11D4_9F32_00C0F02A5F5D__INCLUDED_)

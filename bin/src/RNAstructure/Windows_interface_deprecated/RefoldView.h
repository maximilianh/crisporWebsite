#if !defined(AFX_REFOLDVIEW_H__D1181296_C4FF_47A2_9B77_27437A9A155D__INCLUDED_)
#define AFX_REFOLDVIEW_H__D1181296_C4FF_47A2_9B77_27437A9A155D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// RefoldView.h : header file
//

#include "../src/structure.h"
#include "../src/algorithm.h"


/////////////////////////////////////////////////////////////////////////////
// CRefoldView form view

#ifndef __AFXEXT_H__
#include <afxext.h>
#endif

class CRefoldView : public CFormView
{
protected:
	CRefoldView();           // protected constructor used by dynamic creation
	DECLARE_DYNCREATE(CRefoldView)

// Form Data
public:
	//{{AFX_DATA(CRefoldView)
	enum { IDD = IDD_REFOLDVIEW_FORM };
	CString	m_ctname;
	int		m_number;
	int		m_percent;
	CString	m_savename;
	int		m_window;
	//}}AFX_DATA
	char *startpath;

	
	

// Attributes
public:

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CRefoldView)
	public:
	virtual void OnInitialUpdate();
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	virtual ~CRefoldView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
	//{{AFX_MSG(CRefoldView)
	afx_msg void OnSavebutton();
	afx_msg void OnCtbutton();
	afx_msg void OnStart();
	afx_msg void OnForceCurrent();
	afx_msg void OnForceSaveconstraints();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnEnChangeCtname();
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_REFOLDVIEW_H__D1181296_C4FF_47A2_9B77_27437A9A155D__INCLUDED_)

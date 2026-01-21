#if !defined(AFX_GODIALOG_H__D51E2044_5C53_46A2_973E_7E91955B9534__INCLUDED_)
#define AFX_GODIALOG_H__D51E2044_5C53_46A2_973E_7E91955B9534__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// GoDialog.h : header file
//

#include "oligoobject.h"


/////////////////////////////////////////////////////////////////////////////
// CGoDialog dialog

class CGoDialog : public CDialog
{
// Construction
public:
	CGoDialog(int *Current, COligoObject *Oligoobject, CWnd* pParent /*=NULL*/);   // standard constructor
	BOOL OnInitDialog();
	int *current;
	COligoObject *oligoobject;

	CSpinButtonCtrl m_Spin;
	afx_msg void OnDeltaposSpin(NMHDR* pNMHDR, LRESULT* pResult);
	CWnd* parent;
// Dialog Data
	//{{AFX_DATA(CGoDialog)
	enum { IDD = IDD_GODIALOG };
	int		m_current;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CGoDialog)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CGoDialog)
	virtual void OnOK();
	afx_msg void OnMoststable();
	afx_msg void OnChangeNumber();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_GODIALOG_H__D51E2044_5C53_46A2_973E_7E91955B9534__INCLUDED_)

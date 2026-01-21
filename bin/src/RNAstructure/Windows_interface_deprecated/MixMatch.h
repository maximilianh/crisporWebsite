#if !defined(AFX_MIXMATCH_H__A37E0189_0B3A_496E_BBA4_8F4FEE9539A5__INCLUDED_)
#define AFX_MIXMATCH_H__A37E0189_0B3A_496E_BBA4_8F4FEE9539A5__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// MixMatch.h : header file
//
#include "FoldDoc.h"
/////////////////////////////////////////////////////////////////////////////
// CMixMatch dialog
struct CMixMatchObject
{
	TProgressDialog *progress;
	char *ctinfile;
	datatable *data;
    CFormView *parent;
	char *ctoutfile;
	char *modfile;
    

};

class CMixMatch : public CFormView
{
// Construction
public:
	CMixMatch();   // standard constructor
	CFoldDoc *GetFoldDocument();
	TProgressDialog *progress;
	CMixMatchObject object;
	LRESULT Done(WPARAM wParam, LPARAM lParam);

// Dialog Data
	//{{AFX_DATA(CMixMatch)
	enum { IDD = IDD_MIXMATCH };
	CString	m_ctinname;
	CString	m_ctoutname;
	CString	m_modname;
	//}}AFX_DATA

	void OnInitialUpdate();

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CMixMatch)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	DECLARE_DYNCREATE(CMixMatch)
	// Generated message map functions
	//{{AFX_MSG(CMixMatch)
	afx_msg void OnInct();
	afx_msg void OnMod();
	afx_msg void OnCtout();
	virtual void OnOK();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_MIXMATCH_H__A37E0189_0B3A_496E_BBA4_8F4FEE9539A5__INCLUDED_)

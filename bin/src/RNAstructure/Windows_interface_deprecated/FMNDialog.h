#if !defined(AFX_FMNDIALOG_H__184255A1_97F1_11D4_9F32_00C0F02A5F5D__INCLUDED_)
#define AFX_FMNDIALOG_H__184255A1_97F1_11D4_9F32_00C0F02A5F5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// FMNDialog.h : header file
//
#include "../src/structure.h"

/////////////////////////////////////////////////////////////////////////////
// CFMNDialog dialog

class CFMNDialog : public CDialog
{
// Construction
public:
	CFMNDialog(structure *CT, CWnd* pParent = NULL);   // standard constructor
	structure *ct;

// Dialog Data
	//{{AFX_DATA(CFMNDialog)
	enum { IDD = IDD_FMN_DIALOG };
	short	m_fmn;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CFMNDialog)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CFMNDialog)
	afx_msg void OnChangeFmn();
	
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedOk();
	
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_FMNDIALOG_H__184255A1_97F1_11D4_9F32_00C0F02A5F5D__INCLUDED_)

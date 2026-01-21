#if !defined(AFX_DOUBLEDIALOG_H__184255A0_97F1_11D4_9F32_00C0F02A5F5D__INCLUDED_)
#define AFX_DOUBLEDIALOG_H__184255A0_97F1_11D4_9F32_00C0F02A5F5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// DoubleDialog.h : header file
//
#include "../src/structure.h"

/////////////////////////////////////////////////////////////////////////////
// CDoubleDialog dialog

class CDoubleDialog : public CDialog
{
// Construction
public:
	CDoubleDialog(structure *CT, CWnd* pParent = NULL);   // standard constructor
	structure *ct;
// Dialog Data
	//{{AFX_DATA(CDoubleDialog)
	enum { IDD = IDD_DOUBLE_DIALOG };
	short	m_double;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CDoubleDialog)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CDoubleDialog)
	afx_msg void OnChangeNumberDouble();
	virtual void OnCancel();
	virtual void OnOK();
	afx_msg void OnOk2();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_DOUBLEDIALOG_H__184255A0_97F1_11D4_9F32_00C0F02A5F5D__INCLUDED_)

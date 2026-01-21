#if !defined(AFX_PAIRDIALOG_H__184255A2_97F1_11D4_9F32_00C0F02A5F5D__INCLUDED_)
#define AFX_PAIRDIALOG_H__184255A2_97F1_11D4_9F32_00C0F02A5F5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// PairDialog.h : header file
//

#include "../src/structure.h"

/////////////////////////////////////////////////////////////////////////////
// CPairDialog dialog

class CPairDialog : public CDialog
{
// Construction
public:
	CPairDialog(structure *CT, CWnd* pParent = NULL, bool Isforbid = false);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CPairDialog)
	enum { IDD = IDD_PAIR_DIALOG };
	short	m_base1;
	short	m_base2;
	short	m_length;
	//}}AFX_DATA
	structure *ct;
	bool isforbid; //true if the user is forbidding pairs
	void getinfo();
	BOOL OnInitDialog( ); 
	

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CPairDialog)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CPairDialog)
	afx_msg void OnChangeBase1();
	afx_msg void OnChangeBase2();
	afx_msg void OnChangeLength();
	virtual void OnOK();
	afx_msg void OnOk2();
	
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedOk3();
	
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_PAIRDIALOG_H__184255A2_97F1_11D4_9F32_00C0F02A5F5D__INCLUDED_)

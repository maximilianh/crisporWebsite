#if !defined(AFX_SINGLEDIALOG_H__9AD39FA0_97FC_11D4_9F32_00C0F02A5F5D__INCLUDED_)
#define AFX_SINGLEDIALOG_H__9AD39FA0_97FC_11D4_9F32_00C0F02A5F5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "../src/structure.h"

// SingleDialog.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CSingleDialog dialog

class CSingleDialog : public CDialog
{
// Construction
public:
	CSingleDialog(structure *ct, CWnd* pParent = NULL, int Mode = 1);   // standard constructor
	structure *ct;
// Dialog Data
	//{{AFX_DATA(CSingleDialog)
	enum { IDD = IDD_SINGLE_DIALOG };
	short	m_single;
	//}}AFX_DATA
	int mode;
	BOOL OnInitDialog( );


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CSingleDialog)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CSingleDialog)
	afx_msg void OnChangeNumberSingle();
	virtual void OnOK();
	virtual void OnCancel();
	afx_msg void OnOk2();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_SINGLEDIALOG_H__9AD39FA0_97FC_11D4_9F32_00C0F02A5F5D__INCLUDED_)

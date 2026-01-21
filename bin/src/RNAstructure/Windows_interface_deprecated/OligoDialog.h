#if !defined(AFX_OLIGODIALOG_H__6B71FF76_980E_4B08_AECF_ABC3EDB3E6C8__INCLUDED_)
#define AFX_OLIGODIALOG_H__6B71FF76_980E_4B08_AECF_ABC3EDB3E6C8__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// OligoDialog.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// OligoDialog dialog

class COligoDialog : public CDialog
{
// Construction
public:
	COligoDialog(CWnd* pParent = NULL);   // standard constructor
	
// Dialog Data
	//{{AFX_DATA(OligoDialog)
	enum { IDD = IDD_OLIGODIALOG };
		// NOTE: the ClassWizard will add data members here
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(OligoDialog)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(OligoDialog)
	virtual void OnCancel();
	virtual void OnOK();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()

	
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_OLIGODIALOG_H__6B71FF76_980E_4B08_AECF_ABC3EDB3E6C8__INCLUDED_)

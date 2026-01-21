#if !defined(AFX_TPROGRESSDIALOG_H__2A024480_F2C0_11D2_9F31_00C0F02A5F5D__INCLUDED_)
#define AFX_TPROGRESSDIALOG_H__2A024480_F2C0_11D2_9F31_00C0F02A5F5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// TProgressDialog.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// TProgressDialog frame

class TProgressDialog : public CDialog
{
	
public:
	TProgressDialog(); 
	void update(int frac);
	~TProgressDialog();
//	CProgressCtrl *progressbar;

// Attributes
public:

// Operations
public:

private:
	

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(TProgressDialog)
	public:
	
	//}}AFX_VIRTUAL

// Implementation
protected:
	//virtual ~TProgressDialog();
	//afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);

	// Generated message map functions
	//{{AFX_MSG(TProgressDialog)
	
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_TPROGRESSDIALOG_H__2A024480_F2C0_11D2_9F31_00C0F02A5F5D__INCLUDED_)
